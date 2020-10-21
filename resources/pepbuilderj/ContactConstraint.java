/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author yamule
 */
public class ContactConstraint extends ScoringFunction_Global{
	
	//アップデートしたところだけ計算しようとしてこのような設計なのか
	//ArrayList の index が Atom に対応している。キーがコンタクトする Atom 
	//あんまり使う予定ないので、全部ArrayList の方が描き易くていいのでは
	ArrayList<ArrayList<ContactTargetInfo>> contacts = new ArrayList<>();
	ArrayList<AtomInLattice> atoms = new ArrayList<>();
	
	
	double longestDist;
	int latticeDist;
	
	
	
	public double calcScore(){
		double ret = 0;
		//動いたものだけマークしてそれだけスコア更新するか
		for(int ii = 0;ii < contacts.size();ii++){
			AtomInLattice a = atoms.get(ii);
			ArrayList<ContactTargetInfo> lis = contacts.get(ii);
			for(int jj = 0;jj < lis.size();jj++){
				ret += lis.get(jj).calcScore(atoms.get(ii));
			}
		}
		return ret;
	}
	
			
	public void load(LatticeWorld lw,String filename){
		//<Chain><AA1letter><AAposition><PSBAtomcode>\t<Chain><AA1letter><AAposition><PSBAtomcode>\t<distance>\t<positiveeffectscore0-1.0>
		HashMap<String,Integer> atommap = new HashMap<>();
		
		for(int ii = 0;ii < lw.atoms.size();ii++){
			AtomInLattice a = lw.atoms.get(ii);
			if(a.index != ii){
				throw new RuntimeException("Atom Indices are broken!");
			}
			atommap.put(getAtomCode(a.pdbAtom),ii);
			contacts.add(new ArrayList<ContactTargetInfo>());
		}
		BufferedReader br = null;
		try{
			br = new  BufferedReader(new FileReader(new File(filename)));
			String ss;
			Pattern pat = Pattern.compile("^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\r\n\t]+)");
			HashSet<String> nodataatom = new HashSet<>();
			int okcount = 0;
			while((ss = br.readLine()) != null){
				Matcher mat = pat.matcher(ss);
				if(mat.find()){
					String code1 = mat.group(1);
					String code2 = mat.group(2);
					Double dist = Double.parseDouble(mat.group(3));
					Double score = Double.parseDouble(mat.group(4));
					boolean bflag = true;
					if(!atommap.containsKey(code1)){
						nodataatom.add(code1);
						bflag = false;
					}
					if(!atommap.containsKey(code2)){
						nodataatom.add(code2);
						bflag = false;
					}
					if(!bflag){
						continue;
					}
					contacts.get(atommap.get(code1)).add(new ContactTargetInfo(lw.atoms.get(atommap.get(code2)),dist,score)
					);
					contacts.get(atommap.get(code1)).add(new ContactTargetInfo(lw.atoms.get(atommap.get(code2)),dist,score));
					okcount++;
				}
			}
			for(String s:nodataatom){
				System.err.println(ss+" not found.");
			}
			System.err.println(okcount+" pairs loaded.");
			br.close();
		}catch(IOException exx){
			exx.printStackTrace();
		}catch(Exception exx){
			exx.printStackTrace();
			try{
				br.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}
	}
	
	
	public static String getAtomCode(PDBAtom a){
		StringBuffer ret = new StringBuffer();
		ret.append(a.parent.parent.name);
		if(PDBResidue.aaMap.containsKey(a.parent.getName())){
			ret.append(PDBResidue.aaMap.get(a.parent.getName()));
		}else{
			ret.append("X");
		}
		ret.append(a.parent.getResidueNumber());
		ret.append(a.parent.getInsertionCode());
		ret.append(a.pdb_atom_code);
		return ret.toString();
		
	}
	
	
}
class ContactTargetInfo{
	AtomInLattice target;
	double cutOffDistance;
	double score;
	ContactTargetInfo(AtomInLattice a,double d,double s){
		target = a;
		cutOffDistance = d;
		score = s;
	}
	public double calcScore(AtomInLattice a){//シグモイド関数とかにしてもいいかも
		if(cutOffDistance > a.distance(target)){
			return score;
		}
		return 0.0;
	}
}
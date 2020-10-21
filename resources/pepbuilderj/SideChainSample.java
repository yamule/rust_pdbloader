/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class SideChainSample implements Sampled{
	String residueName = "UNK";
	double prob = 0.0;
	int count = 0;
	public static Pattern labelPat = Pattern.compile("###[\\s]*([^\\s]+)[\\s]+([^\\s]+)");
	HashMap<String,PDBAtom> atoms_all = new HashMap<>();
	HashMap<String,PDBAtom> atoms_scoring = new HashMap<>(); //J48 で使うための
	
	
	public static SideChainSample parseBlock(ArrayList<String> block){
		SideChainSample ret = new SideChainSample();
		for(String ss:block){
			Matcher matt = labelPat.matcher(ss);
			if(matt.find()){
				String g = matt.group(1);
				if(g.equalsIgnoreCase("residue")){
					ret.residueName = matt.group(2);
				}else if(g.equalsIgnoreCase("count")){
					ret.count = Integer.parseInt(matt.group(2));
				}else if(g.indexOf("atom") == 0){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					String atomname = mapp.get("name");
					PDBAtom atom = new PDBAtom();
					ret.atoms_all.put(atomname,atom);
					ret.atoms_scoring.put(atomname,atom);
					atom.loc.set(
							Double.parseDouble(mapp.get("x")),
							Double.parseDouble(mapp.get("y")),
							Double.parseDouble(mapp.get("z"))

					);
					atom.pdb_atom_code = mapp.get("name");
					atom.atom_code = mapp.get("name").substring(0,1);
					
				}
			}
		}
		return ret;
	}
	public static ArrayList<SideChainSample> load(InputStream iss){
		ArrayList<SideChainSample> ret = new ArrayList<>();
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(iss));
			String ln = null;
			ArrayList<String> buff = new ArrayList<>();
			while((ln = br.readLine()) != null){
				if(ln.indexOf("//") == 0){
					ret.add(parseBlock(buff));
					buff.clear();
				}else{
					buff.add(ln);
				}
			}
			br.close();
			int count_all = 0;
			for(SideChainSample bs:ret){
				count_all += bs.count;
			}
			if(count_all > 0){
				for(SideChainSample bs:ret){
					bs.prob = bs.count/(double)count_all;
				}
			}
			
			Collections.sort(ret,new SampleComparator());
			Collections.reverse(ret);
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	public static void main(String[] args){
		InputStream iss = BackBoneSample.class.getResourceAsStream("resources/sampledresidues/ALA.rotamers.dat");
		if(iss == null){
			System.out.println("//");
		}
		load(iss);
	}
	public double getProb(){
		return prob;
	}
	
	
}

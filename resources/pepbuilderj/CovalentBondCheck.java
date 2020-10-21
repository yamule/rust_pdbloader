/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 *
 * @author yamule
 */
public class CovalentBondCheck {
//C:\dummy\vbox_share\crystal_contact_data\biogical\Ponstingl
//というフォルダ以下に合ったもの
//多分何かのテストデータ
//典型的な原子間の距離をとる。
	public static void main(String[] args){
		ArrayList<File> al = FileExtractionTest.getAllFiles("C:\\dummy\\work\\lattice");
		HashMap<String,Integer> closecount = new HashMap<>();
		HashMap<String,Integer> allcount = new HashMap<>();
		HashMap<String,ArrayList<Double>> distancedistribution = new HashMap<>();
		
		for(File f:al){
			if(!Pattern.compile("\\.pdb$").matcher(f.getPath()).find()){
				continue;
			}
			//System.out.println(f.getPath());
			PDBData p = PDBData.loadPDBFile(f.getPath());
			for(String c:p.chains.keySet()){
				PDBChain cc = p.chains.get(c);
				for(PDBResidue r:cc.residues){
					if(r.isValidAA()){
						if(r.alternativeCode.length() == 0 || r.alternativeCode.equals("A")){
							for(int ii = 0;ii < r.atoms.size();ii++){
								PDBAtom a = r.atoms.get(ii);
								if(a.isAlternative()){
									continue;
								}
								if(a.pdb_atom_code.indexOf("H") == 0){
									continue;
								}
								if(a.pdb_atom_code.equals("OXT")){
									continue;
								}
								for(int jj = ii+1;jj < r.atoms.size();jj++){
									PDBAtom b = r.atoms.get(jj);
									if(b.isAlternative()){
										continue;
									}
									
									if(b.pdb_atom_code.indexOf("H") == 0){
										continue;
									}
									if(b.pdb_atom_code.equals("OXT")){
										continue;
									}
									String code = r.getName()+":"+a.pdb_atom_code+"_"+b.pdb_atom_code;
									String code2 = r.getName()+":"+b.pdb_atom_code+"_"+a.pdb_atom_code;
									double dist = a.distance(b);
									
									if(!distancedistribution.containsKey(code)){
										distancedistribution.put(code,new ArrayList<Double>());
									}
									if(!distancedistribution.containsKey(code2)){
										distancedistribution.put(code2,new ArrayList<Double>());
									}
									distancedistribution.get(code).add(dist);
									distancedistribution.get(code2).add(dist);
								}	
							}
						}
					}
				}
			}
		}
		
		for(String kk:distancedistribution.keySet()){
			Collections.sort(distancedistribution.get(kk));
			ArrayList<Double> dal = distancedistribution.get(kk);
			int siz = dal.size();
			int fivepercent = (int)(siz*0.05);
			for(int ii = 0;ii < fivepercent;ii++){
				dal.remove(0);
				dal.remove(dal.size()-1);
			}
		}
		System.out.println("static HashMap<String,Double> lbound = new HashMap<>();");
		System.out.println("static HashMap<String,Double> ubound = new HashMap<>();");
		System.out.println("static{");
		Pattern opat = Pattern.compile("\\:O_");
		Pattern opat2 = Pattern.compile("_O$");
		for(String kk:distancedistribution.keySet()){
			ArrayList<Double> dal = distancedistribution.get(kk);
			System.out.println("lbound.put(\""+kk+"\","+dal.get(0)+");");
			System.out.println("ubound.put(\""+kk+"\","+dal.get(dal.size()-1)+");");
			if(kk.indexOf(":O_") > -1){
				String sk = kk.replaceAll("O_","OXT_");
				System.out.println("lbound.put(\""+sk+"\","+dal.get(0)+");");
				System.out.println("ubound.put(\""+sk+"\","+dal.get(dal.size()-1)+");");
			}
			if((kk+";").indexOf("_O;") > -1){
				String sk = (kk+";").replaceAll("_O;","_OXT");
				System.out.println("lbound.put(\""+sk+"\","+dal.get(0)+");");
				System.out.println("ubound.put(\""+sk+"\","+dal.get(dal.size()-1)+");");
			}
			
		}
		System.out.println("};");
	}
}

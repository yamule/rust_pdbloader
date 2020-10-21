/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class MyWorld {
	public static HashMap<String,String> toOneLetter = new HashMap<>();
	public static HashMap<String,String> toThreeLetter = new HashMap<>();

	static{
		toOneLetter.put("ALA","A");
		toOneLetter.put("ARG","R");
		toOneLetter.put("ASN","N");
		toOneLetter.put("ASP","D");
		toOneLetter.put("CYS","C");
		toOneLetter.put("GLN","Q");
		toOneLetter.put("GLU","E");
		toOneLetter.put("GLY","G");
		toOneLetter.put("HIS","H");
		toOneLetter.put("ILE","I");
		toOneLetter.put("LEU","L");
		toOneLetter.put("LYS","K");
		toOneLetter.put("MET","M");
		toOneLetter.put("PHE","F");
		toOneLetter.put("PRO","P");
		toOneLetter.put("SER","S");
		toOneLetter.put("THR","T");
		toOneLetter.put("TRP","W");
		toOneLetter.put("TYR","Y");
		toOneLetter.put("VAL","V");
		toOneLetter.put("SEC","U");
		toOneLetter.put("UNK","X");
		for(String s:toOneLetter.keySet()){
			toThreeLetter.put(toOneLetter.get(s),s);
		}
	}
	
	
	public static HashSet<String> backbone_atoms = new HashSet<>();
	public static HashMap<String,HashSet<String>> sidechain_atoms = new HashMap<>();
	public static HashMap<String,HashSet<String>> all_atoms = new HashMap<>();
	static{
		backbone_atoms.add("N");
		backbone_atoms.add("C");
		backbone_atoms.add("O");
		backbone_atoms.add("CA");
		String[][] atomsets_ = {
		{"#VAL","C","CG1","N","CA","O","CB","CG2",},
			{"#SER","C","OG","N","CA","O","CB",},
			{"#ILE","CD1","C","CG1","N","CA","O","CB","CG2",},
			{"#LYS","CD","CE","C","CG","NZ","N","CA","O","CB",},
			{"#GLN","CD","C","CG","OE1","NE2","N","CA","O","CB",},
			{"#PRO","CD","C","CG","N","CA","O","CB",},
			{"#PHE","CD2","CD1","CE2","C","CG","CZ","N","CA","CE1","O","CB",},
			{"#TYR","CD2","CD1","CE2","C","CG","CZ","OH","N","CA","CE1","O","CB",},
			{"#GLU","CD","C","CG","OE1","OE2","N","CA","O","CB",},
			{"#TRP","C","CG","CH2","N","O","CD2","CE3","CD1","CE2","CZ2","NE1","CZ3","CA","CB",},
			{"#HIS","CD2","C","CG","ND1","NE2","N","CA","CE1","O","CB",},
			{"#GLY","C","N","CA","O",},
			{"#ARG","CD","C","CG","NH1","NE","CZ","NH2","N","CA","O","CB",},
			{"#ALA","C","N","CA","O","CB",},
			{"#CYS","C","SG","N","CA","O","CB",},
			{"#ASN","C","CG","OD1","ND2","N","CA","O","CB",},
			{"#LEU","CD2","CD1","C","CG","N","CA","O","CB",},
			{"#MET","SD","CE","C","CG","N","CA","O","CB",},
			{"#ASP","C","CG","OD2","OD1","N","CA","O","CB",},
			{"#THR","C","OG1","N","CA","O","CB","CG2"}
		};
		
		for(int ii = 0;ii < atomsets_.length;ii++){
			String code = atomsets_[ii][0].replaceFirst("#","");
			
			sidechain_atoms.put(code,new HashSet<String>());
			all_atoms.put(code,new HashSet<String>());
			for(int jj  = 1;jj < atomsets_[ii].length;jj++){
				if(!backbone_atoms.contains(atomsets_[ii][jj])){
					sidechain_atoms.get(code).add(atomsets_[ii][jj]);
				}
				all_atoms.get(code).add(atomsets_[ii][jj]);
			}
		}
	}
	


	public static HashMap<String,String> lineToHash(String str){
		Pattern pat = Pattern.compile("^([^\\:]+)\\:(.+)");
		String[] ppt = str.split("\t");
		HashMap<String,String> ret = new HashMap<>();
		for(String p:ppt){
			Matcher mat = pat.matcher(p);
			if(mat.find()){
				ret.put(mat.group(1),mat.group(2));
			}
		}
		return ret;
	}
		
}

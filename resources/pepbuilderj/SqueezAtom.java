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
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 *
 * @author yamule
 */
public class SqueezAtom {
//典型的な原子間の距離をとる。
//C:\dummy\vbox_share\crystal_contact_data\biogical\Ponstingl
//というフォルダ以下に合ったもの
//多分何かのテストデータ
	public static void main(String[] args){
		ArrayList<File> al = FileExtractionTest.getAllFiles("C:\\dummy\\work\\lattice");
		HashMap<String,PDBAtom> lines = new HashMap<>();
		
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
								a.loc.x = 0.0;
								a.loc.y = 0.0;
								a.loc.z = 0.0;
								a.bfactor = 25;
								a.occupancy = 1.0;
								a.atom_number = "1";
								r.setAlternativeCode("");
								r.setInsertionCode("");
								r.setResidueNumber(0);
								cc.name = "A";
								lines.put(r.getName()+a.makePDBString(),a);
							}
						}
					}
				}
			}
		}
		ArrayList<String> sal = new ArrayList<>(lines.keySet());
		Collections.sort(sal);
		for(String s:sal){
			PDBAtom a = lines.get(s);
			
			System.out.println("\""+a.makePDBString()+"\",");
			if(a.pdb_atom_code.equals("O")){
				a.pdb_atom_code = "OXT";
				System.out.println("\""+a.makePDBString()+"\",");
				a.pdb_atom_code = "O";
			}
		}
	}
}

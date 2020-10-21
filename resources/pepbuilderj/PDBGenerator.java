/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author yamule
 */
public class PDBGenerator {
	static final String allatomlines[] = {
"ATOM      1    C ALA A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N ALA A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O ALA A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT ALA A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA ALA A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB ALA A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N ARG A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O ARG A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT ARG A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CD ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CZ ARG A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   NE ARG A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1  NH1 ARG A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1  NH2 ARG A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    C ASN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N ASN A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O ASN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT ASN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA ASN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB ASN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG ASN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  ND2 ASN A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1  OD1 ASN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C ASP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N ASP A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O ASP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT ASP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA ASP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB ASP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG ASP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  OD1 ASP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OD2 ASP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C CYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N CYS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O CYS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT CYS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA CYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB CYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   SG CYS A   0       0.000   0.000   0.000  1.00 25.00           S  ",
"ATOM      1    C GLN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N GLN A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O GLN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT GLN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA GLN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB GLN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CD GLN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG GLN A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  NE2 GLN A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1  OE1 GLN A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C GLU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N GLU A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O GLU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT GLU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA GLU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB GLU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CD GLU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG GLU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  OE1 GLU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OE2 GLU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C GLY A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N GLY A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O GLY A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT GLY A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA GLY A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N HIS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O HIS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT HIS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD2 HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE1 HIS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  ND1 HIS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1  NE2 HIS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    C ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N ILE A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O ILE A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT ILE A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD1 ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CG1 ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CG2 ILE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N LEU A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O LEU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT LEU A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD1 LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD2 LEU A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N LYS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O LYS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT LYS A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CD LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CE LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG LYS A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   NZ LYS A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    C MET A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N MET A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O MET A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT MET A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA MET A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB MET A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CE MET A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG MET A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   SD MET A   0       0.000   0.000   0.000  1.00 25.00           S  ",
"ATOM      1    C PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N PHE A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O PHE A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT PHE A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CZ PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD1 PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD2 PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE1 PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE2 PHE A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C PRO A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N PRO A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O PRO A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT PRO A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA PRO A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB PRO A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CD PRO A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG PRO A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C SER A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N SER A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O SER A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT SER A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA SER A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB SER A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   OG SER A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C THR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N THR A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O THR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT THR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA THR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB THR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CG2 THR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  OG1 THR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1    C TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N TRP A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O TRP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT TRP A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD1 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD2 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE2 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE3 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CH2 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CZ2 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CZ3 TRP A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  NE1 TRP A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    C TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N TYR A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O TYR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT TYR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CG TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CZ TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   OH TYR A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  CD1 TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CD2 TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE1 TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CE2 TYR A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    C VAL A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1    N VAL A   0       0.000   0.000   0.000  1.00 25.00           N  ",
"ATOM      1    O VAL A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1  OXT VAL A   0       0.000   0.000   0.000  1.00 25.00           O  ",
"ATOM      1   CA VAL A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1   CB VAL A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CG1 VAL A   0       0.000   0.000   0.000  1.00 25.00           C  ",
"ATOM      1  CG2 VAL A   0       0.000   0.000   0.000  1.00 25.00           C  ",
};
	static final HashMap<String,ArrayList<PDBAtom>> allAtoms = new HashMap<>(); 
	static final HashMap<String,String> aaMap_rev = new HashMap<>();
	static{
		for(String s:allatomlines){
			String residue_name = s.substring(17,20).toUpperCase();
			PDBAtom a = generateAtomFromLine(s);
			if(!allAtoms.containsKey(residue_name)){
				allAtoms.put(residue_name,new ArrayList<PDBAtom>());
			}
			allAtoms.get(residue_name).add(a);
		}
		for(String kk:PDBResidue.aaMap.keySet()){
			if(PDBResidue.validAASet.contains(kk)){
				aaMap_rev.put(PDBResidue.aaMap.get(kk),kk);
			}
		}
	}
	public static PDBChain generateChain(String aacode,String chainid){
		PDBChain ret = new PDBChain(chainid);
		
		String[] ppt = aacode.replaceAll("[\\s]","").split("");
		int rescount = 0;
		int atomcount = 0;
		ArrayList<String> aaz = new ArrayList<>();
		for(String ss:ppt){
			if(ss.length() > 0){
				if(aaMap_rev.containsKey(ss)){
					aaz.add(ss);
				}else{
					System.err.println(ss+" was not parsed.");
				}
			}
		}
		for(int ii = 0;ii < aaz.size();ii++){
			String ss = aaz.get(ii);
			if(ss.length() > 0){
				if(aaMap_rev.containsKey(ss)){
					rescount ++;
					PDBResidue res = new PDBResidue();
					ret.residues.add(res);
					res.setParent(ret);
					
					res.setResidueNumber(rescount);
					res.setName(aaMap_rev.get(ss));
					ArrayList<PDBAtom> a = allAtoms.get(aaMap_rev.get(ss));
					for(PDBAtom atom:a){
						if(atom.pdb_atom_code.equals("OXT")){
							if(ii != aaz.size()-1){
								continue;
							}
						}
						atomcount++;
						
						PDBAtom ac = atom.getCopy();
						res.addAtom(ac);
						ac.atom_number = String.valueOf(atomcount);
					}
				}
			}
		}
		return ret;
	}
	public static PDBAtom generateAtomFromLine(String line){
		String record_name = line.substring(0,6);
		String atom_serial_num = line.substring(6,11);
		String atom_name = line.substring(12,16);
		String alternate_location = line.substring(16,17);
		String residue_name = line.substring(17,20).toUpperCase();
		String chain_id = line.substring(21,22);
		String residue_sequence_number = line.substring(22,26);
		String insertion_code = line.substring(26,27);
		String x_coord = line.substring(30,38);
		String y_coord = line.substring(38,46);
		String z_coord = line.substring(46,54);
		String occupancy = "";
		String temperature_factor = "";
		if(line.length() >= 60){
			occupancy = line.substring(54,60);
		}
		if(line.length() >= 61){
			temperature_factor = line.substring(60,Math.min(line.length(),66));
		}
		String element_symbol = "";
		if(line.length() >= 78){
			 element_symbol = line.substring(76,78).replaceAll(" ","");;
		}
		//String charge = line.substring(78,80);
		if(element_symbol.length() > 0){
		}else{
			element_symbol = atom_name.substring(0,2).replaceAll(" ","");
		}
		PDBAtom newp= new PDBAtom(atom_name
				,element_symbol
				,atom_serial_num
				,new Point3D(Float.parseFloat(x_coord)
				,Float.parseFloat(y_coord)
				,Float.parseFloat(z_coord)));
		newp.setLine(line.replaceAll("[\\r\\n]","")+"\n");
		if(temperature_factor.length() > 0){
			newp.setBFactor(Double.parseDouble(temperature_factor.replaceAll("[\\s]","")));
		}
		if(occupancy.length() > 0){
			newp.setOccupancy(Double.parseDouble(occupancy.replaceAll("[\\s]","")));
		}
		if( line.indexOf("HETATM") == 0){
			newp.setHETATM(true);
		}
		int residue_number = -1000;
		String altflag = "";
		try{
			residue_number  = Integer.parseInt(residue_sequence_number.replaceAll(" ",""));
		}catch(Exception exx){
			residue_number  = Integer.parseInt(residue_sequence_number.replaceAll("[^0-9\\-]",""));
			altflag = residue_sequence_number.replaceAll(" ","");
		}
		if(!alternate_location.equals(" ")){
			newp.setAltCode(alternate_location);
		}
		return newp;
	}
	public static void main(String[] args){
		PDBChain c = generateChain("KHGFMNRTRYTIREWERTYIIKHGFUIYR","A");
		for(PDBResidue r:c.residues){
			for(PDBAtom a:r.atoms){
				System.out.println(a.makePDBString());
			}
		}
	}
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author yamule
 */
public class BasicParamSet {
	static HashSet<String> covalentBondSet = new HashSet<>();
	static String[] covalentBond = {
"ALA:CA_C", "ALA:CA_CB", "ALA:CA_N", "ALA:CB_CA", "ALA:C_CA", "ALA:C_O", "ALA:N_CA", "ALA:O_C", "ARG:CA_C", "ARG:CA_CB", "ARG:CA_N", "ARG:CB_CA", "ARG:CB_CG"
, "ARG:CD_CG", "ARG:CD_NE", "ARG:CG_CB", "ARG:CG_CD", "ARG:CZ_NE", "ARG:CZ_NH1", "ARG:CZ_NH2", "ARG:C_CA", "ARG:C_O", "ARG:NE_CD", "ARG:NE_CZ", "ARG:NH1_CZ"
, "ARG:NH2_CZ", "ARG:N_CA", "ARG:O_C", "ASN:CA_C", "ASN:CA_CB", "ASN:CA_N", "ASN:CB_CA", "ASN:CB_CG", "ASN:CG_CB", "ASN:CG_ND2", "ASN:CG_OD1", "ASN:C_CA", "ASN:C_O"
, "ASN:ND2_CG", "ASN:N_CA", "ASN:OD1_CG", "ASN:O_C", "ASP:CA_C", "ASP:CA_CB", "ASP:CA_N", "ASP:CB_CA", "ASP:CB_CG", "ASP:CG_CB", "ASP:CG_OD1", "ASP:CG_OD2"
, "ASP:C_CA", "ASP:C_O", "ASP:N_CA", "ASP:OD1_CG", "ASP:OD2_CG", "ASP:O_C", "CYS:CA_C", "CYS:CA_CB", "CYS:CA_N", "CYS:CB_CA", "CYS:CB_SG", "CYS:C_CA", "CYS:C_O"
, "CYS:N_CA", "CYS:O_C", "CYS:SG_CB", "GLN:CA_C", "GLN:CA_CB", "GLN:CA_N", "GLN:CB_CA", "GLN:CB_CG", "GLN:CD_CG", "GLN:CD_NE2", "GLN:CD_OE1", "GLN:CG_CB", "GLN:CG_CD"
, "GLN:C_CA", "GLN:C_O", "GLN:NE2_CD", "GLN:N_CA", "GLN:OE1_CD", "GLN:O_C", "GLU:CA_C", "GLU:CA_CB", "GLU:CA_N", "GLU:CB_CA", "GLU:CB_CG", "GLU:CD_CG", "GLU:CD_OE1"
, "GLU:CD_OE2", "GLU:CG_CB", "GLU:CG_CD", "GLU:C_CA", "GLU:C_O", "GLU:N_CA", "GLU:OE1_CD", "GLU:OE2_CD", "GLU:O_C", "GLY:CA_C", "GLY:CA_N", "GLY:C_CA", "GLY:C_O"
, "GLY:N_CA", "GLY:O_C", "HIS:CA_C", "HIS:CA_CB", "HIS:CA_N", "HIS:CB_CA", "HIS:CB_CG", "HIS:CD2_CG", "HIS:CD2_NE2", "HIS:CE1_ND1", "HIS:CE1_NE2", "HIS:CG_CB"
, "HIS:CG_CD2", "HIS:CG_ND1", "HIS:C_CA", "HIS:C_O", "HIS:ND1_CE1", "HIS:ND1_CG", "HIS:NE2_CD2", "HIS:NE2_CE1", "HIS:N_CA", "HIS:O_C", "ILE:CA_C", "ILE:CA_CB"
, "ILE:CA_N", "ILE:CB_CA", "ILE:CB_CG1", "ILE:CB_CG2", "ILE:CD1_CG1", "ILE:CG1_CB", "ILE:CG1_CD1", "ILE:CG2_CB", "ILE:C_CA", "ILE:C_O", "ILE:N_CA", "ILE:O_C", "LEU:CA_C"
, "LEU:CA_CB", "LEU:CA_N", "LEU:CB_CA", "LEU:CB_CG", "LEU:CD1_CG", "LEU:CD2_CG", "LEU:CG_CB", "LEU:CG_CD1", "LEU:CG_CD2", "LEU:C_CA", "LEU:C_O", "LEU:N_CA", "LEU:O_C"
, "LYS:CA_C", "LYS:CA_CB", "LYS:CA_N", "LYS:CB_CA", "LYS:CB_CG", "LYS:CD_CE", "LYS:CD_CG", "LYS:CE_CD", "LYS:CE_NZ", "LYS:CG_CB", "LYS:CG_CD", "LYS:C_CA", "LYS:C_O"
, "LYS:NZ_CE", "LYS:N_CA", "LYS:O_C", "MET:CA_C", "MET:CA_CB", "MET:CA_N", "MET:CB_CA", "MET:CB_CG", "MET:CE_SD", "MET:CG_CB", "MET:CG_SD", "MET:C_CA", "MET:C_O", "MET:N_CA"
, "MET:O_C", "MET:SD_CE", "MET:SD_CG", "PHE:CA_C", "PHE:CA_CB", "PHE:CA_N", "PHE:CB_CA", "PHE:CB_CG", "PHE:CD1_CE1", "PHE:CD1_CG", "PHE:CD2_CE2", "PHE:CD2_CG", "PHE:CE1_CD1"
, "PHE:CE1_CZ", "PHE:CE2_CD2", "PHE:CE2_CZ", "PHE:CG_CB", "PHE:CG_CD1", "PHE:CG_CD2", "PHE:CZ_CE1", "PHE:CZ_CE2", "PHE:C_CA", "PHE:C_O", "PHE:N_CA", "PHE:O_C", "PRO:CA_C"
, "PRO:CA_CB", "PRO:CA_N", "PRO:CB_CA", "PRO:CB_CG", "PRO:CD_CG", "PRO:CD_N", "PRO:CG_CB", "PRO:CG_CD", "PRO:C_CA", "PRO:C_O", "PRO:N_CA", "PRO:N_CD", "PRO:O_C", "SER:CA_C"
, "SER:CA_CB", "SER:CA_N", "SER:CB_CA", "SER:CB_OG", "SER:C_CA", "SER:C_O", "SER:N_CA", "SER:OG_CB", "SER:O_C", "THR:CA_C", "THR:CA_CB", "THR:CA_N", "THR:CB_CA", "THR:CB_CG2"
, "THR:CB_OG1", "THR:CG2_CB", "THR:C_CA", "THR:C_O", "THR:N_CA", "THR:OG1_CB", "THR:O_C", "TRP:CA_C", "TRP:CA_CB", "TRP:CA_N", "TRP:CB_CA", "TRP:CB_CG", "TRP:CD1_CG"
, "TRP:CD1_NE1", "TRP:CD2_CE2", "TRP:CD2_CE3", "TRP:CD2_CG", "TRP:CE2_CD2", "TRP:CE2_CZ2", "TRP:CE2_NE1", "TRP:CE3_CD2", "TRP:CE3_CZ3", "TRP:CG_CB", "TRP:CG_CD1", "TRP:CG_CD2"
, "TRP:CH2_CZ2", "TRP:CH2_CZ3", "TRP:CZ2_CE2", "TRP:CZ2_CH2", "TRP:CZ3_CE3", "TRP:CZ3_CH2", "TRP:C_CA", "TRP:C_O", "TRP:NE1_CD1", "TRP:NE1_CE2", "TRP:N_CA", "TRP:O_C"
, "TYR:CA_C", "TYR:CA_CB", "TYR:CA_N", "TYR:CB_CA", "TYR:CB_CG", "TYR:CD1_CE1", "TYR:CD1_CG", "TYR:CD2_CE2", "TYR:CD2_CG", "TYR:CE1_CD1", "TYR:CE1_CZ", "TYR:CE2_CD2"
, "TYR:CE2_CZ", "TYR:CG_CB", "TYR:CG_CD1", "TYR:CG_CD2", "TYR:CZ_CE1", "TYR:CZ_CE2", "TYR:CZ_OH", "TYR:C_CA", "TYR:C_O", "TYR:N_CA", "TYR:OH_CZ", "TYR:O_C", "VAL:CA_C"
, "VAL:CA_CB", "VAL:CA_N", "VAL:CB_CA", "VAL:CB_CG1", "VAL:CB_CG2", "VAL:CG1_CB", "VAL:CG2_CB", "VAL:C_CA", "VAL:C_O", "VAL:N_CA", "VAL:O_C" 	};
	
	
	
	 
	
	
	static{
		for(String c:covalentBond){
			covalentBondSet.add(c);
		}
	}
	BasicParamSet(){
		
		
	}
	//ここから
	//同一残基内の共有結合をとる。
	//途中を削除されても良いようにConstraintをかける。
	//完全な場合は共有結合ペアだけ考える?
	//public void addIntraResidueConstraint(ArrayList<AtomInLattice> al){
	//	HashMap<PDBResidue,ArrayList<AtomInLattice>> 
///		
//		
//		
//	}
	
	
	public void addPeptideBond(ArrayList<AtomInLattice> al){
		PDBChain pc = al.get(0).pdbAtom.parent.parent;
		HashMap<PDBAtom,AtomInLattice> mapp = new HashMap<>();
		for(AtomInLattice a:al){
			mapp.put(a.pdbAtom, a);
		}
		ArrayList<AtomInLattice> peptide = new ArrayList<>();
		for(int ii = 0;ii < pc.residues.size();ii++){
			PDBResidue cur = pc.residues.get(ii);
			PDBResidue nex = null;
			PDBResidue prev = null;
			if(ii > 0){
				prev = pc.residues.get(ii-1);
				if(mapp.containsKey(prev.getC())){
					if(mapp.containsKey(cur.getN())){
						peptide.add(mapp.get(prev.getC()));
						peptide.add(mapp.get(cur.getN()));
					}
				}
			}
		}
		for(int ii = 0;ii < peptide.size();ii+=2){
			peptide.get(ii).addCovalentBond(peptide.get(ii+1));
			peptide.get(ii+1).addCovalentBond(peptide.get(ii));
		}
	}
	
	
	
	
}

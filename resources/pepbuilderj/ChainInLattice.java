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
public class ChainInLattice {
	PDBChain pdbChain;
	HashMap<PDBResidue,ArrayList<AtomInLattice>> residueAtomMap = new HashMap<>();
	ChainInLattice(LatticeWorld lw,PDBChain chain){
		for(PDBResidue res:chain.residues){
			ArrayList<AtomInLattice> al = new ArrayList<>();
			for(PDBAtom a:res.atoms){
				AtomInLattice atm = new AtomInLattice(lw,a);
				al.add(atm);
				lw.addAtom(atm);
			}
			residueAtomMap.put(res, al);
		}
	}
}

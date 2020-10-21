/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author kimidori
 */
public class SideChainSet {
	HashMap<String,ArrayList<SideChainSample>> sidechains = new HashMap<>();
	String rotamerDir = "resources/sampledresidues/";
	
	SideChainSet(){
		String[] aaname = {
			"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
			,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
		};
		for(String a:aaname){
			sidechains.put(a,SideChainSample.load(
					ChainBuilder.class.getResourceAsStream(rotamerDir+a+".rotamers.dat")));
			
		}
	}
	
	public SideChainSample getRandomly(String sname){
		
		ArrayList<SideChainSample> sal = sidechains.get(sname);
		return sal.get((int)(Math.random()*sal.size()));
		
	}
	public SideChainSample getNext(FloatingResidue f){
		ArrayList<SideChainSample> sal = sidechains.get(f.getName());
		if(f.sidechainIndex >= sal.size()-1){
			f.sidechainIndex = -1;
		}
		f.sidechainIndex++;
		return sal.get(f.sidechainIndex);
	}
	
	
}

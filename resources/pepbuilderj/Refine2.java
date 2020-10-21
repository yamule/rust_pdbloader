/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import static pepbuilderj.ResidueRefine.parseArgs;
import static pepbuilderj.ResidueRefine.parseTargetFile;

/**
 *
 * @author kimidori
 */
public class Refine2 {
	
	public static void main(String[] args){
		//HashMap<String,String> hs = parseArgs(args);
		//String infile = hs.get("-in");
		//String outfile = hs.get("-out");
		String infile = "C:\\dummy\\vbox_share\\casp13\\queries\\T0979\\model0.pfas.pdb";
		String outfile = "C:\\dummy\\vbox_share\\casp13\\queries\\T0979\\model1_re.pdb";
		PDBData pdb = PDBData.loadPDBFile(infile);
		ArrayList<FloatingResidue> fl = new ArrayList<>();
		ArrayList<FloatingResidue> targets = new ArrayList<>();
		
		ChainBuilder cb = new ChainBuilder();
		for(String c:pdb.chains.keySet()){
			ArrayList<PDBResidue> rr = pdb.chains.get(c).residues;
			ArrayList<PDBResidue> ch = new ArrayList<>();
			for(int ii = 0;ii < rr.size();ii++){
				if(rr.get(ii).isMissing() || rr.get(ii).isLigand()){
					ch.add(null);
				}else{
					ch.add(rr.get(ii));
				}
			}
			ArrayList<FloatingResidue> fr = new ArrayList<>();
			for(int ii = 0;ii < rr.size();ii++){
				PDBResidue prev = null;
				PDBResidue next = null;
				if(ii > 0){
					prev = rr.get(ii-1);
				}
				if(ii < rr.size()-1){
					next = rr.get(ii+1);
				}
				if(rr.get(ii) == null){
					fr.add(null);
					continue;
				}
				FloatingResidue ff = new FloatingResidue(rr.get(ii),prev,next);
				fr.add(ff);
				
			}
			for(int ii = 0;ii < fr.size();ii++){
				if(fr.get(ii) == null){
					continue;
				}
				if(ii > 0){
					fr.get(ii).setPrevFloating(fr.get(ii-1));
				}
				if(ii < rr.size()-1){
					fr.get(ii).setNextFloating(fr.get(ii+1));
				}
				fl.add(fr.get(ii));
				if(fr.get(ii).getResidueNumber() > 31 && fr.get(ii).getResidueNumber() < 74){
					targets.add(fr.get(ii));
				}
			}
		}
		ArrayList<PDBResidue> pp = new ArrayList<>();
		for(FloatingResidue f:targets){
			f.changeSideChain(cb.sidechains.getNext(f));
		}
		for(FloatingResidue f:fl){
			f.changeSideChain(cb.sidechains.getNext(f));
			pp.add(f);
		}
		cb.prepare(fl);
		
		for(int ii = 0;ii < 1000;ii++){
			cb.refiner.loopRefile(fl,targets, cb, 10,ii%10 == 9);
			ChainBuilder.saveChains(fl,outfile+"."+ii+".pdb");
		}
	}
}

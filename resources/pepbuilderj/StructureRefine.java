/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import static pepbuilderj.AtomDistancePenalty.TYPE_BACKBONE;
import static pepbuilderj.AtomDistancePenalty.TYPE_SIDECHAIN;

/**
 *
 * @author kimidori
 */
public class StructureRefine {
	ResidueRefine refiner = new ResidueRefine();
	ArrayList<FloatingResidue> residues = new ArrayList<>();
	ArrayList<ArrayList<FloatingResidue>> groups = new ArrayList<>();
	ArrayList<PDBResidue> residues_ = new ArrayList<>();
	ChainBuilder cb = new ChainBuilder();
	StructureRefine(ArrayList<PDBResidue> res){
		residues = ChainBuilder.changeToFloating(res);
		residues_.addAll(residues);
		refiner.prepare(residues_);
		cb.prepare(residues);
		for(FloatingResidue fr:residues){
			fr.changeSideChain(cb.sidechains.getRandomly(fr.getName()));
		}
	}
	public static double windowAverage(double[] d,int index,int windowsize){
		int count = 0;
		double sum = 0;
		for(int ii = -windowsize;ii <=windowsize;ii++){
			int pos = index+ii;
			if(pos*(d.length-pos-1) >= 0){
				sum += d[ii];
				count ++;
			}
		}
		if(count > 0){
			return sum/count;
		}
		return 0;
	}
	public void makeGroup(){
		ArrayList<PDBResidue> pl = new ArrayList<>();
		double pp = 0;
		double csum = 0;
		double[] rres = new double[residues.size()];
		double[] smoothed = new double[residues.size()];
		for(int ii = 0;ii < residues_.size();ii++){
			PDBResidue f = residues_.get(ii);
			double cc = cb.calcAtomCrash_Higher_is_better(residues_,f);
			csum += cc;
			rres[ii] = cc;
		}
		HashMap<PDBResidue,Double> crashscore = new HashMap<>();
		for(int ii = 0;ii < residues_.size();ii++){
			smoothed[ii] = windowAverage(rres,ii,3);
			crashscore.put(residues_.get(ii),smoothed[ii]);
		}
		HashSet<PDBResidue> rm = new HashSet<>();
		
		
		for(int ii = 0;ii < residues_.size();ii++){
			ArrayList<PDBResidue> p = getCrashWith(residues_,residues.get(ii));
			
			for(int jj = 0;jj < residues_.size();jj++){
				if(ii == jj){
					continue;
				}
			}
		}
		//ここから
	}
	
	

	public ArrayList<PDBResidue> getCrashWith(ArrayList<PDBResidue> allres
			,PDBResidue targetresidue){
		ArrayList<PDBResidue> ret = new ArrayList<>();
		for(int ii = 0;ii < allres.size();ii++){
			PDBResidue rr = allres.get(ii);
			if(rr == targetresidue){
				continue;
			}
			if(rr.getCA().distance(targetresidue.getCA()) > 20){
				continue;
			}
			for(PDBAtom a:rr.atoms){
				String acode = a.parent.getName()+":"+a.pdb_atom_code;
				if(a.parent == targetresidue){
					continue;
				}
				for(PDBAtom r:targetresidue.atoms){
					double dist = r.distance(a);
					if(dist > 4.0){
						continue;
					}
					if(cb.distPenalty.backboneNextInteraction_min.containsKey(r) && cb.distPenalty.backboneNextInteraction_min.get(r).containsKey(a)){

						if(cb.distPenalty.backboneNextInteraction_min.get(r).get(a) > (cb.distPenalty.trelance_ratio+1.0)*dist){
							
							ret.add(rr);
							break;
						}
					}else{
						
						if(dist < 1.8){// CYS_CYS が 2.0 くらいなので
							ret.add(rr);
							break;
						}
						String rcode = r.parent.getName()+":"+r.pdb_atom_code;
						if(cb.distPenalty.pairsMinDist.get(rcode+"_"+acode) > (cb.distPenalty.trelance_ratio+1.0)*dist){
							
							ret.add(rr);
							break;
						}
					}
				}
			}
		}
		return ret;
	}
}
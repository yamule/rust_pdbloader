/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static pepbuilderj.AtomDistancePenalty.BACKBONE_PREVC_N;
import static pepbuilderj.ResidueRefine.parseArgs;
import static pepbuilderj.ResidueRefine.parseTargetFile;

/**
 *
 * @author kimidori
 */
public class Refine3 {
	
	public static double calcScore(ChainBuilder cb
			,ArrayList<PDBResidue>residues_for_scoring
			,ArrayList<FloatingResidue> floating_residues_for_crash){
		double scoresum = 0;
		double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(residues_for_scoring, cb.scoring);
		for(double ss:scores){
			scoresum += ss;
		}
		double[] penal = new double[floating_residues_for_crash.size()];
		for(int r = 0;r < floating_residues_for_crash.size();r++){
			penal[r] = cb.distPenalty.calcCrashPenalty(residues_for_scoring,floating_residues_for_crash.get(r),AtomDistancePenalty.TYPE_BACKBONE);
			scoresum += penal[r];
		}
		return scoresum;
	}
	/**
	 * target residues はひとつながりになってないといけない。
	 * @param cb
	 * @param allresidues_in_environment
	 * @param allresidues_in_environment_r
	 * @param targetresidues
	 * @param targetresidues_r
	 * @param save
	 * @param iternum
	 * @return 
	 */
	public static double changeBackBone(ChainBuilder cb
			,ArrayList<FloatingResidue> allresidues_in_environment
			,ArrayList<PDBResidue> allresidues_in_environment_r
			,ArrayList<FloatingResidue> targetresidues,
			ArrayList<PDBResidue> targetresidues_r
			,boolean save
			,int iternum){
		double sssum = 0;
		double prevscore = -1000000;
		double firstscore = 0;
		if(save){
			firstscore = calcScore(cb,allresidues_in_environment_r,allresidues_in_environment);
			prevscore = firstscore - 10000000;
		}else{
			firstscore = calcScore(cb,allresidues_in_environment_r,targetresidues)/targetresidues.size();
			prevscore = -1000000;
		}
		
		for(FloatingResidue fr:targetresidues){
			fr.saveLoc();
		}
		for(int jj = 0;jj < iternum;jj++){
			double scoresum = 0;
			if(save){
				scoresum = calcScore(cb,allresidues_in_environment_r,allresidues_in_environment);
			}else{
				scoresum = calcScore(cb,allresidues_in_environment_r,targetresidues)/targetresidues.size()-firstscore;
			}
			sssum += Math.max(0,scoresum-firstscore);
			if(scoresum > prevscore){
				if(save){
					for(FloatingResidue fr:targetresidues){
						fr.saveLoc();
					}
				}
				prevscore = scoresum;
			}else{
				for(FloatingResidue fr:targetresidues){
					fr.restoreLoc();
				}
				
			}
			
			
			
			FloatingResidue last = targetresidues.get(targetresidues.size()-1);
			FloatingResidue first = targetresidues.get(0);
			
			if(first.prev == null){
				ArrayList<FloatingResidue> fresidue = new ArrayList<>(targetresidues);
				Collections.reverse(fresidue);
				for(int ff = 0;ff < fresidue.size();ff++){
					FloatingResidue fr = fresidue.get(ff);
					fr.nextOmega = Math.PI;
					fr.processBackBone(cb.backbones.getRandomly(fr),fr.nextOmega,false,true);
				}
			}else{
				for(int ff = 0;ff < targetresidues.size();ff++){
					FloatingResidue fr = targetresidues.get(ff);
					fr.prevOmega = Math.PI;
					fr.processBackBone(cb.backbones.getRandomly(fr));
				}
				if(last.next != null){
				
					int tcount = 0;
					int failcount = 0;
					double prevdist = Math.abs(last.getC().distance(last.next.getN())- cb.distPenalty.backbone_average[BACKBONE_PREVC_N]);
					//System.out.println(last.getC().distance(last.next.getN())+"###");
					while(prevdist > 0.3){
						int targetp = (int)(Math.random()*targetresidues.size());
						FloatingResidue tp = targetresidues.get(targetp);
						int previd = tp.backboneIndex;
						double pomega = tp.prevOmega;
						int code = (int)(Math.random()*20);

						if(tp.getName().equals("PRO")){
							if(Math.random() < 0.1){
								tp.prevOmega = 0;
							}
						}else{
							if(Math.random() < 0.03){
								tp.prevOmega = 0;
							}
						}
						if(code == 0){
							tp.prevOmega -= 0.05;
						}
						if(code == 1){
							tp.prevOmega += 0.05;
						}
						tp.processBackBone(cb.backbones.getRandomly(tp));

						for(int ff = targetp+1;ff < targetresidues.size();ff++){
							targetresidues.get(ff).fitBackbone(true);

						}
						double currentdist = Math.abs(last.getC().distance(last.next.getN())- cb.distPenalty.backbone_average[BACKBONE_PREVC_N]);
						if(currentdist < prevdist){
							prevdist = currentdist;
						}else{
							tp.prevOmega = pomega;
							tp.processBackBone(cb.backbones.getOneAt(tp,previd));

							for(int ff = targetp+1;ff < targetresidues.size();ff++){
								targetresidues.get(ff).fitBackbone(true);
							}
						}
						tcount ++;
						if(tcount > iternum*3){
							break;
						}
					}
				}
			}
		}
		
		for(FloatingResidue fr:targetresidues){
			fr.restoreLoc();
		}
		//return sssum/iternum;
		return prevscore;
	}
	
	public static void fixAllUnmapped(ArrayList<FloatingResidue> targets
			,HashSet<FloatingResidue> unmapped
			,int iternum){
		ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(targets,unmapped);
		
		ChainBuilder cb = new ChainBuilder();
		for(FloatingResidue fr:targets){
			if(fr.backboneIndex == -1){
				fr.calcPseudoNeighbour(cb.backbones.getRandomly(fr.getName()));
			}
		}

		FuzzyDecisionTreeScoring_generator j 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
		cb.scoring.clear();
		cb.scoring.add(j);
		cb.prepare(targets);

		for(FloatingResidue fr:targets){
			fr.saveLoc();
			if(fr.sidechainIndex == -1){
				fr.changeSideChain(cb.sidechains.getNext(fr));
			}
		}
		cb.distPenalty.prepare_f(targets);
		ArrayList<PDBResidue> target_r = new ArrayList<>();
		target_r.addAll(targets);
		double[] stepsizez = {5,10};
		for(int ii = 0;ii < iternum;ii++){
			for(double stepsize:stepsizez){
				double firstscore = calcScore(cb,target_r,targets);
				double maxscore = -1000000000.0;
				int maxindex_regionid = -1;
				int maxindex_regionpos = -1;
				ArrayList<VSorter> sal = new ArrayList<>();
				int poscode = 0;//processed segments の id になる
				ArrayList<ArrayList<FloatingResidue>> processed_segments = new ArrayList<>();
				for(int gg = 0;gg < regions.size();gg++){
					ArrayList<FloatingResidue> tregion = regions.get(gg);
							
					for(int rr = 0;rr < tregion.size();rr+=stepsize){
						ArrayList<FloatingResidue> pd = new ArrayList<>();
						ArrayList<PDBResidue> pd_r = new ArrayList<>();
						for(int r = rr;r < rr+stepsize;r++){
							if(r < tregion.size()){
								pd.add(tregion.get(r));
							}
						}
						processed_segments.add(pd);
						pd_r.addAll(pd);
						double csum = 0;
						for(int cc = 0;cc < 10;cc++){
							double cscore = changeBackBone(cb,targets,target_r,pd,pd_r,false,50);
							csum += cscore;
						}
						//double cscore = calcScore(cb,target_r,targets);
						sal.add(new VSorter(csum,poscode));
						poscode++;//processed segments の id になる
						if(csum > maxscore){//使ってない
							maxscore = csum;
							maxindex_regionpos = rr;
							maxindex_regionid = gg;
						}
						for(FloatingResidue fr:pd){
							fr.restoreLoc();
						}		 
					}
				}
				Collections.sort(sal,new VComparator());
				Collections.reverse(sal);
				int tcou = 0;
				for(VSorter vs:sal){
					if(tcou > sal.size() /2.0){
						break;
					}
					tcou++;
					int mm = vs.index;
					ArrayList<FloatingResidue> pd = processed_segments.get(mm);
					ArrayList<PDBResidue> pd_r = new ArrayList<>();
					pd_r.addAll(pd);
					changeBackBone(cb,targets,target_r,pd,pd_r,true,1000);
					double cscore = calcScore(cb,target_r,targets);
					System.out.println(firstscore+"->"+cscore);
				}
			}
		}
	}
	public static ArrayList<ArrayList<FloatingResidue>> getUnmappedRegions(ArrayList<FloatingResidue> allres
			,HashSet<FloatingResidue> unmapped){
		HashSet<FloatingResidue> allhs = new HashSet<>(allres);
		HashSet<FloatingResidue> used = new HashSet<>();
		ArrayList<ArrayList<FloatingResidue>> ret = new ArrayList<>();
		for(int ii = 0;ii < allres.size();ii++){
			FloatingResidue rs = allres.get(ii);
			if(!unmapped.contains(rs) || used.contains(rs) || !allhs.contains(rs)){
				
			}else{
				ArrayList<FloatingResidue> pres = new ArrayList<>();
				pres.add(rs);
				FloatingResidue nex = rs.next;
				while(nex != null){
					if(!unmapped.contains(nex) || used.contains(nex) || !allhs.contains(nex)){
						break;
					}
					pres.add(nex);
					used.add(nex);
					nex = nex.next;
				}
				
				FloatingResidue pre = rs.prev;
				while(pre != null){
					if(!unmapped.contains(pre) || used.contains(pre) || !allhs.contains(pre)){
						break;
					}
					pres.add(0,pre);
					used.add(pre);
					pre = pre.prev;
				}
				ret.add(pres);
			}
		}
		return ret;
	}
	
	
	public static void main(String[] args){
		//HashMap<String,String> hs = parseArgs(args);
		//String infile = hs.get("-in");
		//String outfile = hs.get("-out");
		for(int kk = 1;kk < 6;kk++){
			String infile = "th.T0986s2."+kk+".pdb";
			String outfile = "th.T0986s2."+kk+".ref.pdb";
			PDBData pdb = PDBData.loadPDBFile(infile);
			ArrayList<FloatingResidue> targets = new ArrayList<>();

			ChainBuilder cb = new ChainBuilder();
			for(String c:pdb.chains.keySet()){
				ArrayList<PDBResidue> rr = pdb.chains.get(c).residues;
				targets.addAll(cb.changeToFloating(rr));
			}
			for(FloatingResidue fr:targets){
				fr.calcPseudoNeighbour(cb.backbones.getRandomly(fr.getName()));
			}

			FuzzyDecisionTreeScoring_generator j 
					= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
			cb.scoring.clear();
			cb.scoring.add(j);
			cb.prepare(targets);

			for(FloatingResidue fr:targets){
				fr.saveLoc();
				fr.changeSideChain(cb.sidechains.getNext(fr));
			}
			cb.distPenalty.prepare_f(targets);
			ArrayList<PDBResidue> target_r = new ArrayList<>();
			target_r.addAll(targets);
			double[] stepsizez = {5,10};
			for(int ii = 0;ii < 3;ii++){
				for(double stepsize:stepsizez){
					double firstscore = calcScore(cb,target_r,targets);
					double maxscore = -1000000000.0;
					int maxindex = -1;
					ArrayList<VSorter> sal = new ArrayList<>();
					for(int rr = 0;rr < targets.size();rr+=stepsize){
						ArrayList<FloatingResidue> pd = new ArrayList<>();
						ArrayList<PDBResidue> pd_r = new ArrayList<>();
						for(int r = rr;r < rr+stepsize;r++){
							if(r < targets.size()){
								pd.add(targets.get(r));
							}
						}
						pd_r.addAll(pd);
						double csum = 0;
						for(int cc = 0;cc < 10;cc++){
							double cscore = changeBackBone(cb,targets,target_r,pd,pd_r,false,50);
							csum += cscore;
						}
						//double cscore = calcScore(cb,target_r,targets);
						sal.add(new VSorter(csum,rr));
						if(csum > maxscore){
							maxscore = csum;
							maxindex = rr;
						}
						for(FloatingResidue fr:pd){
							fr.restoreLoc();
						}		 
					}
					Collections.sort(sal,new VComparator());
					Collections.reverse(sal);
					int tcou = 0;
					for(VSorter vs:sal){
						if(tcou > sal.size() /2.0){
							break;
						}
						tcou++;
						int mm = vs.index;
						
						ArrayList<FloatingResidue> pd = new ArrayList<>();
						ArrayList<PDBResidue> pd_r = new ArrayList<>();
						for(int r = mm;r < mm+stepsize;r++){
							if(r < targets.size()){
								pd.add(targets.get(r));
							}
						}
						pd_r.addAll(pd);
						changeBackBone(cb,targets,target_r,pd,pd_r,true,500);
						double cscore = calcScore(cb,target_r,targets);
						System.out.println(firstscore+"->"+cscore);
					}
				}


				double scoresum = 0;
				double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(target_r, cb.scoring);
				double[] penal = new double[target_r.size()];
				for(int r = 0;r < target_r.size();r++){
					penal[r] = cb.distPenalty.calcCrashPenalty(target_r,target_r.get(r),AtomDistancePenalty.TYPE_BACKBONE);
					scores[r] += penal[r];
					scoresum += scores[r];
				}
				System.out.println(scoresum+";;;");
				ChainBuilder.saveChains(targets,outfile+"."+ii+".pdb");

			}
		}
	}
}

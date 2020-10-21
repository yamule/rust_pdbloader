/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Pattern;
import static pepbuilderj.ChainBuilder.changeToFloating;
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;

/**
 *
 * @author kimidori
 */
public class RefineInterFace {
	public static void main__(String[] args){
//		C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p1394
//C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p1
//C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p397
//C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p855
//C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p1038
//C:\dummy\vbox_share\casp13\queries\T0999s\T0999_p1288

		String dirname = "C:\\dummy\\vbox_share\\casp13\\queries\\T1015s1\\singles\\submitter";
		File dir = new File(dirname);
		File[] files = dir.listFiles();
		for(File f:files){
		ArrayList<PDBResidue> allresidue = new ArrayList<>();
			String filename = f.getPath();
			if(!Pattern.compile("\\.pdb$").matcher(filename).find()){
				continue;
			}
			
			if(!Pattern.compile("\\.pdb$").matcher(filename).find()){
				continue;
			}
			PDBData pdb = PDBData.loadPDBFile(filename);
			if(pdb.chains.keySet().size() > 2){
				continue;
			}
			for(String c:pdb.chains.keySet()){
				PDBChain cc = pdb.chains.get(c);
				allresidue.addAll(cc.residues);
			}
			int count = 0;	
			int crashcount = 0;

			for(int ii = 0;ii < allresidue.size();ii++){
				PDBResidue r1 = allresidue.get(ii);
				for(int jj = ii+1;jj < allresidue.size();jj++){
					PDBResidue r2 = allresidue.get(jj);
					if(r1.parent != r2.parent){
						if(r1.getCB().distance(r2.getCB()) < 6.0){
							count++;
						}
						if(r1.getCA().distance(r2.getCA()) < 3.5){//GRY no CA-CA ga 3.7 miman nanode
							crashcount++;
						}
					}
				}
			}
			System.out.println(f.getAbsolutePath()+"\t"+count+"\t"+crashcount);
		}
		
	}
	public static void main_maxrotamer(String[] args){
		String filename = "C:\\dummy\\vbox_share\\casp13\\queries\\T1003\\model_rank1_T1003.5txr.model_0.chain_B.pdb.align.fas.glob.pdb.mul.169.pdb";
		PDBData pdb = PDBData.loadPDBFile(filename);
		ArrayList<ArrayList<PDBAtom>> connected = new ArrayList<>();
		double threshold = 8.0;
		ArrayList<PDBResidue> allresidue = new ArrayList<>();
		
		ChainBuilder cb = new ChainBuilder();
		ArrayList<FloatingResidue> allfloating = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			allfloating.addAll(changeToFloating(cc.residues));
			
		}
		cb.prepare(allfloating);
		for(FloatingResidue r:allfloating){
			BackBoneSample bss = cb.backbones.getRandomly(r.getName());
			r.calcPseudoNeighbour(bss);
		}
		
		for(FloatingResidue r:allfloating){
			cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
		}
		ChainBuilder.saveChains(allfloating, filename+".c2.pdb");
	}
	public static void main(String[] args){
		
		String dirname = "C:\\dummy\\vbox_share\\casp13\\queries\\T1015s1\\singles\\submitter";
		File dir = new File(dirname);
		File[] files = dir.listFiles();
		for(File ff:files){
			if(ff.isDirectory()){
				continue;
			}
						
			if(!Pattern.compile("\\.pdb$").matcher(ff.getAbsolutePath()).find()){
				continue;
			}
			PDBData pdb = PDBData.loadPDBFile(ff.getAbsolutePath());
			ArrayList<ArrayList<PDBAtom>> connected = new ArrayList<>();
			double threshold = 8.0;
			ArrayList<PDBResidue> allresidue = new ArrayList<>();

			ChainBuilder cb = new ChainBuilder();
			ArrayList<FloatingResidue> allfloating = new ArrayList<>();
			for(String c:pdb.chains.keySet()){
				PDBChain cc = pdb.chains.get(c);
				allfloating.addAll(changeToFloating(cc.residues));

			}

			for(FloatingResidue r:allfloating){
				BackBoneSample bss = cb.backbones.getRandomly(r.getName());
				r.calcPseudoNeighbour(bss);
			}
			for(FloatingResidue fr:allfloating){
				fr.saveLoc();
				if(fr.sidechainIndex == -1){
					fr.changeSideChain(cb.sidechains.getNext(fr));
				}
			}

			ArrayList<PDBResidue> allfloating_r = new ArrayList<>();
			allfloating_r.addAll(allfloating);
			ArrayList<ArrayList<FloatingResidue>> regions = new ArrayList<>();

			HashSet<FloatingResidue> unmapped = new HashSet<>();
			for(int ii = 0;ii < allfloating.size();ii++){
				FloatingResidue r1 = allfloating.get(ii);
				for(int jj = ii+1;jj < allfloating.size();jj++){
					FloatingResidue r2 = allfloating.get(jj);
					if(r1.parent != r2.parent){
						ArrayList<PDBAtom> aa = new ArrayList<>();
						aa.add(r1.getC());
						aa.add(r1.getCA());
						aa.add(r1.getO());
						aa.add(r1.getN());
						if(!r1.getName().equals("GLY")){
							aa.add(r1.getCB());
						}
						ArrayList<PDBAtom> aa2 = new ArrayList<>();
						aa2.add(r2.getC());
						aa2.add(r2.getCA());
						aa2.add(r2.getO());
						aa2.add(r2.getN());
						if(!r2.getName().equals("GLY")){
							aa2.add(r2.getCB());
						}
						outer:for(PDBAtom a1:aa){
							for(PDBAtom a2:aa2){
								if(a1.distance(a2) < 1.8){
									unmapped.add(r1);
									unmapped.add(r2);
									if(r1.prev != null){
										unmapped.add(r1.prev);
									}
									if(r1.next != null){
										unmapped.add(r1.next);
									}
									if(r2.prev != null){
										unmapped.add(r2.prev);
									}
									if(r2.next != null){
										unmapped.add(r2.next);
									}
									break outer;
								}
							}
						}
					}
				}	
			}




			regions.addAll(getUnmappedRegions(allfloating,unmapped));
			FuzzyDecisionTreeScoring_generator j 
					= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
			cb.scoring.clear();
			cb.scoring.add(j);
			cb.prepare(allfloating);
			HashSet<FloatingResidue> changed = new HashSet<>();
			for(int ii = 0;ii < 2;ii++){
				for(ArrayList<FloatingResidue> pd:regions){
					ArrayList<PDBResidue> pd_r = new ArrayList<>();	
					pd_r.addAll(pd);
					changed.addAll(pd);
					//System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
					changeBackBone(cb,allfloating,allfloating_r,pd,pd_r,true,5*pd.size());
				}
				ChainBuilder.saveChains(allfloating, ff.getAbsolutePath()+"."+ii+".b.pdb");
			}
			if(true){
				
				for(int ii = 0;ii < 2;ii++){
					for(FloatingResidue r:allfloating){
						cb.refiner.maxRotamer(r, -5,0.3, true, cb.sidechains, false);
						
					}
				}
				ChainBuilder.saveChains(allfloating, ff.getAbsolutePath()+".fin.b.pdb");
				continue;
			}
			HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
			for(int ii = 0;ii < allfloating.size();ii++){
				FloatingResidue r = allfloating.get(ii);
				if(r.next != null && r.getC().distance(r.next.getN()) > 1.4){
					//System.out.println(";;"+modelled.get(ii).getRepresentativeCode());
					unmapped_withneighbour.addAll(cb.getGapFillingResidues(r));
				}

			}
			regions.clear();
			regions.addAll(getUnmappedRegions(allfloating,unmapped_withneighbour));

			for(int ii = 0;ii < 2;ii++){
				for(ArrayList<FloatingResidue> pd:regions){
					ArrayList<PDBResidue> pd_r = new ArrayList<>();	
					pd_r.addAll(pd);
					changed.addAll(pd);
					//System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
					changeBackBone(cb,allfloating,allfloating_r,pd,pd_r,true,5*pd.size());
				}
				ChainBuilder.saveChains(allfloating, ff.getAbsolutePath()+"."+ii+".cc.pdb");
			}

			for(int ii = 0;ii < 2;ii++){
				for(FloatingResidue r:allfloating){
					if(changed.contains(r)){
						cb.refiner.maxRotamer(r, -5,0.3, true, cb.sidechains, false);
					}else{
						r.restoreLoc();
					}
				}
			}
			ChainBuilder.saveChains(allfloating, ff.getAbsolutePath()+".cfin.pdb");
		}
	}
}

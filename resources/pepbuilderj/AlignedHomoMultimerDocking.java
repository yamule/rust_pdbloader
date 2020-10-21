/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static pepbuilderj.ChainBuilder.changeToFloating;
import static pepbuilderj.DockingRefine.center;
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;

/**
 *
 * @author kimidori
 */
public class AlignedHomoMultimerDocking {
	
	
	public static void main_(String[] args){//同じ部分構造からモデリングされたときクラッシュしているので適当に移動させる
		String filename = "C:\\dummy\\vbox_share\\casp13\\queries\\T1004\\try2\\T1004_5m9f_1.pdb";
		PDBData pdb = PDBData.loadPDBFile(filename);
		ArrayList<ArrayList<PDBAtom>> connected = new ArrayList<>();
		double threshold = 8.0;
		ArrayList<PDBResidue> allresidue = new ArrayList<>();
		
		ArrayList<FloatingResidue> allfloating = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			for(PDBResidue r:cc.residues){
				double xshif = 0;
				double yshif = 0;
				double zshif = 0;
				if(r.getResidueNumber() >= 247){
					xshif += 0;
					yshif += 0;
					zshif -= 120;
				}
				
				for(PDBAtom a :r.atoms){
					a.loc.x += xshif;
					a.loc.y += yshif;
					a.loc.z += zshif;
				}
							
			}
			
			allfloating.addAll(changeToFloating(cc.residues));
		}
		
		
		ChainBuilder.saveChains(allfloating, filename+".c.pdb");
	}
	/**
	 * ホモ多量体を作った時に、コンタクトしている奴は一緒に動かしてギャップを埋める。
	 * @param args 
	 */
	public static void main(String[] args){
		String filename = "C:\\dummy\\vbox_share\\casp13\\queries\\T1004\\try2\\th.T1004_nterm.4.pdb.mul.3.pdb.m1.pdb";
		PDBData pdb = PDBData.loadPDBFile(filename);
		ArrayList<ArrayList<PDBAtom>> connected = new ArrayList<>();
		double threshold = 4.0;
		ArrayList<PDBResidue> allresidue = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			for(PDBResidue r:cc.residues){
				allresidue.add(r);
				ArrayList<PDBAtom> closer = null;
				for(ArrayList<PDBAtom> al:connected){
					for(PDBAtom ll:al){
						if(ll.distance(r.getCA()) < threshold){
							closer = al;
						}
					}
				}
				
				if(closer == null){
					ArrayList<PDBAtom> tmp = new ArrayList<PDBAtom>();
					tmp.add(r.getCA());
					connected.add(tmp);
				}else{
					closer.add(r.getCA());
				}
			}
		}
		while(true){
			boolean updated = false;
			for(int ii = 0;ii < connected.size();ii++){
				int mergeid = -1;
				outer:for(int jj = ii+1;jj < connected.size();jj++){
					for(PDBAtom a1:connected.get(ii)){
						for(PDBAtom a2:connected.get(jj)){
							if(a1.distance(a2) < threshold){
								mergeid = jj;
								break outer;
							}
						}	
					}
				}
				if(mergeid > -1){
					connected.get(ii).addAll(connected.get(mergeid));
					connected.get(mergeid).clear();
					updated = true;
					break;
				}
			}
			if(!updated){
				break;
			}
		}
		Iterator<ArrayList<PDBAtom>> ite = connected.iterator();
		while(ite.hasNext()){
			ArrayList<PDBAtom> au = ite.next();
			if(au.size() == 0){
				ite.remove();
			}
		}
		//System.out.println(";;");
		ArrayList<ArrayList<Point3D>> blocks = new ArrayList<>();
		for(int ii = 0;ii < connected.size();ii++){
			ArrayList<Point3D> block = new ArrayList<>();
			for(PDBAtom aa:connected.get(ii)){
				for(PDBAtom a:aa.parent.atoms){
					block.add(a.loc);
				}
			}
			HashSet<Point3D> chk = new HashSet<>(block);
			if(chk.size() != block.size()){
				System.err.println("???? error in code. but anyway it can be processed.");
				block.clear();
				block.addAll(chk);
			}
			blocks.add(block);
		}
		
		ArrayList<ArrayList<PDBAtom>> peptidebonds = new ArrayList<>();
		
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			for(int ii = 0;ii < cc.residues.size();ii++){
				if(ii > 0){
					ArrayList<PDBAtom> p = new ArrayList<>();
					p.add(cc.residues.get(ii).getN());
					p.add(cc.residues.get(ii-1).getC());
					peptidebonds.add(p);
				}
				if(ii < cc.residues.size() -1){
					ArrayList<PDBAtom> p = new ArrayList<>();
					p.add(cc.residues.get(ii+1).getN());
					p.add(cc.residues.get(ii).getC());
					peptidebonds.add(p);
				}
			}
		}
		double distsum = 0;
		HashMap<Point3D,Point3D> prevpoint = new HashMap<>();
		HashMap<Point3D,HashSet<Point3D>> crashedpair = new HashMap<>();
		
		
		for(ArrayList<Point3D> b:blocks){
			for(Point3D a:b){
				crashedpair.put(a,new HashSet<Point3D>());
				prevpoint.put(a,new Point3D(a));
			}
		}
		ArrayList<Point3D> allpoint = new ArrayList<>(crashedpair.keySet());
		for(Point3D p:allpoint){
			for(Point3D p2:allpoint){
				if(p.distance(p2) < 2.0){
					crashedpair.get(p).add(p2);
					crashedpair.get(p2).add(p);
				}
			}
		}
		
		for(ArrayList<PDBAtom> al:peptidebonds){
			distsum += Math.abs(al.get(0).distance(al.get(1)) - 1.30);
		}
		for(int ii = 0;ii < 300;ii++){
			
			boolean cflag = false;
			for(int bb = 0;bb < blocks.size();bb++){
				int randcode = (int)(Math.random()*blocks.size());
				ArrayList<Point3D> a = blocks.get(randcode);
				double fact = 1.0;
				if(ii%2 == 0){
					fact = 5.0;
				}
				Point3D direc = new Point3D(
						(Math.random()-0.5)*fact
						,(Math.random()-0.5)*fact
						,(Math.random()-0.5)*fact
				);
				RandomDocker.moves(a, direc);

				double rotx = (Math.random()*Math.PI-Math.PI/2)*(fact/5.0);
				double roty =  (Math.random()*Math.PI-Math.PI/2)*(fact/5.0);
				double rotz =  (Math.random()*Math.PI-Math.PI/2)*(fact/5.0);
				Point3D c = center(a);
				RandomDocker.rotates(a,c,rotx,roty,rotz);

				outer:for(int bi = 0;bi < blocks.size();bi++){
					for(int bi2 = bi+1;bi2 < blocks.size();bi2++){

						for(Point3D p:blocks.get(bi)){
							for(Point3D p2:blocks.get(bi2)){

								if(!crashedpair.get(p).contains(p2)){
									if(p.distance(p2) < 2.0){
										cflag = true;
										break outer;
									}
								}
							}
						}
					}
				}
			}
			if(cflag){
				for(Point3D p:prevpoint.keySet()){
					p.set(prevpoint.get(p));
				}
				continue;
			}
			double cdistsum = 0;
			for(ArrayList<PDBAtom> al:peptidebonds){
				cdistsum +=  Math.abs(al.get(0).distance(al.get(1)) - 1.30);
			}
			if(cdistsum > distsum){
				for(Point3D p:prevpoint.keySet()){
					p.set(prevpoint.get(p));
				}
			}else{
				distsum = cdistsum;
				
				for(Point3D p:prevpoint.keySet()){
					prevpoint.get(p).set(p);
				}
			}
			if(ii%5 == 0){
				//ChainBuilder.savePDBChains(allresidue, "testout."+ii+".pdb");
			}
		}
		
		
		for(Point3D p:prevpoint.keySet()){
			p.set(prevpoint.get(p));
		}
		ChainBuilder.savePDBChains(allresidue, filename+".m.pdb");
		
		pdb = PDBData.loadPDBFile( filename+".m.pdb");
		ArrayList<FloatingResidue> allfloating = new ArrayList<>();
		ArrayList<ArrayList<FloatingResidue>> regions = new ArrayList<>();
		ChainBuilder cb = new ChainBuilder();
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			ArrayList<FloatingResidue> res = changeToFloating(cc.residues);
			allfloating.addAll(res);
			for(FloatingResidue r:res){
				BackBoneSample bss = cb.backbones.getRandomly(r.getName());
				r.calcPseudoNeighbour(bss);
			}
			for(FloatingResidue fr:res){
				fr.saveLoc();
				if(fr.sidechainIndex == -1){
					fr.changeSideChain(cb.sidechains.getNext(fr));
				}
			}
			
			//fixme 露出残基をとって、露出している方から削っていく
			//インサーションが入った所為で別残基とペプチド結合を作ってしまっているようなのはのける
			//ChainBuilder.saveChains(modelled,outfilename+".p.pdb");
			HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
			for(int ii = 0;ii < res.size()-1;ii++){
				if(res.get(ii).getC().distance(res.get(ii+1).getN()) > 1.4){
					//System.out.println(";;"+modelled.get(ii).getRepresentativeCode());
					unmapped_withneighbour.addAll(cb.getGapFillingResidues(res.get(ii)));
				}

			}
			regions.addAll(getUnmappedRegions(res,unmapped_withneighbour));

		}
		ArrayList<PDBResidue> allfloating_r = new ArrayList<>();
		allfloating_r.addAll(allfloating);
		
		
		FuzzyDecisionTreeScoring_generator j 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
		cb.scoring.clear();
		cb.scoring.add(j);
		cb.prepare(allfloating);
		for(int ii = 0;ii < 5;ii++){
			for(ArrayList<FloatingResidue> pd:regions){
				ArrayList<PDBResidue> pd_r = new ArrayList<>();	
				pd_r.addAll(pd);
				//System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
				changeBackBone(cb,allfloating,allfloating_r,pd,pd_r,true,10*pd.size());
			}
			ChainBuilder.saveChains(allfloating, filename+"."+ii+".c.pdb");
		}
		
		ChainBuilder.saveChains(allfloating, filename+".c.pdb");
	}
}

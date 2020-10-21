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
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;
import static pepbuilderj.TemplateBaseModeller.makeAtomLines;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class ChainBuilder {
	
	ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
	BackBoneSet backbones = new BackBoneSet();
	SideChainSet sidechains = new SideChainSet();
	AtomDistancePenalty distPenalty = new AtomDistancePenalty();
	int sampleCount_Threshold = 3;//sample 数がこれ以下しかない場合削除する
	
	ArrayList<FloatingResidue> res = new ArrayList<>();
	ArrayList<PDBResidue> _res = new ArrayList<>();
	ResidueRefine refiner = new ResidueRefine();
	ScoringDocker docker = new ScoringDocker();
	
	
	ChainBuilder(){
		for(String bb:backbones.backbones.keySet()){
			ArrayList<BackBoneSample> al = backbones.backbones.get(bb);
			Iterator<BackBoneSample> ite = al.iterator();
			while(ite.hasNext()){
				BackBoneSample b = ite.next();
				if(b.count < sampleCount_Threshold){
					ite.remove();
				}
			}
			
			if(al.size() == 0){
				System.err.println("????Did you load right sample file??");
				throw new RuntimeException();
			}
			
		}
		for(String bb:sidechains.sidechains.keySet()){
			ArrayList<SideChainSample> al = sidechains.sidechains.get(bb);
			Iterator<SideChainSample> ite = al.iterator();
			while(ite.hasNext()){
				SideChainSample b = ite.next();
				if(b.count < sampleCount_Threshold){
					ite.remove();
				}
			}
			if(al.size() == 0){
				System.err.println("????Did you load right sample file??");
				throw new RuntimeException();
			}	
		}
		
		
		FuzzyDecisionTreeScoring_generator j 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorAtom_20180503());
		scoring.add(j);
		
	}
	public ArrayList<AtomTriangle> makeTriangles(ArrayList<FloatingResidue> p){
		ArrayList<AtomTriangle> ret = new ArrayList<>();
		for(int ii = 0;ii < p.size();ii++){
			FloatingResidue r = p.get(ii);
			FloatingResidue prev = r.prev;
			FloatingResidue next = r.next;
			
	//public static final int BACKBONE_N_CA = 0;
	//public static final int BACKBONE_CA_C = 1;
	//public static final int BACKBONE_C_O = 2;
	//public static final int BACKBONE_N_C = 3;
	//public static final int BACKBONE_PREVC_CA = 4;
	//public static final int BACKBONE_PREVCA_N = 5;
	//public static final int BACKBONE_PREVC_N = 6;
			if(prev != null){
				AtomTriangle t1 = new AtomTriangle(
				prev.getC(),r.getN(),r.getCA(),
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_N_CA],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_CA]
				);
				ret.add(t1);
				
				if(ii == 0){
					t1.freeze(0);
				}
			}
			AtomTriangle t2 = new AtomTriangle(
				r.getN(),r.getCA(),r.getC(),
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_N_CA],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_CA_C],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_N_C]
				);
			if(next != null){
				AtomTriangle t3 = new AtomTriangle(
				r.getCA(),r.getC(),next.getN(),
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_CA_C],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N],
				distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVCA_N]
				);
				ret.add(t3);
				if(ii == p.size() -1){
					t3.freeze(2);
				}
			}
		}
		return ret;
	}
	
	public void prepare(ArrayList<FloatingResidue> al){
		res.clear();
		_res.clear();
		res.addAll(al);
		_res.addAll(al);
		refiner.prepare(_res);
		distPenalty.prepare_f(res);
		for(FuzzyDecisionTreeScoring_generator f:scoring){
			f.prepare(_res);
		}
	}
	
	
	
	public void waveAllTriangle(ArrayList<AtomTriangle> triangles
	,int itermax,double threshold){
		HashMap<PDBResidue,ArrayList<AtomTriangle>> residues = new HashMap<>();
		for(AtomTriangle a:triangles){
			PDBResidue af = a.atoms[1].parent;
			if(!residues.containsKey(af)){
				residues.put(af,new ArrayList<AtomTriangle>());
			}
			residues.get(af).add(a);
		}
		int iter = 0;
		double dsumprev = 10000;
		while(true){
			for(int ii = 0;ii < 1000;ii++){
				for(AtomTriangle t:triangles){
					t.fitThree(0.8);
				}
			}
			double dsum = 0;
			
			for(AtomTriangle t:triangles){
				dsum += t.diffSum();
			}
			double psum = 0;
			
			for(PDBResidue r:residues.keySet()){
				fitO((FloatingResidue)r);
			}
			for(PDBResidue r:residues.keySet()){
				double penal = distPenalty.calcCrashPenalty(_res, r,AtomDistancePenalty.TYPE_ALL);
				psum += penal;
				ArrayList<AtomTriangle> al = residues.get(r);
				for(AtomTriangle a:al){
					if(penal > a.savedCrash){
						a.savedCrash = penal;
						a.saveLoc();
					}else{
						if(penal == 0){
							if(dsumprev < dsum){
								a.saveLoc();
							}
						}else{
							a.restoreLoc();
						}
					}
				}
			}
			dsumprev = dsum;
			iter++;
			if(iter > itermax){
				break;
			}
			//System.out.println(dsum+";;;"+psum);
			if(dsum < threshold && psum == 0){
				break;
			}
			for(AtomTriangle a:triangles){
				for(int ii = 0;ii < 3;ii++){
					a.atoms[ii].loc.x += Math.random()*2-1.0;
					a.atoms[ii].loc.y += Math.random()*2-1.0;
					a.atoms[ii].loc.z += Math.random()*2-1.0;
				}
			}
		}
		
		for(int ii = 0;ii < 1000;ii++){
			for(AtomTriangle t:triangles){
				t.fitThree(0.8+0.1*ii/(double)1000);
			}
		}
		
		for(PDBResidue r:residues.keySet()){
			fitO((FloatingResidue)r);
		}
	}
	
	/**
	 * Triangle は アラニンかグリシンを想定しており、Triangle を作成した後元の残基を戻す。
	 * 
	 * @param triangleresidue
	 * @param realresidue 
	 */
	public void fitToTriangle(ArrayList<FloatingResidue> triangleresidue
	,ArrayList<FloatingResidue> realresidue
	){
		for(int ii = 0;ii < triangleresidue.size();ii++){
			FloatingResidue fr = triangleresidue.get(ii);
			if(fr.prev != null){
				fr.prevc.loc.set(fr.prev.getC().loc);
			}
			fitO(fr);
			if(fr.next != null){
				fr.nextn.loc.set(fr.next.getN().loc);
			}
			FloatingResidue rr = realresidue.get(ii);
			rr.nextn.loc.set(fr.nextn.loc);
			rr.prevc.loc.set(fr.prevc.loc);
			rr.getC().loc.set(fr.getC().loc);
			rr.getN().loc.set(fr.getN().loc);
			rr.getCA().loc.set(fr.getCA().loc);
			rr.getO().loc.set(fr.getO().loc);
			rr.fitSideChain();
		}
		if(triangleresidue.get(0).prev != null){
			triangleresidue.get(0).prev.setNextFloating(realresidue.get(0));
		}
		
		if(triangleresidue.get(triangleresidue.size()-1).next != null){
			triangleresidue.get(triangleresidue.size()-1).next.setPrevFloating(
			realresidue.get(realresidue.size()-1));
		}
		
	}
	
	
	public ArrayList<FloatingResidue> generateDummyG(
	ArrayList<FloatingResidue> realresidue
	){
		ArrayList<String> al = new ArrayList<>();
		for(int ii = 0;ii < realresidue.size();ii++){
			al.add("GLY");
		}
		ArrayList<FloatingResidue> triangleresidue = generateChain(al);
		
		for(int ii = 0;ii < realresidue.size();ii++){
			FloatingResidue fr = realresidue.get(ii);
			if(fr.prev != null){
				fr.prevc.loc.set(fr.prev.getC().loc);
			}
			fitO(fr);
			if(fr.next != null){
				fr.nextn.loc.set(fr.next.getN().loc);
			}
			FloatingResidue rr = triangleresidue.get(ii);
			rr.nextn.loc.set(fr.nextn.loc);
			rr.prevc.loc.set(fr.prevc.loc);
			rr.getC().loc.set(fr.getC().loc);
			rr.getN().loc.set(fr.getN().loc);
			rr.getCA().loc.set(fr.getCA().loc);
			rr.getO().loc.set(fr.getO().loc);
			rr.fitSideChain();
		}
		if(realresidue.get(0).prev != null){
			realresidue.get(0).prev.setNextFloating(triangleresidue.get(0));
			triangleresidue.get(0).setPrevFloating(realresidue.get(0).prev);
		}
		
		if(realresidue.get(realresidue.size()-1).next != null){
			realresidue.get(realresidue.size()-1).next.setPrevFloating(
			triangleresidue.get(triangleresidue.size()-1));
			triangleresidue.get(triangleresidue.size()-1).setNextFloating(
					realresidue.get(realresidue.size()-1).next);
		}
		return triangleresidue;
	}
	
	public static ArrayList<FloatingResidue> changeToFloating(ArrayList<PDBResidue> rr){
		ArrayList<FloatingResidue> fr = new ArrayList<>();
		PDBChain c = rr.get(0).parent;
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
			ff.setParent(c);
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
		}
		
		
		return fr;
	}
	
	
	
	public ArrayList<FloatingResidue> sequencialDock(ArrayList<ArrayList<FloatingResidue>> groups ){
		ChainBuilder cb = new ChainBuilder();
		PDBChain pc = groups.get(0).get(0).parent;
		PDBChain dummy = new PDBChain("B");
		ArrayList<PDBResidue> template = new ArrayList<>();
		ArrayList<FloatingResidue> current = new ArrayList<>();
		current = groups.remove(0);
		template.addAll(current);
		cb.refiner.prepare(template);
		cb.distPenalty.prepare(template);
		for(FuzzyDecisionTreeScoring_generator f:cb.scoring){
			f.prepare(template);
		}
		while(groups.size() > 0){
			current = groups.remove(0);
			ArrayList<PDBResidue> _current = new ArrayList<>();
			_current.addAll(current);
			ArrayList<ArrayList<PDBResidue>> pl = cb.docker.dock(template
			, _current, 20,3);
			PDBAtom ac = template.get(template.size()-1).getC();
			double cndist = 100000;
			int minindex = 0;
			double okthreshold = 5;
			for(int jj = 0;jj < pl.size();jj++){
				ArrayList<PDBResidue> p = pl.get(jj);
				double ddd = p.get(0).getN().distance(ac);
				if(ddd < okthreshold){
					minindex = jj;
					break;
				}
				if(ddd < cndist){
					minindex = jj;
					cndist = ddd;
				}
				
			}
			template.addAll(pl.get(minindex));
			cb.refiner.prepare(template);
			cb.distPenalty.prepare(template);
			for(FuzzyDecisionTreeScoring_generator f:cb.scoring){
				f.prepare(template);
			}
			
			for(PDBResidue r:template){
				r.setParent(pc);
			}
		}
		
		
		ArrayList<FloatingResidue> ret = changeToFloating(template);
		cb.prepare(ret);
		for(FloatingResidue r:ret){
			cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
		}
		//cb.fillGapsAll(ret,1000);
		ArrayList<FloatingResidue> allresidue = ret;
		ArrayList<ArrayList<FloatingResidue>> regions = new ArrayList<>();
		for(FloatingResidue r:allresidue){
			BackBoneSample bss = cb.backbones.getRandomly(r.getName());
			r.calcPseudoNeighbour(bss);
		}
		for(FloatingResidue fr:allresidue){
			fr.saveLoc();
			if(fr.sidechainIndex == -1){
				fr.changeSideChain(cb.sidechains.getNext(fr));
			}
		}

		HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
		for(int ii = 0;ii < allresidue.size();ii++){
			if(allresidue.get(ii).next != null){
				if(allresidue.get(ii).getC().distance(allresidue.get(ii).next.getN())> 1.4){
				//System.out.println(";;"+modelled.get(ii).getRepresentativeCode());
					unmapped_withneighbour.addAll(cb.getGapFillingResidues(allresidue.get(ii)));
				}
			}

		}
		regions.addAll(getUnmappedRegions(allresidue,unmapped_withneighbour));

		ArrayList<PDBResidue> allfloating_r = new ArrayList<>();
		allfloating_r.addAll(allresidue);
		FuzzyDecisionTreeScoring_generator j 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
		cb.scoring.clear();
		cb.scoring.add(j);
		cb.prepare(allresidue);
		for(int ii = 0;ii < 5;ii++){
			for(ArrayList<FloatingResidue> pd:regions){
				ArrayList<PDBResidue> pd_r = new ArrayList<>();	
				pd_r.addAll(pd);
				//System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
				changeBackBone(cb,allresidue,allfloating_r,pd,pd_r,true,10*pd.size());
			}
		}
			
		
		return ret;
	}
	
	
	
	
	
	public ArrayList<FloatingResidue> halfDock(ArrayList<FloatingResidue> head
			,ArrayList<FloatingResidue> target,ArrayList<FloatingResidue> tail){
		PDBChain pc = target.get(0).parent;
		ArrayList<PDBResidue> template = new ArrayList<>();
		ArrayList<PDBResidue> tar = new ArrayList<>();
		template.addAll(head);
		template.addAll(tail);
		tar.addAll(target);
		ArrayList<ArrayList<PDBResidue>> pl = docker.dock(template
		, tar, 20,3);
		
		double maxscore = -100000;
		ArrayList<PDBResidue> maxlist = null;
		
		ChainBuilder cb = new ChainBuilder();
		for(ArrayList<PDBResidue> rl:pl){
			for(FloatingResidue f:head){
				f.saveLoc();
			}
			for(FloatingResidue f:tail){
				f.saveLoc();
			}
			
			ArrayList<FloatingResidue> fl = new ArrayList<>();
			for(FloatingResidue f:head){
				fl.add(f);
			}
			for(int jj = 0;jj < rl.size();jj++){
				fl.add(new FloatingResidue(rl.get(jj),
						(jj > 0)?(rl.get(jj-1)):(null),(jj < rl.size()-1)?(rl.get(jj+1)):(null)));
			}
			for(FloatingResidue f:tail){
				fl.add(f);
			}
			for(int ii = 1;ii < fl.size();ii++){
				fl.get(ii).setPrevFloating(fl.get(ii-1));
			}
			for(int ii = 0;ii < fl.size()-1;ii++){
				fl.get(ii).setNextFloating(fl.get(ii+1));
			}
				
			for(FloatingResidue f:fl){
				f.changeSideChain(sidechains.sidechains.get(f.getName()).get(0));
			}
			cb.prepare(fl);
			cb.elongTerm(fl.get(0)
			,false,100,false);
			cb.elongTerm(fl.get(fl.size()-1)
			,true,100,false);
			
		
			
			
			prepare(fl);
			ArrayList<PDBResidue> fll = new ArrayList<>();
			fll.addAll(fl);
			for(FloatingResidue r:fl){
				refiner.maxRotamer(r, -5,0.1, true, sidechains, false);
			}
			
			
			fillGapsAll(fl,1000);
			
			
			//ChainBuilder.saveChains(res,"test10.pdb");
			
			FuzzyDecisionTreeScoring_generator.prepareAll(fll,scoring);
			double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(fll, scoring);
			double pp = 0;
			for(PDBResidue f:fll){
				pp += calcAtomCrash_Higher_is_better(fll,f);
			}
			double sc = pp;
			for(double dd:scores){
				sc+=dd;
			}
			if(sc > maxscore){
				maxlist = new ArrayList<>();
				for(FloatingResidue f:fl){
					maxlist.add(f.getCopy());
				}
				maxscore = sc;
			}
			for(FloatingResidue f:fl){
				f.restoreLoc();
			}
		}
		ArrayList<FloatingResidue> fl = new ArrayList<>();
		for(PDBResidue p:maxlist){
			FloatingResidue fff = new FloatingResidue(p,null,null);
			fl.add(fff);
			fff.setParent(pc);
		}
		for(int ii = 1;ii < fl.size();ii++){
			fl.get(ii).setPrevFloating(fl.get(ii-1));
		}
		for(int ii = 0;ii < fl.size()-1;ii++){
			fl.get(ii).setNextFloating(fl.get(ii+1));
		}
		return fl;
	}
	public void fillGapsTarget(ArrayList<FloatingResidue> target,int iternum){

		int ite = 0;
		while(ite < iternum){
			boolean flag = false;
			for(int ii = 0;ii < target.size();ii++){
				//double ddis = fl.get(ii).nextn.distance(fl.get(ii+1).getN());
				if(target.get(ii).next == null){
					elongTerm(target.get(ii)
					,true,100,false);
					continue;
				}
				if(target.get(ii).prev == null){
					elongTerm(target.get(ii)
					,false,100,false);
					continue;
				}
				double ddis = target.get(ii).getC().distance(target.get(ii).next.getN());
				if(Math.abs(ddis - distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N])
						> 0.3){
					
					fillGapBetween(target.get(ii),target.get(ii).next,
					Math.PI,100);
					flag = true;
					ddis = target.get(ii).nextn.distance(target.get(ii).next.getN());
				}
			}
			if(!flag){
				
				break;
			}
			ite++;
		}
		
		
		//ChainBuilder.saveChains(res, "tess1.pdb");
		//ChainBuilder.saveChains(res, "tess2.pdb");
		//fl.get(0).processBackBone(backbones.backbones.get(fl.get(0).getName()).get(0),Math.PI,true,false);
		//fl.get(fl.size()-1).processBackBone(backbones.backbones.get(fl.get(fl.size()-1).getName()).get(0),Math.PI,true,false);
		
		//ChainBuilder.saveChains(res, "tess3.pdb");
	}
	
	public void fillGapsAll(ArrayList<FloatingResidue> fl,int iternum){
		for(int ii = 1;ii < fl.size();ii++){
			fl.get(ii).setPrevFloating(fl.get(ii-1));
		}
		for(int ii = 0;ii < fl.size()-1;ii++){
			fl.get(ii).setNextFloating(fl.get(ii+1));
		}
		
		this.prepare(fl);
		int ite = 0;
		while(ite < iternum){
			boolean flag = false;
			for(int ii = 0;ii < fl.size()-1;ii++){
				//double ddis = fl.get(ii).nextn.distance(fl.get(ii+1).getN());
				
				double ddis = fl.get(ii).getC().distance(fl.get(ii+1).getN());
				if(Math.abs(ddis - distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N])
						> 0.3){
					
					fillGapBetween(fl.get(ii),fl.get(ii+1),
					Math.PI,100);
					flag = true;
					ddis = fl.get(ii).nextn.distance(fl.get(ii+1).getN());
				}
			}
			if(!flag){
				
				break;
			}
			ite++;
		}
		
		for(int ii = 0;ii < fl.size();ii++){
			FloatingResidue fr = fl.get(ii);
			if(fr.prev != null){
				fr.prevc.loc.set(fr.prev.getC().loc);
			}
			if(fr.next != null){
				fr.nextn.loc.set(fr.next.getN().loc);
			}
		}
		if(fl.size() > 1){
			fl.get(1).processBackBone(this.backbones.getOneAt(
					fl.get(1),Math.max(0,fl.get(1).backboneIndex))
					,fl.get(1).nextOmega,false,true);
			fl.get(fl.size()-2).processBackBone(this.backbones.getOneAt(
					fl.get(fl.size()-2),Math.max(0,fl.get(fl.size()-2).backboneIndex)));
		}
		elongTerm(fl.get(0)
		,false,100,false);
		elongTerm(fl.get(fl.size()-1)
		,true,100,false);
		//ChainBuilder.saveChains(res, "tess2.pdb");
		//fl.get(0).processBackBone(backbones.backbones.get(fl.get(0).getName()).get(0),Math.PI,true,false);
		//fl.get(fl.size()-1).processBackBone(backbones.backbones.get(fl.get(fl.size()-1).getName()).get(0),Math.PI,true,false);
		
		//ChainBuilder.saveChains(res, "tess3.pdb");
	}
	
	/**
	 * CA C nextn の位置を基に O の位置を決定する
	 * @param r 
	 */
	public  void fitO(FloatingResidue r){
		if(r.next != null){
			//O の位置を決定
			Point3D point1 = r.getCA().loc;
			Point3D point3 = r.next.getN().loc;
			Point3D veccca = PepProcess.subtracted(point1,point3);
			Point3D c = new Point3D();//面倒くさいので中心点にしている。
			c.x = point3.x+veccca.x/2;
			c.y = point3.y+veccca.y/2;
			c.z = point3.z+veccca.z/2;
			double ddist = r.getC().loc.distance(c);
			double ccdist = distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_C_O];
			if(ddist > 0){
				Point3D lvec = PepProcess.subtracted(r.getC().loc, c);
				lvec.x *= ccdist/ddist;
				lvec.y *= ccdist/ddist;
				lvec.z *= ccdist/ddist;
				
				lvec.x += r.getC().loc.x;
				lvec.y += r.getC().loc.y;
				lvec.z += r.getC().loc.z;
				r.getO().loc.set(lvec);
			}
		}else{
			//最後の残基の場合CAとれないので
			//ついでに一個しかつけない
			//Fixit
			Point3D point1 = r.getC().loc;
			Point3D point2 = r.getCA().loc;
			Point3D veccca = PepProcess.subtracted(point1,point2);
			double ddist = veccca.x*veccca.x+ veccca.y*veccca.y+ veccca.z*veccca.z;
			double ccdist = distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_C_O];
			if(ddist > 0){
				ddist = Math.sqrt(ddist);
				Point3D lvec = PepProcess.subtracted(r.getC().loc, r.getCA().loc);
				
				lvec.x *= ccdist/ddist;
				lvec.y *= ccdist/ddist;
				lvec.z *= ccdist/ddist;
				
				lvec.x += r.getC().loc.x;
				lvec.y += r.getC().loc.y;
				lvec.z += r.getC().loc.z;
				r.getO().loc.set(lvec);
			}
		}
	}
	
	
	
	//crash は負の値が返る
	public double calcBackBoneCrash(ArrayList<PDBResidue> al,int index){
		return distPenalty.calcCrashPenalty(al,  al.get(index), distPenalty.TYPE_BACKBONE);
	}
	
	public double calcSideChainCrash(ArrayList<PDBResidue> al,int index){
		return distPenalty.calcCrashPenalty(al, al.get(index), distPenalty.TYPE_SIDECHAIN);
	}
	
	public double calcAtomCrash_Higher_is_better(ArrayList<PDBResidue> al,PDBResidue target){
		return distPenalty.calcCrashPenalty(al, target, distPenalty.TYPE_ALL);
	}
	
	public double calcBackBoneCrash_Higher_is_better(ArrayList<PDBResidue> al,PDBResidue target){
		return distPenalty.calcCrashPenalty(al, target, distPenalty.TYPE_BACKBONE);
	}
	
	/**
	 * 距離は一つだけの簡易チェック
	 * CYSCYS 以外はこれでいいかも
	 * @param r
	 * @param atoms
	 * @param threshold
	 * @return 
	 */
	public boolean checkSideChainCrash(FloatingResidue r,ArrayList<PDBAtom> atoms,double threshold){
		for(PDBAtom a:r.atoms_sidechain.values()){
			for(PDBAtom aa:atoms){
				if(aa.parent != r){
					if(aa.distance(a) < threshold){
						return true;
					}
				}
			}
		}
		return false;
	}
	
	/**
	 * J48 Scoring を使ってリストに含まれる Residue のスコアを計算して返す。
	 * @param r
	 * @return 
	 */
	public double[] calcChainScore(ArrayList<PDBResidue> r){
		FuzzyDecisionTreeScoring_generator.prepareAll(r,scoring);
		double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(r, scoring);
		double[] backboneprob = new double[r.size()];
		double probsum = 0;
		int probcount = 0;
		for(int ii = 1;ii < r.size()-1;ii++){
			double pp = backbones.getProb(r,ii);
			backboneprob[ii] = pp;
			probsum += pp;
			probcount++;
		}
		if(probcount == 0){
			for(int ii = 0;ii < r.size();ii++){
				backboneprob[ii] = 1;
			}
		}else{
			backboneprob[0] = probsum/probcount;
			backboneprob[r.size()-1] = probsum/probcount;
		}
		double ret = 0;
		for(int ii = 0;ii < r.size();ii++){
			if(scores[ii] > 0){//Backbone の確率と元の確率をかけたもの
			scores[ii] *= backboneprob[ii];
			}else{
				scores[ii] *= 1.0-backboneprob[ii];
			}
		}
		return scores;
	}
	
	public ArrayList<FloatingResidue> buildDummy(ArrayList<String> aa){
		ArrayList<FloatingResidue> ret = new ArrayList<>();
		int cou = 0;
		for(String s:aa){
			if(s.length() == 1){
				FloatingResidue r = new FloatingResidue(MyWorld.toThreeLetter.get(s));
				r.setResidueNumber(cou+1);
				cou++;
				ret.add(r);
			}else{
				FloatingResidue r = new FloatingResidue(s);
				r.setResidueNumber(cou+1);
				cou++;
				ret.add(r);

			}
		}
		for(int ii = 0;ii < ret.size();ii++){
			if(ii > 0){
				ret.get(ii).prev = ret.get(ii-1);
			}
			if(ii < ret.size()-1){
				ret.get(ii).next = ret.get(ii+1);
			}
			
		}
		FloatingResidue fr = ret.get(0);
		ArrayList<BackBoneSample> ssal = backbones.backbones.get(ret.get(0).getName());
		BackBoneSample s = ssal.get(0);
		fr.getC().loc.set(s.atoms[BackBoneSample.C].loc);
		fr.getCA().loc.set(s.atoms[BackBoneSample.CA].loc);
		fr.getN().loc.set(s.atoms[BackBoneSample.N].loc);
		fr.getO().loc.set(s.atoms[BackBoneSample.O].loc);
		for(int ii = 1;ii < ret.size();ii++){
			FloatingResidue r = ret.get(ii);
			r.processBackBone(backbones.backbones.get(r.getName()).get(0),Math.PI,true,false);
			r.changeSideChain(sidechains.sidechains.get(r.getName()).get(0));
		}
		
		fr.processBackBone(backbones.backbones.get(fr.getName()).get(0),Math.PI,false,false);
		fr.changeSideChain(sidechains.sidechains.get(fr.getName()).get(0));
		return ret;
	}
	
	
	/**
	 * 平行移動して最も近い位置でぶつからない場所に移動する
	 */
	public static void moveToClosest(ArrayList<FloatingResidue> al,FloatingResidue f){
		Point3D l = f.getCA().loc;
		ArrayList<Point3D> direc = new ArrayList<>();
		double span = 0.5;
		direc.add(new Point3D(l.x-span,l.y,l.z));
		direc.add(new Point3D(l.x,l.y-span,l.z));
		direc.add(new Point3D(l.x,l.y,l.z-span));
		direc.add(new Point3D(l.x+span,l.y,l.z));
		direc.add(new Point3D(l.x,l.y+span,l.z));
		direc.add(new Point3D(l.x,l.y,l.z+span));
		
		HashMap<Point3D,PDBAtom> maxpoint = new HashMap();
		HashMap<Point3D,Double> maxdist = new HashMap();
		
		for(FloatingResidue ff:al){
			if(ff == f){
				continue;
			}
			for(PDBAtom a:ff.atoms){
				Point3D np = direc.get(0);
				double mindist = direc.get(0).distance(a.loc);
				for(int ii = 1;ii < direc.size();ii++){
					double mdist = direc.get(ii).distance(a.loc);
					if(mdist < mindist){
						mindist = mdist;
						np = direc.get(ii);
					}
				}
				if(!maxdist.containsKey(np)){
					maxdist.put(np,mindist);
					maxpoint.put(np,a);
				}else if(mindist > maxdist.get(np)){
					maxdist.put(np,mindist);
					maxpoint.put(np,a);
				}
			}
		}
		Point3D mindirec = null;
		Point3D basepoint = null;
		double dmindist = Double.MAX_VALUE;
		boolean noatomflag = false;
		for(Point3D p:direc){
			if(!maxdist.containsKey(p)){
				mindirec  = p;
				noatomflag = true;
			}else{
				if(dmindist > maxdist.get(p)){
					mindirec = p;
					dmindist = maxdist.get(p);
					basepoint = maxpoint.get(p).loc;
				}
			}
		}
		mindirec.x -= l.x;
		mindirec.y -= l.y;
		mindirec.z -= l.z;

		
		HashSet<FloatingResidue> ignore = new HashSet<>();
		for(int ii = 0;ii < 20;ii++){
			boolean crashes = false;
			double mindist = Double.MAX_VALUE;
			for(FloatingResidue ff:al){
				if(f == ff){
					continue;
				}
				if(ignore.contains(ff)){
					continue;
				}
				double dd = f.getDistance(ff);
				mindist= Math.min(mindist,dd);
				if(dd < 1.5){
					crashes = true;
					break;
				}
				if(dd > 20){
					ignore.add(ff);
				}
			}
			//System.out.println("++"+mindist);
			if(!crashes){
				break;
			}
			if(ii == 0){
				if(!noatomflag){
					Point3D mp = new Point3D(basepoint);
					mp.x-=l.x;
					mp.y-=l.y;
					mp.z-=l.z;
					f.move(mp);
				}
			}
			f.move(mindirec);
			
		}
		//System.out.println("#########");
	}
	
	
	/**
	 * BackBone と SideChain のリストを Probability が高い順にソートする
	 */
	
	public void sortSamples(){
		for(String s:backbones.backbones.keySet()){
			ArrayList<BackBoneSample> ssal = backbones.backbones.get(s);
			Collections.sort(ssal,new SampleComparator());
			Collections.reverse(ssal);
		}
		
		for(String s:sidechains.sidechains.keySet()){
			ArrayList<SideChainSample> ssal = sidechains.sidechains.get(s);
			Collections.sort(ssal,new SampleComparator());
			Collections.reverse(ssal);
		}
		
		
	}
	/**
	 * Backbone と SideChain をランダムに変える。
	 * @param al
	 * @param index 
	 */
	
	
	public void changeRandom(ArrayList<FloatingResidue> al,int index){
		FloatingResidue r = al.get(index);
		if(Math.random() > 0.5){
			//r.processBackBone(backbones.getRandomly(r.getName()));
			r.processBackBone(backbones.getNext(r));
			r.fitSideChain();
		}else{
			//r.changeSideChain(sidechains.getRandomly(r.getName()));
			r.changeSideChain(sidechains.getNext(r));
		}
	}
	public static void saveOneChain(ArrayList<FloatingResidue> res,String filename){
		ArrayList<PDBResidue> res_ = new ArrayList<>();
		for(FloatingResidue r:res){
			res_.add(r);
		}
		writeToFile(makeAtomLines(res_,1,1,"A"),filename);	
	}
	
	public static void saveChains(ArrayList<FloatingResidue> res,String filename){
		ArrayList<PDBResidue> res_ = new ArrayList<>();
		int acount = 1;
		for(FloatingResidue r:res){
			res_.add(r);
			for(PDBAtom a:r.atoms){
				a.atom_number = String.valueOf(acount);
				acount++;
			}
		}
		writeToFile(makeAtomLines(res_,-1,-1,"A",true),filename);	
	}
	
	public static void savePDBChains(ArrayList<PDBResidue> res_,String filename){
		int acount = 1;
		for(PDBResidue r:res_){
			for(PDBAtom a:r.atoms){
				a.atom_number = String.valueOf(acount);
				acount++;
			}
		}
		writeToFile(makeAtomLines(res_,-1,-1,"A",true),filename);	
	}
	
	public static void main_(String[] args){
		ChainBuilder cb = new ChainBuilder();
		cb.sortSamples();
		ArrayList<String> al = new ArrayList<>();
		String[] str = "GG".split("");
		for(String ss:str){
			if(ss.length() > 0){
				al.add(ss);
			}
		}
		ArrayList<FloatingResidue> res = cb.buildDummy(al);
		for(FloatingResidue r:res){
			r.processBackBone(cb.backbones.getNext(r));
		}
		FloatingResidue chk = res.get(0);
		chk.getN().loc.set(0,0,0);
		chk.getCA().loc.set(1,1,0);
		chk.getC().loc.set(2,0,0);
		chk.nextn.loc.set(3,1,0);
		chk.prevc.loc.set(-1,1,0);
		for(FloatingResidue r:res){
			if(r != chk){
				r.processBackBone(null);
			}
		}
		saveOneChain(res,"C:\\Users\\kimidori\\Desktop\\TBM\\test.pdb");	
	}
	
	
	
	/**
	 * r1 の backbone と Sidechain を変更して最もスコアの高い組み合わせを探す
	 * @param r1
	 * @param distancepenalty
	 */
	public double checkOne(FloatingResidue r1,FloatingResidue r2,
			PDBAtom penaltybaseatom,//これと penaltytarget の距離に distancepenalty が掛けられる
			Point3D penaltytarget,//PrevC nextN どちらかが望ましい
			double distancepenalty){
		double mindist = Double.MAX_VALUE;
		BackBoneSample r1_backboneid = null;
		double r1_omega = 0;
		if(distancepenalty > 0){
			distancepenalty *= -1;
		}
		ArrayList<BackBoneSample> r1_bss = new ArrayList<>(this.backbones.backbones.get(r1.getName()));
		ArrayList<SideChainSample> r1_sc = new ArrayList<>(this.sidechains.sidechains.get(r1.getName()));
		ArrayList<OmegaSample> r1_omegas = null;
		if(r1.prev == r2){
			r1_omegas = new ArrayList<>(r1_bss.get(0).prevOmega.omegas);
		}else if(r1.next == r2){
			r1_omegas = new ArrayList<>(r1_bss.get(0).nextOmega.omegas);
		}else{
			throw new RuntimeException("second arg should be next or previous residue.");
		}
		int itenum = 1000;
		
		SideChainSample currentsc = null;
		//int r1index = r1_bss.size()/2;
		//int r1index_sc = r1_sc.size()/2;
		//int r1index_omega = r1_omegas.size()/2;
		int r1index = r1_bss.size();
		int r1index_sc = r1_sc.size();
		int r1index_omega = r1_omegas.size();
		
		double maxscore = Double.MAX_VALUE/2*-1;
		
		int itecount = 0;
		
		int r1index_h = r1index;
		int r1index_omega_h = r1index_omega;
		int r1scindex_h = r1index_sc;
		ArrayList<PDBResidue> target = new ArrayList<>();
		target.add(r1);
		while(itecount < itenum){
			itecount++;
			
			if(itecount%3 == 0){
				r1index_omega--;
			}
			if(itecount%3 == 1){
				r1index--;
			}
			if(itecount%3 == 2){
				r1index_sc--;
			}
			if(r1index < 0){
				r1index = r1_bss.size()-1;
			}
			if(r1index_omega < 0){
				r1index_omega = r1_omegas.size()-1;
			}
			if(r1index_sc < 0){
				r1index_omega = r1_sc.size()-1;
			}
			
			r1index = Math.min(Math.max(0,r1index),r1_bss.size()-1);
			r1index_omega = Math.min(Math.max(0,r1index_omega),r1_omegas.size()-1);
			r1index_sc = Math.min(Math.max(0,r1index_omega),r1_sc.size()-1);

			OmegaSample oo = r1_omegas.get(r1index_omega);
			r1.processBackBone(r1_backboneid,oo.value,r1.prev == r2,true);
			
			SideChainSample sid = sidechains.sidechains.get(r1.getName()).get(r1index_sc);
			if(sid != currentsc){
				r1.changeSideChain(sid);
				currentsc = sid;
			}
//Prepare ちゃんと入ってないかも
			 FuzzyDecisionTreeScoring_generator.prepareAll(_res,scoring);
			double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(target, scoring);
			double penal = calcAtomCrash_Higher_is_better(_res,r1);
			double dist = 0;
			if(penaltybaseatom != null){
				dist = penaltybaseatom.distance(penaltytarget);
			}
			
			double sc = scores[0]+penal+dist*distancepenalty;
			if(sc > maxscore){
				maxscore = sc;
				
				
				
				r1index_omega_h = r1index_omega;
				r1index_h = r1index;
				r1scindex_h = r1index_sc;
				
				if(itecount%3 != 0){
					r1index_omega = r1index_omega_h-3;
				}else{
				}
				if(itecount%3 != 1){
					r1index = r1index_h+3;
				}else{
				}
				if(itecount%3 != 2){
					r1index_sc = r1scindex_h +3;
				}else{
				}
				
			}
			//System.out.println(maxscore);
		}
		r1.processBackBone(r1_bss.get(r1index_h),r1_omegas.get(r1index_omega_h).value,r1.prev == r2,true);
		//System.out.println(r1.getName());
		SideChainSample sid = r1_sc.get(r1scindex_h);
		if(sid != currentsc){
			r1.changeSideChain(sid);
		}
		
		
		
		FuzzyDecisionTreeScoring_generator.prepareAll(_res,scoring);
		double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(target, scoring);
		double penal = calcAtomCrash_Higher_is_better(_res,r1);
		double dist = 0;
		if(penaltybaseatom != null){
			dist = penaltybaseatom.distance(penaltytarget);
		}

		double sc = scores[0]+penal+dist*distancepenalty;
		return sc;
	}
	
	
	/**
	 * Omega および BackBoneSample をチェックし、最も prevc, currentn と対応する原子の距離が近くなるよう変更する
	 * r1 > 0, r2 < length-1 でなければならない。
	 * @param r1
	 * @param r2 
	 */
	public void connect2(FloatingResidue r1,FloatingResidue r2){
		double mindist = Double.MAX_VALUE;
		BackBoneSample r1_backboneid = null;
		BackBoneSample r2_backboneid = null;
		double r1_omega = 0;
		double r2_omega = 0;
		
		ArrayList<BackBoneSample> r1_bss = new ArrayList<>(this.backbones.backbones.get(r1.getName()));
		ArrayList<BackBoneSample> r2_bss = new ArrayList<>(this.backbones.backbones.get(r2.getName()));
		System.out.println(r1.getCA().loc.x+";"+r2.getCA().loc.x);
		r1_bss.get(0).prevOmega.filt(4);
		r2_bss.get(0).nextOmega.filt(4);
		ArrayList<OmegaSample> r1_prevomegas = new ArrayList<>(r1_bss.get(0).prevOmega.omegas);
		ArrayList<OmegaSample> r2_omegas = new ArrayList<>(r2_bss.get(0).nextOmega.omegas);
		
		
		Point3D r1_prevvec =  PepProcess.subtracted(r1.prev.getC().loc,r1.prev.nextn.loc);
		Point3D r2_nextvec =  PepProcess.subtracted(r2.next.prevc.loc,r2.next.getN().loc);
		Point3D r2nn = r2.next.getN().loc;
		int itenum = 1000;
		
		Point3D r1c = r1.getC().loc;
		Point3D r1n = r1.getC().loc;
		
		Point3D r2c = r2.prevc.loc;
		Point3D r2n = r2.getN().loc;
		/*
		int r1index = 100000;
		int r1index_omega = 100000;
		int r2index = 100000;
		int r2index_omega = 100000;
		*/
		
		int r1index = r1_bss.size()/2;
		int r1index_omega = r1_prevomegas.size()/2;
		int r2index =  r2_bss.size()/2;
		int r2index_omega = r2_omegas.size()/2;
		
		
		int itecount = 0;
		
		int r1index_h = r1index;
		int r1index_omega_h = r1index_omega;
		int r2index_h = r2index;
		int r2index_omega_h = r2index_omega;
		
		while(itecount < itenum){
			itecount++;
			
			if(r1index_omega < 1){
				r1index_omega = r1index_omega_h+5;
			}
			if(r2index_omega < 1){
				r2index_omega = r2index_omega_h+5;
			}
			if(r2index < 1){
				r2index = r2index_h+5;
			}
			if(r1index < 1){
				r1index = r1index_h+5;
			}
			
			
			r1index_omega += 5;
			r1index += 5;
			r2index_omega += 5;
			r2index += 5;
			
			int pitenum = 40;
			for(int ii = 0;ii < pitenum;ii++){
				if(itecount%2 == 0){
					//前方を変える
					if(ii%2 == 0){
						r1index--;
					}else{
						r1index_omega--;
					}
					
					r1index = Math.min(Math.max(0,r1index),r1_bss.size()-1);
					r1index_omega = Math.min(Math.max(0,r1index_omega),r1_prevomegas.size()-1);

					BackBoneSample b1 = r1_bss.get(r1index);
					Point3D nextn = new Point3D(b1.atoms[BackBoneSample.NEXTN].loc);
					Point3D currentc = new Point3D(b1.atoms[BackBoneSample.C].loc);
					ArrayList<Point3D> ploc = new ArrayList<>();
					ploc.add(nextn);
					ploc.add(currentc);
					PepProcess.adjustVector3D(
							b1.atoms[BackBoneSample.N].loc,b1.atoms[BackBoneSample.PREVC].loc,b1.atoms[BackBoneSample.CA].loc
							,r1.prev.nextn.loc,	r1.prev.getC().loc,	r1.prev.getCA().loc,
							ploc
					);
					OmegaSample oo = r1_prevomegas.get(r1index_omega);
					double dd = oo.value;
					for(Point3D p:ploc){
						p.x -= r1.prev.nextn.loc.x;
						p.y -= r1.prev.nextn.loc.y;
						p.z -= r1.prev.nextn.loc.z;
						Point3D ccc = Point3D.rotate(p,r1_prevvec,dd);
						p.set(ccc);
						p.x += r1.prev.nextn.loc.x;
						p.y += r1.prev.nextn.loc.y;
						p.z += r1.prev.nextn.loc.z;
					}

					double ddist = r2c.distance(currentc)+r2n.distance(nextn);
					if(ddist < mindist){
						r1_backboneid = b1;
						r1_omega = dd;
						r1c = currentc;
						r1n = nextn;
						if(r1_backboneid != null && r2_backboneid != null){
							mindist = ddist;
						}
						
						if(ii%2 == 0){
							r1index_h = r1index;
							r1index_omega = r1index_omega_h;
						}else{
							r1index_omega_h = r1index_omega;
							r1index = r1index_h;
						}
						
						r2index_omega = r2index_omega_h+10;
						r2index = r2index_h+10;
						break;
					}
				}else{
					if(ii%2 == 0){
						r2index--;
					}else{
						r2index_omega--;
					}
					
					
					r2index = Math.min(Math.max(0,r2index),r2_bss.size()-1);
					r2index_omega = Math.min(Math.max(0,r2index_omega),r2_omegas.size()-1);


					BackBoneSample b2 = r2_bss.get(r2index);
					Point3D prevc = new Point3D(b2.atoms[BackBoneSample.PREVC].loc);
					Point3D currentn = new Point3D(b2.atoms[BackBoneSample.N].loc);
					ArrayList<Point3D> nloc = new ArrayList<>();
					nloc.add(prevc);
					nloc.add(currentn);

					PepProcess.adjustVector3D(
							b2.atoms[BackBoneSample.NEXTN].loc,
							b2.atoms[BackBoneSample.C].loc,
							b2.atoms[BackBoneSample.CA].loc

							,r2.next.getN().loc,r2.next.prevc.loc,	r2.next.getCA().loc,
							nloc
					);
					double dd2 = r2_omegas.get(r2index_omega).value;
					for(Point3D p:nloc){
						p.x -= r2nn.x;
						p.y -= r2nn.y;
						p.z -= r2nn.z;
						Point3D ccc = Point3D.rotate(p,r2_nextvec,dd2);
						p.set(ccc);
						p.x += r2nn.x;
						p.y += r2nn.y;
						p.z += r2nn.z;
					}
					double ddist = r1n.distance(currentn)+r1c.distance(prevc);
					if(ddist < mindist ){
						r2_backboneid = b2;
						r2_omega = dd2;
						r2c = prevc;
						r2n = currentn;
						if(r1_backboneid != null && r2_backboneid != null){
							mindist = ddist;
						}
						
						if(ii%2 == 0){
							r2index_h = r2index;
							r2index_omega = r2index_omega_h;
						}else{
							r2index_omega_h = r2index_omega;
							r2index = r2index_h;	
						}
						r1index_omega = r1index_omega_h+10;
						r1index = r1index_h+10;
						break;
					}
				}
			}
			System.out.println(r1index+";"+r2index+";"+r1index_omega+";"+r2index_omega+"; "+mindist+"++");
		}
		if(r2_backboneid != null && r1_backboneid != null){
			r1.processBackBone(r1_backboneid,r1_omega,true,true);
			r2.processBackBone(r2_backboneid,r2_omega,false,true);
		}else{
			throw new RuntimeException("debug");
		}
	}
	
	
	
	public static void main___(String[] args){
		ChainBuilder cb = new ChainBuilder();
		cb.sortSamples();
		ArrayList<String> al = new ArrayList<>();
		String[] str = "GAAG".split("");
		for(String ss:str){
			if(ss.length() > 0){
				al.add(ss);
			}
		}
		
		ArrayList<FloatingResidue> res = cb.buildDummy(al);
		int cou = 0;
		
		cb.prepare(res);
		for(FloatingResidue r:res){
			r.processBackBone(cb.backbones.getNext(r));
			r.changeSideChain(cb.sidechains.getNext(r));
		}
		
		for(int ii = 0;ii < res.get(0).atoms_backbone_loc.size();ii++){
			res.get(3).atoms_backbone_loc.get(ii).set(res.get(0).atoms_backbone_loc.get(ii));
		}
		for(int ii = 0;ii < res.get(0).atoms_backbone_loc.size();ii++){
			res.get(3).atoms_backbone_loc.get(ii).x += 2.5;
			res.get(3).atoms_backbone_loc.get(ii).y -= 7.5;
			res.get(3).atoms_backbone_loc.get(ii).z += 1.5;
		}
		
		saveOneChain(res,"test1.pdb");
		cb.checkOne(res.get(2),res.get(3),res.get(2).prevc,res.get(1).getC().loc,-0.5);
		saveOneChain(res,"test2.pdb");
		cb.checkOne(res.get(1),res.get(0),res.get(1).nextn,res.get(2).getN().loc,-0.5);
		saveOneChain(res,"test3.pdb");
		cb.checkOne(res.get(2),res.get(3),res.get(2).prevc,res.get(1).getC().loc,-0.5);
		saveOneChain(res,"test4.pdb");
		cb.checkOne(res.get(1),res.get(0),res.get(1).nextn,res.get(2).getN().loc,-0.5);
		saveOneChain(res,"test5.pdb");
	}
	public static void main__(String[] args){
		ChainBuilder cb = new ChainBuilder();
		cb.sortSamples();
		ArrayList<String> al = new ArrayList<>();
		String[] str = "AAA".split("");
		for(String ss:str){
			if(ss.length() > 0){
				al.add(ss);
			}
		}
		
		ArrayList<FloatingResidue> res = cb.buildDummy(al);
		for(FloatingResidue r:res){
			r.processBackBone(cb.backbones.getNext(r));
		}
		cb.prepare(res);
		ArrayList<PDBAtom> atoms = new ArrayList<>();
		ArrayList<PDBResidue> res_ = new ArrayList<>();
		for(FloatingResidue r:res){
			res_.add(r);
			atoms.addAll(r.atoms);
		}
		
		ArrayList<VSorter> sorter = new ArrayList<>();
		double prevscore = -1000000;
		
		for(int ii = 0;ii < res.size();ii++){
			sorter.add(new VSorter(0,0));
		}
		for(int ii = 0;ii < res.size();ii++){
			FloatingResidue r = res.get(ii);
			r.processBackBone(cb.backbones.getNext(r));
			r.fitSideChain();
			cb.moveToClosest(res,r);
		}
		int maxite = 5;
		for(int jj = 0;jj < maxite;jj++){
			
			for(int ii = 0;ii < res.size();ii++){
				res.get(ii).saveLoc();
			}
			double sc = 0;
			double[] scores = null;

			if(jj == (int)(maxite*0.8)){
				sc = 0;
				scores = cb.calcChainScore(res_);
				for(int ii = 0;ii < res.size();ii++){
					sc += scores[ii];
					if(jj >= (int)(maxite*0.8)){
						double cc = cb.calcAtomCrash_Higher_is_better(res_, res.get(ii));
						sc += cc;
						sorter.get(ii).val = scores[ii]+cc;
					}else{
						sorter.get(ii).val = scores[ii];
					}
					sorter.get(ii).index = ii;
				}
				prevscore = sc;
			}
			if(jj < 2000){
				for(int ii = 0;ii < res.size();ii++){
					cb.changeRandom(res, ii);
				}
			}
			
			sc = 0;
			scores = cb.calcChainScore(res_);
			for(int ii = 0;ii < res.size();ii++){
				sc += scores[ii];
				if(jj >= (int)(maxite*0.8)){
					double cc = cb.calcAtomCrash_Higher_is_better(res_,res_.get(ii));
					sc += cc;
					sorter.get(ii).val = scores[ii]+cc;
				}else{
					sorter.get(ii).val = scores[ii];
				}
				sorter.get(ii).index = ii;
			}
			
			Collections.sort(sorter,new VComparator());;
			for(int ii = 0;ii < res.size();ii++){
				if(ii < res.size()*Math.max(0.5,0)){
					cb.changeRandom(res, sorter.get(ii).index);
				}
			}
			for(int ii = 0;ii < res.size();ii++){
				FloatingResidue f = res.get(ii);
				f.processBackBone(null);
				f.fitSideChain();
				double maxscore = cb.calcAtomCrash_Higher_is_better(res_, res_.get(ii));
				int maxindex = f.sidechainIndex;
				if(maxscore < 0){
					int maxcount =cb.sidechains.sidechains.get(f.getName()).size();
					for(int kk = 0;kk < maxcount;kk++){
						f.changeSideChain(cb.sidechains.getNext(f));
						double scc = cb.calcSideChainCrash(res_, ii);
						if(scc >= 0){
							maxindex = f.sidechainIndex;
							break;
						}
						if(scc > maxscore){
							maxscore = scc;
							maxindex = f.sidechainIndex;
						}
					}
					f.sidechainIndex = Math.max(maxindex -1,-1);
					f.changeSideChain(cb.sidechains.getNext(f));
				}
			}
			
			
			sc = 0;
			scores = cb.calcChainScore(res_);
			
			for(int ii = 0;ii < res.size();ii++){
				sc += scores[ii];
				if(jj >= (int)(maxite*0.8)){
					double cc = cb.calcAtomCrash_Higher_is_better(res_, res_.get(ii));
					sc += cc;
				}
			}
			System.out.println(prevscore+":"+sc);
			if(sc < prevscore){
				for(int ii = 0;ii < res.size();ii++){
					res.get(ii).restoreLoc();
				}
			}else{
				prevscore = sc;
				
			}
		}
		writeToFile(makeAtomLines(res_,1,1,"A"),"C:\\Users\\kimidori\\Desktop\\TBM\\test.pdb");	
	}
	
	
	/**
	 * C 末方向に伸ばす
	 * @param f 
	 */
	public void elongTerm(FloatingResidue target,
			boolean forward,
			int checkindex
			,boolean checkomega0){
		
		ArrayList<PDBResidue> al = new ArrayList<>();
		al.add(target);
		ArrayList<Double> omegas = new ArrayList<>();
		omegas.add(Math.PI);
		if(checkomega0){
			omegas.add(0.0);
		}
		double maxomega = Math.PI;
		double maxscore = -100;
		int maxindex = 0;

		
		
		for(Double omega:omegas){
			for(int ii = 0;ii < checkindex;ii++){
				if(forward){
					
					//前の CN にマップしつつバックボーン変更
					target.processBackBone(backbones.getOneAt(target,ii),
							omega,true,true);
				}else{
					//後ろの CN にマップしつつバックボーン変更
					target.processBackBone(backbones.getOneAt(target,ii),
							omega,false,true);

				}
				double pp = calcAtomCrash_Higher_is_better(_res,target);
				if(pp < 0){
					target.sidechainIndex = -1;
					for(int kk = 0;kk < 20;kk++){
						target.changeSideChain(sidechains.getNext(target));
						pp = calcAtomCrash_Higher_is_better(_res,target);
						if(pp == 0){
							break;
						}
					}
				}
				double scc[] = FuzzyDecisionTreeScoring_generator.calcResidueScores(al,scoring);
				double sc = pp+scc[0];
				if(sc > maxscore){
					maxscore = sc;
					maxomega = omega;
					maxindex = ii;
				}
			}
		}
		
		
		///System.out.println(maxscore+";;;");
		if(forward){
			target.processBackBone(backbones.getOneAt(target,maxindex),
					maxomega,true,true);
		}else{
			target.processBackBone(backbones.getOneAt(target,maxindex),
					maxomega,false,true);
		}
		//System.out.println(target.getC().distance(target.getCA()));
		//System.out.println(target.getC().distance(target.getCA()));
	}
	
	
	
	public  ArrayList<FloatingResidue> generateChain(ArrayList<String> threecode){
		ArrayList<FloatingResidue> fl = new ArrayList<>();
		int ii = 0;
		for(String ss:threecode){
			FloatingResidue fr = new FloatingResidue(ss);
			fr.setResidueNumber(ii+1);
			fl.add(fr);
			if(ii > 0){
				fr.changeSideChain(sidechains.getNext(fr));
				fr.setPrevFloating(fl.get(ii-1));
				fr.setOmega(Math.PI, true);
			}else{
				fr.changeSideChain(sidechains.getNext(fr));
				
			}
			//System.err.println(ii);
			fr.backboneIndex = -1;
			fr.processBackBone(backbones.getNext(fr));
			ii++;
		}
		for(ii = 0;ii < threecode.size()-1;ii++){
			FloatingResidue fr = fl.get(ii);
			fr.setNextFloating(fl.get(ii+1));
		}
		return fl;
	}
	
	
	
	
	public  ArrayList<FloatingResidue> generateChain(int length,String code){
		ArrayList<FloatingResidue> fl = new ArrayList<>();
		for(int ii = 0;ii < length;ii++){
			FloatingResidue fr = new FloatingResidue(code);
			fr.setResidueNumber(ii+1);
			fl.add(fr);
			if(ii > 0){
				fr.changeSideChain(sidechains.getNext(fr));
				fr.setPrevFloating(fl.get(ii-1));
				fr.setOmega(Math.PI, true);
			}else{
				fr.changeSideChain(sidechains.getNext(fr));
				
			}
			//System.err.println(ii);
			fr.backboneIndex = -1;
			fr.processBackBone(backbones.getNext(fr));
		}
		for(int ii = 0;ii < length-1;ii++){
			FloatingResidue fr = fl.get(ii);
			fr.setNextFloating(fl.get(ii+1));
		}
		return fl;
	}
	
	/**
	 * クエリの方にインサーションになっている部分について処理する。
	 * @param al 
	 */
	public void processInsersion(ArrayList<FloatingResidue> al,HashSet<FloatingResidue> mapped){
		if(al.size() < 4){
			return;
		}
		for(FloatingResidue fr:al){
			fr.processBackBone(backbones.getRandomly(fr));
		}
		System.out.println(al.size()+";;;;");
		FloatingResidue first = al.get(0);
		//スコアリングについて Prepare ちゃんとできてないかも
		FloatingResidue last = al.get(al.size()-1);
		if(first.prev != null && last.next != null){
			double ddist = first.prev.getC().distance(last.next.getN());
			ArrayList<ArrayList<FloatingResidue>> fl 
					= generateRandomLoop(al,ddist,5.0,40,1,100);
			if(fl.size() == 0){
				makeDecoy(al,100);
			}else{

				for(int ii = 0;ii < fl.get(0).size();ii++){
					al.get(ii).fitAtoms(fl.get(0).get(ii));
				}
			}
		}else{
			makeDecoy(al,100);
		}
		ArrayList<PDBResidue> template = new ArrayList<>();
		ArrayList<PDBResidue> tar = new ArrayList<>();
		template.addAll(mapped);
		tar.addAll(al);
		ArrayList<ArrayList<PDBResidue>> pl = docker.dock(template
		, tar, 20,5);
		if(pl.size () == 0){
			System.out.println();
		}
		double cndist = 100000;
		int minindex = 0;
		double okthreshold = 5;
		for(int jj = 0;jj < pl.size();jj++){
			ArrayList<PDBResidue> p = pl.get(jj);
			double ddd = 0;
			if(first.prev != null){
				ddd += p.get(0).getN().distance(first.prev.getC());
			}
			if(last.next != null){		
				ddd += p.get(p.size()-1).getC().distance(last.next.getN())+3.0;//FixMe 結合距離ちゃんと追加
			}
			if(ddd < cndist){//FixMe ok threshold 以内でスコアが高いものに変更
				minindex = jj;
				cndist = ddd;
			}
		}

		for(int ii = 0;ii < al.size();ii++){
			al.get(ii).fitAtoms(pl.get(minindex).get(ii));
		}
		
	}
	
	/**
	 * target-next NC 距離からギャップを埋めるのに必要な残基を見積もる
	 * @param target
	 * @param next
	 * @param omega
	 * @param checkindex
	 * @return 
	 */
	public ArrayList<FloatingResidue> getGapFillingResidues(FloatingResidue target){
		
		ArrayList<FloatingResidue> ret = new ArrayList<>();
		ret.add(target);
		if(target.next == null || target.prev == null){
			return ret;
		}
		double ddist = target.prev.getC().distance(target.next.getN());
		double dist_oneres = this.distPenalty.dist_prevc_nextn[AtomDistancePenalty.MEAN5_95];
		FloatingResidue dprev = target.prev;
		FloatingResidue dnext = target.next;
		
		
		while(true){
			if(dprev == null || dnext == null){
				break;
			}
			ret.add(dprev);
			ret.add(dnext);
			if((ret.size()-2)*dist_oneres > dprev.getC().distance(dnext.getN())+2.6){
				break;
			}
			dprev = dprev.prev;
			dnext = dnext.next;
			if(dprev == null || dnext == null){
				break;
			}
		}
		return ret;
	}
	
	
	/**
	 * target および next の間を埋めようとする
	 * @param target
	 * @param omega
	 * @param prevflag
	 * @param checkindex 
	 */
	
	public void fillGapBetween(FloatingResidue target,
			FloatingResidue next,
			double omega,int checkindex){
		
		
		if(target.next == null){
			elongTerm(target,
			true,100,false);
			return;
		}else if(target.prev == null){
			elongTerm(target,
			false,100,false);
			return ;
		}
		double ddist = target.prev.getC().distance(target.next.getN());
		double dist_oneres = (this.distPenalty.dist_prevc_nextn[AtomDistancePenalty.MAX95]*0.5
				+this.distPenalty.dist_prevc_nextn[AtomDistancePenalty.MEAN5_95]*0.5);
		int resnum = (int)(ddist/dist_oneres+0.5);
		ArrayList<FloatingResidue> al = new ArrayList<>();
		al.add(target);
		while(al.size() < resnum){
			FloatingResidue t = al.get(0);
			if(t.prev != null){
				al.add(0,t.prev);
			}
			if(al.get(0).prev == null){	//N 末全部使用
				break;
			}
			ddist = al.get(0).getC().distance(al.get(al.size()-1).next.getN());
			resnum = (int)(ddist/dist_oneres+0.5);
			if(al.size() >= resnum){
				break;
			}	
			FloatingResidue tn = al.get(al.size()-1);
			if(tn.next != null){
				al.add(tn.next);
			}
			if(al.get(al.size()-1).next == null){//C 末全部使用
				break;
			}
			ddist = al.get(0).getC().distance(al.get(al.size()-1).next.getN());
			resnum = (int)(ddist/dist_oneres+0.5);
		}
		ArrayList<FloatingResidue> store = new ArrayList<>(res);
		ArrayList<FloatingResidue> dummy = new ArrayList<>(res);
		ArrayList<FloatingResidue> dal = generateDummyG(al);
		dummy.removeAll(al);
		dummy.addAll(dal);
		this.prepare(dummy);
		
		ArrayList<AtomTriangle> zal = makeTriangles(dal);
		waveAllTriangle(zal,1000,0.5*dal.size());
		fitToTriangle(dal, al);
		this.prepare(store);
		
	}
	
	
	
	
	
	
	/**
	 * target は Next と Prev がある Residue で、最も nextn nextc が近くなる BackBone で停止する
 omega は通常
	 * @param target 
	 */
	
	public void fillGapWith(FloatingResidue target,double omega
			,boolean forward,int checkindex
			,double penaltythreshold){
		
		
		if(target.next == null){
			elongTerm(target,
			true,100,false);
		}else if(target.prev == null){
			elongTerm(target,
			false,100,false);
		}
		if(target.next == null ||  target.prev == null){
			return;
		}
		if(penaltythreshold > 0){
			penaltythreshold *= -1;
		}
		double mindiff = Double.MAX_VALUE;
		double penalty = -10000;
		int minindex = 0;
		ArrayList<PDBResidue> al = new ArrayList<>();
		al.add(target);
		
		for(int ii = 0;ii < checkindex;ii++){
			double dd = 0;
			if(forward){
				target.processBackBone(backbones.getOneAt(target,ii),
						omega,true,true);
				dd = target.nextn.distance(target.next.getN())
					+target.getC().distance(target.next.prevc);
			}else{
				target.processBackBone(backbones.getOneAt(target,ii),
						omega,false,true);
				dd = target.prevc.distance(target.prev.getC())
					+target.getN().distance(target.prev.nextn);
				
			}
			double pp = calcAtomCrash_Higher_is_better(al,target);
			if(dd < mindiff && pp > penaltythreshold ){
				mindiff = dd;
				minindex = ii;
				penalty = pp;
			}
		}
		//逆方向に置こうかと思ったがなぜか距離が全然違った。Fixme
		//double domega = PepProcess.omega(target,target.next)/180*Math.PI;
		//target.processBackBone(backbones.getOneAt(target,minindex),
		//		domega,false,true);
		//ChainBuilder.saveChains(res,"test1a.pdb");
		//System.out.println(domega+";"+(PepProcess.omega(target,target.next)/180*Math.PI));
		
		if(forward){
			target.processBackBone(backbones.getOneAt(target,minindex),
					omega,true,true);
		}else{
			target.processBackBone(backbones.getOneAt(target,minindex),
					omega,false,true);
			
		}
		
		//System.out.println(omega+";"+(PepProcess.omega(target.prev,target)/180*Math.PI));
		//ChainBuilder.saveChains(res,"test1b.pdb");
		//double ddd = target.prevc.distance(target.prev.getC())
		//		+target.getN().distance(target.prev.nextn);
		//System.out.println(ddd);
	}
	
	/**
	 * 前後に残基がある 3 残基までのアミノ酸をフィットさせる。
	 * 無い場合 elongterm で伸ばす。
	 * @param fl
	 * @param distance
	 * @param torelance
	 * @param strnum
	 * @param indexcheck 
	 */
	
	public boolean fitThreeResidue(ArrayList<FloatingResidue> fl
	,ChainBuilder cb//this でもよいと思うが
	,double torelance
	,int iternum){
		FloatingResidue first = fl.get(0);
		FloatingResidue last = fl.get(fl.size()-1);
		
		FloatingResidue prev = first.prev;
		FloatingResidue next = last.next;
		ArrayList<Double> omegas = new ArrayList<>();
		omegas.add(Math.PI);
		
		
		double crashpenal_saved = -100000;
		double dist_saved = 10000;
		boolean forward = prev != null;
		for(FloatingResidue target:fl){
			target.sidechainIndex = -1;
		}
		for(FloatingResidue target:fl){
			if(forward){
				//前の CN にマップしつつバックボーン変更
				target.processBackBone(backbones.getNext(target),
						Math.PI,true,true);
			}else{
				//後ろの CN にマップしつつバックボーン変更
				target.processBackBone(backbones.getNext(target),
						Math.PI,false,true);
			}
		}
		
		for(int ii = 0;ii < iternum;ii++){
			for(FloatingResidue target:fl){
				if(ii != 0 && Math.random() > 0.3){
					target.processBackBone(backbones.getNext(target),
							Math.PI,(target.prev != null),true);
				}
				for(Double omega:omegas){
					//if(forward && target.prev != null){
					//	//前の CN にマップしつつバックボーン変更
					//	target.processBackBone(backbones.getOneAt(target,ii),
					//			omega,true,true);
					//}else{
						//後ろの CN にマップしつつバックボーン変更
					//	target.processBackBone(backbones.getOneAt(target,ii),
					//			omega,false,true);

					//}
					for(FloatingResidue ff:fl){
						ff.fitBackbone(forward);
					}
					
					double pp = calcAtomCrash_Higher_is_better(_res,target);
					if(prev != null && next != null){
						double ddist = prev.nextn.distance(first.getN())
								+ next.prevc.distance(last.getC())
								+ first.prevc.distance(prev.getC())
								+last.nextn.distance(next.getN())
								;
						/*
						System.out.println(prev.nextn.distance(first.getN()));
						System.out.println(next.prevc.distance(last.getC()));
						System.out.println(first.prevc.distance(prev.getC()));
						System.out.println(last.nextn.distance(next.getN()));
						System.out.println(fl.size()+";;;"+ddist);
						*/
						
						if(ddist > torelance && dist_saved > ddist){
							continue;
						}
						
						if(pp < 0){
							target.sidechainIndex = -1;
							for(int kk = 0;kk < 20;kk++){
								target.changeSideChain(sidechains.getNext(target));
								pp = calcAtomCrash_Higher_is_better(_res,target);
								if(pp == 0){
									break;
								}
							}
						}
						if(pp > crashpenal_saved){
							crashpenal_saved = pp;
							dist_saved = ddist;
							for(FloatingResidue ff:fl){
								ff.saveLoc();
							}
						}
						
						if(pp == 0 && ddist <= torelance ){
							return true;
						}
					}else{
						
						if(pp < 0){
							target.sidechainIndex = -1;
							for(int kk = 0;kk < 20;kk++){
								target.changeSideChain(sidechains.getNext(target));
								pp = calcAtomCrash_Higher_is_better(_res,target);
								if(pp == 0){
									break;
								}
							}
						}
						if(pp > crashpenal_saved){
							crashpenal_saved = pp;
							for(FloatingResidue ff:fl){
								ff.saveLoc();
							}
						}
						if(pp == 0){
							return true;
						}
					}
					
				}
			}
		}
		
		for(FloatingResidue ff:fl){
			ff.restoreLoc();
		}
		return false;
	}
	
	
	
	/**
	 * N 末の N と C 末の C が
	 * distance+-dist_torelance の長さのチェーンを作る
	 * @param fl
	 * @param distance
	 * @param dist_torelance
	 * @param crashpenalty_threshold
	 * @param strnum
	 * @param indexcheck 
	 */
	public ArrayList<ArrayList<FloatingResidue>>  generateRandomLoop(ArrayList<FloatingResidue> fl
	,double distance
	,double dist_torelance
	,double crashpenalty_threshold
	,int strnum
	,int indexcheck){
		int iterthreshold = strnum*10;
		int iter = 0;
		if(	crashpenalty_threshold > 0){
			System.err.println("Penalty should be negative. I did *-1");
			crashpenalty_threshold *= -1;
		}
		ArrayList<ArrayList<FloatingResidue>> ret = new ArrayList<>();
		ChainBuilder cb = new ChainBuilder();
		ArrayList<String> pep = new ArrayList<>();
		for(FloatingResidue f:fl){
			pep.add(f.getName());
		}
		PDBChain c = new PDBChain("A");
		while(ret.size() < strnum){
			ArrayList<FloatingResidue> res = cb.generateChain(pep);
			for(FloatingResidue r:res){
				r.setParent(c);
			}
			cb.prepare(res);
			FloatingResidue first = res.get(0);
			FloatingResidue last = res.get(res.size()-1);
			
			ArrayList<FloatingResidue> ures = new ArrayList<>();
			for(int ii = 0;ii < res.size();ii++){
				ures.add(res.get(ii));
				cb.prepare(ures);
				cb.elongTerm(res.get(ii),
				true,
				100
				,false);
					
			}
			for(int jj = 0;jj < res.size();jj++){
				res.get(jj).changeSideChain(sidechains.sidechains.get(res.get(jj).getName()).get(0));
				res.get(jj).saveLoc();
			}
			double prevdiff = 10000;
			double prevpenal  = -100000;
			for(int ii = 0;ii < res.size()*10;ii++){
				int ri = (int)(Math.random()*res.size());
				FloatingResidue target = res.get(ri);
				target.processBackBone(backbones.getOneAt(target,ri),
						Math.PI,true,true);
				
				for(int jj = ii+1;jj < res.size();jj++){
					res.get(jj).fitBackbone(true);
				}
				double ddist = first.getN().distance(last.getC());
				double cdiff = Math.abs(ddist-distance);
				//System.out.println(cdiff+";;"+ddist);
				if(prevdiff < cdiff &&  Math.abs(ddist-distance) > dist_torelance){
					continue;
				}

				double penal = 0;
				for(int jj = 0;jj < res.size();jj++){
					res.get(jj).fitBackbone(true);
					double pp = calcBackBoneCrash_Higher_is_better(_res,res.get(jj));
					penal+=pp;
				}
				if(prevdiff > cdiff && penal >= prevpenal){
					for(int jj = 0;jj < res.size();jj++){
						res.get(jj).saveLoc();
					}
					prevdiff = cdiff;
					prevpenal = penal;
					System.out.println(prevpenal+";;"+prevdiff);
				}else{
					for(int jj = 0;jj < res.size();jj++){
						res.get(jj).restoreLoc();
					}
				}
				if(prevpenal >= crashpenalty_threshold && Math.abs(prevdiff-distance) < dist_torelance){
					ret.add(res);
					break;
				}
			}
			
			if(iter > iterthreshold){
				break;
			}
			iter++;
		}
		return ret;
	}
	
	/**
	 * ある一定距離にわたるペプチド鎖を作る。
	 * @param fl
	 * @param distance
	 * @param torelance 
	 */
	public void makeDistantChain(ArrayList<FloatingResidue> fl,double distance,double torelance){
		FloatingResidue first = fl.get(0);
		FloatingResidue last = fl.get(fl.size()-1);
		
		//for(int ii = 0;ii < fl.size();ii++){
		while(true){
			
			int ri = (int)(Math.random()*fl.size());
			FloatingResidue fr = fl.get(ri);
			fr.processBackBone(backbones.getNext(fr));
			for(int jj = ri+1;jj<fl.size();jj++){
				 fl.get(jj).fitBackbone(true);
				//refiner.waveBackBone(fr,true,true);
			}
			if(Math.abs(first.prevc.distance(last.nextn) -distance) < torelance){
				break;
			}
		}
		
		
	}
	/**
	 * 完全にランダムでチェーンを作成する
	 * @param fl
	 * @param distance
	 * @param torelance 
	 */
	public void makeDecoy(ArrayList<FloatingResidue> fl,int iternum){
		double maxscore = - 100000;
		for(FloatingResidue fr:fl){
			fr.processBackBone(backbones.getRandomly(fr));
		}
		for(int ii = 0;ii < iternum;ii++){
			
			int ri = (int)(Math.random()*fl.size());
			FloatingResidue fr = fl.get(ri);
			fr.processBackBone(backbones.getNext(fr));
			for(int jj = ri+1;jj<fl.size();jj++){
				 fl.get(jj).fitBackbone(true);
				//refiner.waveBackBone(fr,true,true);
			}
			double psum = 0;
			this.prepare(fl);
			for(FloatingResidue r:fl){
				refiner.maxRotamer(r, -5,0.1, true, sidechains, false);
			}
			for(FloatingResidue r:fl){
				double penal = distPenalty.calcCrashPenalty(_res, r,AtomDistancePenalty.TYPE_ALL);
				psum += penal;
			}
			double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(_res, scoring);
			for(double s:scores){
				psum+= s;
			}
			if(psum > maxscore){
				System.out.println("decoy:" +maxscore);
				for(FloatingResidue r:fl){
					r.saveLoc();
				}
				maxscore = psum;
				
			}else{
				
				for(FloatingResidue r:fl){
					r.restoreLoc();
				}
			}
		}
		
		
	}
	
	
	public static void buildUnmapped(ArrayList<FloatingResidue> allseq,HashSet<FloatingResidue> mapped_f,ResidueRefine refiner){
		//ここから
		//MCTS のようなことをさせる
		//2DSTRも入れるか
		//CRASH をのけて、PDP して塊として扱え
		
		
		ChainBuilder cb = new ChainBuilder();
		ArrayList<FloatingResidue> fr = new ArrayList<>();
		for(FloatingResidue f:allseq){
			if(mapped_f.contains(f)){
				fr.add(f);
			}
		}
		
		
		ArrayList<PDBResidue> _fr = new ArrayList<>();
		_fr.addAll(fr);
		refiner.prepare(_fr);
		int strnum = 1;
		for(int si = 0;si < strnum;si++){
			int mapped1 = -1;
			for(int rr = 0;rr < fr.size();rr++){
				if(!mapped_f.contains(fr.get(rr))){
					ArrayList<FloatingResidue> crashed = new ArrayList<>();
					crashed.add(fr.get(rr));
					for(int zz = rr+1;zz < fr.size();zz++){
						if(!mapped_f.contains(fr.get(zz))){
							crashed.add(fr.get(zz));
							rr = zz;
						}else{
							break;
						}
					}
					int mapped2 = -1;
					if(rr < fr.size()-1){
						mapped2 = rr+1;
					}
					ArrayList<FloatingResidue> nocrash = new ArrayList<>(fr);
					nocrash.removeAll(crashed);
					cb.prepare(nocrash);
					for(int ll = 0;ll < crashed.size();ll++){
						if(mapped2 > -1 && crashed.size()/2 < ll){
							break;
						}
						cb.elongTerm(crashed.get(ll), true, 10, true);
					}

					for(int ll = crashed.size()-1;ll >= 0;ll--){
						if(mapped1 > -1 && crashed.size()/2-1 > ll){
							break;
						}
						cb.elongTerm(crashed.get(ll), false, 10, true);
					}


					refiner.loopRefile(fr, crashed, cb, 1000,true);
				}else{
					mapped1 = rr;
				}
			}
		}
		refiner.prepare(_fr);
		for(FloatingResidue r:fr){
			refiner.maxRotamer(r,-5,0.1, true, cb.sidechains, false);
		}
		cb.fillGapsAll(fr,1000);
	}
	
	
	
	
	
	/**
	 * ループを作った後に両端が近くなるように移動する
	 * @param fl
	 * @param prev
	 * @param next 
	 */
	
	public void triangleFit(ArrayList<FloatingResidue> fl ,FloatingResidue prev,FloatingResidue next){
		FloatingResidue first = fl.get(0);
		FloatingResidue last = fl.get(fl.size()-1);
		ArrayList<Point3D> atoms_backbone  = new ArrayList<>();
		for(FloatingResidue f:fl){
			atoms_backbone.addAll(f.atoms_backbone_loc);
		}
		PepProcess.adjustVector3D(
		first.getN().loc,
		first.prevc.loc,
		last.nextn.loc,
		prev.nextn.loc,
		prev.getC().loc,
		next.getN().loc,
		atoms_backbone
		);
		
		for(FloatingResidue f:fl){
			f.fitSideChain();
		}
	}
	
	public static void main_zd(String[] args){
		PDBChain c = new PDBChain("A");
		ChainBuilder cb = new ChainBuilder();
		char[] aax = "GGLAIHLGGGGG".toCharArray();
		ArrayList<String> pep = new ArrayList<>();
		for(char a:aax){
			String th = PepProcess.one_to_three.get(a);
			if(th == null){
				System.out.println(a);
				System.exit(0);
			}
			pep.add(th);
			
		}
		
		ArrayList<FloatingResidue> res = cb.generateChain(pep);
		for(FloatingResidue r:res){
			r.setParent(c);
		}
		ArrayList<FloatingResidue> ures = new ArrayList<>();
		for(int ii = 0;ii < res.size();ii++){
			ures.add(res.get(ii));
			cb.prepare(ures);
			cb.elongTerm(res.get(ii),
			true,
			100
			,false);
			//cb.prepare(ures);
		}
		cb.prepare(ures);
		
		ChainBuilder.saveChains(ures,"test1.pdb");
		ures.remove(0);
		ures.remove(ures.size()-1);
		ArrayList<FloatingResidue> dal = cb.generateDummyG(ures);
		res.removeAll(ures);
		res.addAll(dal);
		cb.prepare(res);
		ArrayList<AtomTriangle> al = cb.makeTriangles(dal);
		cb.waveAllTriangle(al,10000,0.5*ures.size());
		cb.fitToTriangle(dal, ures);
		res.removeAll(dal);
		res.addAll(ures);
		ChainBuilder.saveChains(res,"test2.pdb");
	}
	
	
	
	
	public static void main(String[] args){
		
		for(int kk = 0;kk < 1;kk++){
			PDBChain c = new PDBChain("A");
			ChainBuilder cb = new ChainBuilder();
			char[] aax = "SATMAFHPMLCRLELSVSAAPPASPIDATLLRSLITSVL".toCharArray();
			ArrayList<String> pep = new ArrayList<>();
			for(char a:aax){
				String th = PepProcess.one_to_three.get(a);
				pep.add(th);
			}
			ArrayList<FloatingResidue> res = cb.generateChain(pep);
			for(FloatingResidue r:res){
				BackBoneSample bss = cb.backbones.getRandomly(r.getName());
				r.calcPseudoNeighbour(bss);
				System.out.println(r.nextn.distance(r.getC())-bss.atoms[BackBoneSample.NEXTN].distance(bss.atoms[BackBoneSample.C]));
				System.out.println(r.prevc.distance(r.getN())-bss.atoms[BackBoneSample.PREVC].distance(bss.atoms[BackBoneSample.N]));
				System.out.println();
			}
			
			
		}
	}
	
}

class FloatingResidue extends PDBResidue{
	PDBAtom prevc = null;
	PDBAtom nextn = null;
	FloatingResidue prev;
	FloatingResidue next;
	
	int sidechainIndex = -1;
	int backboneIndex = -1;
	
	ArrayList<Point3D> atoms_loc = new ArrayList<>();
	ArrayList<Point3D> atoms_backbone_loc = new ArrayList<>();
	HashMap<String,PDBAtom> atoms_sidechain = new HashMap<>();
	ArrayList<PDBAtom> atoms_backbone = new ArrayList<>();
	Point3D dummyPoint = new Point3D(0,0,0);
	ArrayList<Point3D> atoms_loc_all = new ArrayList<>();
	ArrayList<Point3D> atoms_loc_all_prev = new ArrayList<>();
	SideChainSample sideChain = null;
	
	
	double prevOmega = Math.PI;
	double nextOmega = Math.PI;
	
	FloatingResidue(){
	}
	
	FloatingResidue(FloatingResidue r){
		this(r.getName());
		this.parent = r.parent;
		this.setResidueNumber(r.getResidueNumber());
		this.setName(r.getName());
		this.setAlternativeCode(r.getAlternativeCode());
		this.setInsertionCode(r.getInsertionCode());
		this.setHETATM(r.isHETATM());
		fitAtoms(r);
		atoms_backbone.get(BackBoneSample.PREVC).loc.set(r.prevc.loc);
		atoms_backbone.get(BackBoneSample.NEXTN).loc.set(r.nextn.loc);
		this.sideChain = r.sideChain;
	
	}
	FloatingResidue(PDBResidue r,PDBResidue prev,PDBResidue next){
		this(r.getName());
		this.parent = r.parent;
		this.setResidueNumber(r.getResidueNumber());
		this.setName(r.getName());
		this.setAlternativeCode(r.getAlternativeCode());
		this.setInsertionCode(r.getInsertionCode());
		this.setHETATM(r.isHETATM());
		fitAtoms(r);
		if(prev != null){

			atoms_backbone.get(BackBoneSample.PREVC).loc.set(prev.getC().loc);
		}else{
			//調整の必要あり
			atoms_backbone.get(BackBoneSample.PREVC).loc.set(this.getN().loc);
		}
		
		if(next != null){
			atoms_backbone.get(BackBoneSample.NEXTN).loc.set(next.getN().loc);
		}else{
			//調整の必要あり
			atoms_backbone.get(BackBoneSample.PREVC).loc.set(this.getC().loc);
		}
		
	}
	
	FloatingResidue(String s){
		this.setName(s);
		for(String a:MyWorld.all_atoms.get(s)){
			PDBAtom atom = new PDBAtom();
			//atom.loc.x = Math.random()*3;
			//atom.loc.y = Math.random()*3;
			//atom.loc.z = Math.random()*3;
			atom.pdb_atom_code = a;
			atom.atom_code = a.substring(0,1);
			this.addAtom(atom);
			if(MyWorld.backbone_atoms.contains(a)){
				//sidechain に含まれる backbone の ATOM については本体を別で管理。sidechain のものと重複するため。
				//Sidechain のアラインにのみ使われる。
				atoms_sidechain.put(a, atom.getCopy());
				atoms_loc.add(atoms_sidechain.get(a).loc);
				atoms_backbone_loc.add(atom.loc);
			}else{
				atoms_sidechain.put(a, atom);
				atoms_loc.add(atom.loc);
			}
		}
		for(int ii = 0;ii < BackBoneSample.aindex.size();ii++){
			atoms_backbone.add(null);
		}
		atoms_backbone.set(BackBoneSample.C,this.getC());
		atoms_backbone.set(BackBoneSample.N,this.getN());
		atoms_backbone.set(BackBoneSample.CA,this.getCA());
		atoms_backbone.set(BackBoneSample.O,this.getO());
		prevc = new PDBAtom();
		nextn = new PDBAtom();
		atoms_backbone.set(BackBoneSample.NEXTN,nextn);
		atoms_backbone.set(BackBoneSample.PREVC,prevc);
		atoms_backbone_loc.add(nextn.loc);
		atoms_backbone_loc.add(prevc.loc);
		
		HashSet<Point3D> att = new HashSet<>(atoms_backbone_loc);
		
		att.addAll(atoms_loc);
		
		atoms_loc_all.addAll(att);
		for(Point3D p:atoms_loc_all){
			atoms_loc_all_prev.add(new Point3D(p));
		}
	}
	
	public void fitAtoms(PDBResidue r){
		
		for(PDBAtom aa:r.atoms){
			PDBAtom bb = this.getAtomByName(aa.pdb_atom_code);
			if(bb == null){
				System.err.println(aa.pdb_atom_code+" is not supported atom.");
			}else{
				bb.loc.set(aa.loc);
			}
		}
		for(String ss:atoms_sidechain.keySet()){
			PDBAtom bb = r.getAtomByName(ss);
			if(bb != null){
				atoms_sidechain.get(ss).loc.set(bb.loc);
			}
		}
		
	}
	public void placeBackBone(PDBResidue template,PDBResidue p,PDBResidue n){
		atoms_backbone.get(BackBoneSample.C).loc.set(template.getC().loc);
		atoms_backbone.get(BackBoneSample.N).loc.set(template.getN().loc);
		atoms_backbone.get(BackBoneSample.CA).loc.set(template.getCA().loc);
		atoms_backbone.get(BackBoneSample.O).loc.set(template.getO().loc);
		if(n != null){
			atoms_backbone.get(BackBoneSample.NEXTN).loc.set(n.getN().loc);
		}
		if(p != null){
			atoms_backbone.get(BackBoneSample.PREVC).loc.set(p.getC().loc);
		}
	}
	
	public void saveLoc(){
		for(int ii = 0;ii < atoms_loc_all.size();ii++){
			atoms_loc_all_prev.get(ii).set(atoms_loc_all.get(ii));
		}
	}
	
	public void restoreLoc(){
		for(int ii = 0;ii < atoms_loc_all.size();ii++){
			atoms_loc_all.get(ii).set(atoms_loc_all_prev.get(ii));
		}
	}
	
	
	public double getDistance(FloatingResidue r){
		double ret = atoms.get(0).distance(r.atoms.get(0));
		for(PDBAtom p:atoms){
			for(PDBAtom pp:r.atoms){
				ret =Math.min(p.distance(pp),ret);
			}
		}
		return ret;
	}
	
	
	public void move(Point3D direc){
		for(Point3D p:atoms_loc_all){
			p.x += direc.x;
			p.y += direc.y;
			p.z += direc.z;
		}
	}
	
	public void setNextFloating(FloatingResidue f){next = f;	}
	
	public void setPrevFloating(FloatingResidue f){prev = f;	}
	
	public void processBackBone(BackBoneSample bs){
		processBackBone(bs,(this.prev != null)?(prevOmega):(nextOmega),this.prev != null,true);
	}
	public void processBackBone(BackBoneSample bs,double radian,boolean prevflag,boolean fitsidechain){
		if(bs != null){
			for(int ii = 0;ii < bs.atoms.length;ii++){
				atoms_backbone.get(ii).loc.set(bs.atoms[ii].loc);
			}
		}
		if(prev != null || next != null){
			setOmega(radian,prevflag);
		}
		if(fitsidechain){
			fitSideChain();
		}
	}
	public void calcPseudoNeighbour(BackBoneSample bs){
		
		PDBAtom n1 = bs.atoms[BackBoneSample.N];
		PDBAtom c1 = bs.atoms[BackBoneSample.C];
		PDBAtom ca1 = bs.atoms[BackBoneSample.CA];
		
		PDBAtom n2 = atoms_backbone.get(BackBoneSample.N);
		PDBAtom c2 = atoms_backbone.get(BackBoneSample.C);
		PDBAtom ca2 = atoms_backbone.get(BackBoneSample.CA);
		ArrayList<Point3D> att = new ArrayList<>();
		att.add(atoms_backbone.get(BackBoneSample.PREVC).loc);
		att.add(atoms_backbone.get(BackBoneSample.NEXTN).loc);
		
		atoms_backbone.get(BackBoneSample.PREVC).loc.set(bs.atoms[BackBoneSample.PREVC].loc);
		atoms_backbone.get(BackBoneSample.NEXTN).loc.set(bs.atoms[BackBoneSample.NEXTN].loc);
		PepProcess.adjustVector3D(
		ca1.loc,n1.loc,c1.loc,
		ca2.loc,n2.loc,c2.loc,
		att);
	}

	public void fitBackbone(boolean forward){
		if(forward){
			setOmega(prevOmega,forward);
		}else{
			setOmega(nextOmega,forward);
		}
		
		fitSideChain();
	}
	public void setOmega(double radian,boolean prevflag){
		
		if(prevflag && this.prev != null){//Prev に何か残基がある時
			prevOmega = radian;
			PepProcess.adjustVector3D(
				this.getN().loc,
				this.prevc.loc,
				this.getCA().loc,
				prev.nextn.loc,
				prev.getC().loc,
				prev.getCA().loc,
				atoms_backbone_loc
				);
			Point3D vec = PepProcess.subtracted(prev.getC().loc,prev.nextn.loc);
			for(Point3D p:atoms_backbone_loc){
				p.x -= prev.nextn.loc.x;
				p.y -= prev.nextn.loc.y;
				p.z -= prev.nextn.loc.z;
				Point3D ccc = Point3D.rotate(p,vec,radian);
				p.set(ccc);
				p.x += prev.nextn.loc.x;
				p.y += prev.nextn.loc.y;
				p.z += prev.nextn.loc.z;
			}
			
		}else{
			if(next == null){
				System.out.println("!!!");
			}
			nextOmega = radian;
			
			PepProcess.adjustVector3D(
				this.nextn.loc,
				this.getC().loc,
				this.getCA().loc,
				
				next.getN().loc,
				next.prevc.loc,
				next.getCA().loc,
				atoms_backbone_loc
				);
			Point3D vec = PepProcess.subtracted(next.prevc.loc,next.getN().loc);
			//vec.standarize();
			for(Point3D p:atoms_backbone_loc){
				p.x -= next.getN().loc.x;
				p.y -= next.getN().loc.y;
				p.z -= next.getN().loc.z;
				Point3D ccc = Point3D.rotate(p,vec,radian);
				p.set(ccc);
				p.x += next.getN().loc.x;
				p.y += next.getN().loc.y;
				p.z += next.getN().loc.z;
			}
		}
	}
	
	public void processBackBone_ignoreOmega(BackBoneSample bs){
		if(bs != null){
			for(int ii = 0;ii < bs.atoms.length;ii++){
				atoms_backbone.get(ii).loc.set(bs.atoms[ii].loc);
			}
		}
		if(prev != null && next != null){
			/*
			PepProcess.adjustVector3D(
				this.prevc.loc,
				this.getN().loc,
				this.nextn.loc,

				prev.getC().loc,
				prev.nextn.loc,
				next.getN().loc,
				atoms_backbone_loc
				);
			*/
			
			PepProcess.adjustVector3D(
				this.prevc.loc,
				this.getN().loc,
				this.nextn.loc,

				prev.getC().loc,
				prev.nextn.loc,
				next.getN().loc,
				atoms_backbone_loc
				);
		}else{
			//回転軸をNCとして回転が必要。
			if(prev == null){
				PepProcess.adjustVector3D(
				this.getC().loc,
				this.nextn.loc,
				next.prevc.loc,
				next.getN().loc,
				atoms_backbone_loc
				);
				//回転軸をNCとして回転が必要。
			}
			if(next == null){
				PepProcess.adjustVector3D(
				this.getN().loc,
				this.prevc.loc,
				prev.nextn.loc,
				prev.getC().loc,
				atoms_backbone_loc
				);
			}
		}
	}
	
	public void changeSideChain(SideChainSample s){
		sideChain = s;
		for(String ss:s.atoms_all.keySet()){
			atoms_sidechain.get(ss).loc.set(s.atoms_all.get(ss).loc);
		}
		
		PDBAtom nn = atoms_sidechain.get("N");//ここに入っているのはダミー
		PDBAtom cc = atoms_sidechain.get("C");
		PDBAtom ca = atoms_sidechain.get("CA");
		
		PDBAtom n2 = atoms_backbone.get(BackBoneSample.N);
		PDBAtom c2 = atoms_backbone.get(BackBoneSample.C);
		PDBAtom ca2 = atoms_backbone.get(BackBoneSample.CA);
		
		PepProcess.adjustVector3D(
		ca.loc,nn.loc,cc.loc,
		ca2.loc,n2.loc,c2.loc,
		this.atoms_loc);
	}
	
	/**
	 * 現在の Side chain を、Backbone に合わせる。
	 */
	public void fitSideChain(){
		if(sideChain == null){
			System.err.println("Side chain data is null.");
		}else{
			changeSideChain(sideChain);
		}
		
	}
	
}



class AtomTriangle{
	PDBAtom[] atoms = new PDBAtom[3];
	Point3D[] orig = new Point3D[3];
	double[] dist = {0,0,0};
	Point3D[] saved = new Point3D[3];
	double savedCrash = -10000;
	HashMap<PDBAtom,Integer> atomFreezed = new HashMap<>();
	AtomTriangle(PDBAtom a1,PDBAtom a2,PDBAtom a3,double a1_a2,double a2_a3,double a3_a1){
		atoms[0] = a1;
		atoms[1] = a2;
		atoms[2] = a3;
		dist[0] = a1_a2;
		dist[1] = a2_a3;
		dist[2] = a3_a1;
		for(int ii = 0;ii < 3;ii ++ ){
			orig[ii] = new Point3D(atoms[ii].loc);
			saved[ii] = new  Point3D(atoms[ii].loc);
		}
	}
	
	public void saveLoc(){
		for(int ii = 0;ii < 3;ii ++ ){
			saved[ii].set(atoms[ii].loc);
		}
	}
	
	
	public void restoreLoc(){
		for(int ii = 0;ii < 3;ii ++ ){
			atoms[ii].loc.set(saved[ii]);
		}
	}
	
	
	public void freeze(int i){
		atomFreezed.put(atoms[i],i);
	}
	
	public void fit_dist(PDBAtom a,PDBAtom b,double tdist,double ratio){
		double ddist = a.distance(b);
		double cx = a.loc.x/2+b.loc.x/2;
		double cy = a.loc.y/2+b.loc.y/2;
		double cz = a.loc.z/2+b.loc.z/2;
		
		
		Point3D vec = new Point3D(a.loc.x-b.loc.x,a.loc.y-b.loc.y,a.loc.z-b.loc.z);
		vec.standarize();
		double kdist = tdist*ratio+(1.0-ratio)*ddist;
		vec.x *= kdist;
		vec.y *= kdist;
		vec.z *= kdist;
		double vcx = cx-vec.x/2;
		double vcy = cy-vec.y/2;
		double vcz = cz-vec.z/2;
		vec.x += vcx;
		vec.y += vcy;
		vec.z += vcz;
		//System.out.println(a.loc.distance(vec));
		//System.out.println(b.loc.distance(new Point3D(vcx,vcy,vcz)));
		a.loc.set(vec);
		b.loc.set(vcx,vcy,vcz);
		//System.out.println(a.distance(b));
		
		
	}
	
	
	public double diffSum(){
		double ret = 0;
		for(int ii = 0;ii < 3;ii++){
			ret += Math.abs(dist[ii]
			-atoms[ii].distance(atoms[(ii+1)%3]));
		}
		return ret;
	}
	public void fitThree(double ratio){
		fit_dist(atoms[0],atoms[1],dist[0],ratio);
		fit_dist(atoms[1],atoms[2],dist[1],ratio);
		fit_dist(atoms[2],atoms[0],dist[2],ratio);
		
		if(atomFreezed.size() > 0){
			
			double px = 0;
			double py = 0;
			double pz = 0;
			for(PDBAtom z:atoms){
				if(atomFreezed.containsKey(z)){
					//固定されるAtomがあった場合処理
					Point3D pos = orig[atomFreezed.get(z)];
					px = pos.x - z.loc.x;
					py = pos.y - z.loc.y;
					pz = pos.z - z.loc.z;
					break;
				}
			}
			
			for(PDBAtom z:atoms){
				z.loc.x += px;
				z.loc.y += py;
				z.loc.z += pz;
			}
		}
	}
	
	
}
	

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static pepbuilderj.AtomDistancePenalty.TYPE_BACKBONE;
import static pepbuilderj.AtomDistancePenalty.TYPE_SIDECHAIN;

/**
 *
 * @author kimidori
 */
public class ResidueRefine {
	AtomDistancePenalty adp = new AtomDistancePenalty();
	FuzzyDecisionTreeScoring_generator scoring 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorAtom_20180503());
	ArrayList<PDBResidue> envresidues = new ArrayList<>();
	AtomDistancePenalty distPenal = new AtomDistancePenalty();
	Random randomizer = null;
	ResidueRefine(){
		randomizer = new Random(System.currentTimeMillis());
	}
	public void prepare(ArrayList<PDBResidue> r){
		scoring.prepare(r);
		envresidues.clear();
		envresidues.addAll(r);
	}
	public double maxRotamer(FloatingResidue r
			,double crashpenalty
			,double crash_tolerance_ratio
			,boolean sidechainonly
			,SideChainSet samples
	,boolean otherchainonly){
		
		ArrayList<PDBAtom> targetatoms = new ArrayList<>();
		if(sidechainonly){
			for(PDBAtom p:r.atoms){
				if(p.pdb_atom_code.equals("CA")
						||p.pdb_atom_code.equals("C")
						||p.pdb_atom_code.equals("O")
						||p.pdb_atom_code.equals("N")
						||p.pdb_atom_code.equals("OXT")
						){
				}else{
					targetatoms.add(p);
				}
			}
		}else{
			targetatoms.addAll(r.atoms);
		}
		ArrayList<SideChainSample> ss =  samples.sidechains.get(r.getName());
		double maxscore = -100000;
		int maxindex = 0;
		
		for(int ii = 0;ii <  ss.size();ii++){
			SideChainSample scc = ss.get(ii);
			if(maxscore > scc.prob*200){//200点以上になるのはないと思う。適当。
				break;
			}
			r.changeSideChain(scc);
			int crash = getCrashNum(r,targetatoms,envresidues,crash_tolerance_ratio,otherchainonly);
			
			
			
			if(crashpenalty > 0){
				crashpenalty *= -1;
			}
			double score = crash*crashpenalty;
			ArrayList<PDBAtom> al = scoring.featureGen.getTargetAtoms(r);
			for(PDBAtom aa:al){
				String code = r.getName()+"_"+aa.pdb_atom_code;
				double df[] = scoring.featureGen.generateFeaturesA(aa);
				if(df == null){
					continue;
				}
				double sc =  scoring.getLabelScore(df, code);
				if(sc == scoring.SCORE_INVALIDRESIDUE){
					System.err.print("cannot get score for "+code);
				}else{			
					score += sc;
				}
			}
			if(score > 0){//スコアについては適当
				score *= scc.prob;
			}
			if(score > maxscore){
				maxscore = score;
				maxindex = ii;
			}
		}
		r.changeSideChain(ss.get(maxindex));
		return maxscore;
	}
	
	
	/**
	 * 他原子と Crash している原子数を返す。
	 * targets_ が null の場合全部の原子をチェックする。
	 */
	public int getCrashNum(PDBResidue target,ArrayList<PDBAtom> targets_
			,ArrayList<PDBResidue> allres //atom にしていいとも思うが、あまりコストはないと思うので
			,double torelance_ratio
	,boolean otherchainonly){
		int ret = 0;
		ArrayList<PDBAtom> targets = new ArrayList<>();
		if(targets_ != null){
			targets.addAll(targets_);
		}else{
			targets.addAll(target.atoms);
		}
		for(int ii = 0;ii < allres.size();ii++){
			PDBResidue aa = allres.get(ii);
			if(otherchainonly && aa.parent == target.parent){
				continue;
			}
			if(aa == target){
				continue;
			}
			if(aa.getCA().distance(target.getCA()) > 20){
				continue;
			}
			for(PDBAtom a:aa.atoms){
				
				String acode = a.parent.getName()+":"+a.pdb_atom_code;
				for(PDBAtom r:targets){
					String rcode = r.parent.getName()+":"+r.pdb_atom_code;
					double dist = a.distance(r);
					if(!adp.pairsMinDist.containsKey(rcode+"_"+acode)){
						if((torelance_ratio+1.0)*dist < 1.8){
							ret += 1;
							break;
						}
					}else{
						if(adp.pairsMinDist.get(rcode+"_"+acode) > (torelance_ratio+1.0)*dist){
							ret += 1;
							break;
						}
					}

				}
			}
		}
		return ret;
	}
	
	public ArrayList<Double> refineSideChains(ArrayList<FloatingResidue> list,SideChainSet sidechains){
		return this.refineSideChains(list,sidechains,false);
	}
	
	/**
	 * backbone は変更せずスコアの高い Rotamer を探す。
	 * List 順に調べるので List への登録順を変えて調整。
	 * @param list
	 * @param sidechains
	 * @param otherchainonly 
	 */
	public ArrayList<Double> refineSideChains(ArrayList<FloatingResidue> list
			,SideChainSet sidechains,boolean otherchainonly){
		ArrayList<Double> ret = new ArrayList<>();
		for(int ii = 0;ii < list.size();ii++){
			ret.add(this.maxRotamer(list.get(ii),
					-5,
					0.2,
					true,
					sidechains
			,otherchainonly));
		}
		return ret;
	}
	
	/**
	 * multimer について考える場合は中の Crash については無視。
	 * これがいいのか悪いのかは不明。
	 * Crash してると Score が低くなると思う。
	 * @param list
	 * @param sidechains 
	 */
	public ArrayList<Double> RefineSideChains_multimer(ArrayList<FloatingResidue> list,SideChainSet sidechains){
		return this.refineSideChains(list,sidechains,true);
	}
	
	/**
	 * point2 を頂点とする三角形が dist12 ~ dist13 を満たすようにする
	 * @param point1
	 * @param point2
	 * @param point3
	 * @param dist12
	 * @param dist23
	 * @param dist13 
	 */
	
	public static void refineTriangle(Point3D point1,Point3D point2, Point3D point3,double dist12,double dist23,double dist13
	
	){
		
		Point3D veccca = PepProcess.subtracted(point1,point3);
		Point3D c = new Point3D();
		c.x = point3.x+veccca.x/2;
		c.y = point3.y+veccca.y/2;
		c.z = point3.z+veccca.z/2;

		//System.out.println(c.distance(new Point3D(point1.x/2+point3.x/2,point1.y/2+point3.y/2,point1.z/2+point3.z/2)));
		
		double d_13 = point1.distance(point3);

		if(d_13 > 0){
			point1.set(c.x+veccca.x/d_13*dist13/2,
					c.y+veccca.y/d_13*dist13/2,
					c.z+veccca.z/d_13*dist13/2
			);
			point3.set(c.x-veccca.x/d_13*dist13/2,
					c.y-veccca.y/d_13*dist13/2,
					c.z-veccca.z/d_13*dist13/2
			);
			
			d_13 = point1.distance(point3);//must be same with dist13
			//System.out.println(d_13-dist13);
		
			double d_12 = point1.distance(point2);
			double d_23 = point2.distance(point3);
			double p = (d_23*d_23+d_13*d_13-d_12*d_12)/2/d_13;
			double pppc = (dist23*dist23+dist13*dist13-dist12*dist12)/2/dist13;
			//実際来るべき垂線の足の位置を計算する。
			if(p == 0){
				return;
			}
			
			Point3D pv = PepProcess.subtracted(point1,point3);
			double ratio = p/d_13;
			pv.x *= ratio;
			pv.y *= ratio;
			pv.z *= ratio;
			
			pv.x += point3.x;
			pv.y += point3.y;
			pv.z += point3.z;
			//point2 から降ろした垂線の足
			
			Point3D pv2 = PepProcess.subtracted(point2,pv);
			
			double dd = pv2.x*pv2.x+pv2.y*pv2.y+pv2.z*pv2.z;
			
			
			double pratio = pppc/dist13;
			double gdist = (dist23)*(dist23)-
					(dist13*pratio)*(dist13*pratio);
			
			if(dd > 0 && gdist > 0){
				c.x = point3.x+veccca.x*pratio;
				c.y = point3.y+veccca.y*pratio;
				c.z = point3.z+veccca.z*pratio;
				
				
				gdist = Math.sqrt(gdist);
				dd = Math.sqrt(dd);
				pv2.x *= gdist/dd;
				pv2.y *= gdist/dd;
				pv2.z *= gdist/dd;

				pv2.x += c.x;
				pv2.y += c.y;
				pv2.z += c.z;
				point2.set(pv2);
			}

		}
		
		
	}
	
	
	
	//public static final int BACKBONE_N_CA = 0;
	//public static final int BACKBONE_CA_C = 1;
	//public static final int BACKBONE_C_O = 2;
	//public static final int BACKBONE_N_C = 3;
	//public static final int BACKBONE_PREVC_CA = 4;
	//public static final int BACKBONE_PREVCA_N = 5;
	/**
	 * backbone が結合角を満たすように変更する。
	 * @param r
	 * @param freezeprev
	 * @param freezenext 
	 */
	public void waveBackBone(FloatingResidue r,boolean freezeprev,boolean freezenext){
		
		
		if(r.prev != null){
			double pxx = r.prev.getC().loc.x;
			double pyy = r.prev.getC().loc.y;
			double pzz = r.prev.getC().loc.z;
			refineTriangle(r.getCA().loc,r.getN().loc,r.prev.getC().loc
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_N_CA]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_CA]
			);
			if(freezeprev){
				double dxx = pxx-r.prev.getC().loc.x;
				double dyy = pyy-r.prev.getC().loc.y;
				double dzz = pzz-r.prev.getC().loc.z;
				r.getCA().loc.x += dxx;
				r.getCA().loc.y += dyy;
				r.getCA().loc.z += dzz;

				r.getN().loc.x += dxx;
				r.getN().loc.y += dyy;
				r.getN().loc.z += dzz;


				r.prev.getC().loc.x += dxx;
				r.prev.getC().loc.y += dyy;
				r.prev.getC().loc.z += dzz;
			}
			
		}
			
		if(r != null){//全部通す
			refineTriangle(r.getN().loc,r.getCA().loc,r.getC().loc
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_N_CA]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_CA_C]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_N_C]
			);
		}
		
		
		if(r.next != null){
			
			double pxx = r.next.getN().loc.x;
			double pyy = r.next.getN().loc.y;
			double pzz = r.next.getN().loc.z;
			
			
			refineTriangle(r.getCA().loc,r.getC().loc,r.next.getN().loc
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_CA_C]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N]
			,distPenal.backbone_average[AtomDistancePenalty.BACKBONE_PREVCA_N]
			);
			
			//O の位置を決定
			Point3D point1 = r.getCA().loc;
			Point3D point3 = r.next.getN().loc;
			Point3D veccca = PepProcess.subtracted(point1,point3);
			Point3D c = new Point3D();//面倒くさいので中心点にしている。
			c.x = point3.x+veccca.x/2;
			c.y = point3.y+veccca.y/2;
			c.z = point3.z+veccca.z/2;
			double ddist = r.getC().loc.distance(c);
			double ccdist = distPenal.backbone_average[AtomDistancePenalty.BACKBONE_C_O];
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
			
			if(freezenext){
				double dxx = pxx-r.next.getN().loc.x;
				double dyy = pyy-r.next.getN().loc.y;
				double dzz = pzz-r.next.getN().loc.z;
				r.getCA().loc.x += dxx;
				r.getCA().loc.y += dyy;
				r.getCA().loc.z += dzz;

				r.getC().loc.x += dxx;
				r.getC().loc.y += dyy;
				r.getC().loc.z += dzz;


				r.next.getN().loc.x += dxx;
				r.next.getN().loc.y += dyy;
				r.next.getN().loc.z += dzz;
			}
		}else{
			//最後の残基の場合CAとれないので
			//ついでに一個しかつけない
			//Fixit
			Point3D point1 = r.getC().loc;
			Point3D point2 = r.getCA().loc;
			Point3D veccca = PepProcess.subtracted(point1,point2);
			double ddist = veccca.x*veccca.x+ veccca.y*veccca.y+ veccca.z*veccca.z;
			double ccdist = distPenal.backbone_average[AtomDistancePenalty.BACKBONE_C_O];
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
		if(r.next != null){
			r.nextn.loc.set(r.next.getN().loc);
		}
		if(r.prev != null){
			r.prevc.loc.set(r.prev.getC().loc);
		}
		
		
	}
	
	
	public void randomBackboneChange(){
		
	}
	
	public void refine(ArrayList<FloatingResidue> target_
			,ArrayList<PDBResidue> allresidues//スコアリングのため
			,int iternum
			,boolean onlysidechain
			,SideChainSet scc,BackBoneSet bs){
		
		ArrayList<FloatingResidue> target = new ArrayList<>(target_);
		for(FloatingResidue tt:target){
			tt.saveLoc();
		}
		ArrayList<FuzzyDecisionTreeScoring_generator> slis = new ArrayList<>();
		slis.add(scoring);
		double[] ss = FuzzyDecisionTreeScoring_generator.calcResidueScores(allresidues, slis);
		double maxscores = 0;
		for(double s:ss){
			maxscores += s;
		}
		for(int ii = 0;ii < iternum;ii++){
			if(!onlysidechain){
				for(int jj = 0;jj < target.size();jj++){
					FloatingResidue r = target.get(jj);
					if(r.getName().equals("PRO")){
						if(randomizer.nextDouble() < 0.1){
							//適当
							r.processBackBone(bs.getNext(r),0.0,false,false);
						}else{
							r.processBackBone(bs.getNext(r),Math.PI,false,false);
						}
					}else{
						r.processBackBone(bs.getNext(r),Math.PI,r.prev != null,false);
					}
					waveBackBone(r,false,false);
				}
			}
			
			ArrayList<Double> scores = refineSideChains(target,scc);
			double[] dd = FuzzyDecisionTreeScoring_generator.calcResidueScores(allresidues, slis);
			double sd = 0;
			for(double d:dd){
				sd += d;
			}
			if(sd < maxscores){
				for(FloatingResidue tt:target){
					tt.restoreLoc();
				}
				for(int jj = 0;jj < target.size()*5;jj++){
					int t1 = (int)(randomizer.nextDouble()*target.size());
					int t2 = (int)(randomizer.nextDouble()*target.size());
					FloatingResidue f1 = target.get(t1);
					FloatingResidue f2 = target.get(t2);
					target.set(t1,f2);
					target.set(t2,f1);
				}
			}else{
				for(FloatingResidue tt:target){
					tt.saveLoc();
				}
				maxscores = sd;
			}
		}
	}
	
	
	public void loopRefile(ArrayList<FloatingResidue> allresidues
			,ArrayList<FloatingResidue> target,ChainBuilder cb
	,int iternum,boolean rotamercheck){//fixme Chainbuilder に移す。
		ArrayList<PDBResidue> ress = new ArrayList<>();
		ArrayList<PDBResidue> resse = new ArrayList<>();
		resse.addAll(allresidues);
		prepare(resse);
		
		
		ress.addAll(target);
		ArrayList<FuzzyDecisionTreeScoring_generator> slis = new ArrayList<>();
		slis.add(scoring);
		double prevscores = -1000000;
		int convnum = 0;
		for(int ii = 0;ii < iternum;ii++){
			double[] ss = FuzzyDecisionTreeScoring_generator.calcResidueScores(ress, slis);
			double maxscores = 0;
			for(double s:ss){
				maxscores += s;
			}
			double sumpenal = 0;
			for(PDBResidue p:ress){
				double penal = distPenal.calcCrashPenalty(envresidues, p,AtomDistancePenalty.TYPE_ALL);
				sumpenal+=penal;
			}
			System.out.println(maxscores+";;;"+sumpenal);
			maxscores += sumpenal;
			double randomfactor = ((ii%2 == 0)?(15):(3))*((iternum-ii)/iternum+1)+1;//Math.min(5, 1+Math.abs(sumpenal/300));
			
			if(prevscores < maxscores){
				for(FloatingResidue r:target){
					r.saveLoc();
				}
				prevscores = maxscores;
				convnum = 0;
			}else{
				for(FloatingResidue r:target){
					r.restoreLoc();
				}

				convnum ++ ;
				if(convnum > 10){
					for(FloatingResidue r:target){
						r.restoreLoc();
					}

					break;
				}
				for(FloatingResidue r:target){
					double mx = Math.random()*randomfactor-randomfactor/2;
					double my = Math.random()*randomfactor-randomfactor/2;
					double mz = Math.random()*randomfactor-randomfactor/2;
					r.move(new Point3D(mx,my,mz));
				}
				
				
				cb.fillGapsTarget(target,1000);
				if(rotamercheck){
					for(FloatingResidue r:target){
						maxRotamer(r,-5,0.1, true, cb.sidechains, false);
					}
				}
			}
		}
	}
	
	public static HashMap<String,String> parseline(String s){
		HashMap<String,String> ret = new HashMap<>();
		String[] pt = s.replaceAll("[\\r\\n]","").split("[\\s]");
		Pattern cpat = Pattern.compile("([^\\:]+)\\:([^\\:]+)");
		for(String pp:pt){
			
			Matcher mat = cpat.matcher(pp);
			if(mat.find()){
				ret.put(mat.group(1),mat.group(2));
			}
		}
		return ret;
	}
	public static ArrayList<ResidueIdentifier> parseTargetFile(String filename){
		ArrayList<ResidueIdentifier> ret = new ArrayList<>();
		
		BufferedReader br = null;
		try{
			br = new  BufferedReader(new FileReader(new File(filename)));
			String line = null;
			while((line = br.readLine()) != null){
				HashMap<String,String> ss = parseline(line);
				if(!ss.containsKey("chain_name")){
					System.err.println("cannnot parse "+line);
					continue;
				}
				ResidueIdentifier ri = new ResidueIdentifier(ss.get("chain_name"),
				ss.get("residue_name"),
				Integer.parseInt(ss.get("residue_number")));
				ret.add(ri);
			}
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	
	public static HashMap<String,String> parseArgs(String args[]){
		HashMap<String,String> ret = new HashMap<>();
		
		for(int ii = 0;ii < args.length;ii++){
			if(args[ii].indexOf("-") == 0){
				if(args.length <= ii+1){
					ret.put(args[ii],"TRUE");
				}else{
					if(args[ii+1].indexOf("-") == 0){
						ret.put(args[ii],"TRUE");
					}else{
						ret.put(args[ii],args[ii+1]);
					}
				}
			}
		}
		return ret;
	}
	
	
	
	public static void main(String[] args){
		HashMap<String,String> hs = parseArgs(args);
		String infile = hs.get("-in");
		String outfile = hs.get("-out");
		String targetFile =  "";
		if(hs.containsKey("-target")){
			targetFile = hs.get("-target"); 
		}
		int iternum = 5;
		
		if(hs.containsKey("-iter")){
			iternum = Integer.parseInt(hs.get("-iter")); 
		}
		boolean onlysidechain = false;
		if(hs.containsKey("-sidechain")){
			onlysidechain = true; 
		}
		
		ArrayList<ResidueIdentifier> targetlist = null;
		if(targetFile.length() != 0){
			targetlist = parseTargetFile(targetFile);
		}
		PDBData pdb = PDBData.loadPDBFile(infile);
		ArrayList<FloatingResidue> fl = new ArrayList<>();
		ArrayList<FloatingResidue> targets = new ArrayList<>();
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
				if(targetlist == null){
					targets.add(fr.get(ii));
				}else{
					Iterator<ResidueIdentifier> ite = targetlist.iterator();
					while(ite.hasNext()){
						ResidueIdentifier ri = ite.next();
						if(ri.hit(fr.get(ii))){
							targets.add(fr.get(ii));
							ite.remove();
						}
					}
				}
			}
		}
		ArrayList<PDBResidue> pp = new ArrayList<>();
		for(FloatingResidue f:fl){
			pp.add(f);
		}
		ResidueRefine rr = new ResidueRefine();
		rr.prepare(pp);
		
		rr.refine(targets, pp,iternum,onlysidechain,
				new SideChainSet(), new BackBoneSet());
		
		
		ChainBuilder.saveChains(fl,outfile);
	}
	
}



class ResidueIdentifier{
	String chainName = "";
	String residueName = "";
	int residueNumber = -1;
	ResidueIdentifier(String c,String n,int i){
		chainName = c;
		residueName = n;
		residueNumber = i;
	}
	public boolean hit(PDBResidue r){
		if(r.getResidueNumber() == residueNumber){
			if(r.getName().equals(residueName)){
				if(r.parent.getName().equals(chainName)){
					return true;
				}
			}
		}
		return false;
	}
	
}
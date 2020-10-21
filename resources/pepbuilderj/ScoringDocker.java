/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static pepbuilderj.TemplateBaseModeller.getRotationAngle;
import static pepbuilderj.TemplateBaseModeller.makeAtomLines;
import static pepbuilderj.TemplateBaseModeller.rotateWith;
import static pepbuilderj.TemplateBaseModeller.rotateXZero;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class ScoringDocker {
	ArrayList<PDBAtom> templateAtoms = new ArrayList<>();
	ArrayList<PseudoAtomPoint> pseudoAtoms = new ArrayList<>();
	static double atomRadius = 1.8;
	static double probeRadius = 1.5;
	static double pseudoatomDistance = 4.0;
	
	//個別計算のため分けているが、同じクラスからのインスタンスである必要がある。
	FuzzyDecisionTreeScoring_generator scoring_template 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorAtom_20180503());
	FuzzyDecisionTreeScoring_generator scoring_target 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorAtom_20180503());
	
	
	ScoringDocker(){
		
	}
	
	
	public static ArrayList<Point3D> mapBase = new ArrayList<>();
	
	static{
		double[] xx = {-1,0.0,1.0};
		double[] yy = {-1,0.0,1.0};
		double[] zz = {-1,0.0,1.0};
		for(int x = 0;x < xx.length;x++){
			for(int y = 0;y < yy.length;y++){
				for(int z = 0;z < zz.length;z++){
					if(xx[x] == 0.0 && yy[y] == 0.0 && zz[z] == 0.0){
						
					}else{
						mapBase.add(new Point3D(xx[x],yy[y],zz[z]));
					}
				}
			}
		}
		for(Point3D p3:mapBase){
			double dlen = Math.sqrt(p3.x*p3.x+p3.y*p3.y+p3.z*p3.z);
			p3.x /= dlen/pseudoatomDistance;
			p3.y /= dlen/pseudoatomDistance;
			p3.z /= dlen/pseudoatomDistance;
		}
	}
	/**
	 * base を平行移動させて原子を囲む点を作成する
	 * @param p
	 * @return 
	 */
	public static ArrayList<Point3D> mapSurroundingPoints(Point3D p){
		ArrayList<Point3D> ret = new ArrayList<>();
		for(Point3D r:mapBase){
			Point3D pp = new Point3D(r);
			pp.x += p.x;
			pp.y += p.y;
			pp.z += p.z;
			ret.add(pp);
		}
		return ret;
	}
	/**
	 * base を基に原子を囲む点を作成し、その点が他原子の基準 Radius 内にあると、
	 * 露出していないとみなし削除される。残った点上に仮原子を作成し返す。
	 * @param atoms
	 * @return 
	 */
	public ArrayList<PseudoAtomPoint> getSurroundingAtoms(ArrayList<PDBAtom> atoms){
		ArrayList<PseudoAtomPoint> ret = new ArrayList<>();
		for(PDBAtom a:atoms){
			//スコアで低いAtomだけにしてもいいかも？
			//周辺残基をとれる保証はない？？
			ArrayList<Point3D> pseudo = new ArrayList<>(mapSurroundingPoints(a.loc));
			for(PDBAtom b:templateAtoms){
				double dd = a.distance(b);
				if(dd < pseudoatomDistance + atomRadius +probeRadius){
					Iterator<Point3D> ite = pseudo.iterator();
					while(ite.hasNext()){
						Point3D p3 = ite.next();
						if(p3.distance(b.loc) < atomRadius+probeRadius){
							ite.remove();
						}
					}
					if(pseudo.size()== 0){
						break;
					}
				}
			}
			if(pseudo.size() > 0){
				for(Point3D p3:pseudo){
					ret.add(new PseudoAtomPoint(p3,a));
				}
			}
		}
		return ret;
	}
	
	public void prepareTemplateProtein(ArrayList<PDBResidue> rr){
		scoring_template.prepare(rr);
		ArrayList<PDBAtom> aa = new ArrayList<>();
		for(PDBResidue r:rr){
			if(r.isLigand() || r.isMissing()){
				continue;
			}
			aa.addAll(r.atoms);
		}
		prepare(aa);
		scoring_template.prepare(rr);
		for(PseudoAtomPoint p:pseudoAtoms){
			p.setFeatures(scoring_template.featureGen.generateFeatures(p.pseudo));
		}
	}
	
	
	public void prepare(ArrayList<PDBAtom> atoms){
		templateAtoms.clear();
		pseudoAtoms.clear();
		for(PDBAtom a:atoms){
			
			templateAtoms.add(a);
		}
		PDBResidue dummy = new PDBResidue();
		dummy.setName("DUM");
		ArrayList<PseudoAtomPoint> pp = getSurroundingAtoms(atoms);
		for(PseudoAtomPoint pz:pp){
			pz.pseudo.parent = dummy;
			pseudoAtoms.add(pz);
		}
	}
	
	/**
	 * xa -> xb ベクトルを ta-> tb ベクトルに合わせる。
	 * 始点は x=tb になる。
	 * @param xa
	 * @param xb
	 * @param ta
	 * @param tb
	 * @param allatom 
	 */
	
	public static void dockWith(PDBAtom xa,PDBAtom xb
	,Point3D ta,Point3D tb,
	ArrayList<PDBAtom> allatom){
		
		Point3D p1 = new Point3D(xa.loc);
		Point3D p2 = new Point3D(xb.loc);
		Point3D p3 = ta;
		Point3D p4 = tb;
		
		double angle[] = getRotationAngle(p4.x-p3.x,p4.y-p3.y,p4.z-p3.z);
		double angle2[] = getRotationAngle(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
			
		for(PDBAtom a:allatom){
			a.loc.x -= p1.x;
			a.loc.y -= p1.y;
			a.loc.z -= p1.z;
			double d1[] = rotateXZero(a.loc.x,a.loc.y,a.loc.z,angle2);
			double dd[] = rotateWith(d1[0],d1[1],d1[2],angle);
			a.loc.set(dd);
			a.loc.x += p3.x;
			a.loc.y += p3.y;
			a.loc.z += p3.z;
		}
	}
	
	
	
	
	
	
	
	/**
	 * Residue に配置された PseudoAtom の数を返す。
	 * @param al
	 * @param threshold
	 * @return 
	 */
	public static HashMap<PDBResidue,Integer> getSurfaceResidues(ArrayList<PDBResidue> al,boolean nobackbone){
		ArrayList<PDBAtom> atoms = new ArrayList<>();
		HashMap<PDBResidue,Integer> surfaces = new HashMap<>();
		for(PDBResidue r:al){
			surfaces.put(r,0);
			boolean dflag = false;
			for(PDBAtom aa:r.atoms){
				if(aa.isAlternative()){
					if(!dflag){
						System.err.println(r.getRepresentativeCode()+" alternative atom is ignored.");
						dflag = true;
					}
					continue;
				}
				atoms.add(aa);
			}
		}
		for(PDBAtom a:atoms){
			if(nobackbone){
				if(a.pdb_atom_code.equals("C") ||
				a.pdb_atom_code.equals("CA") ||
				a.pdb_atom_code.equals("O") ||
				a.pdb_atom_code.equals("N")){
					continue;
				}
			}
			ArrayList<Point3D> pseudo = mapSurroundingPoints(a.loc);
			for(PDBAtom b:atoms){
				if(a == b){
					continue;
				}
				double dd = a.distance(b);
				if(dd < pseudoatomDistance + atomRadius +probeRadius){
					Iterator<Point3D> ite = pseudo.iterator();
					while(ite.hasNext()){
						Point3D p3 = ite.next();
						if(p3.distance(b.loc) < atomRadius+probeRadius){
							ite.remove();
						}
					}
					if(pseudo.size() == 0){
						break;
					}
				}
			}
			
			surfaces.put(a.parent,surfaces.get(a.parent)+pseudo.size());
		}
		
		
		return surfaces;
	}
	
	
	public ArrayList<DockingPoint> calcMergedScores(PDBAtom a,double basescore){
		ArrayList<Double> ad = this.scoring_target.featureGen.generateFeatures(a);
		String label = a.parent.getName()+"_"+a.pdb_atom_code;
		ArrayList<DockingPoint> ret = new ArrayList<>();
		for(PseudoAtomPoint pp:this.pseudoAtoms){
			ArrayList<Double> merged = this.scoring_target.featureGen.mergeFeatures(ad, pp.features);
			double s = this.scoring_target.getLabelScore(FuzzyDecisionTreeScoring_generator.listToArray(merged),label);
			DockingPoint dd = new DockingPoint(a,pp,s,s-basescore);
			ret.add(dd);
		}
		Collections.sort(ret,new DockingScoreComparator());
		Collections.reverse(ret);
		return ret;
	}
	
	
	
	
	public static void main__(String[] args){
		PDBData p  = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\3KB5.pdb");
		for(String cname:p.chains.keySet()){
			PDBChain c = p.chains.get(cname);
			HashMap<PDBResidue,Integer> pi = ScoringDocker.getSurfaceResidues(c.residues,false);
			for(PDBResidue rr:c.residues){
				System.out.println(rr.getRepresentativeCode()+"\t"+pi.get(rr));	
			}
		}
	}
	
	public static void rotate(Point3D start,Point3D end,ArrayList<Point3D> points
			,double radian){
		Point3D basevec = new Point3D(end.x-start.x,end.y-start.y,end.z-start.z);
		
		for(Point3D p:points){
			if(p == start || p == end){
				throw new RuntimeException("list should not contain basepoints.");
			}
			p.x -= start.x;
			p.y -= start.y;
			p.z -= start.z;
		}
		
		for(Point3D p:points){
			Point3D.rotate(p,basevec, radian);
		}
		
		for(Point3D p:points){
			p.x += start.x;
			p.y += start.y;
			p.z += start.z;
		}
		
		
	}
	
	public static void standarize(Point3D p){
		double len = p.x*p.x+ p.y * p.y+ p.z*p.z;
		if(len == 0){
			return;
		}
		len = Math.sqrt(len);
		p.x /= len;
		p.y /= len;
		p.z /= len;
		
	}
	public static Point3D calcNorm(
			Point3D center,Point3D p1,Point3D p2){
		
		return calcNorm(
		p1.x -center.x,
		p1.y -center.y,
		p1.z -center.z,
		p2.x -center.x,
		p2.y -center.y,
		p2.z -center.z				
		);
	
	}
	
	
	
	
	public static Point3D calcNorm(double v1x,double v1y,double v1z,double v2x,double v2y,double v2z){
		Point3D  ret = new Point3D(
		(v1y*v2z-v1z*v2y)*-1,
		(v1z*v2x-v1x*v2z)*-1,
		(v1x*v2y-v1y*v2x)*-1);
		
		standarize(ret);
		return ret;
		
	}
	
	public static void rotateCheck(
			ArrayList<PDBResidue> target
			,ArrayList<PDBResidue> target_ca
			,Point3D p1,Point3D p2
			,ArrayList<PDBResidue> template
			,double resolution
			,AtomDistancePenalty crashpenalty
			
		){
		ArrayList<Point3D> pp = new ArrayList<>();
		for(PDBResidue r:target){
			for(PDBAtom a:r.atoms){
				pp.add(a.loc);
			}
		}
		
		double mincrash = -10000;
		for(double dd = 0.0;dd < Math.PI*2;dd+=resolution){
			rotate(p1,p2,pp,resolution);
			double pt = crashpenalty.dockingCrashPenalty(target_ca, template,crashpenalty.crashPenalty*20);
			if(pt < crashpenalty.crashPenalty*20){
				continue;
			}
			 pt = crashpenalty.dockingCrashPenalty(target, template,crashpenalty.crashPenalty*40);
			if(pt < crashpenalty.crashPenalty*40){
				continue;
			}
			if(mincrash < pt){
				mincrash = pt;
				System.out.println(mincrash);
			}
		}
	}
	public static void moves(ArrayList<Point3D> al,Point3D direc){
		
		for(Point3D pp:al){
			if(pp == direc){
				throw new RuntimeException("exx");
			}
			pp.x += direc.x;
			pp.y += direc.y;
			pp.z += direc.z;
		}
	}
	public static void rotates(ArrayList<Point3D> al,Point3D center,double rx,double ry,double rz){
		
		for(Point3D pp:al){
			pp.x -= center.x;
			pp.y -= center.y;
			pp.z -= center.z;
			double px = pp.x;
			double py = pp.y;
			pp.x = Math.cos(rz)*px-Math.sin(rz)*py;
			pp.y = Math.sin(rz)*px+Math.cos(rz)*py;
			
			px = pp.x;
			py = pp.y;
			double pz = pp.z;
			
			pp.x = Math.cos(ry)*px-Math.sin(ry)*pz;
			pp.z = Math.sin(ry)*px+Math.cos(ry)*pz;
			
			
			
			px = pp.x;
			py = pp.y;
			pz = pp.z;
			
			pp.z = Math.cos(rx)*pz-Math.sin(rx)*py;
			pp.y = Math.sin(rx)*pz+Math.cos(rx)*py;
			
			
			
			pp.x += center.x;
			pp.y += center.y;
			pp.z += center.z;
		}
	}
	public  ArrayList<ArrayList<PDBResidue>> dock(ArrayList<PDBResidue> template
			,ArrayList<PDBResidue> target
	,int docknum
	,int crashtorelance
	){
		return dock(template,target,docknum,crashtorelance,null,null,1000);
	}
	public  ArrayList<ArrayList<PDBResidue>> dock(ArrayList<PDBResidue> template
			,ArrayList<PDBResidue> target
	,int docknum
	,int crashtorelance
	,Point3D restrictA
	,Point3D restrictB
	,double restrictdist){
		ArrayList<PDBResidue> templateres = new ArrayList<>();
		HashMap<PDBResidue,PDBChain> parentchain = new HashMap<>();
		Point3D firstN = target.get(0).getN().loc;
		Point3D lastC = target.get(target.size()-1).getC().loc;
		
		PDBChain c1 = new PDBChain("A");;
		PDBChain c2 = new PDBChain("B");;
		for(PDBResidue t:template){
			parentchain.put(t,t.parent);
			templateres.add(t);
			t.setParent(c1);
		}
		ArrayList<ArrayList<PDBResidue>> ret = new ArrayList<>();
		
		prepareTemplateProtein(templateres);
		AtomDistancePenalty crashpenalty = new  AtomDistancePenalty();
		ArrayList<PDBResidue> targetres = new ArrayList<>();
		
		for(PDBResidue t:target){
			parentchain.put(t,t.parent);
			targetres.add(t);
			t.setParent(c2);
		}
		HashSet<PDBResidue> chk = new HashSet<>(targetres);
		chk.addAll(templateres);
		if(chk.size() != targetres.size()+templateres.size()){
			System.out.println();
		}
		
		
		HashMap<PDBResidue,Integer> pi = ScoringDocker.getSurfaceResidues(targetres,false);
		int 	surfacepoint_threshold = 10;
		scoring_target.prepare(targetres);
		ArrayList<DockingPoint> candidates = new ArrayList<>();
		for(PDBResidue r:targetres){
			ArrayList<PDBAtom> al = scoring_target.featureGen.getTargetAtoms(r);
			if(pi.get(r) < surfacepoint_threshold ){
				continue;
			}
			
			for(PDBAtom a:al){
				ArrayList<Double> d =  scoring_target.featureGen.generateFeatures(a);
				String code = a.parent.getName()+"_"+a.pdb_atom_code;
				double df[] = scoring_target.featureGen.generateFeaturesA(a);
				
				if(df == null){
					continue;
				}
				double sc = scoring_target.getLabelScore(df, code);
				
				ArrayList<DockingPoint> dp = calcMergedScores(a,sc);
				candidates.addAll(dp);
			}
		}
		Collections.sort(candidates,new DockingScoreComparator());
		Collections.reverse(candidates);
		ArrayList<PDBAtom> allatoms = new ArrayList<>();
		ArrayList<Point3D> allpoints = new ArrayList<>();
		ArrayList<PDBResidue> crashchecker_ca = new ArrayList<>();
		
		for(PDBResidue rr:targetres){
			allatoms.addAll(rr.atoms);
			for(PDBAtom a:rr.atoms){
				allpoints.add(a.loc);
			}
			PDBResidue r = new PDBResidue();
			r.setName(rr.getName());
			r.setResidueNumber(rr.getResidueNumber());
			r.atoms.add(rr.getCA());
			crashchecker_ca.add(r);
		}
				
		
		
		
		outer:for(int ii = 0;ii < candidates.size();ii++){
			DockingPoint dp1 = candidates.get(ii);
			for(int jj = ii+1;jj < candidates.size();jj++){
				DockingPoint dp2 = candidates.get(jj);
				if(dp1.targetAtom == dp2.targetAtom
					|| dp1.mappedPosition.pseudo == dp2.mappedPosition.pseudo
						){
					continue;
				}
				for(int kk = jj+1;kk < candidates.size();kk++){
					DockingPoint dp3 = candidates.get(kk);
					if(dp3.targetAtom == dp2.targetAtom
						|| dp3.mappedPosition.pseudo
						== dp2.mappedPosition.pseudo
						|| dp3.targetAtom == dp1.targetAtom
						|| dp3.mappedPosition.pseudo
						== dp1.mappedPosition.pseudo
							){
						continue;
					}
					double dist1 = dp1.targetAtom.distance(dp2.targetAtom);
					double dist2 = dp2.targetAtom.distance(dp3.targetAtom);
					double dist3 = dp3.targetAtom.distance(dp1.targetAtom);
					double pdist1 = dp1.mappedPosition.pseudo.loc.distance(dp2.mappedPosition.pseudo.loc);
					double pdist2 = dp2.mappedPosition.pseudo.loc.distance(dp3.mappedPosition.pseudo.loc);
					double pdist3 = dp3.mappedPosition.pseudo.loc.distance(dp1.mappedPosition.pseudo.loc);
					
					double dthreshold = 2.0;
					if(Math.abs(dist1-pdist1) < dthreshold
					&& Math.abs(dist2-pdist2) < dthreshold
					&&  Math.abs(dist3-pdist3) < dthreshold
					){
						PepProcess.adjustVector3D(
								dp1.targetAtom.loc,
								dp2.targetAtom.loc,
								dp3.targetAtom.loc,
								dp1.mappedPosition.pseudo.loc,
								dp2.mappedPosition.pseudo.loc,
								dp3.mappedPosition.pseudo.loc,
								allpoints);
						Point3D norm = calcNorm(dp1.targetAtom.loc,
								dp2.targetAtom.loc,
								dp3.targetAtom.loc);
						
						double penaltythreshold =crashtorelance;
						double pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
						if(pt < crashpenalty.crashPenalty*penaltythreshold/2){
							moves(allpoints,norm);
							double ppt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							if(ppt < pt){
								norm.x *= -1;
								norm.y *= -1;
								norm.z *= -1;
								moves(allpoints,norm);
								moves(allpoints,norm);
							}
							 pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							while(pt < crashpenalty.crashPenalty*penaltythreshold/2){
								moves(allpoints,norm);
								 pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							}
						}
						
						 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
						if(pt < crashpenalty.crashPenalty*penaltythreshold){
							moves(allpoints,norm);
							double ppt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							if(ppt < pt){
								norm.x *= -1;
								norm.y *= -1;
								norm.z *= -1;
								moves(allpoints,norm);
								moves(allpoints,norm);
							}
							 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							while(pt < crashpenalty.crashPenalty*penaltythreshold){
								moves(allpoints,norm);
								 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							}
							
						}
						
						ArrayList<PDBResidue> pp = new ArrayList<>();
						for(PDBResidue d:targetres){
							pp.add(d.getCopy());
						}
						if(restrictA != null && restrictB != null){
							if(restrictA.distance(firstN)+restrictB.distance(lastC) < restrictdist){
								ret.add(pp);
							}
						}else if(restrictA != null){
							if(restrictA.distance(firstN)< restrictdist){
								ret.add(pp);
							}
						}else if(restrictB != null){
							if(restrictB.distance(lastC)< restrictdist){
								ret.add(pp);
							}
						}else{
							ret.add(pp);
							
						}
						if(ret.size() >= docknum){
							break outer;
						}
					}
				}
			}
			
			if(ret.size() >= docknum){
				break;
			}
		}
		for(PDBResidue p:parentchain.keySet()){
			p.setParent(parentchain.get(p));
		}
		
		return ret;
	}
	public static void main(String[] args){
		ScoringDocker sd = new ScoringDocker();
		PDBData p  = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\casp13\\queries\\T0980s1\\model0.pfas.pdb");
		ArrayList<PDBResidue> templateres = new ArrayList<>();
		for(String cname:p.chains.keySet()){
			PDBChain c = p.chains.get(cname);
			templateres.addAll(c.residues);
		}
		
		PDBData p2  = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\casp13\\queries\\T0980s2\\model0.pfas.pdb");
		ArrayList<PDBResidue> targetres = new ArrayList<>();
		for(String cname:p2.chains.keySet()){
			PDBChain c = p2.chains.get(cname);
			targetres.addAll(c.residues);
		}
		
		ArrayList<FloatingResidue> templateres_f = ChainBuilder.changeToFloating(templateres);
		ArrayList<FloatingResidue> targetres_f = ChainBuilder.changeToFloating(targetres);
		templateres.clear();
		targetres.clear();
		templateres.addAll(templateres_f);
		targetres.addAll(targetres_f);
		
		sd.prepareTemplateProtein(templateres);
		AtomDistancePenalty crashpenalty = new  AtomDistancePenalty();
		
		HashMap<PDBResidue,Integer> pi = ScoringDocker.getSurfaceResidues(targetres,false);
		HashMap<PDBResidue,Integer> pi2 = ScoringDocker.getSurfaceResidues(templateres,false);
		ArrayList<FloatingResidue> flres = new ArrayList<>();
		
		int 	surfacepoint_threshold = 10;
		for(FloatingResidue r:targetres_f){
			if(pi.get(r) < surfacepoint_threshold ){
				continue;
			}
			flres.add(r);
		}
		
		for(FloatingResidue r:templateres_f){
			if(pi2.get(r) < surfacepoint_threshold ){
				continue;
			}
			flres.add(r);
		}
		
		
		sd.scoring_target.prepare(targetres);
		ArrayList<DockingPoint> candidates = new ArrayList<>();
		for(PDBResidue r:targetres){
			if(pi.get(r) < surfacepoint_threshold ){
				continue;
			}
			ArrayList<PDBAtom> al = sd.scoring_target.featureGen.getTargetAtoms(r);
			for(PDBAtom a:al){
				ArrayList<Double> d =  sd.scoring_target.featureGen.generateFeatures(a);
				String code = a.parent.getName()+"_"+a.pdb_atom_code;
				double df[] = sd.scoring_target.featureGen.generateFeaturesA(a);
				
				if(df == null){
					continue;
				}
				double sc = sd.scoring_target.getLabelScore(df, code);
				
				ArrayList<DockingPoint> dp = sd.calcMergedScores(a,sc);
				candidates.addAll(dp);
				//System.out.println(pi.get(r)+"\t"+a.parent.getRepresentativeCode()+"\t"+a.pdb_atom_code+"\t"+sc+"\t"+dp.get(0).scorediff);
			}
		}
		Collections.sort(candidates,new DockingScoreComparator());
		Collections.reverse(candidates);
		
		
		if(candidates.size() > 10000){
			candidates = new ArrayList(candidates.subList(0, 10000));
		}
		ArrayList<PDBAtom> allatoms = new ArrayList<>();
		ArrayList<Point3D> allpoints = new ArrayList<>();
		ArrayList<PDBResidue> crashchecker_ca = new ArrayList<>();
		
		for(PDBResidue rr:targetres){
			allatoms.addAll(rr.atoms);
			for(PDBAtom a:rr.atoms){
				allpoints.add(a.loc);
			}
			PDBResidue r = new PDBResidue();
			r.setName(rr.getName());
			r.setResidueNumber(rr.getResidueNumber());
			r.atoms.add(rr.getCA());
			crashchecker_ca.add(r);
		}
		ChainBuilder cb = new ChainBuilder();
		cb.prepare(flres);
		int dcount = 0;
		for(int ii = 0;ii < candidates.size();ii++){
			DockingPoint dp1 = candidates.get(ii);
			outer:for(int jj = ii+1;jj < candidates.size();jj++){
				DockingPoint dp2 = candidates.get(jj);
				if(dp1.targetAtom == dp2.targetAtom
					|| dp1.mappedPosition.pseudo == dp2.mappedPosition.pseudo
						){
					continue;
				}
				for(int kk = jj+1;kk < candidates.size();kk++){
					DockingPoint dp3 = candidates.get(kk);
					if(dp3.targetAtom == dp2.targetAtom
						|| dp3.mappedPosition.pseudo
						== dp2.mappedPosition.pseudo
						|| dp3.targetAtom == dp1.targetAtom
						|| dp3.mappedPosition.pseudo
						== dp1.mappedPosition.pseudo
							){
						continue;
					}
					double dist1 = dp1.targetAtom.distance(dp2.targetAtom);
					double dist2 = dp2.targetAtom.distance(dp3.targetAtom);
					double dist3 = dp3.targetAtom.distance(dp1.targetAtom);
					double pdist1 = dp1.mappedPosition.pseudo.loc.distance(dp2.mappedPosition.pseudo.loc);
					double pdist2 = dp2.mappedPosition.pseudo.loc.distance(dp3.mappedPosition.pseudo.loc);
					double pdist3 = dp3.mappedPosition.pseudo.loc.distance(dp1.mappedPosition.pseudo.loc);
					
					double dthreshold = 2.0;
					if(Math.abs(dist1-pdist1) < dthreshold
					&& Math.abs(dist2-pdist2) < dthreshold
					&&  Math.abs(dist3-pdist3) < dthreshold
					){
						/*
						
						public static void adjustVector3D(Point3D targetstart,
			Point3D targetend,
			Point3D targetend2,
			Point3D refstart,
			Point3D refend,
			Point3D refend2,
			ArrayList<Point3D> al)
						*/
						PepProcess.adjustVector3D(
								dp1.targetAtom.loc,
								dp2.targetAtom.loc,
								dp3.targetAtom.loc,
								dp1.mappedPosition.pseudo.loc,
								dp2.mappedPosition.pseudo.loc,
								dp3.mappedPosition.pseudo.loc,
								allpoints);
						Point3D norm = calcNorm(dp1.targetAtom.loc,
								dp2.targetAtom.loc,
								dp3.targetAtom.loc);
						
						double penaltythreshold = 10;
						
						double pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
						if(pt < crashpenalty.crashPenalty*penaltythreshold/2){
							moves(allpoints,norm);
							
							double ppt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							
							if(ppt < pt){
								norm.x *= -1;
								norm.y *= -1;
								norm.z *= -1;
								moves(allpoints,norm);
								moves(allpoints,norm);
							}
							

							 pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							while(pt < crashpenalty.crashPenalty*penaltythreshold/2){
								moves(allpoints,norm);
								 pt = crashpenalty.dockingCrashPenalty(crashchecker_ca, templateres,crashpenalty.crashPenalty*penaltythreshold/2);
							}
						}
						
						for(FloatingResidue r:flres){
							cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
						}
						 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
						if(pt < crashpenalty.crashPenalty*penaltythreshold){
							moves(allpoints,norm);

							for(FloatingResidue r:flres){
								cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
							}

							double ppt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							if(ppt < pt){
								norm.x *= -1;
								norm.y *= -1;
								norm.z *= -1;
								moves(allpoints,norm);
								moves(allpoints,norm);
							}
							
								for(FloatingResidue r:flres){
									cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
								}
							 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							while(pt < crashpenalty.crashPenalty*penaltythreshold){
								moves(allpoints,norm);	
								
								for(FloatingResidue r:flres){
									cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
								}
								 pt = crashpenalty.dockingCrashPenalty(targetres, templateres,crashpenalty.crashPenalty*penaltythreshold);
							}
						}
						
			writeToFile(makeAtomLines(targetres,1,1,"B"),
					"C:\\dummy\\vbox_share\\casp13\\queries\\T0980s1\\dock10_model1."+dcount+".pdb");
						dcount++;
						break outer;
					}
				}
			}
			if(dcount > 20){
				break;
			}
			
		}
		
		
		
	}
	
	
	
}

class DockingPoint{
	PDBAtom targetAtom;
	PseudoAtomPoint mappedPosition;
	double score;
	double scorediff;
	DockingPoint(PDBAtom t,PseudoAtomPoint p,double d,double d2){
		targetAtom = t;
		mappedPosition = p;
		score = d;
		scorediff = d2;
	}
	
}


class PseudoAtomPoint{
	PDBAtom parent = null;
	PDBAtom pseudo = null;
	ArrayList<Double> features = null;
	PseudoAtomPoint(Point3D p3,PDBAtom a){
		parent = a;
		PDBAtom ap = new PDBAtom();
		ap.pdb_atom_code = "H";
		ap.loc.set(p3);
		ap.parent = a.parent;
		pseudo = ap;
	}
	public void setFeatures(ArrayList<Double> f ){
		features = f;
	}
}
class DockingScoreComparator implements Comparator<DockingPoint>{
	@SuppressWarnings("unchecked")
	public int compare(DockingPoint arg1, DockingPoint arg2){
		
		if(arg1.scorediff < arg2.scorediff ){
			return -1;
		}
		if(arg1.scorediff == arg2.scorediff ){
			return 0;
		}
			return 1;
	}
	
}
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
import static pepbuilderj.ChainBuilder.changeToFloating;
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;

/**
 *
 * @author kimidori
 */
public class RemakePDB {
	/**
	 * Refine までされなかったPDBについてギャップFilling 等を行う。
	 */
	
	public static ArrayList<ArrayList<FloatingResidue>> getUniqueResidue(ArrayList<ArrayList<FloatingResidue>> alis){
		HashSet<Integer> resused = new HashSet<>();
		ArrayList<ArrayList<FloatingResidue>> ret = new ArrayList<>();
		for(ArrayList<FloatingResidue> rlis:alis){
			ArrayList<FloatingResidue> res = new ArrayList<>();
			for(FloatingResidue r:rlis){
				if(resused.contains(r.getResidueNumber())){
					if(res.size() > 5){
						for(PDBResidue rr:res){
							resused.add(rr.getResidueNumber());
						}
						ret.add(res);
					}
					res = new ArrayList<>();
				}else{

					res.add(r);
				}
			}
			if(res.size() > 5){
				ret.add(res);
				for(PDBResidue r:res){
					resused.add(r.getResidueNumber());
				}
			}

		}
		return ret;
	}
	public static void main(String[] args){
		for(int s = 1;s < 10;s++){
			String filename = "C:\\dummy\\vbox_share\\casp13\\queries\\T0990\\threadingresult\\T0990."+s+".th.pdb";
			PDBData pdb = PDBData.loadPDBFile(filename);
			ArrayList<ArrayList<PDBAtom>> connected = new ArrayList<>();
			double threshold = 8.0;
			ArrayList<FloatingResidue> allresidue = new ArrayList<>();
			ArrayList<ArrayList<PDBAtom>> peptidebonds = new ArrayList<>();

			HashSet<FloatingResidue> terminal = new HashSet<>();//端っこかマップされなかった部分
			for(String c:pdb.chains.keySet()){
				PDBChain cc = pdb.chains.get(c);
				ArrayList<FloatingResidue> res = changeToFloating(cc.residues);
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				for(int ii = 0;ii < res.size();ii++){
					FloatingResidue r = res.get(ii);
					if(r.next != null ){
						if(r.getC().distance(r.next.getN()) > 1.4){
							terminal.add(r);
							terminal.add(r.next);
						}

						ArrayList<PDBAtom> p = new ArrayList<>();
						p.add(r.getC());
						p.add(r.next.getN());
						peptidebonds.add(p);
					}
					if(r.prev != null ){
						if(r.getN().distance(r.prev.getC()) > 1.4){
							terminal.add(r);
							terminal.add(r.prev);
						}
						ArrayList<PDBAtom> p = new ArrayList<>();
						p.add(r.prev.getC());
						p.add(r.getN());
						peptidebonds.add(p);
					}
				}
				for(FloatingResidue r:res){
					allresidue.add(r);
					if(terminal.contains(r)){
						continue;
					}
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

			//CA が一定距離にある場合結合しているとみなす
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
			//マージされて空になっている場合削除
			Iterator<ArrayList<PDBAtom>> ite = connected.iterator();
			while(ite.hasNext()){
				ArrayList<PDBAtom> au = ite.next();
				if(au.size() == 0){
					ite.remove();
				}
			}
			//System.out.println(";;");
			//移動が楽なように Point3D に直してまた同一残基のものも含める
			//fixme blockid でなく block オブジェクトを作って管理する
			HashMap<PDBAtom,Integer> atom_to_blockid = new HashMap<>();
			ArrayList<ArrayList<Point3D>> blocks = new ArrayList<>();
			HashSet<Point3D> crash_ignore = new HashSet<>();//マップされてない原子なのでクラッシュは無視
			for(int ii = 0;ii < connected.size();ii++){
				ArrayList<Point3D> block = new ArrayList<>();
				int blockid = blocks.size();
				for(PDBAtom aa:connected.get(ii)){
					for(PDBAtom a:aa.parent.atoms){
						block.add(a.loc);
						atom_to_blockid.put(a,blockid);
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
			for(FloatingResidue tt:terminal){
				ArrayList<Point3D> block = new ArrayList<>();
				int blockid = blocks.size();
				for(PDBAtom a:tt.atoms){
					block.add(a.loc);
					atom_to_blockid.put(a,blockid);
					crash_ignore.add(a.loc);
				}
				blocks.add(block);
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
				distsum += Math.abs(al.get(0).distance(al.get(1)) - 1.30)
						*(blocks.get(atom_to_blockid.get(al.get(0))).size()
						+blocks.get(atom_to_blockid.get(al.get(1))).size());
			}
			
			int maxid = 0;
			int maxcount = -1;
			for(int rr = 0;rr < blocks.size();rr++){
				if(maxcount < blocks.get(rr).size()){
					maxcount  = blocks.get(rr).size();
					maxid = rr;
				}
			}
			
			
			for(int ii = 0;ii < 100;ii++){
				System.out.println(distsum);
				
				HashMap<Integer,Point3D> blockid_to_direc = new HashMap<>();
				for(ArrayList<PDBAtom> al:peptidebonds){
					if(al.get(0).distance(al.get(1)) - 1.30 > 0){
						ArrayList<Point3D> b0 = blocks.get(atom_to_blockid.get(al.get(0)));
						ArrayList<Point3D> b1 = blocks.get(atom_to_blockid.get(al.get(1)));
						double factor  = 1.0;
						if(b0 == b1){
						}else{
							if(b0.size() < 50  && b1.size() < 50){
								factor = 0.2;
							}
							if(Math.random() < 0.5){
							//if(b0.size() < b1.size()){
								// fixme 他でも指定されていた場合平均にする
								Point3D p = new Point3D(
								al.get(0).loc.x-al.get(1).loc.x*factor,
								al.get(0).loc.y-al.get(1).loc.y*factor,
								al.get(0).loc.z-al.get(1).loc.z*factor
								);
								blockid_to_direc.put(atom_to_blockid.get(al.get(1)),p);
								
							}else{
								Point3D p = new Point3D(
								al.get(1).loc.x-al.get(0).loc.x*factor,
								al.get(1).loc.y-al.get(0).loc.y*factor,
								al.get(1).loc.z-al.get(0).loc.z*factor
								);
								blockid_to_direc.put(atom_to_blockid.get(al.get(0)),p);
							}
						}
					}
				}
				//for(int rr = 0;rr < blocks.size();rr++){
				//	if(maxid == rr){
				//		continue;
				//	}
				for(Integer rr:blockid_to_direc.keySet()){
					//int randcode = (int)(Math.random()*blocks.size());
					int randcode = rr;
					ArrayList<Point3D> a = blocks.get(randcode);
					Point3D direc = blockid_to_direc.get(rr);
					double fact = 0.8;
					double len = direc.x*direc.x+direc.y*direc.y+direc.z*direc.z;
					if(len > 0){
						len = Math.sqrt(len);
					}
					if(len < 1.2){
						direc.x *= -1;
						direc.y *= -1;
						direc.z *= -1;
					}
					if(ii%2 == 0){
						fact = 0.3;
					}
					direc.x *= fact;
					direc.y *= fact;
					direc.z *= fact;
					
					direc.x += (Math.random()-0.5)*5*fact;
					direc.y += (Math.random()-0.5)*5*fact;
					direc.z += (Math.random()-0.5)*5*fact;
					RandomDocker.moves(a, direc);
					
					boolean cflag = false;
					outer:for(Point3D p:a){
						if(crash_ignore.contains(p)){
							continue;
						}
						for(int ll = 0;ll < allpoint.size();ll++){
							Point3D p2 = allpoint.get(ll);
							
							if(crash_ignore.contains(p2)){
								continue;
							}
							if(p == p2){
								continue;
							}
							
							if(p.distance(p2) < 2.0){
								if(!crashedpair.get(p).contains(p2)){
									cflag = true;
									break outer;
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
						cdistsum +=  Math.abs(al.get(0).distance(al.get(1)) - 1.30)
								*(blocks.get(atom_to_blockid.get(al.get(0))).size()
								+blocks.get(atom_to_blockid.get(al.get(1))).size());
					}
					if(cdistsum > distsum){
						for(Point3D p:a){
							p.set(prevpoint.get(p));
						}
					}else{
						distsum = cdistsum;
						for(Point3D p:a){
							prevpoint.get(p).set(p);
						}
					}
				}
			}
			ChainBuilder cb = new ChainBuilder();
			ChainBuilder.saveChains(allresidue, filename+".f.pdb");
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

			//fixme 露出残基をとって、露出している方から削っていく
			//インサーションが入った所為で別残基とペプチド結合を作ってしまっているようなのはのける
			//ChainBuilder.saveChains(modelled,outfilename+".p.pdb");
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
			ChainBuilder.saveChains(allresidue, filename+".c.pdb");
		}
	}		
	

	
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author kimidori
 */
public class AliFine {
	
	ArrayList<RResidue> target = new ArrayList<>();
	ArrayList<RResidue> template = new ArrayList<>();
	HashMap<RResidue,RResidue> target_template = new HashMap<>();
	HashMap<RResidue,RResidue> template_target = new HashMap<>();
	
	HashMap<PDBResidue,RResidue> target_rr = new HashMap();
	HashMap<PDBResidue,RResidue> template_rr = new HashMap();
	
	ArrayList<RResidue> align_target = new ArrayList<>();
	ArrayList<RResidue> align_template = new ArrayList<>();
	
	
	public static TBMResult refined(TBMResult res){
		AliFine ali = new AliFine();
		TBMResult ret = new TBMResult();
		ret.residues = new ArrayList<>(res.residues);
		ret.template_residues = new ArrayList<>(res.template_residues);
		ret.target_template_map = new HashMap<>(res.target_template_map);
		ret.template_target_map = new HashMap<>(res.template_target_map);
		
		int alisize = ret.template_residues.size()+ret.residues.size();
		alisize*=2;
		for(int ii = 0;ii < alisize;ii++){
			ali.align_target.add(null);
			ali.align_template.add(null);
		}
		
		
		
		
		
		
		
		
		ArrayList<ArrayList<PDBResidue>> tals = new ArrayList<>();
		tals.add(res.residues);
		tals.add(res.template_residues);
		ArrayList<HashMap<PDBResidue,RResidue>> tal_rr = new ArrayList<>();
		tal_rr.add(ali.target_rr);
		tal_rr.add(ali.template_rr);
		ArrayList<ArrayList<RResidue>> tal_as = new ArrayList<>();
		tal_as.add(ali.target);
		tal_as.add(ali.template);
		for(int tt = 0;tt < 2;tt++){
			ArrayList<PDBResidue> tal = tals.get(tt);
			ArrayList<RResidue> tal_a = tal_as.get(tt);
			HashMap<PDBResidue,RResidue> to_rr = tal_rr.get(tt);
			for(int ii = 0;ii  < tal.size();ii++){
				PDBResidue r = tal.get(ii);
				RResidue rr = new RResidue();
				to_rr.put(r,rr);
				rr.current = r;
				if(ii > 0){
					rr.prev = tal.get(ii-1);
				}
				if(ii < tal.size() -1){
					rr.next = tal.get(ii+1);
				}
				tal_a.add(rr);
			}
		}
		for(PDBResidue r:res.target_template_map.keySet()){
			ali.target_template.put(ali.target_rr.get(r),ali.template_rr.get(res.target_template_map.get(r)));
		}
		for(PDBResidue r:res.template_target_map.keySet()){
			ali.template_target.put(ali.template_rr.get(r),ali.target_rr.get(res.template_target_map.get(r)));
		}
		
		HashSet<RResidue> template_used = new HashSet<>();
		int spos = 0;
		//ペアワイズアラインメントを再構築する
		for(int ii = 0;ii < ali.target.size();ii++){
			RResidue rr = ali.target.get(ii);
			if(ali.target_template.get(rr) != null){
				RResidue tt = ali.target_template.get(rr);
				ArrayList<RResidue> prevs = new ArrayList<>();
				while(tt.prev != null && !template_used.contains(ali.template_rr.get(tt.prev))){
					prevs.add(ali.template_rr.get(tt.prev));
					tt = ali.template_rr.get(tt.prev);
				}
				for(int jj = 0;jj < prevs.size();jj++){
					ali.align_target.set(spos,null);
					ali.align_template.set(spos,prevs.get(jj));
					spos++;
				}
				template_used.addAll(prevs);
				ali.align_target.set(spos, rr);
				ali.align_template.set(spos,ali.target_template.get(rr));
				template_used.add(ali.target_template.get(rr));
				spos++;
			}else{
				ali.align_target.set(spos, rr);
				ali.align_template.set(spos,null);
				spos++;
			}
		}
		for(int ii = 0;ii < ali.template.size();ii++){
			if(!template_used.contains(ali.template.get(ii))){
				ali.align_target.set(spos, null);
				ali.align_template.set(spos,ali.template.get(ii));
				spos++;
			}
		}
		
		
		
		HashMap<PDBResidue,Integer> surface = RandomDocker.getSurfaceResidues(ret.template_residues, false);
		
		
		HashMap<RResidue,Double> bury_score = new HashMap<>();
		//Fixme  Residue ごとの 平均値かMAXをとって割合にする
		for(PDBResidue s:surface.keySet()){
			bury_score.put(ali.template_rr.get(s),1.0/(surface.get(s)+1.0));
		}
		//Template にギャップが入っており、モデルしようとしている方がインサーション
		
		int tarpos = -1;
		int tempos = -1;
		
		for(int ii = 0;ii < ali.align_target.size();ii++){
			if(ali.align_target.get(ii) != null){
				ali.align_target.get(ii).previndex = ii;
				ali.align_target.get(ii).currentindex = ii;
				ali.align_target.get(ii).maxindex = ii;
			}
			if(ali.align_template.get(ii) != null){
				ali.align_template.get(ii).previndex = ii;
				ali.align_template.get(ii).currentindex = ii;
				ali.align_template.get(ii).maxindex = ii;
			}
		}
		
		ArrayList<RResidue> ptarget = ali.align_target;
		ArrayList<RResidue> ppair = ali.align_template;
		
		for(int ii = 0;ii < ptarget.size();ii++){
			if(ptarget.get(ii) != null){
				tarpos++;
			}
			if(ppair.get(ii) != null){
				tempos++;
			}
			if(tarpos == -1 ||tempos == -1){
				continue;
			}
			if(ptarget.get(ii) != null && ppair.get(ii) == null){
				//テンプレートの方にギャップ
				int st = ii;
				int rcount = 0;
				int maxj = -1;
				int endj = -1;
				double cscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
				for(int jj = st;jj < ppair.size();jj++){
					if(ppair.get(jj) != null){
						int pss = fillgap_Prev(ppair,ptarget,jj);
						double dscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
						if(dscore < cscore){
							maxj = jj;
							cscore = dscore;
							for(int kk = st;kk <= jj;kk++){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}
						}
						rcount++;
						if(rcount  > 4){
							endj = jj;
							break;
						}
					}
				}
				recoverPrevState(ppair);
				int maxj2 = -1;
				rcount = 0;
				
				for(int jj = st;jj >= 0;jj--){
					
					if(ppair.get(jj) != null){
						
						if(ppair.get(jj).next == null){
							//fixme もっと前で判定しているべき
							break;
						}
						int pss = fillgap_Next(ppair,ptarget,jj);
						double dscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
						if(dscore < cscore){
							maxj2 = jj;
							cscore = dscore;
							for(int kk = st;kk >= jj;kk--){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}

							for(int kk = st;kk <= endj;kk++){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}
						}
						rcount++;
						if(rcount  > 4){
							break;
						}
					}
				}
				
				recoverPrevState(ppair);
				if(maxj > -1 || maxj2 > -1){
					int dend = -1;
					for(int jj = 0;jj < ppair.size();jj++){
						if(ppair.get(jj).currentindex != ppair.get(jj).maxindex){
							dend = jj;
						}
						if(ppair.get(jj).prev == null){
							break;
						}
					}
					if(dend > -1){
						ii = dend;
					}
					setMaxState(ppair);
					
					for(int zz = 0;zz < ali.align_target.size();zz++){
						
						if(ali.align_target.get(zz) != null){
							ali.align_target.get(zz).previndex = zz;
							ali.align_target.get(zz).currentindex = zz;
						}
						
						if(ali.align_template.get(zz) != null){
							ali.align_template.get(zz).previndex = zz;
							ali.align_template.get(zz).currentindex = zz;
						}
					}
				}
				
				while(ptarget.get(ii) != null && ppair.get(ii) == null){
					ii++;
				}
			}
		}
		for(int ii = 0;ii < ali.align_target.size();ii++){
			if(ali.align_target.get(ii) != null){
				ali.align_target.get(ii).previndex = ii;
				ali.align_target.get(ii).currentindex = ii;
				ali.align_target.get(ii).maxindex = ii;
			}
			if(ali.align_template.get(ii) != null){
				ali.align_template.get(ii).previndex = ii;
				ali.align_template.get(ii).currentindex = ii;
				ali.align_template.get(ii).maxindex = ii;
			}
		}
		
		tarpos = -1;
		tempos = -1;
		
		ptarget = ali.align_template;
		ppair = ali.align_target;
		
		for(int ii = 0;ii < ptarget.size();ii++){
			if(ptarget.get(ii) != null){
				tarpos++;
			}
			if(ppair.get(ii) != null){
				tempos++;
			}
			if(tarpos == -1 ||tempos == -1){
				continue;
			}
			if(ptarget.get(ii) != null && ppair.get(ii) == null){
				//クエリの方にギャップ
				int st = ii;
				int rcount = 0;
				int maxj = -1;
				int endj = -1;
				double cscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
				for(int jj = st;jj < ppair.size();jj++){
					if(ppair.get(jj) != null){
						int pss = fillgap_Prev(ppair,ptarget,jj);
						double dscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
						if(dscore < cscore){
							maxj = jj;
							cscore = dscore;
							for(int kk = st;kk <= jj;kk++){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}
						}
						rcount++;
						if(rcount  > 4){
							endj = jj;
							break;
						}
					}
				}
				recoverPrevState(ppair);
				int maxj2 = -1;
				rcount = 0;
				
				for(int jj = st;jj >= 0;jj--){
					
					if(ppair.get(jj) != null){
						
						if(ppair.get(jj).next == null){
							//fixme もっと前で判定しているべき
							break;
						}
						int pss = fillgap_Next(ppair,ptarget,jj);
						double dscore = calcNonBuryScore(ali.align_target,ali.align_template,bury_score);
						if(dscore < cscore){
							maxj2 = jj;
							cscore = dscore;
							for(int kk = st;kk >= jj;kk--){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}

							for(int kk = st;kk <= endj;kk++){
								if(ppair.get(kk) != null){
									ppair.get(kk).maxindex= ppair.get(kk).currentindex;
								}
							}
						}
						rcount++;
						if(rcount  > 4){
							break;
						}
					}
				}
				
				recoverPrevState(ppair);
				if(maxj > -1 || maxj2 > -1){
					int dend = -1;
					for(int jj = 0;jj < ppair.size();jj++){
						if(ppair.get(jj).currentindex != ppair.get(jj).maxindex){
							dend = jj;
						}
						if(ppair.get(jj).prev == null){
							break;
						}
					}
					if(dend > -1){
						ii = dend;
					}
					setMaxState(ppair);
					for(int zz = 0;zz < ali.align_target.size();zz++){
						if(ali.align_target.get(zz) != null){
							ali.align_target.get(zz).previndex = zz;
							ali.align_target.get(zz).currentindex = zz;
						}
						if(ali.align_template.get(zz) != null){
							ali.align_template.get(zz).previndex = zz;
							ali.align_template.get(zz).currentindex = zz;
						}
					}
				}
				while(ptarget.get(ii) != null && ppair.get(ii) == null){
					ii++;
				}
				
			}
			
		}
		ret.target_template_map = new HashMap<>();
		ret.template_target_map = new HashMap<>();
		for(int ii = 0;ii < ali.align_target.size();ii++){
			if(ali.align_target.get(ii) != null && ali.align_template.get(ii) != null){
				ret.target_template_map.put(ali.align_target.get(ii).current,ali.align_template.get(ii).current);
				ret.template_target_map.put(ali.align_template.get(ii).current,ali.align_target.get(ii).current);
			}
		}
		ret.mapped.clear();
		for(int ii = 0;ii < ret.residues.size();ii++){
			if(ret.target_template_map.containsKey(ret.residues.get(ii))){
				ret.mapped.add(true);
			}else{
				ret.mapped.add(false);
			}
		}
		
		
		return ret;
	}
	
	public static void setMaxState(ArrayList<RResidue> tar){
		HashMap<Integer,RResidue> prevmap = new HashMap<>();
		for(int ii = 0;ii < tar.size();ii++){
			if(tar.get(ii) != null){
				tar.get(ii).currentindex = ii;
				if(prevmap.containsKey(tar.get(ii))){
					throw new RuntimeException("error in code ???");
				}
				if(tar.get(ii).maxindex != tar.get(ii).currentindex){
					prevmap.put(tar.get(ii).maxindex,tar.get(ii));		}
				}
		}
		for(Integer ii:prevmap.keySet()){
			tar.set(prevmap.get(ii).currentindex,null);
		}
		for(Integer ii:prevmap.keySet()){
			tar.set(ii,prevmap.get(ii));
			prevmap.get(ii).currentindex = ii;
			prevmap.get(ii).maxindex = ii;
		}
	}
	public static void recoverPrevState(ArrayList<RResidue> tar){
		HashMap<Integer,RResidue> prevmap = new HashMap<>();
		for(int ii = 0;ii < tar.size();ii++){
			if(tar.get(ii) != null){
				tar.get(ii).currentindex = ii;
				if(prevmap.containsKey(tar.get(ii))){
					throw new RuntimeException("error in code ???");
				}
				if(tar.get(ii).previndex != tar.get(ii).currentindex){
					prevmap.put(tar.get(ii).previndex,tar.get(ii));		}
				}
		}
		for(Integer ii:prevmap.keySet()){
			tar.set(prevmap.get(ii).currentindex,null);
		}
		for(Integer ii:prevmap.keySet()){
			tar.set(ii,prevmap.get(ii));
			prevmap.get(ii).currentindex = ii;
		}
	}
	
	
	/**
	 * 前方向に詰める。
	 * @param target
	 * @param pair
	 * @param targetindex 
	 */
	public static int fillgap_Prev(ArrayList<RResidue> target,ArrayList<RResidue> pair,int targetindex){
		if(targetindex == 0){
			throw new RuntimeException("??? error in code");
		}
		int nprev = -1;
		for(int ii = targetindex-1;ii >= 0;ii--){
			if(target.get(ii) != null){
				nprev = ii;
				break;
			}
		}
		for(int ii = nprev+1;ii < targetindex;ii++){
			if(pair.get(ii) != null){
				target.get(targetindex).currentindex = ii;
				target.set(ii,target.get(targetindex));
				target.set(targetindex,null);
				return ii;
			}
		}
		return targetindex;
	}
	/**
	 * 後ろ方向に詰める。
	 * @param target
	 * @param pair
	 * @param targetindex 
	 */
	public static int fillgap_Next(ArrayList<RResidue> target,ArrayList<RResidue> pair,int targetindex){
		if(target.get(targetindex).next == null){
			throw new RuntimeException("??? error in code");
		}
		int nnext = -1;
		for(int ii = targetindex+1;ii < target.size();ii++){
			if(target.get(ii) != null){
				nnext = ii;
				break;
			}
		}
		for(int ii = nnext;ii > targetindex;ii--){
			if(pair.get(ii) != null){
				target.get(targetindex).currentindex = ii;
				target.set(ii,target.get(targetindex));
				target.set(targetindex,null);
				return ii;
			}
		}
		return targetindex;
	}
	/**
	 * bury されているべき部分が露出していることを示すスコアを返す。
	 * 低い方が良い。
	 * @param ali_tar
	 * @param ali_tm
	 * @param buryratio
	 * @return 
	 */
	public static double calcNonBuryScore(ArrayList<RResidue> ali_tar,ArrayList<RResidue> ali_tm
			,HashMap<RResidue,Double> buryratio){
		double ret = 0;
		int tstart = -1;
		int tend = 0;
				
		for(int ii = 0;ii < ali_tar.size();ii++){
			if(ali_tar.get(ii) != null){
				if(tstart == -1){
					tstart = ii;
				}
				tend = ii;
			}
		}
		for(int ii = tstart;ii <= tend;ii++){
			if(ali_tar.get(ii) != null && ali_tm.get(ii) == null){
				int st = ii;
				double sc = 0;
				int cou = 0;
				for(int jj = st;jj < ali_tar.size();jj++){
					if( ali_tm.get(jj) == null){
					}else{
						sc += buryratio.get(ali_tm.get(jj));
						cou++;
						break;	
					}
				}
				for(int jj = st;jj >= 0 ;jj--){
					if( ali_tm.get(jj) == null){
					}else{
						sc += buryratio.get(ali_tm.get(jj));
						cou++;
						break;	
					}
				}
				if(cou > 0){
					ret += sc/cou;
				}
			}
			if(ali_tar.get(ii) == null && ali_tm.get(ii) != null){
				ret += buryratio.get(ali_tm.get(ii));
			}
		}
		
		return ret;
	}
	
	
	public ArrayList<GapSegment> getSegment(boolean temp){
		ArrayList<RResidue> tal = target;
		HashMap<RResidue,RResidue> talmap = target_template;
		ArrayList<GapSegment> gaps = new ArrayList<>();
		if(temp){
			tal = template;
			talmap = template_target;
		}
		for(int ii = 0;ii < tal.size();ii++){
			if(talmap.get(tal.get(ii)) == null){
				int gstart = ii;
				int gend = ii;
				for(int jj = ii+1;jj < tal.size();jj++){
					if(talmap.get(tal.get(ii)) == null){
						gend = jj;
					}else{
						break;
					}
				}
				ii = gend;
				GapSegment gs = new GapSegment(gstart,gend,tal,this);
				gaps.add(gs);
			}
		}
		return gaps;
	}
	
}
class RResidue{
	PDBResidue current = null;
	PDBResidue next = null;
	PDBResidue prev = null;
	int previndex = -1;
	int currentindex = -1;
	int maxindex = -1;
}
class GapSegment{
	int segstart = -1;
	int segend = -1;
	ArrayList<RResidue> list = null;
	AliFine ali = null;
	GapSegment(int s,int e,ArrayList<RResidue> l,AliFine a){
		segstart = s;
		segend = e;
		list = l;
		ali = a;
	}
	public RResidue getStart(){
		return list.get(segstart);
	}
	public RResidue getEnd(){
		return list.get(segend);
	}
	/**
	 * 相手側の次の残基を得る
	 * @return 
	 */
	public RResidue getCompNext(){
		
		int bstart = -1;
		HashMap<RResidue,RResidue> formap = null;
		HashMap<RResidue,RResidue> revmap = null;
		ArrayList<RResidue> revlist = null;
		HashMap<PDBResidue,RResidue> formap_rr = null;
		HashMap<PDBResidue,RResidue> revmap_rr = null;
		if(list == ali.target){
			formap = ali.target_template;
			revmap = ali.template_target;
			revlist = ali.template;
			formap_rr = ali.target_rr;
			revmap_rr = ali.template_rr;
		}else{
			formap = ali.template_target;
			revmap= ali.target_template;
			revlist = ali.target;
			formap_rr = ali.template_rr;
			revmap_rr = ali.target_rr;
		}
			
		for(int ii = segstart;ii <= segend;ii++){
			if(formap.get(list.get(ii)) != null){
				bstart = ii;
			}
		}
		if(bstart == -1){
			for(int ii = segstart;ii > -1;ii--){
				if(formap.get(list.get(ii)) != null){
					bstart = ii;
					break;
				}
			}
		}
		if(bstart > -1){
			RResidue ts = formap.get(list.get(bstart));
			return revmap_rr.get(ts.next);
		}else{
			for(int ii = 0;ii < revlist.size();ii++){
				if(revmap.get(revlist.get(ii)) != null){
					return revlist.get(ii);
				}
			}
		}
		return null;
	}
	
	
	
	/**
	 * 相手側の前の残基を得る
	 * @return 
	 */
	public RResidue getCompPrev(){
			
		int bstart = -1;
		HashMap<RResidue,RResidue> formap = null;
		HashMap<RResidue,RResidue> revmap = null;
		ArrayList<RResidue> revlist = null;
		HashMap<PDBResidue,RResidue> formap_rr = null;
		HashMap<PDBResidue,RResidue> revmap_rr = null;
		if(list == ali.target){
			formap = ali.target_template;
			revmap = ali.template_target;
			revlist = ali.template;
			formap_rr = ali.target_rr;
			revmap_rr = ali.template_rr;
		}else{
			formap = ali.template_target;
			revmap= ali.target_template;
			revlist = ali.target;
			formap_rr = ali.template_rr;
			revmap_rr = ali.target_rr;
		}
			
		for(int ii = segstart;ii <= segend;ii++){
			if(formap.get(list.get(ii)) != null){
				bstart = ii;
				break;
			}
		}
		if(bstart == -1){
			for(int ii = segend;ii < list.size();ii++){
				if(formap.get(list.get(ii)) != null){
					bstart = ii;
					break;
				}
			}
		}
		if(bstart > -1){
			RResidue ts = formap.get(list.get(bstart));
			return revmap_rr.get(ts.prev);
		}else{
			for(int ii =revlist.size()-1;ii > -1 ;ii--){
				if(revmap.get(revlist.get(ii)) != null){
					return revlist.get(ii);
				}
			}
		}
		return null;
	}
	
	
	
}

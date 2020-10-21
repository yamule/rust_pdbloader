/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
import static pepbuilderj.ChainBuilder.changeToFloating;
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class Threading {
	FuzzyTreeAlign fz = new FuzzyTreeAlign();
	TemplateBaseModeller tbm = new TemplateBaseModeller();
	
	ChainBuilder cb_frag = new ChainBuilder();
	/**
	 * Scoring alignment を使用してモデリングする
	 * threshold 未満しか並ばなかった場合 null
	 * @param target
	 * @param template
	 * @param threshold
	 * @return 
	 */
	public ThreadingResult threading(ArrayList<FloatingResidue> target
			,ArrayList<PDBResidue> template,int threshold){
		StringBuffer sb = new StringBuffer();
		
		for(PDBResidue p:target){
			sb.append(PepProcess.getAALetter(p.getName()));
		}
		
		SWResult res = fz.plainAlign(sb.toString(),template);
		ArrayList<String> tar = new ArrayList<>();
		ArrayList<String> tem = new ArrayList<>();
		
		for(int ii = 0;ii < res.qseq.size();ii++){
			char r = res.qseq.get(ii);
			char t = res.sseq.get(ii);
			tar.add(String.valueOf(r));
			tem.add(String.valueOf(t));
		}
		TBMResult tres = tbm.model(tar,tem,template);
		tres.score = res.score;
		return tbmToThreading(tres,target,threshold);
	}
	public ThreadingResult tbmToThreading(TBMResult  tres,ArrayList<FloatingResidue> target,int threshold){
		int start = -1;
		int end = -1;
		for(int ii = 0;ii < target.size();ii++){
			tres.residues.get(ii).setResidueNumber(target.get(ii).getResidueNumber());
			if(!tres.residues.get(ii).getName().equals(target.get(ii).getName())){
				throw new RuntimeException("error in code.");
			}
			if(start == -1 && tres.mapped.get(ii)){
				start = ii;
			}
			 if(tres.mapped.get(ii)){
				end = ii;
			}
		}
		if(tres.mapped.size() < threshold){
			return null;
		}
		ThreadingResult ts = new ThreadingResult();
		ts.score = tres.score;
		for(int ii = start;ii <= end;ii++){
			FloatingResidue f = new FloatingResidue(tres.residues.get(ii),null,null);
			ts.fragment.add(f);
			f.calcPseudoNeighbour(cb_frag.backbones.getRandomly(f.getName()));
			if(tres.mapped.get(ii)){
				ts.modelled.add(f);
				ts.used.put(target.get(ii),f);
				//System.out.println(f.getName()+"\t"+target.get(ii).getName());
			}
		}
		for(int ii = 0;ii < ts.fragment.size();ii++){
			if(ii > 0){
				ts.fragment.get(ii).setPrevFloating(ts.fragment.get(ii-1));
			}
			if(ii < ts.fragment.size()-1){
				ts.fragment.get(ii).setNextFloating(ts.fragment.get(ii+1));
			}
		}
		if(target.get(start).prev != null){
			ts.fragment.get(0).setPrevFloating(target.get(start).prev);
		}
		if(target.get(end).next != null){
			ts.fragment.get(ts.fragment.size()-1).setNextFloating(target.get(end).next);
		}
		
		for(int ii = 0;ii < ts.fragment.size();ii++){
			FloatingResidue f = ts.fragment.get(ii);
			f.changeSideChain(cb_frag.sidechains.getNext(f));
		}
		return ts;
	}
	public void refineLoop(ThreadingResult ts,int ite){
		int mapped1 = -1;
		ArrayList<FloatingResidue> nocrash = new ArrayList<>();
		
		for(int rr = 0;rr < ts.fragment.size();rr++){
			FloatingResidue f = ts.fragment.get(rr);
			if(ts.modelled.contains(f)){
				nocrash.add(f);
			}
			if(f.backboneIndex == -1){
				f.calcPseudoNeighbour(cb_frag.backbones.getRandomly(f.getName()));
			}
		}
		for(int rr = 0;rr < ts.fragment.size();rr++){
			FloatingResidue f = ts.fragment.get(rr);
			if(!ts.modelled.contains(f)){
				ArrayList<FloatingResidue> unmapped = new ArrayList<>();
				unmapped.add(f);
				for(int zz = rr+1;zz < ts.fragment.size();zz++){
					FloatingResidue fz = ts.fragment.get(zz);
					if(!ts.modelled.contains(fz)){
						unmapped.add(fz);
						rr = zz;
					}else{
						break;
					}
				}
				int mapped2 = -1;
				if(rr < ts.fragment.size()-1){
					mapped2 = rr+1;
				}
				cb_frag.prepare(nocrash);
				if(mapped1 > -1){
					for(int ll = 0;ll < unmapped.size();ll++){
						if(mapped2 > -1 && unmapped.size()/2 < ll){
							break;
						}
						cb_frag.elongTerm(unmapped.get(ll), true, 10, true);
					}
				}
				if(mapped2 > -1){
					for(int ll = unmapped.size()-1;ll >= 0;ll--){
						if(mapped1 > -1 && unmapped.size()/2-1 > ll){
							break;
						}
						cb_frag.elongTerm(unmapped.get(ll), false, 10, true);
					}
				}

			}else{
				mapped1 = rr;
			}
		}
		ArrayList<FloatingResidue> unmapped = new ArrayList<>();
		HashSet<FloatingResidue> chk = new HashSet<>();
		//ArrayList<FloatingResidue> unmapped = new ArrayList<>();
		for(int ii = 0;ii < ts.fragment.size();ii++){
			FloatingResidue f = ts.fragment.get(ii);
			if(!ts.modelled.contains(f)){
				if(ii > 0 && !chk.contains(f.prev)){
					unmapped.add(f.prev);
					chk.add(f.prev);
				}
				if(!chk.contains(f)){
					unmapped.add(f);
					chk.add(f);
				}
				if(ii > 0 && !chk.contains(f.next)){
					unmapped.add(f.next);
					chk.add(f.next);
				}
			}
		}
		
		
		//最初と最後は別にあるので一時的によける
		FloatingResidue pprev = ts.fragment.get(0).prev;
		FloatingResidue pnex = ts.fragment.get(ts.fragment.size()-1).next;
		ts.fragment.get(0).setPrevFloating(null);
		ts.fragment.get(ts.fragment.size()-1).setNextFloating(null);
		
		cb_frag.prepare(ts.fragment);
		ArrayList<PDBResidue> modelled_r = new ArrayList<>();
		modelled_r.addAll(ts.fragment);
		
		ArrayList<ArrayList<FloatingResidue>> regionsp = getUnmappedRegions(ts.fragment,new HashSet<FloatingResidue> (unmapped));
		for(ArrayList<FloatingResidue> fal:regionsp){
			ArrayList<PDBResidue> pr = new ArrayList<>();
			pr.addAll(fal);
			changeBackBone(cb_frag,ts.fragment,modelled_r,fal,pr,true,10);
		}
		//ChainBuilder.saveChains(ts.fragment,"test.b.tmp.pdb");
		HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>(unmapped);
		for(int ii = 0;ii < ts.fragment.size();ii++){
			if(ts.fragment.get(ii).next != null){
				if(ts.fragment.get(ii).getC().distance(ts.fragment.get(ii).next.getN()) > 1.4){
					unmapped_withneighbour.addAll(cb_frag.getGapFillingResidues(ts.fragment.get(ii)));
				}
			}
		}
		ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(ts.fragment,new HashSet<FloatingResidue> (unmapped_withneighbour));
		for(ArrayList<FloatingResidue> fal:regions){
			ArrayList<PDBResidue> pr = new ArrayList<>();
			pr.addAll(fal);
			changeBackBone(cb_frag,ts.fragment,modelled_r,fal,pr,true,ite);
		}
		
		ts.fragment.get(0).setPrevFloating(pprev);
		ts.fragment.get(ts.fragment.size()-1).setNextFloating(pnex);

	}
	
	
	
	public ArrayList<ThreadingResult> threadingFromPDB(ArrayList<FloatingResidue> fr
			,String pdbfile
	,boolean refineloop){
		PDBData p = PDBData.loadPDBFile(pdbfile);
		ArrayList<ThreadingResult> ret = new ArrayList<>();
		for(String c:p.chains.keySet()){
			ArrayList<PDBResidue> filted = PepProcess.makeFilteredAA(p.chains.get(c).residues, true);
			if(filted.size() < 10){
				continue;
			}
			ThreadingResult res = threading(fr,filted,5);
			res.templateFile = pdbfile;
			res.templateChain = c;
			res.templateResidueNum = p.chains.get(c).residues.size();
			if(res.templateResidueNum == 0){
			}else{
				//あまり良くなさそうだった。
				//res.score /= p.chains.get(c).residues.size();
			}
			if(res != null){
				if(refineloop){
					refineLoop(res,10);
				}
				ret.add(res);
			}
		}
		return ret;
	}
	public ArrayList<ThreadingResult> tbmModel(String infile,ArrayList<FloatingResidue> fr){
		String targetfas = infile;
		ArrayList<ThreadingResult> ret = new ArrayList<>();
		StringBuffer sb = new StringBuffer();
		for(FloatingResidue p:fr){
			sb.append(PepProcess.getAALetter(p.getName()));
		}
		
		ArrayList<TBMEntry> ls = TBMEntry.loadFasta(targetfas);
		HashMap<Integer,TBMEntry> targets = new HashMap<>();
		ArrayList<TBMEntry> templates = new ArrayList<>();
		for(TBMEntry t:ls){
			if(t.name.equals("target")){
				targets.put(t.id, t);
			}else{
				templates.add(t);
			}
		}
		HashMap<Integer,PDBResidue> dmap = new HashMap<>();
		TemplateBaseModeller tbm = new TemplateBaseModeller();
		ArrayList<String> messages = new ArrayList<>();
		ArrayList<FloatingResidue> allres = new ArrayList<>();
		for(TBMEntry t:templates){
			PDBData pdb1 = PDBData.loadPDBFile(t.file);
			TBMEntry target = targets.get(t.targetId);
			TBMResult res = tbm.model(target.sequence
					,t.sequence,pdb1.chains.get(t.chainName));
			res.score = 10000;
			PDBChain pc = new PDBChain(target.chainName);
			ArrayList<FloatingResidue> ffr = new ArrayList<>();
			for(int ii = 0;ii < res.residues.size();ii++){
				ffr.add(fr.get(target.start+ii-1));
				
			}
			ret.add(tbmToThreading(res,ffr,-100));
		}
		
		
		return ret;
	}
	
	public ArrayList<ThreadingResult> fileToThreading(String infile,ArrayList<FloatingResidue> fr){
		ArrayList<ThreadingResult> ret = new ArrayList<>();
		
		PDBData pdb = PDBData.loadPDBFile(infile);
		for(String cname:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(cname);
			TBMResult tres = new TBMResult();
			
//	ArrayList<PDBResidue> residues = new ArrayList<>();
//	ArrayList<Boolean> mapped = new ArrayList<>();
//	ArrayList<PDBResidue> template_residues = new ArrayList<>();
//	HashMap<PDBResidue,PDBResidue> target_template_map = new HashMap<>();
//	HashMap<PDBResidue,PDBResidue> template_target_map = new HashMap<>();
			ArrayList<FloatingResidue> res = changeToFloating(cc.residues);
			for(FloatingResidue t:res){
				tres.residues.add(t);
				tres.mapped.add(true);
			}
			PDBChain pc = new PDBChain(cname);
			ArrayList<FloatingResidue> ffr = new ArrayList<>();
			for(int ii = 0;ii < tres.residues.size();ii++){
				ffr.add(fr.get(tres.residues.get(ii).getResidueNumber()-1));
			}
			ret.add(tbmToThreading(tres,ffr,-100));
		}
		
		
		return ret;
	}
	
	
	
	/**
	 * 前の C 原子が X Ang の範囲に来るよう平衡移動する
	 * @param fl 
	 */
	public void fitToPrev(ArrayList<ArrayList<FloatingResidue>> fl){
		double cndist = cb_frag.distPenalty.backbone_average[AtomDistancePenalty.BACKBONE_PREVC_N];
		for(int ii = 0;ii < fl.size();ii++){
			FloatingResidue fr = fl.get(ii).get(0);
			if(fr.prev != null){
				double ddist = fr.prev.getC().distance(fr.getN());
				if(Math.abs(ddist
						- cndist) > 0.2){
					
					ArrayList<Point3D> p3 = new ArrayList<>();
					for(FloatingResidue ff:fl.get(ii)){
						p3.addAll(ff.atoms_loc_all);
						
					}
					
					
					Point3D vvec = new Point3D(fr.prev.getC().loc);
					vvec.x -= fr.getN().loc.x;
					vvec.y -= fr.getN().loc.y;
					vvec.z -= fr.getN().loc.z;
					vvec.standarize();
					vvec.x *= (ddist-cndist);
					vvec.y *= (ddist-cndist);
					vvec.z *= (ddist-cndist);
					RandomDocker.moves(p3, vvec);
					//if(ddist < cndist){
					//	System.out.println(fr.prev.getC().distance(fr.getN()));
					//}
				}
			}
		}
	}
	
	public static  double sum(double[]s){
		double ret = 0;
		for(double dd:s){
			ret += dd;
		}
		return ret;
	}
	
	public static void main(String args[]){
		for(int si = 0;si < 1;si++){
			String seq = "";
			String filebase = "";
			if(si == 0){
				seq = "MAKFACKCGYVINLIASPGGDEWRLIPEKTLEDIVDLLDGGEAVDGERFYETLRGKEITVYRCPSCGRLHLEEAGRNKFVTYVKECGEL";
				filebase ="T1015s1";
			}else{
				seq = "MAKFACKCGYVINLIASPGGDEWRLIPEKTLEDIVDLLDGGEAVDGERFYETLRGKEITVYRCPSCGRLHLEEAGRNKFVTYVKECGEL";
				filebase ="T1015s1";
			}
			
			
			String[] ss = seq.split("");
			ArrayList<String> cc = new ArrayList<>();
			Threading th = new Threading();
			for(String s:ss){
				if(s.length() > 0){
					cc.add(s);
				}
			}
			ArrayList<FloatingResidue> fr = th.cb_frag.buildDummy(cc);
			ArrayList<ThreadingResult> th_all = new ArrayList<>();
					String[] tbmresults  = {
};//TBM の結果がある時はここに含める
			
			ArrayList<ThreadingResult> tbmresult = new ArrayList<>();
			for(String t:tbmresults){
				ArrayList<ThreadingResult> zres = th.fileToThreading(t,fr);
				ThreadingResult tzres = zres.get(0);
				tbmresult.add(tzres);
			}
			
			
			Pattern pdbfilepat = Pattern.compile("\\.(ent|pdb)");
			ArrayList<String> pdblist = new ArrayList<>();
			String sourcedirpath = "C:\\dummy\\vbox_share\\bioo\\database\\scop40\\pdbstyle-2.07";
			File sourcedir = new File(sourcedirpath);
			File[] pdbfiles = sourcedir.listFiles();

			for(File pp:pdbfiles){
				if(pdbfilepat.matcher(pp.getPath()).find()){
					pdblist.add(pp.getPath());
				}
				if(pp.isDirectory()){//pdb divided と同じように二階層チェックする
					File[] pdbfiles2 = pp.listFiles();
					
					for(File pdd:pdbfiles2){
						if(pdbfilepat.matcher(pdd.getPath()).find()){
							pdblist.add(pdd.getPath());
						}
					}
				}
			}
			int pcou = 0;
			//threading で構造取ってくる
			for(String fname:pdblist){
				pcou++;
				System.out.println(pcou+"/"+pdblist.size());
				ArrayList<ThreadingResult> ttres = th.threadingFromPDB(fr,fname,false);
				if(ttres != null){
					th_all.addAll(ttres);
				}
				//if(Math.random() > 0.995){
				//	break;
				//}
			}
			Collections.sort(th_all,new TResComparator());
			Collections.reverse(th_all);
			ChainBuilder cb2 = new ChainBuilder();
	
			FuzzyDecisionTreeScoring_generator sccc 
					= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180724());
			cb2.scoring.clear();
			cb2.scoring.add(sccc);
			if(true){
				//スコアの高い奴をモデリングする
				for(int ii = 0;ii < Math.min(300,th_all.size());ii++){
					ArrayList<String> messages = new ArrayList<>();
					ThreadingResult tr = th_all.get(ii);
					messages.add("#score:"+tr.score);
					messages.add("#file:"+tr.templateFile);
					messages.add("#chain:"+tr.templateChain);
					messages.add("#residue stats");
					HashSet<FloatingResidue> unmapped = new HashSet<>();
					cb2.prepare(tr.fragment);
					for(int j = 0;j < tr.fragment.size();j++){
						FloatingResidue r = tr.fragment.get(j);
						
						r.calcPseudoNeighbour(cb2.backbones.getRandomly(r.getName()));
						if(tr.modelled.contains(r)){
							messages.add("chain_name:"+r.parent.getName()+"\t"+"residue_name:"+r.getName()
									+"\t"+"residue_number:"+r.getResidueNumber()+"\t"+"ok");
						}else{
							unmapped.add(r);
							messages.add("chain_name:"+r.parent.getName()+"\t"+"residue_name:"+r.getName()
									+"\t"+"residue_number:"+r.getResidueNumber()+"\t"+"ng");
						}
					}


					ArrayList<PDBResidue> modelled_r = new ArrayList<>();
					modelled_r.addAll(tr.fragment);
					ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(tr.fragment,new HashSet<FloatingResidue> (unmapped));
					for(ArrayList<FloatingResidue> pd:regions){
						ArrayList<PDBResidue> pd_r = new ArrayList<>();	
						pd_r.addAll(pd);
						changeBackBone(cb2,tr.fragment,modelled_r,pd,pd_r,true,10);
					} 
					HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
					unmapped_withneighbour.addAll(unmapped);

					//fixme 露出残基をとって、露出している方から削っていく
					//ChainBuilder.saveChains(modelled,outfilename+".p.pdb");
					for(int ff = 0;ff< tr.fragment.size();ff++){
						FloatingResidue rf = tr.fragment.get(ff);
						if(rf.next != null && rf.getC().distance(rf.next.getN()) > 1.4){
							unmapped_withneighbour.addAll(cb2.getGapFillingResidues(rf));
						}
					}
					if(unmapped_withneighbour.size() < tr.fragment.size()/3.0){
						regions = getUnmappedRegions(tr.fragment,new HashSet<FloatingResidue> (unmapped_withneighbour));
						for(int tt = 0;tt < 5;tt++){
							for(ArrayList<FloatingResidue> pd:regions){
								ArrayList<PDBResidue> pd_r = new ArrayList<>();	
								pd_r.addAll(pd);
								changeBackBone(cb2,tr.fragment,modelled_r,pd,pd_r,true,10*pd.size());
							}
						}
						
						for(int kk = 0;kk < 2;kk++){
							for(FloatingResidue r:tr.fragment){
								cb2.refiner.maxRotamer(r, -5,0.3, true, cb2.sidechains, false);
							}
						}
						writeToFile(messages,filebase+"."+ii+".th.res");
						ChainBuilder.saveChains(th_all.get(ii).fragment,filebase+"."+ii+".th.pdb");
					}
				}
				//System.exit(0);
			}
			int strcount = 0;
			int strnum = 1000;
			while(th_all.size() > 0){
				strcount++;
				
				if(strcount > strnum){
					break;
				}
				HashMap<FloatingResidue,ThreadingResult> flags = new HashMap<>();
				/*TBM は親構造分かった方がいいので
				for(int rr = 0;rr < tbmresult.size()*2;rr++){
					int t1 = (int)(Math.random()*tbmresult.size());
					int t2 = (int)(Math.random()*tbmresult.size());
					ThreadingResult h1 = tbmresult.get(t1);
					ThreadingResult h2 = tbmresult.get(t2);
					tbmresult.set(t1, h2);
					tbmresult.set(t2, h1);
				}
				
				HashSet<ThreadingResult> plis = new HashSet<>(tbmresult);
				th_all.removeAll(plis);
				for(ThreadingResult tr:tbmresult){
					th_all.add(0,tr);
				}
				*/
				HashSet<ThreadingResult> plis = new HashSet<>(tbmresult);
				th_all.removeAll(plis);
				if(tbmresult.size() > 0){
					th_all.add(0,tbmresult.get((strcount-1)%tbmresult.size()));
					//うーんもう少し考えるかここは
				}
				
				while(th_all.size() > 0){
					ArrayList<FloatingResidue> ffr = new ArrayList<>();
					outer :for(int ii = 0;ii < fr.size();ii++){
						if(!flags.containsKey(fr.get(ii))){
							ArrayList<FloatingResidue> chk = new ArrayList<>();
							for(int jj = ii;jj < fr.size();jj++){
								if(flags.containsKey(fr.get(jj))){
								}else{
									chk.add(fr.get(jj));
								}
							}
							if(chk.size() > 10){
								ffr.addAll(chk);
								break outer;
							}
						}
					}
					if(ffr.size() == 0){
						break;
					}
					boolean unmappedflag = true;
					for(int ll = 0;ll < th_all.size();ll++){
						ThreadingResult ttres = th_all.get(ll);
						
						int flagcount = 0;
						for(int ii = 0;ii < ffr.size();ii++){
							FloatingResidue f = ffr.get(ii);
							if(ttres.used.containsKey(f) && !flags.containsKey(f)){
								flagcount++;
							}
						}
						if(flagcount > 10){
							th.refineLoop(ttres,10);
							for(int ii = 0;ii < ffr.size();ii++){
								FloatingResidue f = ffr.get(ii);
								if(ttres.used.containsKey(f) && !flags.containsKey(f)){
									f.fitAtoms(ttres.used.get(f));
									flags.put(f,ttres);
								}
							}
							
							unmappedflag = false;
							th_all.remove(ll);
							break;
						}
					}
					if(fr.size() -flags.size()< 10 || unmappedflag){
						break;
					}
				}
				ArrayList<ArrayList<FloatingResidue>> groups  = new ArrayList<>();
				//マップされていない残基がある場合にひとかたまりにする。
				//途中にインサーションが入っているドメインには対応しておらず分割されてしまう
				for(int ii = 0;ii < fr.size();ii++){
					if(flags.containsKey(fr.get(ii))){
						ThreadingResult tr = flags.get(fr.get(ii));
						int start = ii;
						int end  = ii;
						ArrayList<FloatingResidue> fal = new ArrayList<>();
						for(int jj = ii;jj < fr.size();jj++){
							if(flags.get(fr.get(jj)) != null && flags.get(fr.get(jj)) != tr){
								ii = jj -1;
								break;
							}
							if(flags.get(fr.get(jj)) == tr){
								end = jj;
								ii = jj;
							}
						}
						for(int jj = start;jj <= end;jj++){
							fal.add(fr.get(jj));
							flags.put(fr.get(jj),tr);
						}
						
						
						
						groups.add(fal);
					}
				}

				//一個だけの残基およびインサートになっている残基を単独のグループとして加える
				for(int ii = 0;ii < fr.size();ii++){
					if(!flags.containsKey(fr.get(ii))){
						ArrayList<FloatingResidue> fal  = new ArrayList<>();
						fal.add(fr.get(ii));
						groups.add(fal);
					}else if(!flags.get(fr.get(ii)).used.containsKey(fr.get(ii))){
						ArrayList<FloatingResidue> fal  = new ArrayList<>();
						fal.add(fr.get(ii));
						groups.add(fal);
					}
				}


				//N 末から並ぶようにソートする
				for(ArrayList<FloatingResidue> fl:groups){
					Collections.sort(fl,new FRComparator());
				}

				ArrayList<VSorter> sorter = new ArrayList<>();
				for(int ii = 0;ii < groups.size();ii++){
					sorter.add(new VSorter(groups.get(ii).get(0).getResidueNumber(),ii));
				}
				Collections.sort(sorter,new VComparator());
				ArrayList<ArrayList<FloatingResidue>> tmpp = new ArrayList<>();
				for(int ii = 0;ii < sorter.size();ii++){
					tmpp.add(groups.get(sorter.get(ii).index));
				}
				groups = tmpp;
				//--------------- ソート終了
				th.fitToPrev(groups);
				for(ArrayList<FloatingResidue> fl:groups){
					if(fl.size() > 2){
						th.cb_frag.prepare(fl);
						//th.cb_frag.fillGapsAll(fl, 10);
						ArrayList<FloatingResidue> modelled = fl;
						HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
						for(int ii = 0;ii < modelled.size()-1;ii++){
							if(modelled.get(ii).getC().distance(modelled.get(ii+1).getN()) > 1.4){
								System.out.println(";;"+modelled.get(ii).getRepresentativeCode());
								unmapped_withneighbour.addAll(th.cb_frag.getGapFillingResidues(modelled.get(ii)));
							}

						}
						ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(modelled,new HashSet<FloatingResidue> (unmapped_withneighbour));
						ArrayList<PDBResidue> modelled_r = new ArrayList<>();
						modelled_r.addAll(modelled);
						for(ArrayList<FloatingResidue> pd:regions){
							ArrayList<PDBResidue> pd_r = new ArrayList<>();	
							pd_r.addAll(pd);
							System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
							changeBackBone(th.cb_frag,modelled,modelled_r,pd,pd_r,true,50*pd.size());
						}
					}
				}
				for(int ii = 0;ii < fr.size();ii++){
					if(ii > 0){
						fr.get(ii).setPrevFloating(fr.get(ii-1));
					}
					if(ii <  fr.size()-1){
						fr.get(ii).setNextFloating(fr.get(ii+1));
					}
				}
				
				ArrayList<ArrayList<Point3D>> group_points = new ArrayList<>();
				for(ArrayList<FloatingResidue> fl:groups){
					ArrayList<Point3D> p = new ArrayList<>();
					for(FloatingResidue f:fl){
						p.addAll(f.atoms_loc_all);
					}
					group_points.add(p);
				}
				ChainBuilder cb_all = new ChainBuilder();
				FuzzyDecisionTreeScoring_generator j 
						= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180724());
				cb_all.scoring.clear();
				cb_all.scoring.add(j);

				cb_all.prepare(fr);
				ArrayList<PDBResidue> fr_ = new ArrayList<>();
				fr_.addAll(fr);
				th.fitToPrev(groups);
				for(FloatingResidue f:fr){
					f.saveLoc();
				}
				double score_prev = sum(cb_all.calcChainScore(fr_));
				for(int ii = 0;ii < 5;ii++){
					for(ArrayList<Point3D> pp :group_points){
						for(int jj = 0;jj < 10;jj++){
							double xfactor = 1.0;
							if(jj%2 == 0){
								xfactor = 0.01;
							}
							double rotbase = Math.PI*xfactor;
							double movbase = 2.0*xfactor;
							RandomDocker.rotates(pp
								,pp.get(0)
								, Math.random()*rotbase-rotbase/2
								, Math.random()*rotbase-rotbase/2
								, Math.random()*rotbase-rotbase/2);
							Point3D vec = new Point3D(
									Math.random()*movbase-movbase/2
									,Math.random()*movbase-movbase/2
									,Math.random()*movbase-movbase/2);
							RandomDocker.moves(pp, vec);
							th.fitToPrev(groups);
							
							double sc  = sum(cb_all.calcChainScore(fr_));
							if(sc > score_prev){
								for(FloatingResidue f:fr){
									f.saveLoc();
								}
								score_prev = sc;
							}else{
								for(FloatingResidue f:fr){
									f.restoreLoc();
								}
							}
						}
					}
				//	ChainBuilder.saveChains(fr,"test."+strcount+".tmp.pdb");
				//	System.out.println(score_prev);
				}
				
				
				ChainBuilder cbb = new ChainBuilder();
				cbb.prepare(fr);
				HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
				ArrayList<PDBResidue> frr = new ArrayList<>();
				frr.addAll(fr);
				for(FloatingResidue rf:fr){
					rf.calcPseudoNeighbour(cbb.backbones.getRandomly(rf.getName()));
					
					if(rf.next != null && rf.getC().distance(rf.next.getN()) > 1.4){
						unmapped_withneighbour.addAll(cbb.getGapFillingResidues(rf));
						System.out.println(rf.getRepresentativeCode());
					}
				}
				ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(fr,new HashSet<FloatingResidue> (unmapped_withneighbour));
				for(int tt = 0;tt < 5;tt++){
					for(ArrayList<FloatingResidue> pd:regions){
						ArrayList<PDBResidue> pd_r = new ArrayList<>();	
						pd_r.addAll(pd);
						changeBackBone(cbb,fr,frr,pd,pd_r,true,10*pd.size());
					}
				}
				
				for(int kk = 0;kk < 2;kk++){
					for(FloatingResidue r:fr){
						cbb.refiner.maxRotamer(r, -5,0.1, true, cbb.sidechains, false);
					}
				}
				
				/*//デバッグ
				for(int tt = 0;tt < 1;tt++){
					for(ArrayList<FloatingResidue> pd:regions){
						ArrayList<PDBResidue> pd_r = new ArrayList<>();	
						pd_r.addAll(pd);
						changeBackBone(cbb,fr,frr,pd,pd_r,true,5*pd.size());
					}
				}
				*/
				
				System.out.println("======");
				for(FloatingResidue rf:fr){
				
					if(rf.next != null && rf.getC().distance(rf.next.getN()) > 1.4){
						System.out.println(rf.getRepresentativeCode());
					}
				}
				ChainBuilder.saveChains(fr,"th."+filebase+"."+strcount+".pdb");
				System.out.println();
				System.out.println("##"+score_prev+"\t"+"th."+filebase+"."+strcount+".pdb");
			}

		}
	}
	
}


class ThreadingResult{
	ArrayList<FloatingResidue> fragment = new ArrayList<>();//モデリングされた実際の残基。元の構造とは独立している。
	HashMap<PDBResidue,FloatingResidue> used = new HashMap<>();//元の構造の残基でこの Threading によって構築されたもの
	HashSet<PDBResidue> modelled = new HashSet<>();//テンプレートにマップされた残基。これ以外はループモデリングされた
	String templateFile = "";
	String templateChain = "";
	int templateResidueNum = 1;
	double score = 0;
	ThreadingResult(){
	}
}

  class FRComparator implements Comparator<FloatingResidue>{
	@SuppressWarnings("unchecked")
	public int compare(FloatingResidue arg1, FloatingResidue arg2){
		
		if(arg1.getResidueNumber() < arg2.getResidueNumber() ){
			return -1;
		}
		if(arg1.getResidueNumber() == arg2.getResidueNumber() ){
			return 0;
		}
			return 1;
	}
}

  class TResComparator implements Comparator<ThreadingResult>{
	@SuppressWarnings("unchecked")
	public int compare(ThreadingResult arg1, ThreadingResult arg2){
		
		if(arg1.score < arg2.score){
			return -1;
		}
		if(arg1.score == arg2.score ){
			return 0;
		}
			return 1;
	}
}
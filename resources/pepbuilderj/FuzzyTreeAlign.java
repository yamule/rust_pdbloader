/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;
import static pepbuilderj.FuzzyDecisionTreeScoring_generator.calcEachResidueProbabilities;

/**
 *
 * @author kimidori
 */
public class FuzzyTreeAlign {
	ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
	public static BLOSUM62_X0 blosum = new BLOSUM62_X0();
	int iterationNum = 20;
	int fileSearchDepth = 1;
	int  splitSeedLength = 10;
	int splitGapLength = 8;
	
	FuzzyTreeAlign(){
		this(null);
	}
	FuzzyTreeAlign(FuzzyDecisionTreeScoring_generator gen){
		if(gen == null){
			FuzzyDecisionTreeScoring_generator j
			//		= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
			//		= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180511());
					= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180724());
				
			scoring.add(j);
			
		}else{
			scoring.add(gen);
		}
		//for(FuzzyDecisionTreeScoring_generator j:scoring){
		//	j.rescale(4.0/scoring.size());
		//}
	}
	
	
	static char aachars[] = "ARNDCQEGHILKMFPSTWYVX".toCharArray();
	public static void fillWithBLOSUM(PSSMScoring sc,int rownum){
		fillWithBLOSUM(sc,sc.pssm.seq.get(rownum),rownum);
	}
	
	public static void fillWithBLOSUM(PSSMScoring sc,char tchar,int rownum){
		for(char a:aachars){
			if(sc.pssm.char_index_map[a] > -1){
				sc.pssm.scores.get(rownum).set(sc.pssm.char_index_map[a]
						,(double)blosum.scoreMat[tchar][a]);
				
			}
		}
	}
	
	public static void fillWithZERO(PSSMScoring sc,int rownum){
		for(char a:aachars){
			if(sc.pssm.char_index_map[a] > -1){
				sc.pssm.scores.get(rownum).set(sc.pssm.char_index_map[a]
						,0.0);
			}
		}
	}
	
	public void fillNegative(PSSMScoring sc){
		for(int ii = 0;ii < sc.pssm.seq.size()-1;ii++){
			if(sc.getScoreAt(ii,sc.pssm.seq.get(ii)) <= 0){
				fillWithBLOSUM(sc,ii);
			}else{
			}
		}
		System.out.println("");
	}
	
	public PSSMScoring calcPSSM(ArrayList<PDBResidue> target,ArrayList<PDBResidue> allresidue){
		
		ArrayList<HashMap<String,Double>> lis = calcEachResidueProbabilities(
				target,	allresidue,scoring);
		if(lis.size() == 0){
			PSSMData ret = new PSSMData();
			for(int ii = 0;ii < target.size();ii++){
				HashMap<Character,Double> hcc = new HashMap<>();
				for(char c:PSSMCalc.aa_type){
					hcc.put(c,0.0);
				}
				if(PDBResidue.aaMap.get(target.get(ii).getName()) == null){
					System.out.println(target.get(ii).getName());
				}
				ret.addResidue(PDBResidue.aaMap.get(target.get(ii).getName()).charAt(0), hcc);

			}
			PSSMScoring sret = PSSMScoring.prepare(ret);
			return sret;	
		}
		double[][] freq = new double[lis.size()][20];
		for(int ii = 0;ii < lis.size();ii++){
			HashMap<String,Double> d = lis.get(ii);
			for(char cc:PSSMCalc.aa_type){
				freq[ii][PSSMCalc.aa_to_index[cc]] = 0;
			}
			double sum = 0;
			for(String s:d.keySet()){
				int ic = PSSMCalc.aa_to_index[PDBResidue.aaMap.get(s).charAt(0)];
				if(ic > -1){
					freq[ii][ic] = d.get(s);
					sum += freq[ii][ic];
				}
			}
			for(int jj = 0;jj < 20;jj++){
				freq[ii][jj] /= sum;
			}
		}
		int[][] score = PSSMCalc.convertFrequencyToPSSM(freq,PSSMCalc.aa_to_index);
		
		PSSMData ret = new PSSMData();
		for(int ii = 0;ii < target.size();ii++){
			HashMap<String,Double> hss = lis.get(ii);
			HashMap<Character,Double> hcc = new HashMap<>();
			for(char c:PSSMCalc.aa_type){
				hcc.put(c,(double)score[ii][PSSMCalc.aa_to_index[c]]);
			}
			if(PDBResidue.aaMap.get(target.get(ii).getName()) == null){
				System.out.println(target.get(ii).getName());
			}
			ret.addResidue(PDBResidue.aaMap.get(target.get(ii).getName()).charAt(0), hcc);
			
		}
		PSSMScoring sret = PSSMScoring.prepare(ret);
		return sret;
	}
	public  static ArrayList<String> getAllFiles(String dirpath,int depthstop){
		ArrayList<String> ret = new ArrayList<>();
		File dir = new File(dirpath);
		File[] flis = dir.listFiles();
		Pattern filepat = Pattern.compile("\\.(ent|pdb)[0-9]*$");
		Pattern ignorepat = Pattern.compile("^\\.");
		for(File f:flis){
			
			if(ignorepat.matcher(f.getName()).find()){
				continue;
			}
			if(filepat.matcher(f.getName()).find()){
				ret.add(f.getAbsolutePath());
			}
			if(f.isDirectory()){
				if(depthstop > 0){
					ret.addAll(getAllFiles(f.getAbsolutePath(),depthstop-1));
				}
			}
		}
		return ret;
	}
	
	public void processAllPDBInDir(String fastafilepath,String dirpath){
		Pattern filepat = Pattern.compile("\\.(ent|pdb)[0-9]*$");
		ArrayList<Sequence> seqs = Sequence.loadFasta(fastafilepath);
		if(seqs.size() > 1){
			System.err.println("Multi fasta is not supported.");
			throw new RuntimeException("Multi fasta is not supported.");
		}
		ArrayList<String> fpath = getAllFiles(dirpath,this.fileSearchDepth);
		ArrayList<ArrayList<String>> al = new ArrayList<>();
		for(String s:fpath){
			ArrayList<String> rr = new ArrayList<>();
			rr.add(s);
			al.add(rr);
		}
		processAllPDB(seqs.get(0),al);
	}
	
	public static ArrayList<RegionSplit> getAlignedRegion(ArrayList<Character> q,ArrayList<Character> s
	,int gaplen,int seedlen){
		int[] flag = new int[q.size()];
		int[] qpos = new int[q.size()];
		int[] spos = new int[q.size()];
		
		int qp = 0;
		int sp = 0;
		
		for(int ii = 0;ii < q.size();ii++){
			qpos[ii] = qp;
			spos[ii] = sp;
			if(s.get(ii) != '-' && q.get(ii) != '-'){
				flag[ii] = 1;
			}else{
				flag[ii] = 0;
			}
			if(s.get(ii) != '-'){
				sp++;
			}
			if(q.get(ii) != '-'){
				qp++;
			}
		}
		
		ArrayList<Integer> apos = new ArrayList<>();
		ArrayList<RegionSplit> ret = new ArrayList<>();
		int mstart = -1;
		int mend = -1;
		int gstart = -1;
		
		for(int ii = 0;ii < q.size();ii++){
			if(flag[ii] == 1){
				if(mstart < 0){
					mstart = ii;
				}
				mend = ii;
				gstart = -1;
			}
			if(flag[ii] == 0){
				if(gstart < 0){
					gstart = ii;
				}
				if(ii-gstart+1 >= gaplen){
					if(mend-mstart+1 >= seedlen){
						apos.add(mstart);
						apos.add(mend+1);
						
						RegionSplit rs = new RegionSplit();
						rs.alignmentstart = mstart;
						rs.alignmentend = mend;
						rs.qstart = qpos[mstart];
						rs.qend = qpos[mend];
						rs.sstart = spos[mstart];
						rs.send = spos[mend];
						
						
						
						ret.add(rs);
					}
					mstart = -1;
					mend = -1;
				}
			}
		}
		
		if(mend-mstart+1 >= seedlen){
			apos.add(mstart);
			apos.add(mend+1);
			RegionSplit rs = new RegionSplit();
			rs.alignmentstart = mstart;
			rs.alignmentend = mend;
			rs.qstart = qpos[mstart];
			rs.qend = qpos[mend];
			rs.sstart = spos[mstart];
			rs.send = spos[mend];
			ret.add(rs);
			mstart = -1;
			mend = -1;
		}
		return ret;
	}
	public ArrayList<ChainIntermediateState> joinChains(ArrayList<ChainIntermediateState> r){
		HashMap<String,ArrayList<ChainIntermediateState>> hs = new HashMap<>();
		for(ChainIntermediateState cc:r){
			if(!hs.containsKey(cc.name)){
				hs.put(cc.name,new ArrayList<ChainIntermediateState>());
			}
			hs.get(cc.name).add(cc);
		}
		ArrayList<ChainIntermediateState> ret = new ArrayList<>();
		for(String n:hs.keySet()){
			Collections.sort(hs.get(n),new IntermediateStateComparator());
			ret.add(ChainIntermediateState.merge(hs.get(n)));
		}
		return ret;
		
		
	}
	public void freezeHighscoreSegment(ArrayList<ChainIntermediateState> cres){
		if(cres.size() == 0){
			return;
		}
		ChainIntermediateState maxchain = cres.get(0);
		double maxalign = SmithWaterman.calcPSSMScore_nogap(cres.get(0).target_aligned,cres.get(0).template_aligned,cres.get(0).templatePSSM);
		for(int ii = 0;ii < cres.size();ii++){
			double sc = SmithWaterman.calcPSSMScore_nogap(cres.get(ii).target_aligned,cres.get(ii).template_aligned,cres.get(ii).templatePSSM);
			if(maxalign < sc){
				maxchain = cres.get(ii);
				maxalign = sc;
			}
		}
		maxchain.freeze(true);
	}
	
	public ArrayList<ChainIntermediateState> splitWithAlignment(SWResult res
			,ChainIntermediateState cs,
			boolean keeplongest){
		ArrayList<RegionSplit> spp = getAlignedRegion(res.qseq,res.sseq,this.splitGapLength,this.splitSeedLength);
		
		if(keeplongest && spp.size() > 0){
			RegionSplit maxchain = spp.get(0);
			int maxalign = maxchain.countMatches();
			for(RegionSplit rs:spp){
				if(maxalign < rs.countMatches()){
					maxchain = rs;
					maxalign = rs.countMatches();
				}
			}
			spp.clear();
			spp.add(maxchain);
		}
		
		int offset = 0;
		int len = res.qseq.size();
		int qlen = 0;
		int slen = 0;
		for(int ii = 0;ii < res.qseq.size();ii++){
			if(res.qseq.get(ii) != '-'){
				qlen++;
			}
			if(res.sseq.get(ii) != '-'){
				slen++;
			}
		}
		
		ChainIntermediateState remained = cs;
		
		ArrayList<RegionSplit> tpp = new ArrayList<>();
		for(int ii = 0;ii < spp.size();ii++){
			RegionSplit r = spp.get(ii);
			int qnext = -1000;
			int snext = -1000;
			int qlast = -1000;
			int slast = -1000;
			
			int alast = -1000;
			int anext = -1000;
			if(ii < spp.size()-1){
				qnext = spp.get(ii+1).qstart;
				snext = spp.get(ii+1).sstart;
				anext = spp.get(ii+1).alignmentstart;
			}else{
				qnext = qlen;
				snext = slen;
				anext = res.qseq.size();
			}
			
			if(ii > 0){
				qlast = spp.get(ii-1).qend;
				slast = spp.get(ii-1).send;
				alast = spp.get(ii-1).alignmentend;
			}else{
				qlast = -1;
				slast = -1;
				alast = -1;
			}
			if(alast+1 <=  r.alignmentstart-1){
				RegionSplit rx = new RegionSplit();
				rx.qstart = qlast+1;
				rx.qend = r.qstart-1;
				rx.sstart = slast+1;
				rx.send = r.sstart-1;
				rx.alignmentstart = alast+1;
				rx.alignmentend = r.alignmentstart-1;
				
				tpp.add(rx);
			}
			tpp.add(r);
			if(ii == spp.size()-1){
				
				RegionSplit rx = new RegionSplit();
				rx.qstart = r.qend+1;
				rx.qend = qnext-1;
				rx.sstart = r.send+1;
				rx.send = snext-1;
				rx.alignmentstart = r.alignmentend+1;
				rx.alignmentend = anext-1;
				
				tpp.add(rx);
			}
		}
		ArrayList<ChainIntermediateState> ret = new ArrayList<>(); 
		for(int ii = 0;ii < tpp.size()-1;ii++){
			
			ArrayList<ChainIntermediateState> rs = remained.splitAt(tpp.get(ii).qend+1-tpp.get(ii).qstart
					,tpp.get(ii).send+1-tpp.get(ii).sstart);
			ret.add(rs.get(0));
			remained = rs.get(1);
		}
		ret.add(remained);
		if(tpp.size() == 0){
			
			
			RegionSplit rx = new RegionSplit();
			rx.qstart =0;
			rx.qend = qlen-1;
			rx.sstart = 0;
			rx.send = slen-1;
			rx.alignmentstart = 0;
			rx.alignmentend = len-1;
			tpp.add(rx);
		}
		
		for(int ii = 0;ii < ret.size();ii++){
			ret.get(ii).region = tpp.get(ii);
			List<Character> ta = res.qseq.subList(tpp.get(ii).alignmentstart, tpp.get(ii).alignmentend+1);
			List<Character> sa = res.sseq.subList(tpp.get(ii).alignmentstart, tpp.get(ii).alignmentend+1);
			ret.get(ii).setAlignment(new ArrayList<Character>(ta), new ArrayList<Character>(sa));
		}
		return ret;
		
	}
	
	public SWResult plainAlign(String seq_,ArrayList<PDBResidue> template){
		ArrayList<Character> seq = SmithWaterman.toCharList(seq_);
		ArrayList<PDBResidue> ress = TemplateBaseModeller.prepareTemplateResidues(
				PepProcess.makeFilteredAA(template));
		PSSMScoring testpssm = calcPSSM(ress,ress);
		ChainIntermediateState ci = new ChainIntermediateState("dummy",ress,testpssm
				,seq);
		SmithWaterman sw = new SmithWaterman();
		sw.penalO = 10;
		sw.penalE = 1;
		SWResult res = sw.align_with_PSSM(seq_,ci.templatePSSM);
		return res;
	}
	
	public SWResult refineGap(ArrayList<Character> target,ArrayList<Character> template,
			ArrayList<PDBResidue> templatepdb,int realignwindow){
		int shiftsiz = realignwindow;//ギャップ領域の前後何残基を動かすか
		ArrayList<PDBResidue> ress = TemplateBaseModeller.prepareTemplateResidues(
			PepProcess.makeFilteredAA(templatepdb));
		PSSMScoring testpssm = calcPSSM(ress,ress);
		ArrayList<Character> chkseq = new ArrayList<>();
		StringBuffer targetseq = new StringBuffer();
		for(int ii = 0;ii < target.size();ii++){
			if(target.get(ii) != '-'){
				targetseq.append(target.get(ii));
			}
		}
		int cou = 0;
		for(int ii = 0;ii < template.size();ii++){
			if(template.get(ii) != '-'){
				chkseq.add(template.get(ii));
				if(template.get(ii) != testpssm.pssm.seq.get(cou)){
					System.err.println("Sequence discrepancy."+template.get(ii)+";"+testpssm.pssm.seq.get(cou));
				}
				cou++;
			}
		}
		
		SmithWaterman sw = new SmithWaterman();
		sw.penalO = 10;
		sw.penalE = 1;
		sw.prepareFreezing(targetseq.length(),testpssm.pssm.seq.size(),10000);
		
		int tstart = -1;
		int tend = -1;
		
		int qcou = 0;
		int scou = 0;
		for(int ii = 0;ii < template.size();ii++){
			if(target.get(ii) != '-' && template.get(ii) != '-'){
				if(tstart == -1){
					tstart = ii;
				}
				tend = ii;
				sw.freeze(qcou, scou);
			}
			if(target.get(ii) != '-'){
				qcou++;
			}
			if(template.get(ii) != '-'){
				scou++;
			}
		}
		if(tstart > -1){
			
			qcou = 0;
			scou = 0;
			for(int ii = 0;ii < template.size();ii++){
				
				if(target.get(ii) != '-' && template.get(ii) == '-'){
					for(int jj = qcou-shiftsiz;jj <= qcou+shiftsiz;jj++){
					//ギャップが多い場合冗長になってしまうがまあ許容範囲だろう
						int pos = jj;
						if(pos > -1 && pos < targetseq.length()){
							sw.melt(pos);
						}
					}
				}
				
				if(target.get(ii) == '-' && template.get(ii) != '-'){
					for(int jj = qcou-shiftsiz;jj <= qcou+shiftsiz-1;jj++){
						int pos = jj;
						if(pos > -1 && pos < targetseq.length()){
							sw.melt(pos);
						}
					}
				}
				
				if(target.get(ii) != '-'){
					qcou++;
				}
				if(template.get(ii) != '-'){
					scou++;
				}
			}
		}
		
		SWResult res = sw.align_with_PSSM(targetseq.toString(),testpssm);
		return res;
	}
	
	
	
	public void processAllPDB(Sequence seq,ArrayList<ArrayList<String>> pathlist){
		
		ChainBuilder cb = new ChainBuilder();
		for(ArrayList<String> ss:pathlist){
			System.out.println(ss);
			PDBData p  = PDBData.loadPDBFile(ss.get(0));
			ArrayList<String> cnames;
			if(ss.size() == 1){
				cnames = new ArrayList<>(p.chains.keySet());
				Collections.sort(cnames);
			}else{
				cnames = new ArrayList<>();
				for(int ii = 1;ii < ss.size();ii++){
					cnames.add(ss.get(ii));
				}
			}
			
			for(String cname:cnames){
				ArrayList<PDBResidue> ress = TemplateBaseModeller.prepareTemplateResidues(PepProcess.makeFilteredAA(p.chains.get(cname).residues));
				PSSMScoring testpssm = calcPSSM(ress,ress);
				//>sp|Q8I2J4.1|PROF_PLAF7 RecName: Full=Profilin
				
				ChainIntermediateState ci = new ChainIntermediateState(seq.name,ress,testpssm
						,seq.seq);
				SmithWaterman sw = new SmithWaterman();
				sw.penalO = 10;
				sw.penalE = 1;
				SWResult res = sw.align_with_PSSM(seq.getSeq(),ci.templatePSSM);
				
				//SWResult res = sw.align(seq.getSeq(),SmithWaterman.listToString(ci.template_raw));
				System.out.println(res.score);
				System.out.println(SmithWaterman.listToString(res.qseq));
				System.out.println(SmithWaterman.listToString(res.sseq));
				
				sw.penalO = 10;
				sw.penalE = 1;
				boolean splitFlag = false;
				iterationNum = 1;
				for(int ll = 0;ll < iterationNum;ll++){
					//printPSSM(ci);
					IntermediateSet iss = new IntermediateSet();
					if(splitFlag){
						ArrayList<ChainIntermediateState> cres = splitWithAlignment(res,ci,true);
						ChainIntermediateState freezed = null;
						if(cres.size() > 0){
							ChainIntermediateState maxchain = cres.get(0);

							this.freezeHighscoreSegment(cres);
							ArrayList<ChainIntermediateState> head = new ArrayList<>();
							ArrayList<ChainIntermediateState> tail = new ArrayList<>();
							for(ChainIntermediateState cc:cres){
								if(cc.freezed){
									maxchain = cc;
								}
							}
							int fin = cres.indexOf(maxchain);
							for(int ii = 0;ii < cres.size();ii++){
								if(ii < fin){
									head.add(cres.get(ii));
								}
								if(ii > fin){
									tail.add(cres.get(ii));
								}
							}
							if(head.size() > 0){
								iss.addChain(ChainIntermediateState.merge(head));
							}
							iss.addChain(maxchain);
							if(tail.size() > 0){
								iss.addChain(ChainIntermediateState.merge(tail));
							}
							freezed = maxchain;
							freezed.updateWith(cb,freezed.target_aligned,freezed.template_aligned);
						}else{
							iss.addChain(ci);
						}

						for(ChainIntermediateState css:iss.chains){
							if(css == freezed){
							}else{
								if(css.target.size() < 10 || css.template.size() < 10){
								}else{
									SWResult rres = sw.align_with_PSSM(SmithWaterman.listToString(css.target),css.templatePSSM);
									css.updateWith(cb, rres.qseq, rres.sseq);
									css.setAlignment(rres.qseq,rres.sseq);
								}
							}
						}
						iss.updateEnvironment(this);

						ci = ChainIntermediateState.merge(iss.chains);
						res.qseq = ci.target_aligned;
						res.sseq = ci.template_aligned;
					}else{
						SWResult rres = sw.align_with_PSSM(SmithWaterman.listToString(ci.target),ci.templatePSSM);
						ci.updateWith(cb, rres.qseq, rres.sseq);
						ci.setAlignment(rres.qseq,rres.sseq);
						
						iss.addChain(ci);
						iss.updateEnvironment(this);

					}
				
				}
				//ChainIntermediateState conc = ChainIntermediateState.merge(iss.chains);
				ChainIntermediateState conc =  ci;
				//res = sw.align_with_PSSM(seq.getSeq(),conc.templatePSSM);
				//conc.setAlignment(res.qseq,res.sseq);
				double sc = sw.calcPSSMScore_nogap(conc.target_aligned,conc.template_aligned,conc.templatePSSM);
				
				System.out.println("file:\t"+ss.get(0));
				System.out.println("cnain:\t"+cname);
				System.out.println("score:\t"+sc);
				System.out.println("query:\t"+SmithWaterman.listToString(conc.target_aligned));
				System.out.println("mapped:\t"+SmithWaterman.listToString(conc.template_aligned));
			}
			
		}
	}
	public static ArrayList<PDBResidue> fToP(ArrayList<FloatingResidue> f){
		ArrayList<PDBResidue> r = new ArrayList<>();
		for(FloatingResidue ff:f){
			r.add(ff);
		}
		return r;
	}
	
	
	
	public static void main__(String[] args){
		ArrayList<Character> q = new ArrayList<>();
		ArrayList<Character> s = new ArrayList<>();
		for(char c : "-AA-A-A---A-AA-".toCharArray()) {
			q.add(c);
		}
		
		for(char c : "BBBBB-A---CCCAA-".toCharArray()) {
			s.add(c);
		}
		
		ArrayList<RegionSplit> i = getAlignedRegion(q,s,3,3);
		for(RegionSplit r:i){
			System.out.println(r.alignmentstart+"\t"+r.alignmentend);
			System.out.println(r.qstart+"\t"+r.qend);
			System.out.println(r.sstart+"\t"+r.send);
		}
	}
	public void printPSSM(ChainIntermediateState cs){
		for(int ii = 0;ii < cs.templatePSSM.pssm.seq.size();ii++){
			System.out.print(ii+"\t");
			
			if(cs.flag.get(ii)){
				
				System.out.print("*");
			}
			
			for(Double d:cs.templatePSSM.pssm.scores.get(ii)){
				System.out.print("\t"+d);
			}
			System.out.println("");
			
		}	
		System.out.println("=====");
	}
	
	
	public static void main(String[] args){
		FuzzyTreeAlign ft = new FuzzyTreeAlign();
		ft.processAllPDBInDir( "C:\\dummy\\vbox_share\\bioo\\database\\scop40\\testfas.fas","C:\\dummy\\vbox_share\\bioo\\database\\scop40\\test");
	}
}

class RegionSplit{
	int qstart = -1;
	int qend = -1;
	int sstart = -1;
	int send = -1;
	int alignmentstart = -1;
	int alignmentend = -1;
	
	public int countMatches(){
		int qlen = qend-qstart+1;
		int slen = send-sstart+1;
		int alen = alignmentend-alignmentstart+1;
		int shorter = Math.min(slen,qlen);
		int longer = Math.max(slen,qlen);
		
		return shorter-(alen-longer);
	}
}



class ChainIntermediateState{
	//テンプレート用のリスト。FloatingResidue 用の関数を使うことがあるので変換している
	ArrayList<FloatingResidue> template = new ArrayList<>();
	ArrayList<PDBResidue> template_original = new ArrayList<>();
	
	//間違いチェック用のリスト
	ArrayList<Character> template_raw= new ArrayList<>();
	
	//アラインメントごとに更新される Residue
	ArrayList<FloatingResidue> alignedDummy = new ArrayList<>();
	
	RegionSplit region = null;
	//ターゲットになる配列ギャップ無し
	ArrayList<Character> target = new ArrayList<>();
	
	ArrayList<Character> target_aligned = new ArrayList<>();
	ArrayList<Character> template_aligned = new ArrayList<>();
	
	//ターゲットがマップされた配列。ターゲット側がインサーションだと存在しない
	//テンプレートを変更するのに使う
	ArrayList<Character> alignedResidues = new ArrayList<>();
	
	//テンプレート上でスコアが高く、スコアリングが適応される部分
	ArrayList<Boolean> flag = new ArrayList<>();
	
	PSSMScoring templatePSSM = null;
	
	boolean noTemplate = false;
	boolean noQuery = false;
	
	
	//分割されたものを繋ぎなおす用
	String name = "";
	int qStart = -1;
	int qEnd = -100;
	
	int tStart = -1;
	int tEnd = -100;
	
	boolean freezed = false;
	public void freeze(boolean b){
		freezed = b;
	}
	
	public static ChainIntermediateState merge(ArrayList<ChainIntermediateState> cal){
		ArrayList<PDBResidue> t1 = new ArrayList<>();
		ArrayList<Character> c1 = new ArrayList<>();
		ArrayList<PSSMScoring> ps = new ArrayList<>();
		//ChainIntermediateState(String n,ArrayList<PDBResidue> temp
		//,PSSMScoring tpssm,ArrayList<Character> targ){
		ArrayList<Character> taral = new ArrayList<>();
		ArrayList<Character> temal = new ArrayList<>();
		for(ChainIntermediateState cc:cal){
			t1.addAll(cc.template_original);
			c1.addAll(cc.target);
			ps.add(cc.templatePSSM);
			taral.addAll(cc.target_aligned);
			temal.addAll(cc.template_aligned);
		}
		PSSMScoring mm = PSSMScoring.merge(ps);
		ChainIntermediateState ret =  new ChainIntermediateState(cal.get(0).name,t1,mm,c1);
		ret.region = new RegionSplit();
		RegionSplit las = cal.get(cal.size()-1).region;
		RegionSplit fir = cal.get(0).region;
		ret.region.alignmentend = las.alignmentend;
		ret.region.qend = las.qend;
		ret.region.send = las.send;
		
		ret.region.alignmentstart = fir.alignmentstart;
		ret.region.qstart = fir.qstart;
		ret.region.sstart = fir.sstart;
		ret.name = cal.get(0).name;
		ret.qStart = cal.get(0).qStart;
		ret.tStart = cal.get(0).tStart;
		
		ret.qEnd = cal.get(cal.size()-1).qEnd;
		ret.tEnd = cal.get(cal.size()-1).tEnd;
		
		ret.setAlignment(taral, temal);
		return ret;
	}
	public void setAlignment(ArrayList<Character> targ,ArrayList<Character> temp){
		target_aligned = targ;
		template_aligned = temp;
		
	}
	public ArrayList<ChainIntermediateState> splitAt(int qpos,int tpos){
		ArrayList<PDBResidue> t1 = new ArrayList<>();
		ArrayList<PDBResidue> t2 = new ArrayList<>();
		for(int ii = 0;ii < template_original.size();ii++){
			if(ii < tpos){
				t1.add(template_original.get(ii));
			}else{
				t2.add(template_original.get(ii));
			}
		}
		ArrayList<Character> c1 = new ArrayList<>();
		ArrayList<Character> c2 = new ArrayList<>();

		for(int ii = 0;ii < target.size();ii++){
			if(ii < qpos){
				c1.add(target.get(ii));
			}else{
				c2.add(target.get(ii));
			}
		}
		ArrayList<PSSMScoring> ps = templatePSSM.splitAt(tpos);

		ChainIntermediateState ci1 = new ChainIntermediateState(name,qStart,qStart+qpos-1,tStart,tStart+tpos-1,t1,ps.get(0),c1,true);
		if(qpos  < 0){
			qpos = 0;
		}
		if(tpos < 0){
			tpos = 0;
		}
		ChainIntermediateState ci2 = new ChainIntermediateState(name,qStart+qpos,qEnd,tStart+tpos,tEnd,t2,ps.get(1),c2,true);
		ArrayList<ChainIntermediateState> ret = new ArrayList<>();
		ret.add(ci1);
		ret.add(ci2);
		return ret;
	}
	
	ChainIntermediateState(String n,ArrayList<PDBResidue> temp,PSSMScoring tpssm,ArrayList<Character> targ){
		this(n,0,targ.size()-1,0,temp.size()-1,temp,tpssm,targ,true);
	}
	ChainIntermediateState(
			String n,int s,int e,int ts,int te,ArrayList<PDBResidue> temp
			,PSSMScoring tpssm,ArrayList<Character> targ
	,boolean negativefill){
		name = n;
		for(Character c:targ){
			if(c == '-'){
			}else{
				target.add(Character.toUpperCase(c));
			}
		}
		qStart = s;
		qEnd = e;
		if(qEnd > target.size() -1){
			qEnd = target.size() -1;
		}
		tStart = ts;
		tEnd = te;
		
		if(te < ts){
			noTemplate = true;
		}
		if(e < s){
			noQuery = true;
		}
		
		
		
		templatePSSM = tpssm;
		for(PDBResidue d:temp){//u-n
			if(d.isLigand() || d.isMissing()){
				continue;
			}
			template.add(new FloatingResidue(d,null,null));
			template_original.add(d);
			template_raw.add(PDBResidue.aaMap.get(d.getName()).toCharArray()[0]);
		}
		for(int ii = 0;ii < templatePSSM.pssm.seq.size();ii++){
			if(templatePSSM.pssm.seq.get(ii) != template_raw.get(ii)){
				System.out.println("???"+templatePSSM.pssm.seq.get(ii)+";"+template_raw.get(ii));
				
		//		throw new RuntimeException("incompatibility between PSSM and template");
			}
			//マイナスの場合スコアリングが働かないので BLOSUM を充てる
			
			if(templatePSSM.getScoreAt(ii,templatePSSM.pssm.seq.get(ii)) < 0 && negativefill){
				
				FuzzyTreeAlign.fillWithBLOSUM(templatePSSM,ii);
				//FuzzyTreeAlign.fillWithZERO(templatePSSM,ii);
				flag.add(false);
			}else{
				flag.add(true);
			}
		}
	}
	
	
	public void updateWith(ChainBuilder builder,ArrayList<Character> targetalignment,ArrayList<Character> templatealignment){
		alignedResidues.clear();
		alignedDummy.clear();
		
		int tmpcount = 0;
		for(int ii = 0;ii < targetalignment.size();ii++){
			if(templatealignment.get(ii) == '-'){
			}else{
				alignedResidues.add(targetalignment.get(ii));
				if(targetalignment.get(ii) == '-'){
					//alignedDummy.add(new FloatingResidue("ALA"));
					alignedDummy.add(new FloatingResidue(PepProcess.getAACode(template_raw.get(tmpcount))));
					
					//FloatingResidue fr = new FloatingResidue();	
					//fr.missing = true;
					//alignedDummy.add(fr);
				}else{
					alignedDummy.add(new FloatingResidue(
							PepProcess.getAACode(targetalignment.get(ii))));
				}
				tmpcount++;
			}
		}
		if(alignedDummy.size() == 0){
			System.out.println(";;");
		}
		for(int ii = 0;ii < template.size();ii++){
			FloatingResidue p = null;
			FloatingResidue c = alignedDummy.get(ii);
			FloatingResidue n = null;
			
			if(c.atoms.size() == 0){
				continue;
			}
			
			FloatingResidue pt = null;
			FloatingResidue ct = template.get(ii);
			FloatingResidue nt = null;
			
			if(ii > 0){
				p = alignedDummy.get(ii-1);
				pt = template.get(ii-1);
			}
			if(ii < template.size()-1){
				n = alignedDummy.get(ii+1);
				nt = template.get(ii+1);
			}
			c.prev = p;
			c.next = n;
			c.placeBackBone(ct,pt,nt);
			
			c.sidechainIndex = -1;
			
			c.changeSideChain(builder.sidechains.getNext(c));
			
		}
		/*
		Crash しているかと思ったがあまりしていなかった。
		builder.prepare(alignedDummy);
		for(int ii = 0;ii < template.size();ii++){
			FloatingResidue c = alignedDummy.get(ii);
			double d = builder.calcSideChainCrash(builder._res,ii);
			
			if(d >=0){
				continue;
			}
			
			for(int jj = 0;jj < 10;jj++){
				c.changeSideChain(builder.sidechains.getNext(c));
				d = builder.calcSideChainCrash(builder._res,ii);
				if(d >=0){
					break;
				}
			}
		}
		*/
	}
}
class IntermediateSet{
	
	ArrayList<PDBResidue> environment = new ArrayList<>();
	ArrayList<ChainIntermediateState> chains = new ArrayList<>();
	
	
	
	
	public void addChain(ChainIntermediateState ci){
		chains.add(ci);
	}
	
	//複数チェーンがある場合のために全部を更新する
	public void updateEnvironment(FuzzyTreeAlign ft){
		environment.clear();
		for(ChainIntermediateState cs:chains){
			environment.addAll(cs.alignedDummy);
		}
	
		for(ChainIntermediateState cs:chains){
			
			ArrayList<PDBResidue> al = new ArrayList<>();
			//ダミーに変更する場合
			
			for(FloatingResidue f:cs.alignedDummy){
				al.add(f);
			}
			/*for(FloatingResidue f:cs.template){
				al.add(f);
			}*/
			cs.templatePSSM = ft.calcPSSM(al,environment);
			for(int ii = 0;ii < cs.templatePSSM.pssm.seq.size();ii++){
				if(!cs.flag.get(ii)){
					cs.templatePSSM.pssm.seq.set(ii,cs.template_raw.get(ii));
					FuzzyTreeAlign.fillWithBLOSUM(cs.templatePSSM,ii);
					
					//FuzzyTreeAlign.fillWithZERO(cs.templatePSSM,ii);
				}
			}	
		}
	}
}

class IntermediateStateComparator implements Comparator<ChainIntermediateState>{
	@SuppressWarnings("unchecked")
	public int compare(ChainIntermediateState arg1, ChainIntermediateState arg2){
		int cn = arg1.name.compareTo(arg2.name);
		if(cn != 0){
			return cn;
		}
		if(arg1.qStart == arg2.qStart){
			if(arg1.tStart < arg2.tStart){
				return -1;
			}
			if(arg1.tStart > arg2.tStart){
				return 1;
			}
		}else{
			if(arg1.qStart < arg2.qStart){
				return -1;
			}
			
			if(arg1.qStart > arg2.qStart){
				return 1;
			}
		}
		return 0;
	}
}

class FTPSSMScoring extends PSSMScoring{
	HashSet<Integer> freezed = new HashSet<>();
	public void freezeScore(int rownum){
		freezed.add(rownum);
		for(char a:FuzzyTreeAlign.aachars){
			if(pssm.char_index_map[a] > -1){
				pssm.scores.get(rownum).set(pssm.char_index_map[a]
						,(double)FuzzyTreeAlign.blosum.scoreMat[pssm.seq.get(rownum)][a]);
				
			}
		}
		
	}
}

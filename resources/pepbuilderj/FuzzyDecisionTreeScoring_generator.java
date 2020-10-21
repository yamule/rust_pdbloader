
package pepbuilderj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class FuzzyDecisionTreeScoring_generator{
	/*
	ある空間上の点から見て、何の原子が何Å内に何個あるかを見て、その点に何の原子が来るか予測する。
	*/
	HashMap<String,HashMap<String,Double>> score_with_pred = new HashMap<>();
	HashMap<String,HashMap<String,Double>> prob_with_pred = new HashMap<>();
	HashMap<String,HashMap<String,Double>> count_with_pred = new HashMap<>();
	int count_sum = 0;
	HashMap<String,Integer> atom_count = new HashMap<>();
	public static final double SCORE_INVALIDRESIDUE = -1000;
	RandomForestProcess tree = null;
	double scaleFactor = 3.0;
	FeatureGenerator featureGen = null;
	double log_baseline_factor = 0.00001;
	public static double score_min = -8.0;
	//FeatureGenerator featureGen = new FeatureGeneratorCB_20180501();
	//CB のみの方が大体成績がいい。
	//FeatureGenerator featureGen = new FeatureGeneratorAtom_20180503();
	FuzzyDecisionTreeScoring_generator(FeatureGenerator fgen,String treefile,String scorefile){
		try{
			FileInputStream tf = new FileInputStream(new File(treefile));
			FileInputStream sf = new FileInputStream(new File(scorefile));
			init(fgen,tf,sf);
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	FuzzyDecisionTreeScoring_generator(FeatureGenerator fgen,InputStream tf,InputStream sf){
		try{
			init(fgen,tf,sf);
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	FuzzyDecisionTreeScoring_generator(FeatureGenerator fgen){
		
		ArrayList<String> treefiles = fgen.getResourceFilePath();
		InputStream treefile = 	FuzzyDecisionTreeScoring_generator.class.getResourceAsStream(treefiles.get(0));
		InputStream scorefile = FuzzyDecisionTreeScoring_generator.class.getResourceAsStream(treefiles.get(1));
		init(fgen,treefile,scorefile);
	}
	public void init(FeatureGenerator fgen,InputStream treefile,InputStream scorefile){
		this.featureGen = fgen;
		this.tree = new RandomForestProcess();
		this.tree.loadForestFromStream_SimpleFormat(treefile);
		
		try{
			//予測正答率から計算したスコア。
			BufferedReader br = new BufferedReader(new InputStreamReader(
					scorefile ));
			String ln = null;
			Pattern pat = Pattern.compile("^trueatom_predictedas[\\s]*([^\\s\\:]+)[\\s\\:]+([^\\s\\:].+)");
			Pattern countpat = Pattern.compile("^#count[\\s]+([^\\s\\:]+)[\\s\\:]+([^\\s\\:].+)");
			while((ln = br.readLine()) != null){
				
				Matcher cmat = countpat.matcher(ln);
				if(cmat.find()){
					String k = cmat.group(1).toUpperCase();
					if(k.equals("ALL")){
						continue;
					}
					this.atom_count.put(k,Integer.parseInt(cmat.group(2)));
					this.count_sum += Integer.parseInt(cmat.group(2));
				}
				Matcher mat = pat.matcher(ln);
				if(mat.find()){
					String k = mat.group(1);
					String[] ppt = mat.group(2).replaceAll("[\\r\\n]","").split("[\\s]");
					HashMap<String,Double> hm = new HashMap<String,Double>();
					this.count_with_pred.put(k,hm);
					for(int ii = 0;ii < ppt.length;ii+=2){
						hm.put(ppt[ii],Double.parseDouble(ppt[ii+1]));
					}
				}
			}
			this.calcFreq(this.scaleFactor);
			
			
			
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
	}
	
	public static FuzzyDecisionTreeScoring_generator __generate(InputStream treefile,InputStream scorefile){
		FuzzyDecisionTreeScoring_generator ret = new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180511());
		
		ret.tree = new RandomForestProcess();
		ret.tree.loadForestFromStream_SimpleFormat(treefile);
		
		
		try{
			//予測正答率から計算したスコア。
			BufferedReader br = new BufferedReader(new InputStreamReader(
					scorefile ));
			String ln = null;
			Pattern pat = Pattern.compile("^trueatom_predictedas[\\s]*([^\\s\\:]+)[\\s\\:]+([^\\s\\:].+)");
			Pattern countpat = Pattern.compile("^#count[\\s]+([^\\s\\:]+)[\\s\\:]+([^\\s\\:].+)");
			while((ln = br.readLine()) != null){
				
				Matcher cmat = countpat.matcher(ln);
				if(cmat.find()){
					String k = cmat.group(1).toUpperCase();
					if(k.equals("ALL")){
						continue;
					}
					ret.atom_count.put(k,Integer.parseInt(cmat.group(2)));
					ret.count_sum += Integer.parseInt(cmat.group(2));
				}
				Matcher mat = pat.matcher(ln);
				if(mat.find()){
					String k = mat.group(1);
					String[] ppt = mat.group(2).replaceAll("[\\r\\n]","").split("[\\s]");
					HashMap<String,Double> hm = new HashMap<String,Double>();
					ret.count_with_pred.put(k,hm);
					for(int ii = 0;ii < ppt.length;ii+=2){
						hm.put(ppt[ii],Double.parseDouble(ppt[ii+1]));
					}
				}
			}
			ret.calcFreq(ret.scaleFactor);
			
			
			
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
		
		return ret;
		
	}
	
	/**
	 * Frequency を出す
	 * @param sfact 
	 */
	public void calcFreq(double sfact){
		int allcount = 0;
		int allcount_cb = 0;
		HashSet<String> cbs = new HashSet<>();

		for(String s:count_with_pred.keySet()){
			HashMap<String,Double> hm = count_with_pred.get(s);
			for(String ss:hm.keySet()){
				allcount += hm.get(ss);
				if(ss.indexOf("CB") > -1){
					cbs.add(ss);
					allcount_cb += hm.get(ss);
				}
			}
		}
		scaleFactor = sfact;
		for(String s:count_with_pred.keySet()){
			//System.out.println(s);
			HashMap<String,Double> hm = count_with_pred.get(s);
			HashMap<String,Double> sc = new HashMap<>();
			prob_with_pred.put(s,sc);
			HashMap<String,Double> ssc = new HashMap<>();
			score_with_pred.put(s,ssc);


			int predcount = 0;
			int predcount_cb = 0;
			for(String ss:atom_count.keySet()){
				double cc = 1;
				if(hm.containsKey(ss)){
					cc = hm.get(ss);
				}
				predcount += cc;
				if(cbs.contains(ss)){
					predcount_cb += cc;
				}
			}
			for(String ss:atom_count.keySet()){
				double cc = 1;
				if(hm.containsKey(ss)){
					cc = hm.get(ss);
				}
				//back freq を今入れると邪魔になる
				ssc.put(ss,Math.log(cc/((double)predcount)/(atom_count.get(ss)
						/(double)allcount)+log_baseline_factor));
				sc.put(ss,cc/((double)predcount));
			}
		}
	}
	public void prepare(ArrayList<PDBResidue> al){
		featureGen.prepare(al);
	}
	
	/**
	 * 至らない関数。 Fixme
	 * @param df
	 * @param label
	 * @return 
	 */
	public Double getLabelScore(double[] df,String label){
		String res = this.tree.predict_Majority(df);
		HashMap<String,Double> hm = this.score_with_pred.get(res);
		String rcode = this.featureGen.transCode(label);
		for(String cc:hm.keySet()){
			String ccode = this.featureGen.transCode(cc);
			if(ccode.equals(rcode)){
				return hm.get(cc);
			}
		}
		return null;
	}
	public static double[] listToArray(ArrayList<Double> d){
		if(d == null){
			return null;
		}
		double[] ret = new double[d.size()];
		for(int ii = 0;ii < d.size();ii++){
			ret[ii] = d.get(ii);
		}
		return ret;
	}
	public double getScoreOf(PDBData pdb,PDBResidue target){
		ArrayList<PDBResidue> residues = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			residues.addAll(pdb.chains.get(c).residues);
		}
		this.prepare(residues);
		double ret = 0;
		ArrayList<PDBAtom> tatoms = featureGen.getTargetAtoms(target);
		for(PDBAtom a:tatoms){
			String code = featureGen.transCode(a.parent.getName()+"_"+a.pdb_atom_code);
			
			double df[] = featureGen.generateFeaturesA(a);
			if(df == null){
				continue;
			}
			String res = this.tree.predict_Majority(df);
			HashMap<String,Double> hm = this.prob_with_pred.get(res);
			for(String ss:hm.keySet()){
				
				String rescode = featureGen.transCode(ss);
				if(rescode.equals(code)){
					ret+=hm.get(ss);
				}
			}
		}
		return ret;
		
		
		
		
		//return this.getResidueScoresOf(target).get(target.getName());
	}
	
	
	
	public HashMap<String,Double> getEachResidueLogScoreOn(PDBResidue target){
		return getEachResidueProbabilityOn(target,false);
	}
	/**
	 * 特定の残基の位置に 20 種の残基それぞれが来る確率を返す。
	 * CB 使っている奴だけ用。キーは三文字表記。
	 * @param atoms
	 * @param target
	 * @return 
	 */
	public HashMap<String,Double> getEachResidueProbabilityOn(PDBResidue target){
		return getEachResidueProbabilityOn(target,true);
	}
	public HashMap<String,Double> getEachResidueProbabilityOn(PDBResidue target,boolean prob){
		HashMap<String,Double> ret = new HashMap<>();
		ArrayList<PDBAtom> tatoms = featureGen.getTargetAtoms(target);
		for(PDBAtom a:tatoms){
			String code = featureGen.transCode(a.parent.getName()+"_"+a.pdb_atom_code);
			
			double df[] = featureGen.generateFeaturesA(a);
			if(df == null){
				continue;
			}
			String res = this.tree.predict_Majority(df);
			HashMap<String,Double> hm = null;
			if(prob){
				hm = this.prob_with_pred.get(res);
			}else{
				hm = this.score_with_pred.get(res);
			}
			for(String ss:hm.keySet()){
				String rescode = featureGen.transCode(ss);
				if(!ret.containsKey(rescode)){
					ret.put(rescode,0.0);
				}
				ret.put(rescode,hm.get(ss)+ret.get(rescode));
			}
		}
		if(ret.size() == 0){
			return null;
		}
		return ret;
	}
	public static double max(ArrayList<Double> a){
		double ret = a.get(0);
		for(Double d:a){
			ret = Math.max(ret,d);
		}
		return ret;
	}
	public static double min(ArrayList<Double> a){
		double ret = a.get(0);
		for(Double d:a){
			ret = Math.min(ret,d);
		}
		return ret;
	}
	
	
	public static double average(ArrayList<Double> a){
		double sum = 0;
		for(Double d:a){
			sum += d;
		}
		return sum/a.size();
	}
	
	public static double var(ArrayList<Double> a){
		double ave = average(a);
		double sum = 0;
		for(Double d:a){
			sum += (ave-d)*(ave-d);
		}
		return sum/a.size();
	}
	
	
	/**
	 * 20 種の残基がある残基の位置に来る確率を返す。
	 * @param target
	 * @param back
	 * @param scoring
	 * @return 
	 */
	public static ArrayList<HashMap<String,Double>> calcEachResidueProbabilities(ArrayList<PDBResidue> target,ArrayList<PDBResidue> back,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		ArrayList<HashMap<String,Double>> ret = new ArrayList<>();
		
		for(int ii = 0;ii < target.size();ii++){
			PDBResidue r = target.get(ii);
			ret.add(new HashMap<String,Double>());
		}
		
		for(FuzzyDecisionTreeScoring_generator j:scoring){
			j.prepare(back);
			for(int rr = 0;rr < target.size();rr++){
				PDBResidue r = target.get(rr);
				if(r.isMissing() || r.isLigand()){
					//ret.add(null);
				}else{
					if(r.getCB() == null){
						System.out.println("");
					}
					HashMap<String,Double> scc = j.getEachResidueProbabilityOn(r);
					if(scc != null){
						for(String s:scc.keySet()){
							if(ret.get(rr).containsKey(s)){
								ret.get(rr).put(s,scc.get(s)/scoring.size()+ret.get(rr).get(s));
							}else{
								ret.get(rr).put(s,scc.get(s)/scoring.size());
							}
						}
					}
				}
			}
		}
		
		
		return ret;
	}
	
	/**
	 * 渡されたリストの残基が来る確率を返す。使わないかも？
	 * @param target
	 * @param back
	 * @param scoring
	 * @return 
	 */
	public static double[] calcResidueProbability(ArrayList<PDBResidue> target,ArrayList<PDBResidue> back,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		double[] ret = new double[target.size()];
		ArrayList<HashMap<String,Double>> ps = calcEachResidueProbabilities(target,back,scoring);
		for(int ii = 0;ii < target.size();ii++){
			HashMap<String,Double> p = ps.get(ii);
			if(p.containsKey(target.get(ii).getName())){
				ret[ii] = p.get(target.get(ii).getName());
			}else{
				System.err.println("Invalid residue scores ara added. Sum of scores can not be used.");
				ret[ii] = SCORE_INVALIDRESIDUE;
			}
		}
		return ret;
	}
	
	
	/**
	 * smithwaterman で並べ直すが、あまり意味が無かった。Low complexity region があったりするとちがうかも？
	 * @param filepath
	 * @return 
	 */
	public double[] calcPDBScore_SW(String filepath){
		PDBData p  = PDBData.loadPDBFile(filepath);
		ArrayList<Sequence> seqs = Sequence.PDBToSeq(p);
		HashMap<String,Sequence> cseq = new HashMap<>();
		for(Sequence s:seqs){
			cseq.put(s.name,s);
		}
		FuzzyTreeAlign aligner = new FuzzyTreeAlign(this);
		ArrayList<String> cnames;
		cnames = new ArrayList<>(p.chains.keySet());
		Collections.sort(cnames);
		double score = 0;
		for(String cname:cnames){
			ArrayList<PDBResidue> ress = 
					TemplateBaseModeller
					.prepareTemplateResidues
					(PepProcess.makeFilteredAA(p.chains.get(cname).residues));
			PSSMScoring testpssm = aligner.calcPSSM(ress,ress);
			PSSMData sdd = testpssm.pssm;
			/*
			PSIBLAST の ASCII のようなファイル形式で返す。
			for(int ii = 0;ii <  sdd.scores.size();ii++){
				System.out.print(ii+"\t"+sdd.seq.get(ii)+"\t");
				for(char cc:testpssm.colhead){
					System.out.print(cc+"\t"+sdd.scores.get(ii).get(sdd.char_index_map[cc])+"\t");
				}
				System.out.print("\n");
			}
			*/
			
			
			//>sp|Q8I2J4.1|PROF_PLAF7 RecName: Full=Profilin
//this(n,0,targ.size()-1,0,temp.size()-1,temp,tpssm,targ,true);
			ChainIntermediateState ci = new ChainIntermediateState(
					cseq.get(cname).name,
					0,cseq.get(cname).seq.size()-1,
					0,ress.size()-1
					,ress
					,testpssm
					,cseq.get(cname).seq,false);
			SmithWaterman sw = new SmithWaterman();
			sw.penalO = 10;
			SWResult res = sw.align_with_PSSM(cseq.get(cname).getSeq(),ci.templatePSSM);
			score += res.score;
		}
		double[] ret = new double[4];
		ret[0] = score;
		return ret;
	}
	
	
	
	
	
	
	
	/**
	 * 渡されたファイルのスコアを計算する。
	 * @param filepath
	 * @return 
	 */
	
	public double[] calcPDBScore(String filepath){
		PDBData p  = PDBData.loadPDBFile(filepath);
		double[] ret = {0.0};
		ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
		scoring.add(this);
		ArrayList<PDBResidue> ress = new ArrayList<>();
		for(String cname:p.chains.keySet()){
			ress.addAll(p.chains.get(cname).residues);
		}
		prepareAll(ress,scoring);
		double dd[] = calcResidueScores(ress,scoring);
		for(double d:dd){
			if(d == SCORE_INVALIDRESIDUE){
			}else{
			ret[0] += d;
			}
		}
		
		return ret;
		
	}
	public static PSSMScoring calcPSSM(ArrayList<PDBResidue> target,ArrayList<PDBResidue> allresidue,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		
		ArrayList<HashMap<String,Double>> lis = calcEachResidueProbabilities(
				target,
				allresidue,scoring);
		if(lis.size() == 0){
			PSSMData ret = new PSSMData();
			for(int ii = 0;ii < target.size();ii++){
				HashMap<String,Double> hss = lis.get(ii);
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
	
	public double[] calcBackgroundFrequency(){
		double[] ret = new double[PSSMCalc.background_frequency.length];
		double allcount = 0;
		for(String s:atom_count.keySet()){
			String[] k = s.split("_");
			int di = PSSMCalc.aa_to_index[PDBResidue.aaMap.get(k[0]).charAt(0)];
			if(di < -1){
			}else{
				ret[di] += atom_count.get(s);
				allcount += atom_count.get(s);
			}
		}
		for(int ii = 0;ii < ret.length;ii++){
			ret[ii] /= allcount;
		}
		return ret;
		
	}
	
	/**
	 * ラムダ計算とかしない
	 * 純粋にスコアを計算するだけ
	 * @param target
	 * @param allresidue
	 * @param scoring
	 * @return 
	 */
	public static PSSMScoring calcPSSM_noRescale(ArrayList<PDBResidue> target,ArrayList<PDBResidue> allresidue,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		
		ArrayList<HashMap<String,Double>> lis = calcEachResidueProbabilities(
				target,
				allresidue,scoring);
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
			//大体 1.0 にノーマライズされているはずだが。。
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
		
		
		PSSMData ret = new PSSMData();
		//double[] backfreq = scoring.get(0).calcBackgroundFrequency();
		double[] backfreq = PSSMCalc.background_frequency;
		for(int ii = 0;ii < target.size();ii++){
			HashMap<String,Double> hss = lis.get(ii);
			HashMap<Character,Double> hcc = new HashMap<>();
			for(char c:PSSMCalc.aa_type){
				hcc.put(c,
						Math.log(						
						(double)freq[ii][PSSMCalc.aa_to_index[c]]
								/backfreq[PSSMCalc.aa_to_index[c]]
						+0.00001));
			}
			if(PDBResidue.aaMap.get(target.get(ii).getName()) == null){
				System.out.println(target.get(ii).getName());
			}
			ret.addResidue(PDBResidue.aaMap.get(target.get(ii).getName()).charAt(0), hcc);
			
		}
		PSSMScoring sret = PSSMScoring.prepare(ret);
		return sret;
	}
	
	public static void prepareAll(ArrayList<PDBResidue> template,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		for(FuzzyDecisionTreeScoring_generator s:scoring){
			s.prepare(template);
		}
	}
	/**
	 * prepareAll してないと不具合が起きる。
	 * @param rss
	 * @param scoring
	 * @return 
	 */
	public static double[] calcResidueScores(ArrayList<PDBResidue> rss
			,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		HashMap<PDBResidue,Double> ret = new HashMap<>();
		ArrayList<PDBResidue> target = new ArrayList<>();
		for(PDBResidue ee:rss){
			if(ee.isLigand() || ee.isMissing()){
			}else{
				target.add(ee);
				ret.put(ee,0.0);
			}
			
		}
		for(FuzzyDecisionTreeScoring_generator j:scoring){
			//j.prepare(target);
			for(int rr = 0;rr < target.size();rr++){
				PDBResidue targetres = target.get(rr);
				ArrayList<PDBAtom> tatoms = j.featureGen.getTargetAtoms(targetres);
				for(PDBAtom a:tatoms){
					String code = j.featureGen.transCode(a.parent.getName()+"_"+a.pdb_atom_code);
					double df[] = j.featureGen.generateFeaturesA(a);
					if(df == null){
						continue;
					}
					/*
					double s = j.getLabelScore(df, a.parent.getName()+"_"+a.pdb_atom_code);
					ret.put(targetres,s+ret.get(targetres));
					*/
					String res = j.tree.predict_Majority(df);
					HashMap<String,Double> hm = j.score_with_pred.get(res);
					boolean hitflag = false;
					for(String ss:hm.keySet()){

						String rescode = j.featureGen.transCode(ss);
						if(rescode == null){
							continue;
						}
						score_min = Math.min(hm.get(ss), score_min);
						if(!ret.containsKey(targetres)){
							ret.put(targetres,0.0);
						}
						if(rescode.equals(code)){
							ret.put(targetres,hm.get(ss)+ret.get(targetres));
							hitflag = true;
						}
					}
					if(!hitflag){
						if(!ret.containsKey(targetres)){
							ret.put(targetres,0.0);
						}
						ret.put(targetres,score_min+ret.get(targetres));
					}
					
				}
			}
		}
		double[] dret = new double[rss.size()];
		for(int ii = 0;ii < rss.size();ii++){
			if(ret.containsKey(rss.get(ii))){
				dret[ii] = ret.get(rss.get(ii));
			}else{
				dret[ii] = SCORE_INVALIDRESIDUE;
			}
		}

		return dret;
	}
	
	
	
	
	public static void main_(String[] args){
		//double[] allscore=  calcPDBScore("C:\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain\\1A0K.pdb");
		//System.out.println(allscore[0]);
		
		PDBData p  = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\3KB5.pdb");
		ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
		for(int ii = 0;ii < 1;ii++){
			FuzzyDecisionTreeScoring_generator j = new  FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180511());
			scoring.add(j);
		}
		String cname = p.chains.keySet().iterator().next();
		ArrayList<PDBResidue> ress = p.chains.get(cname).residues;
		FuzzyDecisionTreeScoring_generator.prepareAll(ress,scoring);
		double[] res = FuzzyDecisionTreeScoring_generator.calcResidueScores(ress,scoring);
		
	}
	
	public static void main(String[] args){
		//double[] allscore=  calcPDBScore("C:\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain\\1A0K.pdb");
		//System.out.println(allscore[0]);
		
		ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
		for(int ii = 0;ii < 1;ii++){
			FuzzyDecisionTreeScoring_generator j = new  FuzzyDecisionTreeScoring_generator(new FeatureGeneratorMergeCB_20180511());
			scoring.add(j);
		}
		
		String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\";
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		int posi = 0;
		int nega = 0;
		int zero = 0;
		for(File f:lis){
			if(Pattern.compile("\\.(pdb|ent)$").matcher(f.getName()).find()){
				PDBData p = PDBData.loadPDBFile(f.getPath());
				String cname = p.chains.keySet().iterator().next();
				ArrayList<PDBResidue> ress = p.chains.get(cname).residues;
				prepareAll(ress,scoring);
				double[] res = calcResidueScores(ress,scoring);
				for(double r:res){
					if(r >  Math.sqrt(2)){
						posi++;
					}else if(r <  1.0/Math.sqrt(2)){
						nega++;
					}else{
						zero++;
					}
					
				}
			}
		}
		System.out.println("posi\t"+posi+"\t"+"zero\t"+zero+"\t"+"nega\t"+nega);
	}
}

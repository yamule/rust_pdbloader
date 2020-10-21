
package pepbuilderj;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static pepbuilderj.FuzzyDecisionTreeScoring_generator.calcEachResidueProbabilities;


public class FuzzyDecisionTreeScoring{
	/*
	ある空間上の点から見て、何の原子が何Å内に何個あるかを見て、その点に何の原子が来るか予測する。
	*/
	ArrayList<InputVars> varList = new ArrayList<>();
	HashMap<String,HashMap<String,Double>> score_with_pred = new HashMap<>();
	HashMap<String,HashMap<String,Double>> score_with_pred_cb = new HashMap<>();
	HashMap<String,HashMap<String,Integer>> count_with_pred = new HashMap<>();
	HashMap<String,Integer> atom_count = new HashMap<>();
	
	HashMap<String,Integer> fullname_to_colindex = new HashMap<>();//dist まで入っている
	HashMap<String,HashMap<String,ArrayList<Integer>>> names_to_colindex = new HashMap<>();//res_pdbatomcode
	HashMap<String,ArrayList<Integer>> backboneatom_to_colindex = new HashMap<>();
	RandomForestProcess tree = null;
	double crashPenalty = -5;
	double crashThreshold = 2.0;
	double[] varCassette;
	double scaleFactor = 3.0;
	public static FuzzyDecisionTreeScoring generate(InputStream treefile,InputStream colnamefile,InputStream scorefile){
		FuzzyDecisionTreeScoring ret = new FuzzyDecisionTreeScoring();
		
		ret.tree = new RandomForestProcess();
		ret.tree.loadForestFromStream_SimpleFormat(treefile);
		
		
		try{
			//変数の順番と名前の表記されたファイル。
			BufferedReader br = new BufferedReader(new InputStreamReader(
					colnamefile));
			String ln = null;
			Pattern pat = Pattern.compile("^([^#=]+)=([0-9]+)");
			while((ln = br.readLine()) != null){
				Matcher mat = pat.matcher(ln);
				if(mat.find()){
					if(mat.group(1).equals("target")){
						continue;
					}
					ret.fullname_to_colindex.put(mat.group(1),Integer.parseInt(mat.group(2)));
				}
			}
			br.close();
			ret.varList = new ArrayList<>();
			ret.varCassette = new double[ret.fullname_to_colindex.size()];
			for(int ii = 0;ii < ret.fullname_to_colindex.size();ii++){
				ret.varList.add(new InputVars(null,null,ii));
			}
			
			for(String s:ret.fullname_to_colindex.keySet()){
				int pos = ret.fullname_to_colindex.get(s);
				String[] ppt = s.split("_");
				
				if(ppt.length > 2){
							
					ret.varList.get(pos).residueCode = ppt[0];
					ret.varList.get(pos).atomCode = ppt[1];
					ret.varList.get(pos).setThreshold(Double.parseDouble(ppt[2]));
					if(!ret.names_to_colindex.containsKey(ppt[0])){
						ret.names_to_colindex.put(ppt[0],new HashMap<String,ArrayList<Integer>>());
					}
					if(!ret.names_to_colindex.get(ppt[0]).containsKey(ppt[1])){
						ret.names_to_colindex.get(ppt[0]).put(ppt[1],new ArrayList<Integer>());
					}
					ret.names_to_colindex.get(ppt[0]).get(ppt[1]).add(pos);
					
				}else{
					ret.varList.get(pos).atomCode = ppt[0];
					if(ppt[0].equals("CBALL")){
						ret.varList.get(pos).atomCode = "CB";
					}
					ret.varList.get(pos).residueCode = null;
					ret.varList.get(pos).threshold = Double.parseDouble(ppt[1]);
					
					if(!ret.backboneatom_to_colindex.containsKey(ppt[0])){
						ret.backboneatom_to_colindex.put(ppt[0],new ArrayList<Integer>());
					}
					ret.backboneatom_to_colindex.get(ppt[0]).add(pos);
					
				}
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
		
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
				}
				Matcher mat = pat.matcher(ln);
				if(mat.find()){
					String k = mat.group(1);
					String[] ppt = mat.group(2).replaceAll("[\\r\\n]","").split("[\\s]");
					HashMap<String,Integer> hm = new HashMap<String,Integer>();
					ret.count_with_pred.put(k,hm);
					for(int ii = 0;ii < ppt.length;ii+=2){
						hm.put(ppt[ii],Integer.parseInt(ppt[ii+1]));
					}
				}
			}
			ret.rescale(ret.scaleFactor);
			
			
			
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
		
		return ret;
		
	}
	public void rescale(double sfact){int allcount = 0;
		int allcount_cb = 0;
		HashSet<String> cbs = new HashSet<>();

		for(String s:count_with_pred.keySet()){
			HashMap<String,Integer> hm = count_with_pred.get(s);
			for(String ss:hm.keySet()){
				allcount += hm.get(ss);
				if(ss.indexOf("_CB") > -1 || ss.indexOf("_PCB") > -1){
					cbs.add(ss);
					allcount_cb += hm.get(ss);
				}
			}
		}
		scaleFactor = sfact;
		for(String s:count_with_pred.keySet()){
			//System.out.println(s);
			HashMap<String,Integer> hm = count_with_pred.get(s);
			HashMap<String,Double> sc = new HashMap<>();
			score_with_pred.put(s,sc);

			HashMap<String,Double> sc_cb = new HashMap<>();
			score_with_pred_cb.put(s,sc_cb);

			int predcount = 0;
			int predcount_cb = 0;
			for(String ss:atom_count.keySet()){
				int cc = 1;
				if(hm.containsKey(ss)){
					cc = hm.get(ss);
				}
				predcount += cc;
				if(cbs.contains(ss)){
					predcount_cb += cc;
				}
			}
			for(String ss:atom_count.keySet()){
				int cc = 1;
				if(hm.containsKey(ss)){
					cc = hm.get(ss);
				}
				//LOGとったらなんか性能が下がったので取らない
				sc.put(ss,cc/((double)predcount)/(atom_count.get(ss)/(double)allcount));

				if(cbs.contains(ss)){
					if(predcount_cb == 0){
						sc_cb.put(ss,0.0);
						continue;
					}
					sc_cb.put(ss,
						cc/((double)predcount_cb)/(atom_count.get(ss)/(double)allcount_cb)	
					);
				}
			}
		}
	}
	public double getScoreOf(PDBData pdb,PDBAtom target){
		ArrayList<PDBAtom> atoms = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			for(PDBResidue rr:pdb.chains.get(c).residues){
				for(PDBAtom a:rr.atoms){
					if(!a.isAlternative()){
						atoms.add(a);
					}
				}
			}
		}
		return getScoreOf(atoms,target);
	}
	
	/**
	 * 特定の残基の CB が来る確率を返す。キーは三文字表記。
	 * @param atoms
	 * @param target
	 * @return 
	 */
	public HashMap<String,Double> getResidueCBScoreOf(ArrayList<PDBAtom> atoms,PDBAtom target){
		double penalty = 0;
		if(crashPenalty > 0){
			crashPenalty *= -1;
		}
		for(InputVars i:this.varList){
			i.count = 0;
			i.setTarget(target);
		}
		for(PDBAtom a:atoms){
			if(a.parent == target.parent){
				continue;
			}
			if(a.isAlternative()){
				continue;
			}else{
				if(this.names_to_colindex.containsKey(a.parent.getName())){
					if(this.names_to_colindex.get(a.parent.getName()).containsKey(a.pdb_atom_code)){
						ArrayList<Integer> al  =this.names_to_colindex.get(a.parent.getName()).get(a.pdb_atom_code);
						for(Integer ll:al){
							this.varList.get(ll).checkAndCount(a,true);
						}
					}else if(this.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
						ArrayList<Integer> al  =this.backboneatom_to_colindex.
								get(a.pdb_atom_code);
						for(Integer ll:al){
							this.varList.get(ll).checkAndCount(a,true);
						}
					}else{
						//System.out.println(a.pdb_atom_code);
					}
				}else if(this.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
					ArrayList<Integer> al  =this.backboneatom_to_colindex.
							get(a.pdb_atom_code);
					for(Integer ll:al){
						this.varList.get(ll).checkAndCount(a,true);
					}
				}
			}
		}
		for(int ii = 0;ii < this.varCassette.length;ii++){
			this.varCassette[ii] = 0;
			this.varCassette[ii] = this.varList.get(ii).count;
		}
		String res = this.tree.predict_Majority(varCassette);
		HashMap<String,Double> hm = this.score_with_pred_cb.get(res);
		HashMap<String,Double> ret = new HashMap<>();
		for(String ss:hm.keySet()){
			if(ss.indexOf("_CB") > 0 || ss.indexOf("_PCB") > 0){
				ret.put(ss.split("_")[0],hm.get(ss));
			}
		}
		return ret;
	}
	public double getScoreOf(ArrayList<PDBAtom> atoms,PDBAtom target){
		double penalty = 0;
		if(crashPenalty > 0){
			crashPenalty *= -1;
		}
		for(InputVars i:this.varList){
			i.count = 0;
			i.setTarget(target);
		}
		for(PDBAtom a:atoms){
			if(a.parent == target.parent){
				continue;
			}
			if(a.isAlternative()){
				continue;
			}else{
				if(this.names_to_colindex.containsKey(a.parent.getName())){
					if(this.names_to_colindex.get(a.parent.getName()).containsKey(a.pdb_atom_code)){
						ArrayList<Integer> al  =this.names_to_colindex.get(a.parent.getName()).get(a.pdb_atom_code);
						for(Integer ll:al){
							this.varList.get(ll).checkAndCount(a,false);
						}
					}else if(this.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
						ArrayList<Integer> al  =this.backboneatom_to_colindex.
								get(a.pdb_atom_code);
						for(Integer ll:al){
							this.varList.get(ll).checkAndCount(a,false);
						}
					}else{
						//System.out.println(a.pdb_atom_code);
					}
				}else if(this.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
					ArrayList<Integer> al  =this.backboneatom_to_colindex.
							get(a.pdb_atom_code);
					for(Integer ll:al){
						this.varList.get(ll).checkAndCount(a,false);
					}
				}
			}
		}
		for(int ii = 0;ii < this.varCassette.length;ii++){
			this.varCassette[ii] = 0;
			this.varCassette[ii] = this.varList.get(ii).count;
		}
		String res = this.tree.predict_Majority(varCassette);
		HashMap<String,Double> hm = this.score_with_pred.get(res);
		if(hm.containsKey(target.pdb_atom_code)){//backbone 
			return hm.get(target.pdb_atom_code)+penalty;
		}
		if(hm.containsKey(target.parent.getName()+"_"+target.pdb_atom_code)){
			return hm.get(target.parent.getName()+"_"+target.pdb_atom_code)+penalty;
		}
		return 0+penalty;
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
	
	public static ArrayList<HashMap<String,Double>> ___calcAllResidueProbability
		(ArrayList<PDBResidue> target,ArrayList<PDBResidue> back
				,ArrayList<FuzzyDecisionTreeScoring_generator> scoring){
		ArrayList<HashMap<String,Double>> ret = new ArrayList<>();
		ArrayList<PDBAtom> atoms = new ArrayList<>();
		for(int ii = 0;ii < target.size();ii++){
			PDBResidue r = target.get(ii);
			ret.add(new HashMap<String,Double>());
		}
		
		for(int ii = 0;ii < back.size();ii++){
			atoms.addAll(back.get(ii).atoms);
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
					for(String s:scc.keySet()){
						if(ret.get(rr).containsKey(s)){
							ret.get(rr).put(s,scc.get(s)+ret.get(rr).get(s));
						}else{
							ret.get(rr).put(s,scc.get(s));
						}
					}
				}
			}
		}
		return ret;
	}
	
	public static double[] calcResidueScore(ArrayList<PDBResidue> target,ArrayList<PDBResidue> back,ArrayList<FuzzyDecisionTreeScoring> scoring){
		ArrayList<PDBAtom> atoms = new ArrayList<>();
		double[] ret = new double[target.size()];
		for(int ii = 0;ii < target.size();ii++){
			PDBResidue r = target.get(ii);
			ret[ii] = 0;
		}
		for(int ii = 0;ii < back.size();ii++){
			atoms.addAll(back.get(ii).atoms);
		}
		for(FuzzyDecisionTreeScoring j:scoring){
			for(int rr = 0;rr < target.size();rr++){
				PDBResidue r = target.get(rr);
				if(r.isValidAA()){
					if(r.isMissing()){
						continue;
					}
					for(PDBAtom a:r.atoms){
						if(a.isAlternative()){
							continue;
						}
						//登録されているアミノ酸であるかを調べている
						if(j.names_to_colindex.containsKey(r.getName())){
							if(j.names_to_colindex.get(r.getName()).containsKey(a.pdb_atom_code)){
								double sc =  j.getScoreOf(atoms, a);
								ret[rr] += sc;
							}else if(j.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
								double sc =  j.getScoreOf(atoms, a);
								ret[rr] += sc;
							}
						}
					}
				}
			}
		}
		return ret;
	}
	
	
	
	
	public static double[] calcPDBScore(String filepath){
		
		//String resourcesPath="resources/test.txt";
		//InputStream stream = J48Scoring.class.getResourceAsStream(resourcesPath);
		
		double [] ret = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		
		for(int ii = 0;ii < 3;ii++){
			String treefile="resources/tree_simple."+ii+".dat";
			String varfile="resources/varmap."+ii+".dat";
			String scorefile="resources/scores."+ii+".dat";

			FuzzyDecisionTreeScoring j = FuzzyDecisionTreeScoring.generate(FuzzyDecisionTreeScoring.class.getResourceAsStream(treefile)
					,FuzzyDecisionTreeScoring.class.getResourceAsStream(varfile), 
					FuzzyDecisionTreeScoring.class.getResourceAsStream(scorefile));
			//System.out.println("");
			PDBData p  = PDBData.loadPDBFile(filepath);
			double allscore = 0;
			double allscore_with_back = 0;
			double normalizedscore = 0;
			double normalizedscore_with_back = 0;
			ArrayList<Double> scores = new ArrayList<>();
			ArrayList<Double> rawscores = new ArrayList<>();
			for(String k:p.chains.keySet()){
				PDBChain c = p.chains.get(k);
				ArrayList<PDBResidue> back = new ArrayList<>(c.residues);
				
				//c.residues.clear();
				for(int rr = 0;rr < back.size();rr++){
				
				//	c.residues.add(back.get(rr));
					PDBResidue r = c.residues.get(rr);
					int acou = 0;
					double aval = 0;
					int acouback = 0;
					double avalback = 0;
					if(r.isValidAA()){
						if(r.isMissing()){
							continue;
						}
						for(PDBAtom a:r.atoms){
							if(a.isAlternative()){
								continue;
							}
							//
							if(j.names_to_colindex.containsKey(r.getName())){
								if(j.names_to_colindex.get(r.getName()).containsKey(a.pdb_atom_code)){
									double sc =  j.getScoreOf(p, a);
									sc = Math.max(0,sc);
									allscore += sc;
									allscore_with_back += sc;
									aval+=sc;
									acou++;
									avalback +=sc;
									acouback++;
								}else if(j.backboneatom_to_colindex.containsKey(a.pdb_atom_code)){
									double sc =  j.getScoreOf(p, a);
									sc = Math.max(0,sc);
									allscore_with_back += sc;
									avalback +=sc;
									acouback++;
								}
							}
						}
					}
					if(acou > 0){
						normalizedscore += aval/acou;
					}
					if(acouback > 0){
						normalizedscore_with_back += avalback/acouback;
					}
					//scores.add(normalizedscore/(rr+1));
					scores.add(allscore/(rr+1));
					rawscores.add(aval);
					//rawscores.add(avalback);
				}
			}
			//まあいろいろやっているけど [0] くらいしかとらない
			ret[0] += allscore;
			ret[1] += normalizedscore;
			ret[2] += allscore_with_back;
			ret[3] += normalizedscore_with_back;
			
			//maxpooling
			int window = 10;
			ArrayList<Double> dscores = new ArrayList<>();
			for(int kk = 0;kk < scores.size();kk++){
				double m = scores.get(kk);
				for(int ss = -window;ss <= window;ss++){
					int pos = kk+ss;
					if(pos < scores.size() && pos >=0){
						m = Math.max(m,scores.get(pos));
					}
				}
				dscores.add(m);
			}
			
			ArrayList<Double> mscores = new ArrayList<>();
			for(int kk = 0;kk < rawscores.size();kk++){
				double m = rawscores.get(kk);
				for(int ss = -window;ss <= window;ss++){
					int pos = kk+ss;
					if(pos < rawscores.size() && pos >=0){
						m = Math.max(m,rawscores.get(pos));
					}
				}
				mscores.add(m);
			}
			
			
			ArrayList<Double> sscores = new ArrayList<>();
			for(int kk = 0;kk < rawscores.size();kk++){
				double s = 0;
				int cou = 0;
				for(int ss = -window;ss <= window;ss++){
					int pos = kk+ss;
					if(pos < rawscores.size() && pos >=0){
						s += rawscores.get(pos);
						cou++;
					}
				}
				sscores.add(s/cou);
			}
			
			ArrayList<Double> smscores = new ArrayList<>();
			for(int kk = 0;kk < rawscores.size();kk++){
				double m = sscores.get(kk);
				for(int ss = -1;ss <= 1;ss++){
					int pos = kk+ss;
					if(pos < rawscores.size() && pos >=0){
						m = Math.max(m,sscores.get(pos));
					}
				}
				smscores.add(m);
			}
			
			
			
			ret[4] += average(scores);
			ret[5] += average(dscores);
			ret[6] += average(mscores);
			ret[7] += average(smscores);
			
			ret[8] += scores.get((int)(scores.size()*0.75));
			
		}
		return ret;
		
	}
	
	
	
	public static void main(String[] args){
		//double[] allscore=  calcPDBScore("C:\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain\\1A0K.pdb");
		//System.out.println(allscore[0]);
		
		PDBData p  = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain\\1A0K.pdb");
		ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
		for(int ii = 0;ii < 1;ii++){
			FuzzyDecisionTreeScoring_generator j = new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
			scoring.add(j);
		}
		String cname = p.chains.keySet().iterator().next();
		ArrayList<PDBResidue> ress = p.chains.get(cname).residues;
		ArrayList<HashMap<String,Double>> lis = calcEachResidueProbabilities(
				p.chains.get(cname).residues,
				p.chains.get(cname).residues,scoring);
		for(int ii = 0;ii < ress.size();ii++){
			HashMap<String,Double> hss = lis.get(ii);
			if(hss.containsKey(ress.get(ii).getName())){
				System.out.println(hss.get(ress.get(ii).getName()));
			}else{
				System.out.println(0);
			}
			
		}
	}
}

class InputVars{
	//RESIDUE_ATOM の文字列をいちいち作っているとリソースの無駄なので
	String residueCode = null;
	String atomCode = null;
	double threshold = 4.0;
	PDBAtom target = null;
	int count = 0;
	int colIndex = -1;
	
	InputVars(String rcode,String acode,int i){
		atomCode = acode;
		colIndex = i;
		residueCode = rcode;
	}
	//予測したい原子をセットする
	public void setTarget(PDBAtom t){
		target = t;
	}
	public void setThreshold(double d){
		threshold = d;
	}
	//与えられた Atom が条件を満たす場合にカウントしているようだが
	//List で持っていた方が楽かも？
	//check_ca とすると、CA の方が近い場合無視される
	public void checkAndCount(PDBAtom a,boolean check_ca){
		if(residueCode != null){
			if(!a.parent.getName().equals(residueCode)){
				return;
			}
		}
		if(!a.pdb_atom_code.equals(atomCode)){
			return;
		}
		if(target.parent == a.parent){
			return;
		}
		double dd = a.distance(target);
		if(dd < threshold){
			if(target.parent.getCA().distance(a) < dd){
				return;
			}
			count++;
		}
	}
}

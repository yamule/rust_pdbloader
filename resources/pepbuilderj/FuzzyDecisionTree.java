package pepbuilderj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import static pepbuilderj.RandomForestProcess.getParsedArray_String;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author kimidori
 */
public class FuzzyDecisionTree {
	int classNum = 3;
	ArrayList<DirtySample> samples = new ArrayList<>();
	HashMap<Integer,String> index_to_label = new HashMap<>();
	HashMap<String,Integer> label_to_index = new HashMap<>();
	HashMap<String,Integer> columnname_to_index = new HashMap<>();
	
	ArrayList<DirtySample> testDataset = null;
	double cosineSimilarityThreshold = 1.0/Math.sqrt(2);
	double negativeRatioThreshold = 0.1;
	double[] backgroundFreq = null;
	double[] sampleCount = null;//クラスごとのサンプル数
	double relativeSampleCount = 0;//1クラスにつき1。クラス数に対応
	
	double[] backgroundFreq_test = null;
	
	boolean useRelativeFreq = false;
	double identityThreshold = 1.0;
	double remainedSampleThreshold = 0.01;//sum(remained_samplenum/sampleCount/relativeSampleCount) がこれより小さくならないようにする
	EXRRFNode root = null;
	
	boolean useLogScore=true;
	
	
	public static double logbottom = 0.0001;
	public static int SPLIT_TYPE_ENTROPY = 0;
	public static int SPLIT_TYPE_ENTROPY_ITERATIVE = 1;
	public static int SPLIT_TYPE_LOGSCORE = 2;
	public static int SPLIT_TYPE_GINI = 3;
	int septype = SPLIT_TYPE_GINI;
//ratio/backgroundratio が、log2(x) > 0.5 もしくは < -0.5 であるかで Positive Negative を判断する
	// Smat にしたときに 1 になるかどうかをイメージ
	
	
	public static FuzzyDecisionTree loadTable(String filename,boolean header,boolean hasnamecol){
		FuzzyDecisionTree ret = new FuzzyDecisionTree();
		try{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			
			String line = null;
			if(header){
				line = br.readLine();
				ArrayList<String> al  = getParsedArray_String(line);
				if(hasnamecol){
					al.remove(0);
					for(int ii = 0;ii < al.size();ii++){
						ret.columnname_to_index.put(al.get(ii), ii);
					}
				}
			}
			HashSet<String> answers = new HashSet<>();
			while((line = br.readLine()) != null){
				if(line.replaceAll("[\\s]","").length() == 0){
					continue;
				}
				
				ArrayList<String> al  = getParsedArray_String(line);
				
				int slen = al.size();
				
				String sname ="";
				String answer = "";
				if(hasnamecol){
					slen--;
					sname = al.remove(0);
				}
				slen--;
				answer = al.remove(al.size()-1);
				
				DirtySample ds = new DirtySample();
				ds.classLabel = answer;
				//ds.name = sname;
				ds.name = "";
				ds.values = new ArrayList<>();
				answers.add(answer);
				for(int ii = 0;ii < al.size();ii++){
					ds.values.add(Double.parseDouble(al.get(ii)));
				}
				ret.samples.add(ds);
			}
			ArrayList<String> al = new ArrayList<>(answers);
			for(int ii = 0;ii < al.size();ii++){
				ret.index_to_label.put(ii,al.get(ii));
				ret.label_to_index.put(al.get(ii),ii);
			}
			ret.classNum = al.size();
			for(DirtySample d:ret.samples){
				d.classIndex = ret.label_to_index.get(d.classLabel);
			}
			ret.root = new EXRRFNode();
			ret.root.xsamples = ret.samples;
			ret.calcWeights();
			
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	/**
	 * サンプルの一部をテストデータセットとして分ける
	 * @param ratio 
	 */
	public void extractTestDataset(double ratio){
		double coss = -1;
		int count = 0;
		if(ratio <=0){
			return;
		}
		while(coss < 1.0-cosineSimilarityThreshold/2){
			count ++;
			if(count > 1000){
				System.err.println("Cannot make valiable test dataset.");
				testDataset = null;
				break;
			}
			if(testDataset != null){
				samples.addAll(testDataset);
				testDataset.clear();
			}
			_extractTestDataset(ratio);
			double[] d1 = new double[classNum];
			double[] d2 = new double[classNum];
			for(int ii = 0;ii < classNum;ii++){
				d1[ii] = 0;
				d2[ii] = 0;
			}
			for(DirtySample d:samples){
				d1[d.classIndex] += 1.0;
			}
			for(DirtySample d:samples){
				d2[d.classIndex] += 1.0;
			}
			coss = cosineSimilarity(d1,d2);
		}
		this.root.xsamples = this.samples;
		this.root.testsamples = this.testDataset;
	}
	public void _extractTestDataset(double ratio){
		int tsiz = (int)(samples.size()*ratio);
		int[] tindex = shuffled(samples.size());
		HashSet<Integer> rem = new HashSet<>();
		for(int ii = 0;ii < tsiz;ii++){
			rem.add(tindex[ii]);
		}
		testDataset = new ArrayList<>();
		Iterator<DirtySample> ite = samples.iterator();
		int cou = 0;
		while(ite.hasNext()){
			DirtySample dd = ite.next();
			if(rem.contains(cou)){
				testDataset.add(dd);
				ite.remove();
			}
			cou++;
		}
		calcWeights();
	}
	public void calcWeights(){
		sampleCount = new double[classNum];
		for(int ii = 0;ii < classNum;ii++){
			sampleCount[ii] = 0;
		}
		for(DirtySample d:samples){
			sampleCount[d.classIndex]+=1;
		}
		relativeSampleCount = 0;
		for(int ii = 0;ii < classNum;ii++){
			if(sampleCount[ii] > 0){
				relativeSampleCount+=1;
			}
		}
		
		this.backgroundFreq = calcFrequency_D(samples);
		if(testDataset != null){
			this.backgroundFreq_test = calcFrequency_D(testDataset);
		}
	}
	public double[] calcFrequency_D(Collection<DirtySample> labels){
		double[] ret = new double[this.classNum];
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] = 0;
		}
		for(DirtySample s:labels){
			ret[s.classIndex] += 1;
		}
		
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] /= labels.size();
		}
		return ret;
	}
	public double[] calcFrequency(ArrayList<String> labels){
		double[] ret = new double[this.classNum];
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] = 0;
		}
		for(String s:labels){
			ret[this.label_to_index.get(s)] += 1;
		}
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] /= labels.size();
		}
		return ret;
	}
	
	public void save(String filename){
		RRFTree tree = new RRFTree();
		int nodenum = 0;
		ArrayList<RRFNode> al = new ArrayList<>();
		ArrayList<RRFNode> updated = new ArrayList<>();
		updated.add(root);
		while(updated.size() > 0){
			RRFNode a = updated.remove(0);
			al.add(a);
			if(a.childNodes != null){
				for(int ii = 0;ii < a.childNodes.length;ii++){
					updated.add(a.childNodes[ii]);
				}
			}
		}
		
		tree.nodes = al.toArray(new RRFNode[al.size()]);
		RandomForestProcess p = new RandomForestProcess();
		p.trees = new ArrayList<RRFTree>();
		p.trees.add(tree);
		StringBuffer mapline = new StringBuffer();
		for(String ss:this.columnname_to_index.keySet()){
			mapline.append(ss+"="+String.valueOf(columnname_to_index.get(ss))+"\n");
		}
		if(mapline.length() > 0){
			p.exLabel.put("map",mapline.toString());
		}
		p.saveForestToFile_SimpleFormat(filename);
		
	}
	/**
	 * テーブル形式のデータをパースする用
	 * @param str
	 * @return 
	 */
	public static ArrayList<String> getParsedArray_String(String str){
		String[] ar = ("#"+str).replaceAll("[\\s]+$","").split("[\\s]+");
		ar[0] = ar[0].replaceFirst("^#","");
		ArrayList<String> ret = new ArrayList<String>();
		for(int ii = 0;ii < ar.length;ii++){
			try{
				//if(ar[ii].length() > 0){
					ret.add(ar[ii]);
				//}
			}catch(Exception exx){
				exx.printStackTrace();
			}
		}
		return ret;
	}
	public static double entropy(ArrayList<Integer> al,int classnum){
		int[] count = new int[classnum];
		for(int ii = 0;ii < classnum;ii++){
			count[ii] = 0;
		}
		for(Integer a:al){
			count[a]++;
		}
		double ret = 0;
		for(int ii = 0;ii < classnum;ii++){
			if(count[ii] > 0){
				double ratio = count[ii]/(double)al.size();
				ret -= ratio*Math.log(ratio)/Math.log(2);
			}
		}
		return ret;
	}
	
	public double[] calcThreshold_iterativemerge(ArrayList<DirtySample> d,int targetcolumn,double[] frequency){

		int cdev = classNum;
		cdev = 3;
		
		int[] classcount = new int[cdev];
		int classcount_lower[] = new int[cdev];
		
		
		double[] ff = null;
		if(frequency == null){
			ff = calcFrequency_D(d);
		}else{
			ff = frequency;
		}
		int[] groupid = new int[classNum];
		double basicscore = 0;
		
		
		
		for(int ii = 0;ii < classNum;ii++){
			if(this.useLogScore){
				if(ff[ii]/backgroundFreq[ii] > Math.sqrt(2)){
					groupid[ii] = 2;
				}else if(ff[ii]/backgroundFreq[ii] < 1/Math.sqrt(2)){
					groupid[ii] = 0;
					basicscore += ff[ii];
				}else{
					groupid[ii] = 0;
					basicscore += ff[ii];
				}
			}else{
				
				if(ff[ii]/backgroundFreq[ii] > 1.00){
					groupid[ii] = 2;
				}else if(ff[ii]/backgroundFreq[ii] < 1.00){
					groupid[ii] = 0;
					basicscore += ff[ii];
				}else{
					groupid[ii] = 0;
					basicscore += ff[ii];
				}
			}
		}
		
		if(basicscore == 0){
			
			return null;
		}
		
		int itecount = 3;
		
		double minentropy = Double.MAX_VALUE;
		double minlowratio = Double.MAX_VALUE;
		double minthreshold = Double.MAX_VALUE;
		for(int kk = 0;kk < itecount;kk++){
			double minthreshold_start = 0;
			double minthreshold_end = 0;
			for(DirtySample ds:d){
				ds.classGroupIndex = groupid[ds.classIndex];
			}
			
			for(int ii = 0;ii < cdev;ii++){
				classcount[ii] = 0;
				classcount_lower[ii] = 0;
			}

			for(DirtySample ds:d){
				ds.setTarget(targetcolumn);
				classcount[ds.classGroupIndex]++;
			}
			Collections.sort(d,new DirtySampleComparator());

			for(int jj = 0; jj < d.size()-1;jj++){
				classcount_lower[d.get(jj).classGroupIndex]++;
				if(d.get(jj).getTargetValue() != d.get(jj+1).getTargetValue()){

					double e1 = 0;
					double e2 = 0;
					//分割後のエントロピーの計算
					for(int ii = 0;ii < cdev;ii++){
						if(true){
							double ratio = Math.max(classcount_lower[ii]/(double)(jj+1),logbottom);
							e1 -= ratio*Math.log(ratio)/Math.log(2);
						}
						if(true){
							double ratio = Math.max((classcount[ii]-classcount_lower[ii])/(double)(d.size()-jj-1),logbottom);
							e2 -= ratio*Math.log(ratio)/Math.log(2);
						}
					}
					double ze =  e1*((jj+1)/(double)d.size())+e2*((d.size()-jj-1)/(double)d.size());
					if(ze == minentropy){
						minthreshold_end = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
					}else if(ze < minentropy){
						minentropy = ze;
						minthreshold_start = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
						minthreshold_end = minthreshold_start;
					}
				}
			}
			ArrayList<DirtySample> lower = new ArrayList<>();
			ArrayList<DirtySample> higher = new ArrayList<>();
			for(DirtySample ds:d){
				if(ds.getTargetValue() < minthreshold_start/2+minthreshold_end/2){
					lower.add(ds);
				}else{
					higher.add(ds);
				}
			}
			double[][] yyratio = calcFrequencies(d,targetcolumn,minthreshold_start/2+minthreshold_end/2);
			double lowratio_a = 0;
			double lowratio_b = 0;
			boolean updatedflag = false;
			
			for(int ii = 0;ii < classNum;ii++){
				int classcode = 1;
				if(this.backgroundFreq[ii] > 0){
					if(useLogScore){
						if(yyratio[0][ii]/this.backgroundFreq[ii] < 1.0/Math.sqrt(2)){
							lowratio_a += yyratio[0][ii];
						}
						if(yyratio[1][ii]/this.backgroundFreq[ii] < 1.0/Math.sqrt(2)){
							lowratio_b += yyratio[1][ii];
						}
					}else{
						if(yyratio[0][ii]/this.backgroundFreq[ii] < 1.0){
							lowratio_a += yyratio[0][ii];
						}
						if(yyratio[1][ii]/this.backgroundFreq[ii] < 1.0){
							lowratio_b += yyratio[1][ii];
						}
					}
					boolean highflaga = false;
					boolean highflagb = false;
					if(useLogScore){
						highflaga = yyratio[0][ii]/this.backgroundFreq[ii] > Math.sqrt(2);
						highflagb = yyratio[1][ii]/this.backgroundFreq[ii] > Math.sqrt(2);
					}else{
						
						highflaga = yyratio[0][ii]/this.backgroundFreq[ii] > 1.0;
						highflagb = yyratio[1][ii]/this.backgroundFreq[ii] > 1.0;
					}
					
					
					if(highflaga != highflagb){
						if(highflaga){//低い方に偏っている
							classcode = 0;
						}else{
							classcode = 2;
						}
					}
				}
				if(groupid[ii] != classcode){
					updatedflag = true;
				}
				groupid[ii] = classcode;
			}
			
			
			if(lowratio_b + lowratio_a < minlowratio){
				minlowratio = lowratio_b + lowratio_a;
				minthreshold = minthreshold_start/2+minthreshold_end/2;
			}
			if(!updatedflag){
				break;
			}
			//System.out.println(lowratio_b + lowratio_a);
			
			
		}
		//System.out.println("============");
			
		double[] ret = new double[2];
		ret[0] = minlowratio-basicscore;
		ret[1] = minthreshold;
		//if(ret[0] == 0){
		//	System.out.println(";");
		//}
		return ret;
	}
	
	public double[] calcThreshold(ArrayList<DirtySample> d,int targetcolumn){
		
		int cdev = classNum;
		int[] classcount = new int[cdev];
		int classcount_lower[] = new int[cdev];
		for(int ii = 0;ii < cdev;ii++){
			classcount[ii] = 0;
			classcount_lower[ii] = 0;
		}
		
		
		for(DirtySample ds:d){
			ds.classGroupIndex = ds.classIndex;
		}
		
		for(DirtySample ds:d){
			ds.setTarget(targetcolumn);
			classcount[ds.classGroupIndex]++;
		}
		Collections.sort(d,new DirtySampleComparator());
		double minthreshold_start = 0;
		double minthreshold_end = 0;
		
		double minentropy = Double.MAX_VALUE;
		for(int jj = 0; jj < d.size()-1;jj++){
			classcount_lower[d.get(jj).classGroupIndex]++;
			if(d.get(jj).getTargetValue() != d.get(jj+1).getTargetValue()){
				
				double e1 = 0;
				double e2 = 0;
				//分割後のエントロピーの計算。ブロック重複している
				for(int ii = 0;ii < cdev;ii++){
					if(true){
						double ratio = Math.max(classcount_lower[ii]/(double)(jj+1),logbottom);
						e1 -= ratio*Math.log(ratio)/Math.log(2);
					}
					if(true){
						double ratio = Math.max((classcount[ii]-classcount_lower[ii])/(double)(d.size()-jj-1),logbottom);
						e2 -= ratio*Math.log(ratio)/Math.log(2);
					}
				}
				double ze =  e1*((jj+1)/(double)d.size())+e2*((d.size()-jj-1)/(double)d.size());
				if(ze == minentropy){
					minthreshold_end = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
				}else if(ze < minentropy){
					minentropy = ze;
					minthreshold_start = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
					minthreshold_end = minthreshold_start;
				}
			}
		}
		double[] ret = new double[2];
		ret[0] = minentropy;
		ret[1] = minthreshold_start/2+minthreshold_end/2;
		return ret;
	}
	
	public double[] calcThreshold_Gini(ArrayList<DirtySample> d,int targetcolumn){
		
		int cdev = classNum;
		int[] classcount = new int[cdev];
		int classcount_lower[] = new int[cdev];
		for(int ii = 0;ii < cdev;ii++){
			classcount[ii] = 0;
			classcount_lower[ii] = 0;
		}
		
		
		for(DirtySample ds:d){
			ds.classGroupIndex = ds.classIndex;
		}
		
		for(DirtySample ds:d){
			ds.setTarget(targetcolumn);
			classcount[ds.classGroupIndex]++;
		}
		Collections.sort(d,new DirtySampleComparator());
		double minthreshold_start = 0;
		double minthreshold_end = 0;
		
		double mingini = Double.MAX_VALUE;
		for(int jj = 0; jj < d.size()-1;jj++){
			classcount_lower[d.get(jj).classGroupIndex]++;
			if(d.get(jj).getTargetValue() != d.get(jj+1).getTargetValue()){
				
				double e1 = 0;
				double e2 = 0;
				//分割後のエントロピーの計算。ブロック重複している
				for(int ii = 0;ii < cdev;ii++){
					if(true){
						//https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity
						//double ratio = classcount_lower[ii]/(double)(jj+1)*classcount[ii]/d.size();
						double ratio = classcount_lower[ii]/(double)(jj+1);
						e1 += ratio*ratio;
					}
					if(true){
						//double ratio = (classcount[ii]-classcount_lower[ii])/(double)(d.size()-jj-1)*classcount[ii]/d.size();
						double ratio = (classcount[ii]-classcount_lower[ii])/(double)(d.size()-jj-1);
						e2 += ratio*ratio;
					}
				}
				double ze =  (1.0-e1)*((jj+1)/(double)d.size())+(1.0-e2)*((d.size()-jj-1)/(double)d.size());
				if(ze == mingini){
					minthreshold_end = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
				}else if(ze < mingini){
					mingini = ze;
					minthreshold_start = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
					minthreshold_end = minthreshold_start;
				}
			}
		}
		double[] ret = new double[2];
		ret[0] = mingini;
		ret[1] = minthreshold_start/2+minthreshold_end/2;
		return ret;
	}
	public double[] calcThreshold_maxscoregain(ArrayList<DirtySample> d,int targetcolumn){
		
		int cdev = classNum;
		int[] classcount = new int[cdev];
		int classcount_lower[] = new int[cdev];
		for(int ii = 0;ii < cdev;ii++){
			classcount[ii] = 0;
			classcount_lower[ii] = 0;
		}
		
		
		for(DirtySample ds:d){
			ds.classGroupIndex = ds.classIndex;
		}
		
		for(DirtySample ds:d){
			ds.setTarget(targetcolumn);
			classcount[ds.classGroupIndex]++;
		}
		Collections.sort(d,new DirtySampleComparator());
		double minthreshold_start = 0;
		double minthreshold_end = 0;
		
		double minentropy = Double.MAX_VALUE;
		for(int jj = 0; jj < d.size()-2;jj++){
			classcount_lower[d.get(jj).classGroupIndex]++;
			int samplecount = jj+1;
			if(d.get(jj).getTargetValue() != d.get(jj+1).getTargetValue()){
				
				double e1 = 0;
				double e2 = 0;
				//分割後のスコアの計算
				for(int ii = 0;ii < cdev;ii++){
					
					double ratio1 = Math.max(classcount_lower[ii]/(double)samplecount/this.backgroundFreq[ii],logbottom);
					e1 -= Math.max(-2,Math.min(Math.log(ratio1),2))*classcount_lower[ii];
					
					double ratio2 = Math.max((classcount[ii]-classcount_lower[ii])/((double)d.size()-samplecount)/this.backgroundFreq[ii],logbottom);
					e2 -= Math.max(-2,Math.min(Math.log(ratio2),2))*(classcount[ii]-classcount_lower[ii]);
					
					
					
					//debug
					/*
					if(classcount_lower[ii]/(double)samplecount < this.backgroundFreq[ii]){
						e1 += classcount_lower[ii];
					}else{
						e1 -= classcount_lower[ii];
						
					}
					if((classcount[ii]-classcount_lower[ii])/((double)d.size()-samplecount) < this.backgroundFreq[ii]){
						e2 += classcount[ii]-classcount_lower[ii];
					}else{
						e2 -= classcount[ii]-classcount_lower[ii];
						
					}
					*/
					
					
				}
				double ze =  Math.max(e1,e2);
				if(ze == minentropy){
					minthreshold_end = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
				}else if(ze < minentropy){
					minentropy = ze;
					minthreshold_start = d.get(jj).getTargetValue()/2+d.get(jj+1).getTargetValue()/2;
					minthreshold_end = minthreshold_start;
				}
			}
		}
		double[] ret = new double[2];
		ret[0] = minentropy;
		ret[1] = minthreshold_start/2+minthreshold_end/2;
		return ret;
	}
	
	/**
	 * クラスに含まれるサンプルのうち THreshold 未満になる割合を返す。
	 * @param ds
	 * @param targetindex
	 * @param threshold
	 * @return 
	 */
	public double[] countLowerSampleRatio(ArrayList<DirtySample> ds,int targetindex,double threshold){
		double[] ret = new double[classNum];
		int[] count = new int[classNum];
		for(int ii = 0;ii < classNum;ii++){
			ret[ii] = 0;
			count[ii] = 0;
		}
		for(DirtySample d:ds){
			d.setTarget(targetindex);
			if(d.getTargetValue() < threshold){
				ret[d.classIndex]+=1;
			}
			count[d.classIndex]+=1;
		}
		for(int ii = 0;ii < classNum;ii++){
			if(count[ii] > 0){
				ret[ii] /= count[ii];
			}
		}
		return ret;
	}
	
	public static double cosineSimilarity(double a[],double b[]){
		double alen = 0;
		double blen = 0;
		double xsum = 0;
		for(int ii = 0;ii < a.length;ii++){
			xsum += a[ii]*b[ii];
			alen += a[ii]*a[ii];
			blen += b[ii]*b[ii];
		}
		if(alen == 0|| blen == 0){
			return 1;
		}
		return xsum/Math.sqrt(alen)/Math.sqrt(blen);
		
	}
	
	
	public double[][] calcFrequencies(ArrayList<DirtySample> ds,int targetindex,double threshold){
		ArrayList<DirtySample> upper = new ArrayList<>();
		ArrayList<DirtySample> lower = new ArrayList<>();
		for(DirtySample d:ds){
			d.setTarget(targetindex);
			if(d.getTargetValue() < threshold){
				lower.add(d);
			}else{
				upper.add(d);
			}
		}
		double[][] ret = new double[2][0];
		ret[0] = this.calcFrequency_D(lower);
		ret[1] = this.calcFrequency_D(upper);
		return ret;
	}
	
	
	public double calcCosineSimilarityWithTestSet(ArrayList<DirtySample> train
			,ArrayList<DirtySample> test,int targetindex,double threshold){
		/*
		double[] aratio = countLowerSampleRatio(train,targetindex,threshold);
		double[] bratio = countLowerSampleRatio(test,targetindex,threshold);
		double[] aratio2 = new double[classNum];
		double[] bratio2 = new double[classNum];
		for(int ii = 0;ii < classNum;ii++){
			aratio2[ii] = 1.0-aratio[ii];
			bratio2[ii] = 1.0-bratio[ii];
		}
		*/
		double[][] xxratio = calcFrequencies(train,targetindex,threshold);
		double[][] yyratio = calcFrequencies(test,targetindex,threshold);
		for(int ii = 0;ii < classNum;ii++){
			if(this.backgroundFreq[ii] > 0){
				xxratio[0][ii] /= this.backgroundFreq[ii];
				xxratio[1][ii] /= this.backgroundFreq[ii];
			}
			
			if(this.backgroundFreq_test[ii] > 0){
				yyratio[0][ii] /= this.backgroundFreq_test[ii];
				yyratio[1][ii] /= this.backgroundFreq_test[ii];
			}
		}
		
		/*
		double base = logbottom;
		for(int ii = 0;ii < classNum;ii++){
			xxratio[0][ii] = Math.log(Math.max(xxratio[0][ii],logbottom));
			xxratio[1][ii] = Math.log(Math.max(xxratio[1][ii],logbottom));
			yyratio[0][ii] = Math.log(Math.max(yyratio[0][ii],logbottom));
			yyratio[1][ii] = Math.log(Math.max(yyratio[1][ii],logbottom));
			
		}
		*/
		//まあ大体同じになる
		//System.out.println(cosineSimilarity(xxratio[0],yyratio[0])-cosineSimilarity(aratio,bratio));
		//return Math.min(cosineSimilarity(aratio,bratio),cosineSimilarity(aratio2,bratio2));
		return Math.min(cosineSimilarity(xxratio[0],yyratio[0]),cosineSimilarity(xxratio[1],yyratio[1]));
		
	}
	
	
	public double calcNegativeRatioWithTestSet(ArrayList<DirtySample> test,int targetindex,double threshold){

		double[][] yyratio = calcFrequencies(test,targetindex,threshold);
		double negative_a = 0;
		double negative_b = 0;
		for(int ii = 0;ii < classNum;ii++){
			if(this.backgroundFreq_test[ii] > 0){
				if(useLogScore){
					if(yyratio[0][ii]/this.backgroundFreq_test[ii] < 1/Math.sqrt(2)){
						negative_a += yyratio[0][ii];
					}

					if(yyratio[1][ii]/this.backgroundFreq_test[ii] < 1/Math.sqrt(2)){
						negative_b += yyratio[1][ii];
					}
				}else{
					if(yyratio[0][ii]/this.backgroundFreq_test[ii] < 1.0){
						negative_a += yyratio[0][ii];
					}

					if(yyratio[1][ii]/this.backgroundFreq_test[ii] < 1.0){
						negative_b += yyratio[1][ii];
					}
				}
			}
		}
		return Math.max(negative_a,negative_b);
		
	}
	/**
	 * クラスごとに重みづけされたサンプルの割合の合計
	 * 多いサンプルほど低く見積もられ全体合計は 1 
	 * @param ds
	 * @return 
	 */
	public double calcFreqSum(Collection<DirtySample> ds){
		if(useRelativeFreq){
		int[] count = new int[classNum];
		for(int ii = 0;ii < classNum;ii++){
			count[ii] = 0;
		}
		
		double ret = 0;
		for(DirtySample d:ds){
			count[d.classIndex]++;
		}
		for(int ii = 0;ii < classNum;ii++){
			ret += count[ii]/(double)this.sampleCount[ii]/this.relativeSampleCount;
		}
		
		return ret;
		}else{
			return ds.size()/(double)this.samples.size();
		}
	}
	
	
	public void calcOptimalSplit(EXRRFNode parentnode,int calctype){
		ArrayList<DirtySample> d = parentnode.xsamples;
		int[] si = shuffled(d.get(0).values.size());
		//System.out.println(d.get(0).values.size());
		double minentropy = Double.MAX_VALUE;
		double minthreshold = Double.MAX_VALUE;
		int minindex = -1;
		
		double[] d1 = new double[classNum];
		for(int ii = 0;ii < classNum;ii++){
			d1[ii] = 0;
		}
		double maxcount = 0;
		double sumcount = 0;
		for(DirtySample dd:d){
			d1[dd.classIndex] += 1.0;
			if(maxcount < d1[dd.classIndex]){
				maxcount  = d1[dd.classIndex];
			}
		}
		if(maxcount/d.size() >= identityThreshold || (remainedSampleThreshold > 0 && calcFreqSum(d) < remainedSampleThreshold )){
			
			parentnode.asLeaf(true);
			//parentnode.predValue = index_to_label.get(parentnode.getMajority(classNum));
			//System.out.println(d.size()+" belongs same class.");
			return;
		}
		parentnode.asLeaf(false);
		double[] freq = null;
		if(calctype == SPLIT_TYPE_ENTROPY_ITERATIVE){
			freq = calcFrequency_D(d);
		}
		for(int ii = 0;ii < si.length;ii++){//今のところ全部調べている
			double res[] = null;
			if(calctype == SPLIT_TYPE_ENTROPY_ITERATIVE){
				res = calcThreshold_iterativemerge(d,si[ii],freq);
				if(res == null){
					minindex = si[ii];
					minentropy = 0;
					minthreshold = 0;
					break;
				
				}
			}else if(calctype == SPLIT_TYPE_ENTROPY){
				res = calcThreshold(d,si[ii]);
			}else  if(calctype == SPLIT_TYPE_GINI){
				res = calcThreshold_Gini(d,si[ii]);
			}else  if(calctype == SPLIT_TYPE_LOGSCORE){
				res = this.calcThreshold_maxscoregain(d,si[ii]);
			}
			if(res[0] < minentropy){
				if(parentnode.testsamples != null && parentnode.testsamples.size() > 0){
					double cres = this.calcCosineSimilarityWithTestSet(parentnode.xsamples,
							parentnode.testsamples, si[ii],res[1]);
					
					if(cres < cosineSimilarityThreshold){
						continue;
					}
					
					
					
				}
				
				ArrayList<ArrayList<DirtySample>> ddal = split(d,si[ii],res[1]);
				if(this.calcFreqSum(ddal.get(0)) < this.remainedSampleThreshold
					||	this.calcFreqSum(ddal.get(1)) < this.remainedSampleThreshold
					){
					continue;
				}
				minindex = si[ii];
				minentropy = res[0];
				minthreshold = res[1];
				if(minentropy == 0){
					break;
				}
			}
		}
		
		if(calctype == SPLIT_TYPE_ENTROPY_ITERATIVE && minentropy == 0){
			calcOptimalSplit(parentnode,SPLIT_TYPE_ENTROPY);
			return;
		}
		if(minentropy == Double.MAX_VALUE){
			parentnode.leafFlag = true;
			
			//System.out.println(d.size()+" exceeds cosine similarity.");
			return;
		}
		ArrayList<ArrayList<DirtySample>> al = split(d,minindex,minthreshold);
		
		
		parentnode.splitThreshold = minthreshold;
		parentnode.splitVarNum = minindex;
		
		EXRRFNode r1 = new EXRRFNode();
		r1.parentNode = parentnode;
		r1.xsamples = al.get(0);
		r1.lastScore = minentropy;
		EXRRFNode r2 = new EXRRFNode();
		r2.parentNode = parentnode;
		r2.xsamples = al.get(1);
		r2.lastScore = minentropy;
		parentnode.childNodes = new RRFNode[2];
		parentnode.childNodes[0] = r1;
		parentnode.childNodes[1] = r2;
		if(parentnode.testsamples != null){
			ArrayList<ArrayList<DirtySample>> dal = split(parentnode.testsamples,minindex,minthreshold);
			r1.testsamples = dal.get(0);
			r2.testsamples = dal.get(1);
			
		}
	}
	
	
	/**
	 * 予測クラスに対する正解クラスの数を返す。
	 * @param arghash
	 * @param rg
	 * @return 
	 */
	public static HashMap<String,HashMap<String,Integer>> countTrueClass(HashMap<String,String> arghash,RandomForestProcess rg,boolean shufflefordebug){
		boolean header = arghash.containsKey("-hasheader");
		boolean hasnamecol = arghash.containsKey("-hasnamecol");
		boolean hasanswercol = arghash.containsKey("-hasanswercol");
		HashMap<String,HashMap<String,Integer>> ret = new HashMap<>();
		ArrayList<String> answer_debug = new ArrayList<>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(arghash.get("-infile"))));
			
			String line = null;
			if(header){
				line = br.readLine();
			}
			int okcount = 0;
			int allcount = 0;
			while((line = br.readLine()) != null){
				if(line.replaceAll("[\\s]","").length() == 0){
					continue;
				}
				
				ArrayList<String> al  = getParsedArray_String(line);
				
				int slen = al.size();
				
				String sname ="";
				String answer = "";
				if(hasnamecol){
					slen--;
					sname = al.remove(0);
				}
				if(hasanswercol){
					slen--;
					answer = al.remove(al.size()-1);
				}
				double[] d = new double[slen];
				for(int ii = 0;ii < al.size();ii++){
					d[ii] = Double.parseDouble(al.get(ii));
				}
				if(shufflefordebug){
					answer_debug.add(answer);
					if(answer_debug.size() > 1000){
						answer = answer_debug.remove((int)(answer_debug.size()*Math.random()));
					}
					/*
					for(int ii = 0;ii < al.size()*2;ii++){
						int ks = (int)(Math.random()*al.size());
						int ks1 = (int)(Math.random()*al.size());
						double dkk = d[ks];
						d[ks] = d[ks1];
						d[ks1] = dkk;
					}*/
				}
				String predres = rg.predict_Majority(d);
				
				if(hasanswercol){
					if(predres.equals(answer)){
						okcount++;
					}
					allcount++;
				}
				
				
				if(hasanswercol){
					if(!ret.containsKey(predres)){
						ret.put(predres,new HashMap<String,Integer>());
					}
					
					if(!ret.get(predres).containsKey(answer)){
						ret.get(predres).put(answer,0);
					}
					ret.get(predres).put(answer,ret.get(predres).get(answer)+1);
				}
			}
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
			
		}
		return ret;
	}
	public void makeTree(EXRRFNode start){
		ArrayList<EXRRFNode> updated  =new ArrayList<>();
		updated.add(start);
		
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
			if(r.xsamples.size() == 0){
				System.out.println("?+");
			}

			int npcount[] = countNegaPosiSamples(r.xsamples);
			/*
			if(npcount[1]/(double)r.xsamples.size() > 0.9){
				r.leafFlag = true;
				int za = r.getMajority(classNum);
				r.predValue = index_to_label.get(za)+"%";
				System.out.println(npcount[1]+";"+r.xsamples.size()+";");
			}
			*/
			
			if(r.leafFlag){

			}else{

				if(r.parentNode != null){
					r.depth = r.parentNode.depth+1;
				}
				if(septype == SPLIT_TYPE_ENTROPY_ITERATIVE){
					calcOptimalSplit(r,(r.depth > 3)?(SPLIT_TYPE_ENTROPY_ITERATIVE):(SPLIT_TYPE_ENTROPY));
				}else{
					calcOptimalSplit(r,septype);
				}
				//calcOptimalSplit(r,SPLIT_TYPE_ENTROPY);
				
				boolean retrysplit = false;//Negative が多い場合 leaf でももう一度検討する
				if(!r.leafFlag){
					updated.add((EXRRFNode)r.childNodes[0]);
					updated.add((EXRRFNode)r.childNodes[1]);
				}else if(retrysplit){
					npcount = countNegaPosiSamples(r.xsamples);
					if(npcount[0]/(double)r.xsamples.size() > 0.1){
						r.leafFlag = false;
						double remthresh = remainedSampleThreshold;
						remainedSampleThreshold *= 0.5;
						//calcOptimalSplit(r,SPLIT_TYPE_ENTROPY_ITERATIVE);
						calcOptimalSplit(r,septype);
						//calcOptimalSplit(r,SPLIT_TYPE_LOGSCORE);//LOGSCOREは結果悪い
						remainedSampleThreshold = remthresh;
						
						if(r.leafFlag){
							int zb = r.getMajority(classNum);
							if(index_to_label.get(zb) == null){
								System.out.println("???");
							}
							r.predValue = index_to_label.get(zb);

						}else{
							updated.add((EXRRFNode)r.childNodes[0]);
							updated.add((EXRRFNode)r.childNodes[1]);
						}
					}else{
						int zb = r.getMajority(classNum);
						r.predValue =index_to_label.get(zb);
					}
				}else{
					int zb = r.getMajority(classNum);
					r.predValue =index_to_label.get(zb);
				}
			}
		}
	}
	public void splitWith(EXRRFNode parentnode,int varIndex,double threshold){
		ArrayList<ArrayList<DirtySample>> al = split(parentnode.xsamples,varIndex,threshold);
		parentnode.splitThreshold = threshold;
		parentnode.splitVarNum = varIndex;

		EXRRFNode r1 = new EXRRFNode();
		r1.parentNode = parentnode;
		r1.xsamples = al.get(0);
		EXRRFNode r2 = new EXRRFNode();
		r2.parentNode = parentnode;
		r2.xsamples = al.get(1);
		parentnode.childNodes = new RRFNode[2];
		parentnode.childNodes[0] = r1;
		parentnode.childNodes[1] = r2;
		if(parentnode.testsamples != null){
			ArrayList<ArrayList<DirtySample>> dal = split(parentnode.testsamples,
					varIndex,threshold);
			r1.testsamples = dal.get(0);
			r2.testsamples = dal.get(1);
		}
		makeTree(r1);
		makeTree(r2);
	}
	/**
	 * Negative なサンプルの割合とか計算して、ノードの悪さを測ろうとしている。
	 * @param parentnode
	 * @return 
	 */
	public ArrayList<EXRRFNode> calcLoss(EXRRFNode parentnode){
		ArrayList<EXRRFNode> checked = new ArrayList<>();
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(parentnode);
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
			checked.add(r);
			double[] ff = null;
			ff = calcFrequency_D(r.testsamples);//どちらがいいか、
			double closs = 0;

			for(int ii = 0;ii < classNum;ii++){
				if(useLogScore){
					if(ff[ii]/backgroundFreq[ii] > Math.sqrt(2)){
						closs -= ff[ii];
					}else if(ff[ii]/backgroundFreq[ii] < 1/Math.sqrt(2)){
						closs += ff[ii];
					}else{
					}
				}else{
					
					if(ff[ii]/backgroundFreq[ii] > 1.0){
						closs -= ff[ii];
					}else if(ff[ii]/backgroundFreq[ii] < 1.0){
						closs += ff[ii];
					}else{
					}
				}
			}
			
			r.loss = closs;
			if(r.leafFlag){
				
			}else{
				updated.add((EXRRFNode)r.childNodes[0]);
				updated.add((EXRRFNode)r.childNodes[1]);

			}
		}
		return checked;
	}
	
	
	
	public void tryRollBack(EXRRFNode parentnode,int depth){
		calcLoss(parentnode);
		parentnode.calcPlayoutLoss(depth);
		double firstloss = parentnode.childLoss;
		double bestloss = firstloss;
		SplittingVar bestsplit = new SplittingVar(parentnode.splitVarNum,parentnode.splitThreshold);
		parentnode.collectChildSplit(4);
		for(SplittingVar sv:parentnode.childSplit){
			splitWith(parentnode,sv.varIndex,sv.threshold);
			
			calcLoss(parentnode);
			parentnode.calcPlayoutLoss(depth);
			double currentloss = parentnode.childLoss;
			if(currentloss < bestloss){
				bestloss = currentloss;
				bestsplit = sv;
			}
			//System.out.println(firstloss+";"+bestloss+";"+currentloss);
		}
		splitWith(parentnode,bestsplit.varIndex,bestsplit.threshold);
		//System.out.println("===================");
	}
	
	
	public int[] shuffled(int inum){
		int index[] = new int[inum];
		for(int ii = 0;ii < index.length;ii++){
			index[ii] = ii;
		}
		for(int ii = 0;ii < index.length;ii++){
			int i1 = (int)(Math.random()*index.length);
			int i2 = (int)(Math.random()*index.length);
			int d1 = index[i2];
			index[i2] = index[i1];
			index[i1] = d1;
		}
		return index;
	}
	/**
	 * Threshold より小さいサンプルを 0、大きいサンプルを 1 のリストに入れて返す。
	 * @param d
	 * @param columnnum
	 * @param threshold
	 * @return 
	 */
	public ArrayList<ArrayList<DirtySample>> split(ArrayList<DirtySample> d,int columnnum,double threshold){
		 ArrayList<ArrayList<DirtySample>> ret = new  ArrayList<>();
		 ArrayList<DirtySample> r1 = new  ArrayList<>();
		 ArrayList<DirtySample> r2 = new  ArrayList<>();
		 for(DirtySample dd:d){
			dd.setTarget(columnnum);
			if(dd.getTargetValue() < threshold){
				r1.add(dd);
			}else{
				r2.add(dd);
			}
		 }
		 ret.add(r1);
		 ret.add(r2);
		 return ret;
	}
	
	public int[] countNegaPosiSamples(Collection<DirtySample> cc){
		
		double[] ff = calcFrequency_D(cc);
		int[] groupid = new int[classNum];
		for(int ii = 0;ii < classNum;ii++){
			if(useLogScore){
				if(ff[ii]/backgroundFreq[ii] < 1/Math.sqrt(2)){
					groupid[ii] = 0;
				}else if(ff[ii]/backgroundFreq[ii] > Math.sqrt(2)){
					groupid[ii] = 1;
				}else{
					groupid[ii] = -1;
				}
			}else{
				if(ff[ii]/backgroundFreq[ii] < 1.0){
					groupid[ii] = 0;
				}else if(ff[ii]/backgroundFreq[ii] > 1.0){
					groupid[ii] = 1;
				}else{
					groupid[ii] = -1;
				}
			}
		}
		int negacount = 0;
		int posicount = 0;
		for(DirtySample xs:cc){
			if(groupid[xs.classIndex] == 0){
				negacount++;
			}
			
			if(groupid[xs.classIndex] == 1){
				posicount++;
			}
		}
		int[] ret = new int[2];
		ret[0] = negacount;
		ret[1] = posicount;
		return ret;
	}
	
	
	
	public static double[] getScores(HashMap<String,HashMap<String,Integer>> testcount){
		HashSet<String> predclass_ = new HashSet<>();
		predclass_.addAll(testcount.keySet());
		ArrayList<String> predclass = new ArrayList<>(predclass_);
		HashSet<String> trueclass_ = new HashSet<>();
		for(String p:predclass){
			if(testcount.containsKey(p)){
				trueclass_.addAll(testcount.get(p).keySet());
			}
		}
		ArrayList<String> trueclass = new ArrayList<>(trueclass_);
		HashMap<String,Integer> trueclass_index = new HashMap<>();
		HashMap<String,Integer> all_training_class_count = new HashMap<>();
		for(int ii = 0;ii < trueclass.size();ii++){
			trueclass_index.put(trueclass.get(ii),ii);
			all_training_class_count.put(trueclass.get(ii),0);
		}
		int all_train_count = 0;
		for(String pp:testcount.keySet()){
			for(String tt:testcount.get(pp).keySet()){
				if(tt.equals("HOH")){
					continue;
				}
				all_training_class_count.put(tt,all_training_class_count.get(tt)+
						testcount.get(pp).get(tt));
				all_train_count += testcount.get(pp).get(tt);
			}
		}

		HashMap<String,Double> all_training_class_ratio = new HashMap<>();
		for(String tt:all_training_class_count.keySet()){
			if(tt.equals("HOH")){
				continue;
			}
			all_training_class_ratio.put(tt,all_training_class_count.get(tt)/(double)all_train_count);
		}
		int all_posi = 0;
		int all_nega = 0;
		for(String pp:testcount.keySet()){
			int[] att = new int[trueclass_index.size()];
			for(int ii = 0;ii < att.length;ii++){
				att[ii] = 0;
			}
			double score = 0;
			int predtraincount = 0;
			for(String ttt:testcount.get(pp).keySet()){
				if(ttt.equals("HOH")){
					continue;
				}
				predtraincount += testcount.get(pp).get(ttt);
			}
			int negacount = 0;
			int posicount = 0;
			for(String tt:testcount.get(pp).keySet()){
				if(tt.equals("HOH")){
					continue;
				}
				att[trueclass_index.get(tt)] = testcount.get(pp).get(tt);
				double d = 0;
				if(testcount.get(pp).containsKey(tt)){
					//double ratio = (testcount.get(pp).get(tt))/(double)predtraincount/all_training_class_ratio.get(tt);
					//if(testcount.get(pp).get(tt)/(double)predtraincount < all_training_class_ratio.get(tt)){
					//	d = -5;
					//}else{
					//	d = Math.log(ratio);
					//}
					double ratio = Math.max(testcount.get(pp).get(tt)/(double)predtraincount/all_training_class_ratio.get(tt)
							,logbottom);
					d = Math.log(ratio);
				}else{
					d = -5;
				}
				d = d/Math.log(2);
				if(d < -0.5){
					negacount += att[trueclass_index.get(tt)];
				}else if(d > 0.5){
					posicount += att[trueclass_index.get(tt)];
				}
				score += att[trueclass_index.get(tt)]*d;
			}
			all_posi += posicount;
			all_nega += negacount;
			//System.out.println(pp+"\t"+score+"\t"+posicount+"\t"+negacount);


		}
		double[] ret = new double[3];
		ret[0] = all_posi;
		ret[1] = all_nega;
		ret[2] = all_train_count;
		return ret;
	}
	public static void main(String[] args){
		try{
			//double[] nodesample = {0.001,0.00025,0.0005};//0.0025 は少し下がった
			//double[] cosdist = {0.9,0.8,0.7,0.5,0.25,};
			//double[] testsize = {0.1,0.25,0.5};
			//String parentdir = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sparse\\";
			
			//String parentdir = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\";
			//String parentdir = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_onlycb_jumpcb45\\";
			String parentdir = "_C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\";
			//double[] nodesample = {0.001,0.0005};//0.0025 は少し下がった
			double[] nodesample = {0.0005};//0.0025 は少し下がった
			double[] cosdist = {0.75};
			double[] testsize = {0.2};
		int treecount = 0;
		try{
		//Thread.sleep(60l*10*1000);	
		}catch(Exception exx){
			
		}
		for(double coss:cosdist){
			for(double sizz:testsize){
					
				for(double nsampp:nodesample){
					for(int ll = 0;ll < 20
							;ll++){
						//String infile = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\3fold\\input."+ll+".train.dat";
						//String intestfile ="C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\3fold\\input."+ll+".test.dat"; 
						//System.out.println("=============================");
						//System.out.println( "cosinedist threshold:"+coss);
						//System.out.println("testsize:"+sizz);
						
						String infile = parentdir+"input."+ll+".train.dat";
						String intestfile = parentdir+"input."+ll+".test.dat";
						
						String treefile = parentdir+"tree_simple."+ll+".dat"+"_"+treecount;
							treecount++;
						boolean maketree = true;
						boolean debug = false;
						double all_posi = 0;
						double all_nega = 0;
						if(maketree){
							FuzzyDecisionTree co = FuzzyDecisionTree.loadTable(
								infile
								, true,true);
							ArrayList<EXRRFNode> updated = new ArrayList<>();
							co.extractTestDataset(sizz);
							co.cosineSimilarityThreshold = coss;

							co.identityThreshold = 0.62;
							co.remainedSampleThreshold = nsampp;
							
							//木の作成
							co.makeTree(co.root);
							
							HashSet<EXRRFNode> checked = new HashSet<>();
							for(int ii = 0; ii < 10;ii++){
								ArrayList<EXRRFNode> al = co.calcLoss(co.root);
								double maxloss = -100;
								EXRRFNode maxnode = null;
								ArrayList<EXRRFNode> leafs = new ArrayList<>();
								for(int kk = 0;kk < al.size();kk++){
									EXRRFNode r = al.get(kk);
									if(r.isLeaf()){
										leafs.add(r);
										if(r.loss > maxloss || maxnode == null){
											if(r.depth > 2){
												if(!checked.contains(r.parentNode.parentNode)){
													maxloss = r.loss;
													maxnode = r;
												}
											}
										}
									}
								}
								
								/*if(maxnode.depth > 3 && !checked.contains(maxnode.parentNode.parentNode.parentNode)){
									EXRRFNode p = (EXRRFNode)maxnode.parentNode.parentNode.parentNode;
									checked.add(p);
									co.tryRollBack(p,4);
								}else 
									*/
								if(maxnode != null && maxnode.depth > 2 && !checked.contains(maxnode.parentNode.parentNode)){
									//EXRRFNode p = (EXRRFNode)maxnode.parentNode.parentNode;
									//checked.add(p);
									//co.tryRollBack(p,3);
									break;
								}else{
									break;
								}
							}
							ArrayList<EXRRFNode> ns = co.root.getChildNodes();
							int ncou = 0;
							for(EXRRFNode r:ns){
								if(r.isLeaf()){
									//r.predValue += r.lastScore+"@"+ncou;
									r.predValue += "@"+ncou;
									ncou++;
								}
							}
							
							co.save(treefile);
						}
						//String treefile = "tree_simple."+ll+".dat";

						RandomForestProcess rg = new RandomForestProcess();
						rg.loadForestFromFile_SimpleFormat(treefile);
						
						
						HashMap<String,String> ss = new HashMap<>();
						ss.put("-hasheader", "true");
						ss.put("-hasnamecol", "true");
						ss.put("-hasanswercol", "true");
						ss.put("-infile", infile);
						
						HashMap<String,HashMap<String,Integer>> traincount = FuzzyDecisionTree.countTrueClass(ss,rg,false);
						HashMap<String,HashMap<String,Double>> trainratio = new HashMap<>();
						
						HashMap<String,Integer> countsum = new HashMap<>();//background ratio を計算するためだけに必要
						HashMap<String,Double> backgroundratio = new HashMap<>();//計算に使用したサンプルの中に占めるアミノ酸の割合
						int countall = 0;
						for(String s:traincount.keySet()){
							HashMap<String,Integer> mpp = traincount.get(s);
							for(String kk:mpp.keySet()){
								if(!countsum.containsKey(kk)){
									countsum.put(kk,0);
								}
								countsum.put(kk,countsum.get(kk)+mpp.get(kk));
								countall += mpp.get(kk);
							}
						}
						
						for(String s:countsum.keySet()){
							backgroundratio.put(s,countsum.get(s)/(double)countall);
						}
						
						for(String s:traincount.keySet()){
							HashMap<String,Integer> mpp = traincount.get(s);
							int gcount  = 0;
							for(String kk:mpp.keySet()){
								gcount += mpp.get(kk);
							}
							
							trainratio.put(s,new HashMap<String,Double>());
							//for(String kk:mpp.keySet()){
							for(String kk:countsum.keySet()){
								if(mpp.containsKey(kk)){
									double gratio = mpp.get(kk)/(double)gcount/backgroundratio.get(kk);
									trainratio.get(s).put(kk,gratio);
								}else{
									trainratio.get(s).put(kk,0.0);
								}
							}
							
							
						}
						
						
						
						
						
						double[] scc = getScores( traincount);
						System.out.println(treefile+"\t"+intestfile);
						ss.put("-infile", intestfile);
						//HashMap<String,HashMap<String,Integer>> testcount = CorDecisionTree.countTrueClass(ss,rg,true);
						HashMap<String,HashMap<String,Integer>> testcount = FuzzyDecisionTree.countTrueClass(ss,rg,false);
						
						double[] scc2 = {0.0,0.0,0.0};
						for(String s:testcount.keySet()){
							HashMap<String,Integer> mpp = testcount.get(s);
							for(String kk:mpp.keySet()){
								if(trainratio.get(s).get(kk) 
									< 1.0/Math.sqrt(2)){
									scc2[1] += mpp.get(kk);
								}
								if(trainratio.get(s).get(kk) 
									> Math.sqrt(2)){
									scc2[0] += mpp.get(kk);
								}
								scc2[2] +=  mpp.get(kk);
							}
						}
						
						all_posi = scc[0]+scc2[0];
						all_nega = scc[1]+scc2[1];
						double sscore = all_posi/(scc[2]+scc2[2]);
						sscore = scc2[0]/scc2[2];
						
						System.out.println("##cosinedist threshold\t"+coss
								+"\ttestsize\t"+sizz
								+"\tnodesample\t"+nsampp
								+"\tallposi\t"+all_posi
								+"\tallnega\t"+all_nega
								+"\ttrain\t"+String.valueOf( all_posi/(scc[2]+scc2[2]))+"\ttest\t"+String.valueOf(sscore)
						);
						//System.exit(0);
					}
				}
			}
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	public static void main__(String[] args){
		FuzzyDecisionTree co = FuzzyDecisionTree.loadTable("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\treegencheck\\iris_source_xx.dat", false,false);
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		co.extractTestDataset(0.1);
		updated.add(co.root);
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
			if(r.leafFlag){
				
			}else{
				co.calcOptimalSplit(r,SPLIT_TYPE_ENTROPY);
				if(!r.leafFlag){
					updated.add((EXRRFNode)r.childNodes[0]);
					updated.add((EXRRFNode)r.childNodes[1]);
				}else{
					r.predValue = co.index_to_label.get(r.getMajority(co.classNum));
				//	r.predValue ="";
				//	for(DirtySample ds:r.xsamples){
				//		r.predValue += co.index_to_label.get(ds.classIndex)+";";
				//	}
				}
			}
		}
		co.save("testout.tree.dat");
		
		RandomForestProcess rg = new RandomForestProcess();
		rg.loadForestFromFile_SimpleFormat("testout.tree.dat");
		
		HashMap<String,String> ss = new HashMap<>();
		ss.put("-hasanswercol", "true");
		ss.put("-infile", "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\treegencheck\\iris_chk.dat");
		ss.put("-outfile", "testres.dat");
		
		RandomForestProcess.normalPrediction(ss,rg);
		
		System.out.println(";;");
	}
}
class DirtySample{
	ArrayList<Double> values = new ArrayList<>();
	int targetColumn = 0;
	int classIndex = 0;
	int classGroupIndex = 0;//他のクラスト同じように扱った方が良い場合にまとめられる
	String name = "";
	String classLabel = "";
	boolean ignore = false;
	boolean prevignore =false;
	public void setTarget(int i){
		targetColumn = i;
	}
	public double getTargetValue(){
		return values.get(targetColumn);
	}
}
 class DirtySampleComparator implements Comparator<DirtySample>{
	@SuppressWarnings("unchecked")
	public int compare(DirtySample arg1, DirtySample arg2){
		
		if(arg1.getTargetValue() < arg2.getTargetValue()){
			return -1;
		}
		if(arg1.getTargetValue() == arg2.getTargetValue()){
			return 0;
		}
			return 1;
	}
	
}



class EXRRFNode extends RRFNode{//デフォルトでイコールである場合は（高いとみなして） 1 に進む
	/*RRFNode[] childNodes;
	RRFNode parentNode;
	int splitVarNum;
	boolean leafFlag;
	double splitThreshold;
	String predValue = null;
	boolean includeEquals_Up = true;
	*/
	ArrayList<DirtySample> xsamples = new ArrayList<>();
	ArrayList<DirtySample> testsamples = new ArrayList<>();
	double[] frequency = null;
	double[] frequency_test = null;
	double lastScore= 0;
	double loss = 0;
	double childLoss = 0;//指定 Depth まで下りた時点もしくは葉であるノードのスコアの合計
	ArrayList<SplittingVar> childSplit = null;
	
	
	
	public void calcPlayoutLoss(int dep){
		childLoss = 0;
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(this);
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
			if(r.isLeaf()){
				childLoss += r.loss;
			}else{
				if(r.depth == this.depth+dep){
					childLoss += r.loss;
				}else{
					updated.add((EXRRFNode)r.childNodes[0]);
					updated.add((EXRRFNode)r.childNodes[1]);
				}
			}
		}
	}
	
	
	public void collectChildSplit(int dep){
		childSplit = new ArrayList<SplittingVar>();
		
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(this);
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
			if(r.isLeaf()){
			}else{
				if(r != this){
					childSplit.add(new SplittingVar(r.splitVarNum,r.splitThreshold));
				}
				if(r.depth == this.depth+dep){
					
				}else{
					updated.add((EXRRFNode)r.childNodes[0]);
					updated.add((EXRRFNode)r.childNodes[1]);
				}
			}
		}
	}
	
	public ArrayList<EXRRFNode> getChildNodes(){
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(this);
		ArrayList<EXRRFNode> ret = new ArrayList<>();
		while(updated.size() > 0){
			EXRRFNode r = updated.remove(0);
				ret.add(r);
			if(r.isLeaf()){
			}else{
				updated.add((EXRRFNode)r.childNodes[0]);
				updated.add((EXRRFNode)r.childNodes[1]);

			}
		}
		ret.remove(0);
		return ret;
	}
		
	
	
	public void setScore(double d){
		loss = d;
	}
	public double getScore(){
		return loss;
	}
	public int getMajority(int classNum){
		
		int[] classcount = new int[classNum];
		int maxcount = 0;
		int maxindex = 0;
		for(int ii = 0;ii < classNum;ii++){
			classcount[ii] = 0;
		}
		for(DirtySample ds:xsamples){
			classcount[ds.classIndex]++;
			if(classcount[ds.classIndex] > maxcount){
				maxindex = ds.classIndex;
				maxcount = classcount[ds.classIndex];
			}
		}
		return maxindex;
	}
}


class SplittingVar{
	int varIndex = -1;
	double threshold = -1;
	SplittingVar(int i,double t){
		varIndex = i;
		threshold = t;
	}
}
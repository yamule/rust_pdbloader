/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author kimidori
 */
public class WeightedKMeansPSSM extends WeightedKMeansClustering{
	int classNum = 0;
	HashMap<Integer,String> index_to_label = new HashMap<>();
	HashMap<String,Integer> label_to_index = new HashMap<>();
	double[] backgroundFreq;
	double[] backgroundFreq_test;
	ArrayList<WKSampleData> testData = null;
	public static double logbottom = 0.0001;
	
	public static WeightedKMeansPSSM load(String filename,boolean header,boolean hasnamecol,boolean hasanswercol){
		WeightedKMeansClustering sc = WeightedKMeansClustering.loadTable(filename,header,hasnamecol,hasanswercol);
		WeightedKMeansPSSM ret = new WeightedKMeansPSSM();
		ret.clusters = sc.clusters;
		ret.columnname_to_index = sc.columnname_to_index;
		ret.samples = sc.samples;
		ret.vlist = sc.vlist;
		
		ret.countClassNum();
		return ret;
	}
	
	
	/**
	 * 二つに分けようとしたときに最も低くなる閾値とその際のエントロピーを返す。
	 * エントロピーが高い＝重視すべきでない特徴
	 * @param d
	 * @param targetcolumn
	 * @param frequency
	 * @return 
	 */
	
	
	public double[] calcMinEntropy(ArrayList<WKSampleData> d,int targetcolumn,double[] frequency){

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
		
		if(basicscore == 0){
			
			return null;
		}
		
		double minentropy = Double.MAX_VALUE;
		double minlowratio = Double.MAX_VALUE;
		double minthreshold = Double.MAX_VALUE;
		double minthreshold_start = 0;
		double minthreshold_end = 0;
		for(WKSampleData ds:d){
			ds.classGroupIndex = groupid[label_to_index.get(ds.classLabel)];
		}

		for(int ii = 0;ii < cdev;ii++){
			classcount[ii] = 0;
			classcount_lower[ii] = 0;
		}

		ArrayList<VSorter> sorter = new ArrayList<>();

		for(int ii = 0;ii < d.size();ii++){
			sorter.add(new VSorter(d.get(ii).pos[targetcolumn],ii));
		}
		Collections.sort(sorter,new VComparator());
		for(int jj = 0; jj < sorter.size()-1;jj++){
			WKSampleData sd = d.get(sorter.get(jj).index);
			WKSampleData sd_next = d.get(sorter.get(jj+1).index);
			classcount_lower[sd.classGroupIndex]++;
			if(sd.pos[targetcolumn] != sd_next.pos[targetcolumn]){

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
					minthreshold_end = sd.pos[targetcolumn]/2+sd_next.pos[targetcolumn]/2;
				}else if(ze < minentropy){
					minentropy = ze;
					minthreshold_start = sd.pos[targetcolumn]/2+sd_next.pos[targetcolumn]/2;
					minthreshold_end = minthreshold_start;
				}
			}
		}
		
		if(minentropy == Double.MAX_VALUE){
			return null;
		}
		double[] ret = new double[2];
		ret[0] = minentropy;
		ret[1] = minthreshold_end/2+minthreshold_start/2;
		return ret;
	}
	
	
	
	
	
	/**
	 * クラスの種類数を計算する。
	 * ロードするたびに呼び出す必要がある。
	 */
	public void countClassNum(){
		HashSet<String> cn = new HashSet<>();
		for(WKSampleData sd:samples){
			cn.add(sd.classLabel);
		}
		if(testData != null){
			for(WKSampleData sd:testData){
				cn.add(sd.classLabel);
			}
		}
		ArrayList<String> al = new ArrayList<>(cn);
		for(int ii = 0;ii < al.size();ii++){
			label_to_index.put(al.get(ii),ii);
			index_to_label.put(ii,al.get(ii));
		}
		classNum = cn.size();
		calcBackgroundFreq();
	}
	/**
	 * テスト用と、トレーニング用のサンプル内のクラスの割合を計算する。
	 * テストデータセットを分けたりしたとき呼び出す必要がある。
	 */
	public void calcBackgroundFreq(){
		backgroundFreq = calcFrequency_D(samples);
		if(testData != null){
			backgroundFreq_test = calcFrequency_D(testData);
		}
	}
	
	/**
	 * リスト内のサンプルをカウントして、バックグラウンドより低い割合であるサンプルの割合を返す。
	 * @param sal
	 * @param back
	 * @return 
	 */
	public double calcNegativeRatio(ArrayList<WKSampleData> sal,double[] back){
		double[] ra = calcFrequency_D(sal);
		double ret = 0;
		for(int ii = 0;ii < ra.length;ii++){
			if(ra[ii] < back[ii]){
				ret += ra[ii];
			}
		}
		return ret;
	}
	
	public double[] calcFrequency_D(Collection<WKSampleData> labels){
		double[] ret = new double[this.classNum];
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] = 0;
		}
		for(WKSampleData s:labels){
			ret[label_to_index.get(s.classLabel)] += 1;
		}
		
		for(int ii = 0;ii < this.classNum;ii++){
			ret[ii] /= labels.size();
		}
		return ret;
	}
	
	public void calcNegativeRatio(){
		double negativemax = 0;
		WKCluster maxcluster = null;
		double[] esum = new double[this.vlist.size()];
		int[] nullcount = new int[this.vlist.size()];
		for(int ii = 0;ii < esum.length;ii++){
			esum[ii] = 0;
			nullcount[ii] = 0;
		}
		double entropyall_sum = 0;
		double negativeratio_sum = 0;
		double maxentropy = 0;
		for(WKCluster c:clusters){
			double d = calcNegativeRatio(c.samples,backgroundFreq);
			negativeratio_sum += d;
			if(d > negativemax){
				maxcluster = c;
				negativemax = d;
			}
			if(d > 0){
				for(int ii = 0;ii < esum.length;ii++){
					double[] res = this.calcMinEntropy(c.samples, ii,this.backgroundFreq);
					if(res == null){
						esum[ii] += 1*d;
						entropyall_sum += 1*d;
						maxentropy = Math.max(maxentropy,esum[ii]);
					}else{
						esum[ii] += res[0]*d;
						entropyall_sum += res[0]*d;
						maxentropy = Math.max(maxentropy,esum[ii]);
					}
				}
			}
		}
		double me = 0.0;
		int mindex = 0;
		double wmax = 0;
		double xfactor  = 0.5;
		for(int ii = 0;ii < esum.length;ii++){
			vlist.get(ii).weight *= (1.0+xfactor/2-esum[ii]/maxentropy*xfactor);
			wmax = Math.max(wmax, vlist.get(ii).weight);
		}
		System.out.println("max weight:"+wmax);
		for(int ii = 0;ii < esum.length;ii++){
			vlist.get(ii).weight /= wmax;
		}
		
		for(WKCluster c:clusters){
			if(c.samples.size() == 0){
				continue;
			}
			
			
			double[] ra = calcFrequency_D(c.samples);
			int maxindex = 0;
			double maxratio = ra[0]/backgroundFreq[0];
			for(int ii = 1;ii < ra.length;ii++){
				if(ra[ii]/backgroundFreq[ii] > maxratio){
					maxindex = ii;
					maxratio = ra[ii]/backgroundFreq[ii];
				}
			}
			
			boolean[] dflag = new boolean[classNum];
			
			for(int ii = 0;ii < ra.length;ii++){
				if(ii == maxindex){
					dflag[ii] = true;
				}else{
					dflag[ii] = false;
				}
			}
			/* プラスになるやつは全部取った。
			for(int ii = 0;ii < ra.length;ii++){
				if(ra[ii] < backgroundFreq[ii]){
					dflag[ii] = true;
				}else{
					dflag[ii] = false;
				}
			}
			*/
			Iterator<WKSampleData> ite = c.samples.iterator();
			while(ite.hasNext()){
				WKSampleData w = ite.next();
				if(dflag[this.label_to_index.get(w.classLabel)]){
					ite.remove();
					w.setParent(null);
				}
			}
			c.updateMean();
		}
		/*
		int count = 0;
		
		while(update()){
			count++;
			System.out.println(count);
		}
		*/
		for(int jj = 0;jj < 100;jj++){
			for(WKCluster c:clusters){
				if(c.samples.size() == 0){
					continue;
				}
				double[] ra = calcFrequency_D(c.samples);
				int maxindex = 0;
				double maxratio = ra[0]/backgroundFreq[0];
				for(int ii = 1;ii < ra.length;ii++){
					if(ra[ii]/backgroundFreq[ii] > maxratio){
						maxindex = ii;
						maxratio = ra[ii]/backgroundFreq[ii];
					}
				}
				boolean[] dflag = new boolean[classNum];
				for(int ii = 0;ii < ra.length;ii++){
					if(ii == maxindex){
						dflag[ii] = true;
					}else{
						dflag[ii] = false;
					}
				}
				Iterator<WKSampleData> ite = c.samples.iterator();
				while(ite.hasNext()){
					WKSampleData w = ite.next();
					if(dflag[this.label_to_index.get(w.classLabel)]){
						ite.remove();
						w.setParent(null);
					}
				}
				c.updateMean();
			}
			
			if(!update()){
				break;
			}
			System.out.print("=");
		}
		System.out.println("");
		double negativeratio_sum_new = 0;
		for(WKCluster c:clusters){
			double d = calcNegativeRatio(c.samples,backgroundFreq);
			negativeratio_sum_new += d;
		}
		System.out.println("updated:"+negativeratio_sum+"->"+negativeratio_sum_new);
	}
	public static void main(String[] args){
		String infile = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_onlycb_jumpcb45\\input.0.train.dat";
		WeightedKMeansPSSM wc = WeightedKMeansPSSM.load(infile, true,true,true);
		wc.setClusterNum(60);
		for(int ii = 0;ii < 60;ii++){
			wc.clusters.get(ii).addSample(wc.samples.get(ii));
		}
		wc.prepare();
		wc.update();
		for(int jj = 0;jj < 100;jj++){
			wc.calcNegativeRatio();
			for(int ii = 0;ii < wc.clusters.size();ii++){
				System.out.println(ii+";"+wc.clusters.get(ii).samples.size()+";"+wc.calcNegativeRatio(wc.clusters.get(ii).samples,wc.backgroundFreq));
			}
		}
	}
	
	
}
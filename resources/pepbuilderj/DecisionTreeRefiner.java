/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Pattern;
import static pepbuilderj.FeatureGeneratorCB_20180501.listToArray;

/**
 *
 * @author kimidori
 */
public class DecisionTreeRefiner {
	public static final int POSI = 0;
	public static final int NEGA = 1;
	public static final int ZERO = 2;
	public static final int LOSS_GINI = 0;
	public static final int LOSS_POSINEGA = 1;
	public static final int LOSS_MINUSPOSI = 2;
	public static final int LOSS_POSINEGAPLUS = 3;//POSITIVE-NEGATIVE が最大化するように
	public static final int LOSS_POSINEGAPLUS_WITHFLAG = 4;//POSITIVE-NEGATIVE が最大化するように
	int lossType = LOSS_GINI;
	FuzzyDecisionTree tree;
	HashMap<String,Double> backgroundFreq = new HashMap<>();

	HashMap<EXRRFNode,Integer> original_samplenum = new HashMap<>();
	DecisionTreeRefiner(String tablefile,boolean header,boolean hasnamecol){
		//特徴変数等のマップがあるのでダミーでも必要。
		tree = FuzzyDecisionTree.loadTable(tablefile,header,hasnamecol);
		calcBackgroundFreq();
	}
	public void loadTreeStructure(String treefile){
		RandomForestProcess rg = new RandomForestProcess();
			rg.loadForestFromFile_SimpleFormat(treefile);
		RRFTree rftree = rg.trees.get(0);
		
		ArrayList<RRFNode> rfnode = new ArrayList<>();	
		ArrayList<RRFNode> updated = new ArrayList<>();
		updated.add(rftree.nodes[0]);
		
		HashMap<RRFNode,EXRRFNode> nnmap = new HashMap<>();
		while(updated.size() > 0){
			RRFNode n = updated.remove(0);
			EXRRFNode e = new EXRRFNode();
			nnmap.put(n,e);
			e.splitVarNum = n.splitVarNum;
			e.leafFlag = n.leafFlag;
			e.splitThreshold = n.splitThreshold;
			e.predValue = n.predValue;
			e.includeEquals_Up = n.includeEquals_Up;
	
			rfnode.add(n);
			if(n.childNodes == null){
				continue;
			}
			for(int ii = 0;ii < n.childNodes.length;ii++){
				if(n.childNodes[ii] != null){
					updated.add(n.childNodes[ii]);
				}
			}
		}
		tree.root = nnmap.get(rftree.nodes[0]);
		for(RRFNode n:rfnode){
			EXRRFNode e = nnmap.get(n);
			if(n.parentNode != null){
				e.parentNode = nnmap.get(n.parentNode);
			}
			if(n.childNodes == null){
				continue;
			}else{
				e.childNodes = new EXRRFNode[2];
			}
			
			for(int ii = 0;ii < n.childNodes.length;ii++){
				if(n.childNodes[ii] != null){
					e.childNodes[ii] = nnmap.get(n.childNodes[ii]);
				}
			}
		}
	}
	
	public void reloadSamples(String tablefile,boolean header,boolean hasnamecol){
		FuzzyDecisionTree dtree = FuzzyDecisionTree.loadTable(tablefile,header,hasnamecol);
		tree.samples = dtree.samples;
		tree.root.xsamples = tree.samples;
		calcBackgroundFreq();
	}
	
	public void formatSamples(){
		tree.root.xsamples = tree.samples;
		remapSamples(tree.root);
		
		ArrayList<EXRRFNode> leafs = getAllLeafNode(tree.root);
		for(EXRRFNode l:leafs){
			original_samplenum.put(l,l.xsamples.size());
		}
		calcBackgroundFreq();
		flushStructure();
	}
	
	public void remapSamples(){
		tree.root.xsamples = tree.samples;
		remapSamples(tree.root);
	}
	
	public void calcBackgroundFreq(){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:tree.samples){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		this.backgroundFreq = new HashMap<>();
		for(String l:count.keySet()){
			this.backgroundFreq.put(l,count.get(l)/(double)this.tree.samples.size());
		}
	}
	/**
	 * startnode のサンプルを用いて startnode より下流のノードについてサンプルを更新する
	 * @param startnode 
	 */
	public void remapSamples(EXRRFNode startnode){
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(startnode);
		while(updated.size() > 0){
			EXRRFNode parentnode = updated.remove(0);
			if(parentnode.childNodes != null){
				EXRRFNode r1 = (EXRRFNode)parentnode.childNodes[0];
				EXRRFNode r2 = (EXRRFNode)parentnode.childNodes[1];
				if(parentnode.xsamples != null){
					ArrayList<ArrayList<DirtySample>> dal = tree.split(parentnode.xsamples,parentnode.splitVarNum
							,parentnode.splitThreshold);
					r1.xsamples = dal.get(0);
					r2.xsamples = dal.get(1);
				}
				updated.add(r1);
				updated.add(r2);
			}
		}
	}
	
	
	public String predict(double[] feature){
		RRFNode nex = tree.root;
		while(nex != null){
			if(nex.childNodes != null){
				nex = nex.getNextNode(feature);
			}else{
				return nex.getPredValue();
			}
		}
		return null;
	}
	/**
	 * Child がない Node のリストを返す。
	 * @param startnode
	 * @return 
	 */
	public ArrayList<EXRRFNode> getAllLeafNode(EXRRFNode startnode){
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		ArrayList<EXRRFNode> ret = new ArrayList<>();
		updated.add(startnode);
		while(updated.size() > 0){
			EXRRFNode parentnode = updated.remove(0);
			if(parentnode.childNodes != null){
				EXRRFNode r1 = (EXRRFNode)parentnode.childNodes[0];
				EXRRFNode r2 = (EXRRFNode)parentnode.childNodes[1];
				if(parentnode.xsamples != null){
					ArrayList<ArrayList<DirtySample>> dal = tree.split(parentnode.xsamples,parentnode.splitVarNum
							,parentnode.splitThreshold);
					r1.xsamples = dal.get(0);
					r2.xsamples = dal.get(1);
				}
				updated.add(r1);
				updated.add(r2);
			}else{
				ret.add(parentnode);
			}
		}
		return ret;
	}
	
	public void debug_predsample(){
		HashMap<String,Integer> allanswers = new HashMap<>();//回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		for(DirtySample d:tree.root.xsamples){
			String ans = d.classLabel;
			String pred = predict(listToArray(d.values));
			if(!allanswers.containsKey(ans)){
				allanswers.put(ans,0);
			}
			if(!pred_answers.containsKey(pred)){
				pred_answers.put(pred,new HashMap<String,Integer>());
			}
			if(!pred_answers.get(pred).containsKey(ans)){
				pred_answers.get(pred).put(ans,0);
			}
			allanswers.put(ans,1+allanswers.get(ans));
			pred_answers.get(pred).put(ans,1+pred_answers.get(pred).get(ans));
		}

		ArrayList<String> lines = new ArrayList<>();
		int allcount = 0;
		for(String s:allanswers.keySet()){
			lines.add("#count "+s+":"+allanswers.get(s));
			allcount += allanswers.get(s);
		}
		lines.add("#count all:"+allcount);
		
		for(String s:pred_answers.keySet()){
			HashMap<String,Integer> hs = pred_answers.get(s);
			ArrayList<VSorter> al = new ArrayList<>();
			ArrayList<String> labels = new ArrayList<>(hs.keySet());
			for(int ii = 0;ii < labels.size();ii++){
				al.add(new VSorter(hs.get(labels.get(ii)),ii));
			}
			Collections.sort(al,new VComparator());
			Collections.reverse(al);
			StringBuffer sb = new StringBuffer();
			sb.append("trueatom_predictedas "+s+":");
			for(VSorter v:al){
				sb.append("\t"+labels.get(v.index)+"\t"+v.val);
			}
			lines.add(sb.toString());
		}	
	}
	
	

	public int calcPositiveSampleNum(ArrayList<DirtySample> al){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:al){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		int ret = 0;
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)al.size();
			//どこまでをネガティブとするか・・・
			if(freq > this.backgroundFreq.get(l)){
				ret++;
			}
		}
		return ret;
	}
	
	public int calcNegativeSampleNum(ArrayList<DirtySample> al){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:al){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		int ret = 0;
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)al.size();
			//どこまでをネガティブとするか・・・
			if(freq < this.backgroundFreq.get(l)){
				ret++;
			}else{
				System.out.println(freq);
			}
		}
		return ret;
	}
	
	public double[] calcLoss(ArrayList<DirtySample> lis,int typ){
		if(typ == LOSS_POSINEGAPLUS){
			double sc[] = calcPosiNegaSampleNum(lis);
			double ret[] = new double[3];
			ret[POSI] = sc[POSI];
			ret[NEGA] = sc[POSI]*-1+sc[NEGA];
			ret[ZERO] = sc[ZERO];
			return ret;
		}else if(typ == LOSS_POSINEGAPLUS_WITHFLAG){
			double sc[] = calcPosiNegaSampleNum_WithFlag(lis);
			double ret[] = new double[3];
			ret[POSI] = sc[POSI];
			ret[NEGA] = sc[POSI]*-1+sc[NEGA];
			ret[ZERO] = sc[ZERO];
			return ret;
		}else if(typ == LOSS_MINUSPOSI){
			double sc[] = calcPosiNegaSampleNum(lis);
			double ret[] = new double[3];
			ret[POSI] = sc[NEGA]*-1;
			ret[NEGA] = sc[POSI]*-1;
			ret[ZERO] = sc[ZERO];
			return ret;
		}else if(typ == LOSS_POSINEGA){
			return calcPosiNegaSampleNum(lis);
		}else if(typ == LOSS_GINI){
			double sc[] = calcGiniInpurity(lis);
			
			return sc;
		}
		return null;
	}
	
	
	
	/**
	 * この関数から返る値を最大化しようとする
	 * @param al
	 * @return 
	 */
	public double[] calcPosiNegaSampleNum(ArrayList<DirtySample> al){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:al){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		int ret_posi = 0;
		int ret_nega = 0;
		int ret_zero = 0;
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)al.size();
			//どこまでをネガティブとするか・・・
			/*
			if(freq < this.backgroundFreq.get(l)){
				ret_nega += count.get(l);
			}else if(freq > this.backgroundFreq.get(l)){
				ret_posi += count.get(l);
			}else{
				ret_zero += count.get(l);
			}
			*/
			//四捨五入でマイナス1になるかどうかで判断
			if((freq/this.backgroundFreq.get(l)) < 1.0/Math.sqrt(2)){
				ret_nega += count.get(l);
			}else if(freq/this.backgroundFreq.get(l) > Math.sqrt(2)){
				ret_posi += count.get(l);
			}else{
				ret_zero += count.get(l);
			}
		}
		double[] ret = {ret_posi,ret_nega,ret_zero};
		return ret;
	}
	/**
	 * この関数から返る値を最大化しようとする
	 * @param al
	 * @return 
	 */
	public double[] calcPosiNegaSampleNum_WithFlag(ArrayList<DirtySample> al){
		HashMap<String,Integer> count = new HashMap<>();
		HashMap<String,Integer> count_include = new HashMap<>();
		for(DirtySample ds:al){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
				count_include.put(z,0);
			}
			count.put(z,count.get(z)+1);
			if(!ds.ignore){
				count_include.put(z,count_include.get(z)+1);
			}
		}
		int ret_posi = 0;
		int ret_nega = 0;
		int ret_zero = 0;
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)al.size();
			//四捨五入でマイナス1になるかどうかで判断
			if((freq/this.backgroundFreq.get(l)) < 1.0/Math.sqrt(2)){
				ret_nega += count_include.get(l);
			}else if(freq/this.backgroundFreq.get(l) > Math.sqrt(2)){
				ret_posi += count_include.get(l);
			}else{
				ret_zero += count_include.get(l);
			}
		}
		double[] ret = {ret_posi,ret_nega,ret_zero};
		return ret;
	}
	/**
	 * Positive もしくは Zero の奴を無視するフラグを立てる
	 * @param trainingdataset 
	 */
	public void setPositiveZeroIgnore(ArrayList<DirtySample> trainingdataset){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:trainingdataset){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		HashSet<String> ispositive = new HashSet<>();
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)trainingdataset.size();
			//四捨五入でマイナス1になるかどうかで判断
			if((freq/this.backgroundFreq.get(l)) < 1.0/Math.sqrt(2)){
			}else if(freq/this.backgroundFreq.get(l) > Math.sqrt(2)){
				ispositive.add(l);
			}else{
				ispositive.add(l);
			}
		}
		for(DirtySample ds:trainingdataset){
			String z = ds.classLabel;
			if(!ispositive.contains(z)){
				ds.ignore = false;
			}else{
				ds.ignore = true;
			}
		}
	}
	/**
	 * Positive もしくは Zero と判断されるラベルを返す
	 * @param trainingdataset 
	 */
	public HashMap<String,Double> getLabelPosiNega(ArrayList<DirtySample> trainingdataset){
		HashMap<String,Integer> count = new HashMap<>();
		HashMap<String,Double> ret = new HashMap<>();
		
		for(DirtySample ds:trainingdataset){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)trainingdataset.size();
			//四捨五入でマイナス1になるかどうかで判断
			if((freq/this.backgroundFreq.get(l)) < 1.0/Math.sqrt(2)){
				ret.put(l,-1.0);
			}else if(freq/this.backgroundFreq.get(l) > Math.sqrt(2)){
				ret.put(l,1.0);
			}else{
				ret.put(l,0.0);
			}
		}
		return ret;
	}
	
	/**
	 * leaf の xsamples に入っている sample について、予測されたラベルと正答ラベルの数を数えて String 化し返す。
	 * @param leafs
	 * @return 
	 */
	public static ArrayList<String> countClasses(ArrayList<EXRRFNode> leafs){
		HashMap<String,Integer> allanswers = new HashMap<>();//回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		for(EXRRFNode l:leafs){
			if(pred_answers.containsKey(l.predValue)){
				throw new RuntimeException("?? answer class duplication is not supported.");
			}
			HashMap<String,Integer> pcou = new HashMap<>();
			
			for(DirtySample d:l.xsamples){
				if(!pcou.containsKey(d.classLabel)){
					pcou.put(d.classLabel,0);
				}
				if(!allanswers.containsKey(d.classLabel)){
					allanswers.put(d.classLabel,0);
				}
				allanswers.put(d.classLabel,allanswers.get(d.classLabel)+1);
				pcou.put(d.classLabel,pcou.get(d.classLabel)+1);
			}
			pred_answers.put(l.predValue,pcou);
		}
		
	
		ArrayList<String> lines = new ArrayList<>();
		int allcount = 0;
		for(String s:allanswers.keySet()){
			lines.add("#count "+s+":"+allanswers.get(s));
			allcount += allanswers.get(s);
		}
		lines.add("#count all:"+allcount);
		for(String s:pred_answers.keySet()){
			HashMap<String,Integer> hs = pred_answers.get(s);
			ArrayList<VSorter> al = new ArrayList<>();
			ArrayList<String> labels = new ArrayList<>(hs.keySet());
			for(int ii = 0;ii < labels.size();ii++){
				al.add(new VSorter(hs.get(labels.get(ii)),ii));
			}
			Collections.sort(al,new VComparator());
			Collections.reverse(al);
			StringBuffer sb = new StringBuffer();
			sb.append("trueatom_predictedas "+s+":");
			for(VSorter v:al){
				sb.append("\t"+labels.get(v.index)+"\t"+v.val);
			}
			lines.add(sb.toString());
		}
		return lines;
	}
	public double[] calcGiniInpurity(ArrayList<DirtySample> al){
		HashMap<String,Integer> count = new HashMap<>();
		for(DirtySample ds:al){
			String z = ds.classLabel;
			if(!count.containsKey(z)){
				count.put(z,0);
			}
			count.put(z,count.get(z)+1);
		}
		double ret_posi = 0;
		double ret_nega = 0;
		double ret_zero = 0;
		for(String l:count.keySet()){
			double freq = count.get(l)/(double)al.size();
			ret_posi += freq*freq;
			ret_nega += 1.0-freq*freq;
		}
		double[] ret = {ret_posi,ret_nega,ret_zero};
		return ret;
	}
	public void thresholdShift(EXRRFNode e){
		int itarget = e.splitVarNum;
		int startvar = e.splitVarNum;
		
		
		ArrayList<DirtySample> ds = e.xsamples;
		for(DirtySample d:ds){
			d.setTarget(itarget);
		}
		
		ArrayList<EXRRFNode> leafs = getAllLeafNode(e);
		int negativenum = 0;
		int positivenum = 0;
		int zeronum = 0;
		//元の状態のスコア
		for(EXRRFNode l:leafs){
			double[] rs = calcLoss(l.xsamples,lossType);
			negativenum += rs[NEGA];
			positivenum += rs[POSI];
			zeronum += rs[ZERO];
		}
		boolean randomvarchange = false;
//System.out.println("start:"+negativenum+";"+positivenum+";"+zeronum);
		if(Math.random() < 0.1){
			itarget = (int)(e.xsamples.get(0).values.size()*Math.random());
			e.setSplitVarNum(itarget);
			randomvarchange = true;
		}
		
		int uindex = 0;
		HashSet<Double> unique_hs = new HashSet<>();
		for(DirtySample d:ds){
			d.setTarget(itarget);
			unique_hs.add(d.getTargetValue());
		}
		Collections.sort(ds,new DirtySampleComparator());
		ArrayList<Double> unique = new ArrayList<>(unique_hs);
		Collections.sort(unique);
		
		
		for(int ii = 0;ii < unique.size()-1;ii++){//equals up によって違いそうだが面倒なので。。。
			if(unique.get(ii+1) >= e.splitThreshold 
				&& unique.get(ii) <= e.splitThreshold){
				uindex = ii;//Threshold の上下で別れる Index の下側
				break;
			}
		}
		double start_thresh = e.splitThreshold;
		double minthreshold =  e.splitThreshold;
		int minnegative = negativenum;
		int drand = Math.max((int)(unique.size()*0.5),1);
		
		//プラス方向へのシフト
		for(int ii = uindex;ii < unique.size()-1;ii+=drand){
			if(ii > unique.size() - 1){
				break;
			}
			double tthresh = unique.get(ii)/2+unique.get(ii+1)/2;//+1 からで良い気はするが・・・。
			e.setSplitThreshold(tthresh);
			remapSamples(e);
			int dnegativenum = 0;
			int dpositivenum = 0;
			
			boolean underlimitflag = false;//サンプル数小さくなりすぎ
			for(EXRRFNode l:leafs){
				double[] rs = calcLoss(l.xsamples,lossType);
				dnegativenum += rs[NEGA];
				dpositivenum += rs[POSI];
				if(l.xsamples.size() < this.original_samplenum.get(l)/2){
					underlimitflag = true;
				}
			}
			//System.out.println("shift:"+negativenum+"->"+dnegativenum);
			//System.out.println("shift:"+positivenum+"->"+dpositivenum);
			
			if(!underlimitflag && dnegativenum <= minnegative){
				minnegative = dnegativenum;
				minthreshold = tthresh;
			}else{
				if(!randomvarchange){
					break;
				}
			}
		}
		double minthreshold2 =  e.splitThreshold;
		int minnegative2 = negativenum;
		//マイナス方向へのシフト
		for(int ii = uindex;ii > 0;ii-=drand){
			if(ii < 1){
				break;
			}
			double tthresh = unique.get(ii)/2+unique.get(ii-1)/2;
			e.setSplitThreshold(tthresh);
			remapSamples(e);
			int dnegativenum = 0;
			int dpositivenum = 0;
			boolean underlimitflag = false;//サンプル数小さくなりすぎ
			for(EXRRFNode l:leafs){
				double[] rs = calcLoss(l.xsamples,lossType);
				dnegativenum += rs[NEGA];
				dpositivenum += rs[POSI];
				if(l.xsamples.size() < this.original_samplenum.get(l)/2){
					underlimitflag = true;
				}
			}
			//System.out.println("shift:"+negativenum+"->"+dnegativenum);
			//System.out.println("shift:"+positivenum+"->"+dpositivenum);
			
			if(!underlimitflag && dnegativenum <= minnegative2){
				minnegative2 = dnegativenum;
				minthreshold2 = tthresh;
			}else{
				if(!randomvarchange){
					break;
				}
			}
		}
		if(minnegative2 < negativenum || minnegative < negativenum){
			e.setSplitVarNum(itarget);
			if(minnegative2 < minnegative){
				e.setSplitThreshold(minthreshold2);
			}else{
				e.setSplitThreshold(minthreshold);
			}
			
		}else{
			e.setSplitVarNum(startvar);
			e.setSplitThreshold(start_thresh);
		}
		remapSamples(e);
		
	}
	
	
	public ArrayList<EXRRFNode> getAllNodes(EXRRFNode startnode){
		ArrayList<EXRRFNode> ret = new ArrayList<>();
		ArrayList<EXRRFNode> updated = new ArrayList<>();
		updated.add(startnode);
		while(updated.size() > 0){
			EXRRFNode parentnode = updated.remove(0);
			ret.add(parentnode);
			if(parentnode.childNodes != null){
				EXRRFNode r1 = (EXRRFNode)parentnode.childNodes[0];
				EXRRFNode r2 = (EXRRFNode)parentnode.childNodes[1];
				
				updated.add(r1);
				updated.add(r2);
			}
		}
		return ret;
	}
	
	public static void printStrings(String outfilename,ArrayList<String> al){
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
			new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 
			for(String s:al){
				pw.write(s+"\n");
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public void flushStructure(){
		ArrayList<EXRRFNode> al = getAllNodes(tree.root); 
		for(EXRRFNode l:al){
			l.asLeaf(false);
		}
		ArrayList<EXRRFNode> leafs = getAllLeafNode(tree.root);
		for(EXRRFNode l:leafs){
			l.asLeaf(true);
		}
	}
	public static void main(String args[]){
		int samplenum_min_split = 200;
		int samplenum_min = 200;
		
		for(int tt = 0;tt < 20;tt++){
			HashSet<DirtySample> negative_prev = new HashSet<>();
			//DecisionTreeRefiner dtr = new DecisionTreeRefiner("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_cb\\input."+tt+".train.dat",true,true);
			//dtr.loadTreeStructure("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_cb\\tree_simple."+tt+".dat_"+String.valueOf(tt));
			//DecisionTreeRefiner dtr = new DecisionTreeRefiner("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\input."+tt+".train.dat",true,true);
			//dtr.loadTreeStructure("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\tree_simple."+tt+".dat_"+String.valueOf(tt));
			//DecisionTreeRefiner dtr = new DecisionTreeRefiner("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_cb\\input."+tt+".test.dat",true,true);
			//dtr.loadTreeStructure("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_cb\\tree_simple."+tt+".dat_"+String.valueOf(tt));
			

			DecisionTreeRefiner dtr = new DecisionTreeRefiner("input."+tt+".train.dat",true,true);
			dtr.loadTreeStructure("tree_simple."+tt+".dat_"+String.valueOf(tt));
			dtr.formatSamples();
			ArrayList<DirtySample> original_samples = new ArrayList<>();
			original_samples.addAll(dtr.tree.samples);
			for(int booster = 0;booster < 10;booster++){
				
				ArrayList<EXRRFNode> al = dtr.getAllNodes(dtr.tree.root); 
				ArrayList<EXRRFNode> leafs = dtr.getAllLeafNode(dtr.tree.root);
				if(booster > 0){
					Iterator<DirtySample> ite = dtr.tree.samples.iterator();
					while(ite.hasNext()){
						DirtySample ds = ite.next();
						if(ds.ignore){
							ite.remove();
						}
					}
					dtr.formatSamples();
					dtr.tree.cosineSimilarityThreshold = 0;
					dtr.tree.identityThreshold = 0.62;
					dtr.tree.remainedSampleThreshold = samplenum_min_split*2.0/(double)dtr.tree.samples.size();
					dtr.tree.makeTree(dtr.tree.root);
					dtr.flushStructure();
					al = dtr.getAllNodes(dtr.tree.root); 
					leafs = dtr.getAllLeafNode(dtr.tree.root);
					int lcou = 0;
					for(EXRRFNode l:leafs){
						int zb = l.getMajority(dtr.tree.classNum);
						if(dtr.tree.index_to_label.get(zb) == null){
							System.out.println("???");
						}
						l.asLeaf(true);
						l.predValue = dtr.tree.index_to_label.get(zb)+"@"+lcou;
						lcou++;
					}
					//dtr.tree.samples.clear();
					//dtr.tree.samples.addAll(original_samples);
					//dtr.formatSamples();
				}
				double negativenum = 0;
				double positivenum = 0;
				double zeronum = 0;
				//dtr.lossType = LOSS_POSINEGA;
				//dtr.lossType = LOSS_MINUSPOSI;
				dtr.lossType = LOSS_POSINEGAPLUS_WITHFLAG;
				//dtr.lossType = LOSS_GINI;
				for(EXRRFNode l:leafs){
					double[] rs = dtr.calcLoss(l.xsamples,dtr.lossType);
					negativenum += rs[NEGA];
					positivenum += rs[POSI];
					zeronum += rs[ZERO];
				}
				System.out.println(tt+";root_start:"+negativenum+";"+positivenum+";"+zeronum);
				printStrings("reftree."+booster+"."+tt+".start.dat.scores.source",countClasses(leafs));
				
				dtr.tree.remainedSampleThreshold = -1;
				ArrayList<EXRRFNode> updated = new ArrayList<>();
				updated.addAll(leafs);

				while(updated.size() > 0){
					EXRRFNode l = updated.remove(0);
					if(l.xsamples.size() > samplenum_min_split){
						dtr.tree.calcOptimalSplit(l,FuzzyDecisionTree.SPLIT_TYPE_GINI);
						if(l.childNodes != null){
							if(((EXRRFNode)l.childNodes[1]).xsamples.size() > samplenum_min && 
									((EXRRFNode)l.childNodes[0]).xsamples.size() > samplenum_min){
								
								l.asLeaf(false);
								l.predValue = "";
								updated.add((EXRRFNode)l.childNodes[0]);
								updated.add((EXRRFNode)l.childNodes[1]);
							}else{
								l.asLeaf(true);
								l.childNodes = null;
							}
						}
					}
				}
				dtr.flushStructure();
				al = dtr.getAllNodes(dtr.tree.root); 
				leafs = dtr.getAllLeafNode(dtr.tree.root);
				
				int lcou = 0;
				for(EXRRFNode l:leafs){
					int zb = l.getMajority(dtr.tree.classNum);
					if(dtr.tree.index_to_label.get(zb) == null){
						System.out.println("???");
					}
					l.asLeaf(true);
					l.predValue = dtr.tree.index_to_label.get(zb)+"@"+lcou;
					lcou++;
				}
				
				dtr.tree.save("reftree."+booster+"."+tt+".start.dat");
				printStrings("reftree."+booster+"."+tt+".start.dat.scores",countClasses(leafs));
				dtr.formatSamples();
				al = dtr.getAllNodes(dtr.tree.root); 
				leafs = dtr.getAllLeafNode(dtr.tree.root);
				printStrings("reftree."+booster+"."+tt+".start.dat.scores.chk",countClasses(leafs));
				
				dtr.loadTreeStructure("reftree."+booster+"."+tt+".start.dat");
				dtr.formatSamples();
				al = dtr.getAllNodes(dtr.tree.root); 
				leafs = dtr.getAllLeafNode(dtr.tree.root);
				printStrings("reftree."+booster+"."+tt+".start.dat.scores.chk2",countClasses(leafs));
				
				
				int iternum = 50;
				if(booster > 0){
					iternum = 500;
				}
				
				for(int jj = 0;jj < iternum;jj++){
					for(int ii = 0;ii < 1000;ii++){
						int id = (int)(Math.random()*al.size());
						EXRRFNode ex = al.get(id);
						if(ex.childNodes != null){
							dtr.thresholdShift(ex);
						}
					}
					System.out.print(jj+"========;\t");

					double negativenum2 = 0;
					double positivenum2 = 0;
					double zeronum2 = 0;
					for(EXRRFNode l:leafs){
						double[] rs = dtr.calcLoss(l.xsamples,dtr.lossType);
						negativenum2 += rs[NEGA];
						positivenum2 += rs[POSI];
						zeronum2 += rs[ZERO];
					}
					System.out.println("root_end:"+negativenum2+";"+positivenum2+";"+zeronum2);
					if(jj%100 == 0){
						dtr.tree.save("reftree."+tt+"."+booster+"."+jj+".dat");
						printStrings("reftree."+tt+"."+booster+"."+jj+".dat.scores",countClasses(leafs));
					}
					//LOSS_POSINEGAPLUS の時だけ
					//if(dtr.lossType == LOSS_POSINEGAPLUS_WITHFLAG &&  positivenum2/(negativenum2+positivenum2*2+zeronum) > 0.8){
					//	break;
					//}
				}
				
				for(DirtySample d:original_samples){
					d.classLabel = d.classLabel.replaceAll("^;;+","");
				}
				dtr.formatSamples();
				
				
				dtr.tree.save("reftree."+tt+"."+booster+".fin.dat");
				printStrings("reftree."+tt+"."+booster+".fin.dat.scores",countClasses(leafs));
	
				HashMap<String,HashMap<String,Double>> predres = new HashMap<>();
				for(EXRRFNode l:leafs){//この木でポジティブが出るラベルでないラベルを記録
					predres.put(l.getPredValue(),dtr.getLabelPosiNega(l.xsamples));
				}
				
				dtr.tree.samples.clear();//pred 関数を作るのが面倒なので全部マップしなおす
				dtr.tree.samples.addAll(original_samples);
				dtr.formatSamples();
				
				int negacount = 0;
				int prevnega_nega = 0;//前回も negative で 今回も negative になったもの
				int prevnega_posi = 0;//前回も negative で 今回 positive になったもの
				int posicount = 0;
				int zerocount = 0;
				HashSet<DirtySample> currentnega = new HashSet<>();
				HashSet<DirtySample> currentzero = new HashSet<>();
				for(EXRRFNode l:leafs){
					for(DirtySample d:l.xsamples){
						double res = -1;
						if(predres.get(l.getPredValue()).containsKey(d.classLabel)){
						//ここで nullpo が出ると、コードをアップデートしたせいで使えなくなっている。	
							res = predres.get(l.getPredValue()).get(d.classLabel);
						}
						if(res < -0.5){//1.0, 0.0, -1.0 しか入らないようにしているが。。。
							negacount++;
							if(negative_prev.contains(d)){
								prevnega_nega++;
							}
							currentnega.add(d);
							
						}else if(res > 0.5){
							posicount++;
							if(negative_prev.contains(d)){
								prevnega_posi++;
							}
						}else{
							currentzero.add(d);
							zerocount++;
						}
						if(res >= -0.000001){
							d.ignore = true;
						}else{
							d.ignore = false;
						}
					}
				}
				negative_prev = currentnega;
				for(DirtySample d:original_samples){
					if(!d.ignore){
						//付けると GiNI 計算の時も別物として扱われて上手いこと別れるかと思ったが
						//成績は悪かった。
					//	d.classLabel = ";;"+d.classLabel;
					}else{
						if(Math.random() < currentnega.size()/(double)original_samples.size()/2.0){
							d.ignore = false;
						}
					}
					//if(currentzero.contains(d)){//zero は判断つかないのでのける
					//	d.ignore = true;
					//}else{
					//	d.ignore = true;
					//	d.ignore = false;
					//}
				}
				
				System.out.print(booster+"\tprevnega_negatives:"+prevnega_nega+"\t");
				System.out.print("prevnega_positives:"+prevnega_posi+"\t");
				System.out.print("negatives:"+negacount+"\t");
				System.out.print("positives:"+posicount+"\t");
				System.out.print("zero:"+zerocount+"\n");
				if(negacount < 2000){//tekitou
					break;
				}
			}
		}
	}
	
}

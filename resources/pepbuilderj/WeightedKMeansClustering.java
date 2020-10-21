/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author kimidori
 */
public class WeightedKMeansClustering {
	ArrayList<WKVariable> vlist = new ArrayList<>();
	ArrayList<WKSampleData> samples = new ArrayList<>();
	ArrayList<WKCluster> clusters = new ArrayList<>();
	HashMap<String,Integer> columnname_to_index = new HashMap<>();
	
	public static WeightedKMeansClustering loadTable(String filename,boolean header,boolean hasnamecol,boolean hasanswercol){
		WeightedKMeansClustering ret = new WeightedKMeansClustering();
		try{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			ArrayList<String> vlistname = new ArrayList<>();
			String line = null;
			if(header){
				line = br.readLine();
				ArrayList<String> al  = getParsedArray_String(line);
				if(hasnamecol){
					al.remove(0);
				}
				for(int ii = 0;ii < al.size();ii++){
					ret.columnname_to_index.put(al.get(ii), ii);
					vlistname.add(al.get(ii));
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
				if(hasanswercol){
					answer = al.remove(al.size()-1);
				}
				ArrayList<Double> values = new ArrayList<>();
				answers.add(answer);
				for(int ii = 0;ii < al.size();ii++){
					values.add(Double.parseDouble(al.get(ii)));
				}
				double d[] = new double[values.size()];
				for(int ii = 0;ii < values.size();ii++){
					d[ii] = values.get(ii);
				}
				
				WKSampleData ds = new WKSampleData(sname,answer,d);
				ret.samples.add(ds);
			}
			if(vlistname.size() == 0){
				for(int ii = 0;ii < ret.samples.get(0).pos.length;ii++){
					vlistname.add("V"+ii);
				}
			}
			for(int ii = 0;ii < ret.samples.get(0).pos.length;ii++){
				WKVariable v = new WKVariable();
				v.columnName = vlistname.get(ii);
				v.index = ii;
				ret.vlist.add(v);
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	
	public void setClusterNum(int c){
		clusters.clear();
		for(int ii = 0;ii < c;ii++){
			clusters.add(new WKCluster());
		}
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
	public WKCluster calcNearestCluster(double[] a){
		double ndist = Double.MAX_VALUE;
		WKCluster ret = null;
		for(int ii = 0;ii < clusters.size();ii++){
			WKCluster c = clusters.get(ii);
			double d = distance(c.meanpos,a);
			if(d < ndist){
				ndist = d;
				ret = c;
			}
		}
		return ret;
	}
	public double distance(double[] a,double[] b){
		double sum = 0;
		for(int ii = 0;ii < a.length;ii++){
			sum += vlist.get(ii).weight*(a[ii]-b[ii])*(a[ii]-b[ii]);
		}
		if(sum > 0){
			return Math.sqrt(sum);
		}
		return 0;
	}
	public void prepare(){
		for(WKCluster c:clusters){
			c.meanpos = new double[samples.get(0).pos.length];
			c.updateMean();
		}
	}
	public boolean update(){
		HashSet<WKCluster> updatedc = new HashSet<>();
		for(WKSampleData ss:samples){
			WKCluster c = calcNearestCluster(ss.pos);
			if(c != ss.parent || ss.parent == null){
				if(ss.parent != null){
					updatedc.add(ss.parent);
				}
				updatedc.add(c);
				ss.setParent(c);
			}
			
				if(c == null){
				System.out.println("l;l");	
				}
		}
		for(WKCluster c:updatedc){
			if(c != null){
				c.clearSample(false);
			}
		}
		for(WKSampleData ss:samples){
			if(updatedc.contains(ss.parent)){
				ss.parent.addSample(ss);
			}
		}
		for(WKCluster c:updatedc){
			c.updateMean();
		}
		
		return updatedc.size() > 0;
	}
	public static void main(String[] args){
		WeightedKMeansClustering wc = WeightedKMeansClustering.loadTable("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\kmeanscluster\\testinput.txt", true, true,true);
		wc.setClusterNum(2);
		wc.clusters.get(0).addSample(wc.samples.get(0));
		wc.clusters.get(1).addSample(wc.samples.get(wc.samples.size()-1));
		wc.prepare();
		while(wc.update()){
			System.out.println("updated");
		}
		for(WKCluster cc:wc.clusters){
			System.out.println(cc.meanpos[0]);
			for(WKSampleData s:cc.samples){
				System.out.print(s.name);
			}
			System.out.println("");
		}
	}
	
}



class WKVariable{
	String columnName ="V0";
	int index = 0;
	double weight = 1.0;
}
class WKCluster{
	double meanpos[];
	ArrayList<WKSampleData> samples = new ArrayList<>();
	
	
	public void clearSample(boolean changechildstate){
		if(changechildstate){
			for(WKSampleData ws:samples){
				ws.parent = null;
			}
		}
		samples.clear();
	}
	public void addSample(WKSampleData s){
		s.setParent(this);
		samples.add(s);
	}
	/*
	public void updateSample(Collection<WKSampleData> c){
		samples.clear();
		samples.addAll(c);
		
	}
	*/
	public void updateMean(){
		if(samples.size() == 0){
			return;
		}
		double[] sum = new double[samples.get(0).pos.length];
		for(int ii = 0;ii < sum.length;ii++){
			sum[ii] = 0;
		}
		for(WKSampleData d:samples){
			for(int ii = 0;ii < sum.length;ii++){
				sum[ii] += d.pos[ii];
			}	
		}
		for(int ii = 0;ii < sum.length;ii++){
			meanpos[ii] = sum[ii]/samples.size();
		}
	}
}

class WKSampleData{
	double pos[];
	String classLabel = "";
	String name = "";
	WKCluster parent = null;
	int classGroupIndex = -1;
	WKSampleData(String n,String c,double[] d){
		name = n;
		classLabel = c;
		pos = Arrays.copyOf(d,d.length);
	}
	public void setParent(WKCluster c){
		parent = c;
	}
}



/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class BackBoneSample implements Sampled{
	public static final int N = 0;
	public static final int CA = 1;
	public static final int C = 2;
	public static final int O = 3;
	public static final int PREVC = 4;
	public static final int NEXTN = 5;
	
	String residueName = "UNK";
	double prob = 0.0;
	int count = 0;
	double phi = 0.0;
	double psi = 0.0;
	PDBAtom[] atoms = new PDBAtom[6];
	public static Pattern labelPat = Pattern.compile("###[\\s]*([^\\s]+)[\\s]+([^\\s]+)");
	public static HashMap<String,Integer> aindex = new HashMap<>();
	OmegaSet prevOmega = null;
	OmegaSet nextOmega = null;
	static{
		aindex.put("N", N);
		aindex.put("C", C);
		aindex.put("CA", CA);
		aindex.put("O", O);
		aindex.put("PREVC", PREVC);
		aindex.put("NEXTN", NEXTN);
	}
	
	public static int getIndex(String s){
		String ss = s.toUpperCase();
		if(!aindex.containsKey(ss)){
			return -1;
		}
		return aindex.get(ss);
	}
	
	
	public static BackBoneSample parseBlock(ArrayList<String> block){
		BackBoneSample ret = new BackBoneSample();
		for(String ss:block){
			Matcher matt = labelPat.matcher(ss);
			if(matt.find()){
				String g = matt.group(1);
				if(g.equalsIgnoreCase("residue")){
					ret.residueName = matt.group(2);
				}else if(g.equalsIgnoreCase("code")){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					ret.phi = Double.parseDouble(mapp.get("phi"));
					ret.psi = Double.parseDouble(mapp.get("psi"));
					
				}else if(g.equalsIgnoreCase("count")){
					ret.count = Integer.parseInt(matt.group(2));
				}else if(g.equalsIgnoreCase("prevc")){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					int index = getIndex("PREVC");
					
					ret.atoms[index] = new PDBAtom();
					ret.atoms[index].loc.set(
							Double.parseDouble(mapp.get("x")),
							Double.parseDouble(mapp.get("y")),
							Double.parseDouble(mapp.get("z"))

					);
					ret.atoms[index].pdb_atom_code = "C";
					ret.atoms[index].atom_code = "C";
				}else if(g.equalsIgnoreCase("nextn")){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					int index = getIndex("NEXTN");
					
					ret.atoms[index] = new PDBAtom();
					ret.atoms[index].loc.set(
							Double.parseDouble(mapp.get("x")),
							Double.parseDouble(mapp.get("y")),
							Double.parseDouble(mapp.get("z"))

					);
					ret.atoms[index].pdb_atom_code = "N";
					ret.atoms[index].atom_code = "N";
				}else if(g.indexOf("atom") == 0){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					int index = getIndex(mapp.get("name"));
					if(index == -1){
						System.err.println(mapp.get("name")+" was not considered as backbone atom.");
					}else{
						ret.atoms[index] = new PDBAtom();
						ret.atoms[index].loc.set(
								Double.parseDouble(mapp.get("x")),
								Double.parseDouble(mapp.get("y")),
								Double.parseDouble(mapp.get("z"))
						
						);
						ret.atoms[index].pdb_atom_code = mapp.get("name");
						ret.atoms[index].atom_code = mapp.get("name").substring(0,1);
					}
				}
			}
		}
		return ret;
	}
	
	public double getProb(){
		return prob;
	}
	
	public static ArrayList<OmegaSet> parseOmegaBlock(ArrayList<String> block){
		ArrayList<OmegaSet> ret = new ArrayList<>();
		ArrayList<Double> val = new ArrayList<>();
		ArrayList<Integer> count = new ArrayList<>();
		ArrayList<Double> val_prev = new ArrayList<>();
		ArrayList<Integer> count_prev = new ArrayList<>();
		for(String ss:block){
			Matcher matt = labelPat.matcher(ss);
			if(matt.find()){
				String g = matt.group(1);
				if(g.equalsIgnoreCase("prevomega")){
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					val_prev.add(Double.parseDouble(mapp.get("prevomega"))/180.0*Math.PI);
					count_prev.add(Integer.parseInt(mapp.get("count")));
				}else{
					HashMap<String,String> mapp = MyWorld.lineToHash(ss);
					val.add(Double.parseDouble(mapp.get("omega"))/180.0*Math.PI);
					count.add(Integer.parseInt(mapp.get("count")));
				}
			}
		}
		if(val_prev.size() > 0){
			OmegaSet os = new OmegaSet(val_prev,count_prev,true);
			ret.add(os);
		}
		
		if(val.size() > 0){
			OmegaSet os = new OmegaSet(val,count,false);
			ret.add(os);
		}
		return ret;
		
	}
	
	
	public static ArrayList<BackBoneSample> load(InputStream iss){
		ArrayList<BackBoneSample> ret = new ArrayList<>();
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(iss));
			String ln = null;
			ArrayList<String> buff = new ArrayList<>();
			OmegaSet prevomega = null;
			OmegaSet nextomega = null;
			while((ln = br.readLine()) != null){
				if(ln.indexOf("//") == 0){
					if(buff.get(0).indexOf("###prevomega") == 0 || buff.get(0).indexOf("###nextomega") == 0){
						ArrayList<OmegaSet> al = parseOmegaBlock(buff);
						for(OmegaSet os :al ){
							if(os.prev){
								prevomega = os;
							}else{
								nextomega = os;
							}
						}
					}else{
						ret.add(parseBlock(buff));
					}
					buff.clear();
				}else{
					buff.add(ln);
				}
			}
			br.close();
			int count_all = 0;
			for(BackBoneSample bs:ret){
				count_all += bs.count;
				bs.prevOmega = prevomega;
				bs.nextOmega = nextomega;
			}
			if(count_all > 0){
				for(BackBoneSample bs:ret){
					bs.prob = bs.count/(double)count_all;
				}
			}
			
			Collections.sort(ret,new SampleComparator());
			Collections.reverse(ret);
			prevomega.sort();
			nextomega.sort();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	public static void main(String[] args){
		InputStream iss = BackBoneSample.class.getResourceAsStream("resources/sampledresidues/ALA.backbones.dat");
		if(iss == null){
			System.out.println("//");
		}
		load(iss);
	}
}

class OmegaSample implements Sampled{
	double value = 0;
	int count = 0;
	double prob = 0;
	
	public double getProb(){
		return prob;
	}
}

class OmegaSet{
	boolean prev = false;
	ArrayList<OmegaSample> omegas = new ArrayList<>();
	
	OmegaSet(ArrayList<Double> v,ArrayList<Integer> c,boolean isprev){
		prev = isprev;
		if(v.size() != c.size()){
			throw new RuntimeException();
		}
		for(int ii = 0;ii < v.size();ii++){
			OmegaSample os = new OmegaSample();
			os.value = v.get(ii);
			os.count = c.get(ii);
			omegas.add(os);
		}
		calcProb();
		sort();
	}
	public void sort(){
		Collections.sort(omegas,new SampleComparator());
		Collections.reverse(omegas);
	}
	
	
	
	/**
	 * 数が少なくエラーであると思われるものを削除する
	 */
	public void filt(int threshold){
		Iterator<OmegaSample> ite = omegas.iterator();
		while(ite.hasNext()){
			OmegaSample os = ite.next();
			if(os.count < threshold){
				ite.remove();
			}
		}
		calcProb();
	}
	public void calcProb(){
		int sum = 0;
		for(OmegaSample os:omegas){
			sum += os.count;
			
		}
		if(sum == 0){
			for(OmegaSample os:omegas){
				os.prob = 1.0/omegas.size();
			}
		}else{
			for(OmegaSample os:omegas){
				os.prob = os.count/(double)sum;
			}
		}
	}
}

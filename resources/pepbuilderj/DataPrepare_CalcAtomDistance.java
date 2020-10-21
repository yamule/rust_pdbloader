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
import java.util.regex.Pattern;
import static pepbuilderj.DataPrepare_MakeResidueSample.average;
import static pepbuilderj.DataPrepare_MakeResidueSample.outdirname;
import static pepbuilderj.DataPrepare_MakeResidueSample.sampledirname;

/**
 *
 * @author kimidori
 */
public class DataPrepare_CalcAtomDistance {
	public static String outfilename = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\atom_mindist.dat";


	public static void calcAtomDistance(){
		File dir = new File(sampledirname);
		File outdir = new File(outdirname);
		if(!outdir.exists()){
			outdir.mkdir();
		}
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 


			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			//最短距離トップ 5% を除いた最初の値
			HashMap<String,ArrayList<Double>> mindist= new HashMap<>();
			HashMap<String,ArrayList<Double>> mindist_backbone = new HashMap<>();
			ArrayList<Double> dist_between_prevc_nextn = new ArrayList<>();
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					ArrayList<PDBResidue> rr = new ArrayList<>(PepProcess.makeFilteredAA(c.residues,true));
					boolean[] cbreak = PepProcess.checkChainBreak(rr);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					
					for(int ii = 0;ii < rr.size();ii++){
						allatoms.addAll(rr.get(ii).atoms);
					}
					for(int ii = 0;ii < rr.size();ii++){
						HashSet<PDBAtom> ignore = new HashSet<>();
						if(ii > 0 && ii < rr.size()-1){
							if(!cbreak[ii-1] && !cbreak[ii]){
								PDBAtom pc = rr.get(ii-1).getC();
								PDBAtom nn = rr.get(ii+1).getN();
								if(pc != null && nn != null){
									dist_between_prevc_nextn.add(pc.distance(nn));
								}
								
							}
						}
						
						
						if(ii > 0){
							if(!cbreak[ii-1]){
								ignore.add(rr.get(ii-1).getC());
								ignore.add(rr.get(ii-1).getCA());
								ignore.add(rr.get(ii-1).getN());
								ignore.add(rr.get(ii-1).getO());
								ArrayList<PDBAtom> backbone_prev = new ArrayList<>();
								ArrayList<PDBAtom> backbone_current = new ArrayList<>();
								backbone_prev.add(rr.get(ii-1).getC());
								backbone_prev.add(rr.get(ii-1).getCA());
								backbone_prev.add(rr.get(ii-1).getN());
								backbone_prev.add(rr.get(ii-1).getO());
								
								backbone_current.add(rr.get(ii).getC());
								backbone_current.add(rr.get(ii).getCA());
								backbone_current.add(rr.get(ii).getN());
								backbone_current.add(rr.get(ii).getO());
								for(PDBAtom p2:backbone_prev){
									for(PDBAtom p1:backbone_current){
										String code = "!PREV:"+p2.pdb_atom_code+"\t!CURRENT:"+p1.pdb_atom_code;
										if(!mindist_backbone.containsKey(code)){
											mindist_backbone.put(code, new ArrayList<Double>());
										}
										mindist_backbone.get(code).add(p1.distance(p2));
									}
								}
								for(PDBAtom p2:backbone_current){
									for(PDBAtom p1:backbone_current){
										String code = "!CURRENT:"+p2.pdb_atom_code+"\t!CURRENT:"+p1.pdb_atom_code;
										if(!mindist_backbone.containsKey(code)){
											mindist_backbone.put(code, new ArrayList<Double>());
										}
										mindist_backbone.get(code).add(p1.distance(p2));
									}
								}
							}
						}
						if(ii < rr.size()-1){
							if(!cbreak[ii]){
								ignore.add(rr.get(ii+1).getC());
								ignore.add(rr.get(ii+1).getCA());
								ignore.add(rr.get(ii+1).getN());
								ignore.add(rr.get(ii+1).getO());
							}
						}
						ignore.addAll(rr.get(ii).atoms);
						for(PDBAtom a1:rr.get(ii).atoms){
							String astr = rr.get(ii).getName()+":"+a1.pdb_atom_code;
							for(PDBAtom a2:allatoms){
								if(ignore.contains(a2)){
									continue;
								}
								String bstr = a2.parent.getName()+":"+a2.pdb_atom_code;
								String code = astr+"\t"+bstr;
								double ddist = a1.distance(a2);
								//5.0 より遠いものは interaction しないとして無視する。
								if(ddist > 5.0){
									continue;
								}
								if(astr.compareTo(bstr) > 0){
									code = bstr+"\t"+astr;
								}
								if(!mindist.containsKey(code)){
									mindist.put(code,new ArrayList<Double>());
								}
								mindist.get(code).add(ddist);
							}
						}
					}
				}
			}
			double pcnn[] = get95(dist_between_prevc_nextn);
			pw.write("!DIST_PREVC_NEXTN\t"+"min95:"+pcnn[0]+"\tmax95:"+pcnn[1]+"\tmean5_95:"+average5_95(dist_between_prevc_nextn)+"\n");
			
			for(String ss:mindist_backbone.keySet()){
				double res[] = get95(mindist_backbone.get(ss));
				pw.write(ss+"\t"+"min95:"+res[0]+"\tmax95:"+res[1]+"\tmean5_95:"+average5_95(mindist_backbone.get(ss))+"\n");
			}
			for(String ss:mindist.keySet()){
				double res[] = get95(mindist.get(ss));
				pw.write(ss+"\tmin95:"+res[0]+"\n");
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static double[] get95(ArrayList<Double> al){
		Collections.sort(al);
		double[] ret = new double[2];
		ret[0] = al.get((int)(al.size()*0.05));
		ret[1] = al.get((int)(al.size()*0.95));
		return ret;
	}
	
	/**
	 * 前 5％ 後ろ 5% を除いた平均
	 * @param al
	 * @return 
	 */
	public static double average5_95(ArrayList<Double> al){
		Collections.sort(al);
		int cou = 0;
		double sum = 0;
		int st = (int)(al.size()*0.05);
		int en = (int)(al.size()*0.95);
		for(int ii = st;ii <= en;ii++){
			sum += al.get(ii);
			cou++;
		}
		if(cou == 0){
			return 0;
		}
		return sum/cou;
	}
	
	public static void main(String[] args){
		calcAtomDistance();
	}
}

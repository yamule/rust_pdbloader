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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Pattern;
import static pepbuilderj.PepProcess.calcCrossCenterRatio;
import static pepbuilderj.TemplateBaseModeller.makeAtomLines;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class DataPrepare_MakeResidueSample {
	public static String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\";
	public static String outdirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\rotamer\\";
		
	
	/**
	 * target を、ref の N-CA-C の三角形に合わせる。
	 * 角度を合わせた後は CA の位置を一致させる。
	 * @param target
	 * @param ref 
	 */
	public static void snapTo(PDBResidue target,PDBResidue ref){
		ArrayList<Point3D> al = new ArrayList<>();
		for(PDBAtom a:target.atoms){
			al.add(a.loc);
		}
		PepProcess.adjustVector3D(target.getCA().loc,
		target.getN().loc,
		target.getC().loc,
		ref.getCA().loc,
		ref.getN().loc,
		ref.getC().loc,al);
	}
	
	
	public static double average(ArrayList<Double> al){
		double sum = 0;
		for(Double d:al){
			sum+=d;
		}
		if(al.size() == 0){
			return 0;
		}
		return sum/al.size();
	}
	
	public static Point3D averagePoint(ArrayList<Location3D> al){
		Point3D ret = new Point3D(0,0,0);
		for(Location3D l:al){
			ret.x += l.getX();
			ret.y += l.getY();
			ret.z += l.getZ();
		}
		if(al.size() == 0){
			return ret;
		}
		ret.x /= al.size();
		ret.y /= al.size();
		ret.z /= al.size();
		
		return ret;
	}
	
	
	public static void makeBackboneFragments(){
		File dir = new File(sampledirname);
		File outdir = new File(outdirname);
		if(!outdir.exists()){
			outdir.mkdir();
		}
		File[] lis = dir.listFiles();
		try{

			//td <- read.table("distout.dat",header=F);
			//pca = prcomp(td, scale=T)
			//plot(pca$x[,"PC1"],pca$x[,"PC2"])
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			
			HashMap<String,HashMap<String,ArrayList<PDBResidue>>> rescode_code_list = new HashMap<>();
			HashMap<String,HashMap<String,ArrayList<PDBResidue>>> rescode_pomega_list = new HashMap<>();
			HashMap<String,HashMap<String,ArrayList<PDBResidue>>> rescode_omega_list = new HashMap<>();
			
			HashMap<PDBResidue,PDBResidue> prevr = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> nextr = new HashMap<>();
			
			HashMap<String,ArrayList<Double>> distlist = new HashMap<>();
			distlist.put("NCA",new ArrayList<Double>());
			distlist.put("CAC",new ArrayList<Double>());
			distlist.put("CnexN",new ArrayList<Double>());
			//distlist.put("C_N",new ArrayList<Double>());
			distlist.put("CA_nexN",new ArrayList<Double>());//omega 180 の位置を計算するために保持しておく距離
			distlist.put("C_nexCA",new ArrayList<Double>());//omega 180 の位置を計算するために保持しておく距離
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
					for(int ii = 1;ii < rr.size()-1;ii++){
						if(!cbreak[ii] && !cbreak[ii-1]){
							PDBResidue r = rr.get(ii);
							
							if(!rescode_code_list.containsKey(r.getName())){
								rescode_code_list.put(r.getName(),new HashMap<String,ArrayList<PDBResidue>>());
								rescode_pomega_list.put(r.getName(),new HashMap<String,ArrayList<PDBResidue>>());
								rescode_omega_list.put(r.getName(),new HashMap<String,ArrayList<PDBResidue>>());
							}
							
							PDBResidue pre = rr.get(ii-1);
							PDBResidue nex = rr.get(ii+1);
							prevr.put(r,pre);
							nextr.put(r,nex);

							double phi = PepProcess.phi(r, pre);
							double psi = PepProcess.psi(r, nex);
							double omega = PepProcess.omega(r, nex);
							double prevomega = PepProcess.omega(pre,r);
							
							StringBuffer rcode = new StringBuffer();
							rcode.append("phi:"+((int)((phi+180)/10.0+0.5)*10-180)
									+"\tpsi:"+((int)((psi+180)/10.0+0.5)*10-180));
							
							String rc = rcode.toString();
							String pcode = "prevomega:"+((int)((prevomega+180)/5.0+0.5)*5.0-180);
							String ocode = "omega:"+((int)((omega+180)/5.0+0.5)*5.0-180);
							String resname = r.getName();
							
							distlist.get("NCA").add(r.getCA().distance(r.getN()));
							distlist.get("CAC").add(r.getCA().distance(r.getC()));
							//distlist.get("C_N").add(r.getC().distance(r.getN()));
							distlist.get("CnexN").add(r.getC().distance(nex.getN()));
							distlist.get("CA_nexN").add(r.getCA().distance(nex.getN()));
							distlist.get("C_nexCA").add(r.getC().distance(nex.getCA()));
							
							/*debug 用
							omega が 0 周辺のやつは次が PRO なのが多かった
							prevomega が 0 周辺の奴は PRO でなければ
							自分かその隣にまあなんか強い結合しているようなのが多かった
							if(!resname.equals("PRO") && Math.abs(180-Math.abs(omega)) > 100) {
								System.out.println(f.getName());
								System.out.println(r.getRepresentativeCode());
								System.out.println(pre.getName()+"\t"+nex.getName());
							}
							*/
							
							
							if(!rescode_code_list.get(resname).containsKey(rc)){
								rescode_code_list.get(resname).put(rc,new ArrayList<PDBResidue>());
							}
							if(!rescode_pomega_list.get(resname).containsKey(pcode)){
								rescode_pomega_list.get(resname).put(pcode,new ArrayList<PDBResidue>());
							}
							if(!rescode_omega_list.get(resname).containsKey(ocode)){
								rescode_omega_list.get(resname).put(ocode,new ArrayList<PDBResidue>());
							}
							
							rescode_code_list.get(resname).get(rc).add(r);
							rescode_pomega_list.get(resname).get(pcode).add(r);
							
							rescode_omega_list.get(resname).get(ocode).add(r);
						}
					}
				}
			}
			
			
			for(String s:rescode_code_list.keySet()){
				
				//File xdir = new File(outdir.getPath()+"/"+s);
				//if(!xdir.exists()){
				//	xdir.mkdir();
				//}
				
				File xdir = outdir;//Resource に入れるとサブフォルダにならないため面倒くさい
				String resfilename = xdir.getPath()+"/"+s+".backbones.dat";
				PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(resfilename,false),"UTF-8"))); 
				
				ArrayList<String> pks = new ArrayList<>(rescode_pomega_list.get(s).keySet());
				Collections.sort(pks);
				System.out.println(s);
				for(String pp:pks){
					String lln = "###prevomega\t"+pp+"\tcount:"+rescode_pomega_list.get(s).get(pp).size()+"\n";
					//System.out.print(s+"\t"+lln);
					pw.write(lln);
				}
				pw.write("//\n");
				ArrayList<String> pks2 = new ArrayList<>(rescode_omega_list.get(s).keySet());
				Collections.sort(pks2);
				for(String pp:pks2){
					String lln = "###nextomega\t"+pp+"\tcount:"+rescode_omega_list.get(s).get(pp).size()+"\n";
					//System.out.print(s+"\t"+lln);
					pw.write(lln);
				}
				pw.write("//\n");
				
				for(String code:rescode_code_list.get(s).keySet()){
					ArrayList<PDBResidue> al = rescode_code_list.get(s).get(code);
					PDBResidue ref;//平均的な構成の持った Residue
					
					if(al.size() == 1){
						ref = al.get(0);
					}else{
						
						ArrayList<Double> phi = new ArrayList<>();
						ArrayList<Double> psi = new ArrayList<>();
						ArrayList<Double> distnca = new ArrayList<>();
						ArrayList<Double> distcac = new ArrayList<>();
						ArrayList<Double> distnextn = new ArrayList<>();
						ArrayList<Double> distprevc = new ArrayList<>();
						for(PDBResidue r:al){
							PDBResidue pre = prevr.get(r);
							PDBResidue nex = nextr.get(r);
							double dphi = PepProcess.phi(r, pre);
							double dpsi = PepProcess.psi(r, nex);
							
							phi.add(dphi);
							psi.add(dpsi);
							distnca.add(r.getN().distance(r.getCA()));
							distcac.add(r.getCA().distance(r.getC()));
							
							distnextn.add(r.getC().distance(nextr.get(r).getN()));
							distprevc.add(r.getN().distance(prevr.get(r).getC()));
						}
						double avph = average(phi);
						double avps = average(psi);
						double avnca = average(distnca);
						double avcac = average(distcac);
						double avpc = average(distprevc);
						double avnn = average(distnextn);
						ArrayList<VSorter> s1 = new ArrayList<>();
						ArrayList<VSorter> s2 = new ArrayList<>();
						ArrayList<VSorter> s3 = new ArrayList<>();
						ArrayList<VSorter> s4 = new ArrayList<>();
						ArrayList<VSorter> s5 = new ArrayList<>();
						ArrayList<VSorter> s6 = new ArrayList<>();
						for(int ii = 0;ii < phi.size();ii++){
							s1.add(new VSorter(Math.abs(phi.get(ii)-avph),ii));
							s2.add(new VSorter(Math.abs(psi.get(ii)-avps),ii));
							s3.add(new VSorter(Math.abs(distnca.get(ii)-avnca),ii));
							s4.add(new VSorter(Math.abs(distcac.get(ii)-avcac),ii));
							s5.add(new VSorter(Math.abs(distprevc.get(ii)-avpc),ii));
							s6.add(new VSorter(Math.abs(distnextn.get(ii)-avnn),ii));
						}
						Collections.sort(s1,new VComparator());
						Collections.sort(s2,new VComparator());
						Collections.sort(s3,new VComparator());
						Collections.sort(s4,new VComparator());
						Collections.sort(s5,new VComparator());
						Collections.sort(s6,new VComparator());
						int minrank = s1.size()*4;
						PDBResidue minres = null;
						int[] sumrank = new int[s1.size()];
						
						for(int ii = 0;ii < s1.size();ii++){
							sumrank[ii] = 0;
						}
						
						
						
						for(int ii = 0;ii < s1.size();ii++){
							sumrank[s1.get(ii).index] += ii;
							sumrank[s2.get(ii).index] += ii;
							sumrank[s3.get(ii).index] += ii*2;//キョリには大きめのウェイト
							sumrank[s4.get(ii).index] += ii*2;
							sumrank[s5.get(ii).index] += ii*2;
							sumrank[s6.get(ii).index] += ii*2;
						}
						for(int ii = 0;ii < s1.size();ii++){
							if(sumrank[ii] < minrank){
								minrank = sumrank[ii];
								minres = al.get(ii);
							}
						}
						/*
						String[] debug = new String[s1.size()];
						for(int ii = 0;ii < s1.size();ii++){
							debug[ii] = "";
						}
						for(int ii = 0;ii < s1.size();ii++){
							debug[s1.get(ii).index] += ii;
						}
						
						for(int ii = 0;ii < s2.size();ii++){
							debug[s2.get(ii).index] += ii;
						}
						
						for(int ii = 0;ii < s3.size();ii++){
							debug[s3.get(ii).index] += ii;
						}
						for(int ii = 0;ii < s4.size();ii++){
							debug[s4.get(ii).index] += ii;
						}
						for(int ii = 0;ii < s5.size();ii++){
							debug[s5.get(ii).index] += ii;
						}
						for(int ii = 0;ii < s6.size();ii++){
							debug[s6.get(ii).index] += ii;
						}
						
						
						System.out.println("=================");
						System.out.println(avnca);
						System.out.println(avcac);
						System.out.println(avpc);
						System.out.println(avnn);
						System.out.println("+++");
						System.out.println(minres.getN().distance(minres.getCA()));
						System.out.println(minres.getCA().distance(minres.getC()));
						System.out.println(minres.getN().distance(prevr.get(minres).getC()));
						System.out.println(minres.getC().distance(nextr.get(minres).getN()));
						PDBResidue dres = al.get(0);
						System.out.println("+++");
						System.out.println(dres.getN().distance(dres.getCA()));
						System.out.println(dres.getCA().distance(dres.getC()));
						System.out.println(dres.getN().distance(prevr.get(dres).getC()));
						System.out.println(dres.getC().distance(nextr.get(dres).getN()));
						System.out.println("=================");
						*/
						ref = minres;
					}
					pw.write("###===========================\n");
					pw.write("###code\t"+code+"\n");
					pw.write("###count\t"+al.size()+"\n");
					pw.write("###debug\t"+al.size()+"\t"+code.replaceAll("\\:","\t")+"\n");
					pw.write("###residue\t"+s+"\n");
					Location3D ploc = prevr.get(ref).getC().getLocation();
					Location3D nloc = nextr.get(ref).getN().getLocation();
					Location3D plocca = prevr.get(ref).getCA().getLocation();
					Location3D nlocca = nextr.get(ref).getCA().getLocation();
					pw.write("###prevc\tx:"+ploc.getX()+"\ty:"+ploc.getY()+"\tz:"+ploc.getZ()+"\n");
					pw.write("###nextn\tx:"+nloc.getX()+"\ty:"+nloc.getY()+"\tz:"+nloc.getZ()+"\n");
					
					pw.write("###prevca\tx:"+plocca.getX()+"\ty:"+plocca.getY()+"\tz:"+plocca.getZ()+"\n");
					pw.write("###nextca\tx:"+nlocca.getX()+"\ty:"+nlocca.getY()+"\tz:"+nlocca.getZ()+"\n");
					for(PDBAtom aa:ref.atoms){
						String sss = aa.pdb_atom_code;
						if(sss.equals("O") || sss.equals("C") || sss.equals("CA") || sss.equals("N")){
							Location3D aloc = aa.loc;
							pw.write("###atom\tname:"+sss+"\tx:"+aloc.getX()+"\ty:"+aloc.getY()+"\tz:"+aloc.getZ()+"\n");
						}else{
						}
					}
					pw.write("//\n");
				}
				pw.close();
			}
			for(String s:distlist.keySet()){
				System.out.println(s+"\t"+average(distlist.get(s)));
			}
			/*
			2018/02/10
			CnexN	1.3303962296321674
			CAC	1.525136543237837
			CA_nexN	2.4318112083177663
			C_nexCA	2.4342773976940855
			NCA	1.459856507426044
			*/
			
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
	}
	/**
	 * Omega が 180 になるような CA の位置を計算するための
	 * prevCA->x->CA となる n->prevc を結ぶ線分上にある x の、n からキョリを
	 * n->prevc の距離で割った値
	 * x->CA の距離を計算する
	 * @return 
	 */
	public static double[] calcOmega180Ratio(){
		double[] ret = new double[3];
		
			/*
			2018/02/10
			CnexN	1.3303962296321674
			CAC	1.525136543237837
			CA_nexN	2.4318112083177663
			C_nexCA	2.4342773976940855
			NCA	1.459856507426044
			*/
			
	
		double canexn = 2.4318112083177663;
		double nca = 1.459856507426044;
		double cnexca = 2.4342773976940855;
		double cac = 1.525136543237837;
		double cnexn = 1.3303962296321674;
		ret[0] = calcCrossCenterRatio(canexn,nca,cnexca,cac,cnexn);
		double u = (nca*nca+cnexn*cnexn-cnexca*cnexca)/(2*cnexn);
		double r = Math.sqrt(nca*nca-u*u);
		ret[1] = Math.sqrt(r*r+(u-ret[0]*cnexn)*(u-ret[0]*cnexn));
		
		
		double u2 = (canexn*canexn+cnexn*cnexn-cac*cac)/(2*cnexn);
		double r2 = Math.sqrt(canexn*canexn-u2*u2);
		
		ret[2] = Math.sqrt(r2*r2+(u2-ret[0]*cnexn)*(u2-ret[0]*cnexn));
		
		return ret;
	}
	public static double[] calcOmega180Ratio_debug(){
		double[] ret = new double[2];
			
			double canexn = 1;
			double nca = 1;
			double cnexca = 1;
			double cac = 1;
			double cnexn = 1;
		ret[0] = calcCrossCenterRatio(canexn,nca,cnexca,cac,cnexn);
		double u = (nca*nca+cnexn*cnexn-cnexca*cnexca)/(2*cnexn);
		double r = Math.sqrt(nca*nca-u*u);
		ret[1] = Math.sqrt(r*r+(u-ret[0]*cnexn)*(u-ret[0]*cnexn));
		return ret;
	}
	public static void makeRotamer(){
		
		HashMap<String,HashSet<String>> atomsets = MyWorld.sidechain_atoms;
		for(String tss:atomsets.keySet()){
			atomsets.get(tss).addAll(MyWorld.backbone_atoms);
		}
		
		File dir = new File(sampledirname);
		File outdir = new File(outdirname);
		File[] list = dir.listFiles();
		if(!outdir.exists()){
			outdir.mkdir();
		}
		HashMap<String,ArrayList<PDBResidue>> allaa = new HashMap<>();
		Pattern pdbpat = Pattern.compile("\\.(pdb|ent)$");
		for(File ff:list){
			if(!pdbpat.matcher(ff.getName()).find()){
				continue;
			}
			PDBData pdb = PDBData.loadPDBFile(ff.getPath());
			for(String s:pdb.chains.keySet()){
				PDBChain c = pdb.chains.get(s);
				ArrayList<PDBResidue> filtered = PepProcess.makeFilteredAA(c.residues,true);
				for(PDBResidue rr:filtered){
					if(!allaa.containsKey(rr.getName())){
						allaa.put(rr.getName(),new ArrayList<PDBResidue>());
					}
					allaa.get(rr.getName()).add(rr);
				}
			}
		}
		for(String ss:allaa.keySet()){
			ArrayList<PDBResidue> al = allaa.get(ss);
			
			int scou = al.size();
			Iterator<PDBResidue> ite = al.iterator();
			HashSet<String> atomset = atomsets.get(ss);
			
			while(ite.hasNext()){
				PDBResidue r = ite.next();
				HashSet<String> att = new HashSet<>();
				for(PDBAtom a:r.atoms){
					att.add(a.pdb_atom_code);
				}
				if(att.size() != atomset.size()){
					ite.remove();
				}
			}
			if(scou > al.size()*2){
				System.err.println("???? is there some unusual atoms? "+ss);
				for(String s:atomset){
					System.err.print(s+",");
				}
				throw new RuntimeException();
			}else{
				//System.out.print("{\"#"+ss+"\"");
				//for(String s:atomset){
				//	System.out.print(",\""+s+"\"");
				//}
				//System.out.println("}");
			}
			System.out.println(al.size()+"/"+scou+";"+ss);
			
			
			
			//Collections.reverse(al);
			//ちゃんと想定通りに動いているか。ランダムにやってみたが、残った Residue の数は同じだった。
			//中身が同じかは見てない
			//for(int ii = 0;ii < al.size();ii++){
			//	int sst = (int)(Math.random()*al.size());
			//	int sst2 = (int)(Math.random()*al.size());
			//	PDBResidue aaa = al.get(sst);
			//	al.set(sst,al.get(sst2));
			//	al.set(sst2,aaa);
			//}
			
			PDBResidue base = al.remove(0);
			for(int ii = 0;ii < al.size();ii++){
				snapTo(al.get(ii),base);
			}
			al.add(base);
			
			
			HashMap<String,ArrayList<Location3D>> allpos = new HashMap<>();
			for(int ii = 0;ii < al.size();ii++){
				PDBResidue r = al.get(ii);
				for(PDBAtom a:r.atoms){
					if(!allpos.containsKey(a.pdb_atom_code)){
						allpos.put(a.pdb_atom_code,new ArrayList<Location3D>());
					}
					allpos.get(a.pdb_atom_code).add(a.loc);
				}
			}
			
			HashMap<String,Point3D> avepos = new HashMap<>();
			for(String s:allpos.keySet()){
				avepos.put(s,averagePoint(allpos.get(s)));
			}
			ArrayList<VSorter> vss = new ArrayList<>();
			for(int ii = 0;ii < al.size();ii++){
				double maxdist = 0;
				for(PDBAtom a:al.get(ii).atoms){
					maxdist = Math.max(maxdist,a.distance(avepos.get(a.pdb_atom_code)));
				}
				vss.add(new VSorter(maxdist,ii));
			}
			
			Collections.sort(vss,new VComparator());
			Collections.reverse(vss);
			ArrayList<PDBResidue> dal = new ArrayList<>();
			for(int ii = 0;ii < vss.size();ii++){
				dal.add(al.get(vss.get(ii).index));
			}
			al.clear();
			ArrayList<PDBResidue> remained = new ArrayList<>();
			ArrayList<Integer> similarcount = new ArrayList<>();
			while(dal.size() > 0){
				PDBResidue candidate = dal.remove(0);
				remained.add(candidate);
				int similar = 1;
				Iterator<PDBResidue> pr = dal.iterator();
				while(pr.hasNext()){
					PDBResidue c = pr.next();
					double maxdist = 0;
					for(String sss:atomset){
						if(sss.equals("O") || sss.equals("C") || sss.equals("CA") || sss.equals("N")){
							continue;
						}
						PDBAtom a = candidate.getAtomByName(sss);
						PDBAtom b = c.getAtomByName(sss);
						maxdist = Math.max(maxdist,a.distance(b));
					}
					if(maxdist < 0.5){
						pr.remove();
						similar++;
					}
				}
				similarcount.add(similar);
			}
			//File xdir = new File(outdir.getPath()+"/"+ss);
			//if(!xdir.exists()){
			//	xdir.mkdir();
			//}
			File xdir = outdir;
			String fname = xdir.getPath()+"/"+ss+".rotamers.dat";
			//String simname = xdir.getPath()+"/"+ss+".count.dat";
			ArrayList<PDBResidue> res = new ArrayList<>();
			try{
				PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fname,false),"UTF-8"))); 
				
				for(int jj = 0;jj < remained.size();jj++){
					pw.write("###===========================\n");
					pw.write("###index\t"+jj+"\n");
					pw.write("###residue\t"+ss+"\n");
					
					pw.write("###count\t"+similarcount.get(jj)+"\n");
					PDBResidue rr = remained.get(jj);
					for(PDBAtom a:rr.atoms){
						pw.write("###atom\tname:"+a.pdb_atom_code);
						pw.write("\tx:"+a.loc.x);
						pw.write("\ty:"+a.loc.y);
						pw.write("\tz:"+a.loc.z+"\n");
					}
					pw.write("//\n");
				}
				System.out.println("remained "+ss+":"+remained.size());
				pw.close();
			}catch(Exception exx){
				exx.printStackTrace();
			}
			//writeToFile(makeAtomLines(remained,1,1,"A"),fname);
		}
	}
	public static void main(String[] args){
		makeRotamer();
		makeBackboneFragments();
		//double[] ret = calcOmega180Ratio();
		//System.out.println(ret[0]);
		//System.out.println(ret[1]);
		//System.out.println(ret[2]);
		/*
		2018/02/10
			0.4237513761543342
			1.8184933638016645
			1.9902838487766326
		*/
		
		
		
		//System.out.println(calcOmega180Ratio_debug());
		
		/*
		ArrayList<PDBResidue> al = new ArrayList<>();
		for(int ii = 0;ii < 3;ii++){
			PDBResidue r = new PDBResidue();
			r.setName("GLY");
			PDBAtom a = new PDBAtom();
			a.pdb_atom_code = "N";
			PDBAtom b = new PDBAtom();
			b.pdb_atom_code = "CA";
			PDBAtom c = new PDBAtom();
			c.pdb_atom_code = "C";
			PDBAtom d = new PDBAtom();
			d.pdb_atom_code = "O";
			r.addAtom(a);
			r.addAtom(b);
			r.addAtom(c);
			r.addAtom(d);
			al.add(r);
		}
		al.get(0).getC().loc.set(-2,0,0);
		al.get(1).getN().loc.set(-1,-1,0);
		al.get(1).getCA().loc.set(0,0,0);
		al.get(1).getC().loc.set(1,-1,0);
		al.get(2).getN().loc.set(2,0,0);
		
		
		double phi = PepProcess.phi(al.get(1),al.get(0));
		double psi = PepProcess.psi(al.get(1),al.get(2));
		System.out.println(phi+"\t"+psi);
		*/
	}
	
}


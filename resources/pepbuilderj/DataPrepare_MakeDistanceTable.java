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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
/**
 *PSSM ALIGNMENT を作るためのスコアリングするための Decision Tree を作るためのテーブルを作成する
 * @author kimidori
 */
public class DataPrepare_MakeDistanceTable {
	public static String outfilename = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_onlycb_jumpcb45\\table_merged.dat";
	public static String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\";
	
	public static boolean hasBackbone(PDBResidue r){
		return r.getC() != null && r.getCA() != null && r.getN() != null;
	}
	public static double[] calcangles(PDBResidue target,PDBResidue prev,PDBResidue next){
		if(!hasBackbone(target) || !hasBackbone(prev) || !hasBackbone(next)){
			return null;
		}
		double x = PepProcess.phi(target,prev);
		double y = PepProcess.psi(target,next);	
		double o = Math.abs(180-Math.abs(PepProcess.omega(prev,target)));
		double o2 = Math.abs(180-Math.abs(PepProcess.omega(target,next)));
		double[] ret = {x,y,o,o2};
		return ret;
	}
	public static Point3D calcReverseCB(PDBResidue r){
		Point3D ca = r.getCA().loc;
		Point3D oc = r.getC().loc;
		Point3D n = r.getN().loc;
		
		
		if(ca == null || oc == null || n == null){
			return null;
		}
		Point3D vec = new Point3D();
		vec.x = oc.x-n.x;
		vec.y = oc.y-n.y;
		vec.z = oc.z-n.z;
		
		
		Point3D pb = new Point3D();
		pb.x = n.x/2+oc.x/2 -ca.x;
		pb.y = n.y/2+oc.y/2 -ca.y;
		pb.z = n.z/2+oc.z/2 -ca.z;
		
		//Point3D ret = Point3D.rotate(pb,vec,-2.186025);
		Point3D ret = Point3D.rotate(pb,vec,2.186025);
		
		double rlen = Math.sqrt(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);
		if(rlen == 0){
			return null;
			
		}
		rlen = 1.53/rlen;
		ret.x *= rlen;
		ret.y *= rlen;
		ret.z *= rlen;
		
		
		ret.x += ca.x;
		ret.y += ca.y;
		ret.z += ca.z;
		
		
		return ret;
	}
	public static Point3D calcCenterCB(PDBResidue r){
		
		if(r.getCA() == null || r.getC() == null || r.getN() == null){
			return null;
		}
		//if(r.getCB() != null){
		//	return r.getCB().loc;
		//}
		Point3D ca = r.getCA().loc;
		Point3D oc = r.getC().loc;
		Point3D n = r.getN().loc;
		
		
		if(ca == null || oc == null || n == null){
			return null;
		}
		Point3D vec = new Point3D();
		vec.x = oc.x-n.x;
		vec.y = oc.y-n.y;
		vec.z = oc.z-n.z;
		
		
		Point3D ret = new Point3D();
		ret.x = n.x/2+oc.x/2 -ca.x;
		ret.y = n.y/2+oc.y/2 -ca.y;
		ret.z = n.z/2+oc.z/2 -ca.z;
		
		
		double rlen = Math.sqrt(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);
		if(rlen == 0){
			return null;
			
		}
		rlen = 1.53/rlen;
		ret.x *= rlen;
		ret.y *= rlen;
		ret.z *= rlen;
		
		ret.x *= -1;
		ret.y *= -1;
		ret.z *= -1;
		
		
		ret.x += ca.x;
		ret.y += ca.y;
		ret.z += ca.z;
		
		
		return ret;
	}
		
	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_phipsi_removing(){
		//特徴量を減らすのを試みる
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {3,6,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_neighbor = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_next = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_prev = new HashMap<>();
			
			HashMap<PDBResidue,Double> phiangles = new HashMap<>();
			HashMap<PDBResidue,Double> psiangles = new HashMap<>();
			HashMap<PDBResidue,Double> omegaangles = new HashMap<>();
			HashMap<PDBResidue,Double> omega2angles = new HashMap<>(); 
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			for(String a:atomnames){
				if(a.equals("CA")){
					for(double d:distthreshold){
						headertext.add(a+"_"+d);
					}
				}
			}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			//for(double d:distthreshold){
			//	headertext.add(jcbdist+"CB_ALL_"+d);
			//}

			//headertext.add("neighborCA");
			headertext.add("phi");
			headertext.add("psi");
			headertext.add("omega1");
			headertext.add("omega2");
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_neighbor.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					for(int ii = 1;ii < c.residues.size();ii++){
						PDBResidue prev = null;
						PDBResidue curr = c.residues.get(ii);
						if(ii > 0){
							prev = c.residues.get(ii-1);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_a){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}


					for(int ii = 0;ii < c.residues.size();ii++){
						for(int jj = 0;jj < c.residues.size();jj++){
							if(ii == jj){//無くてもいい気がする
								continue;
							}
							PDBResidue curr = c.residues.get(ii);
							PDBResidue prev = c.residues.get(jj);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_b){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}
					
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					//前後残基の CA 距離
					HashMap<PDBResidue,Double> cadists = new HashMap<>();
					
					
					for(PDBResidue rr:c.residues){
						if(res_neighbor.containsKey(rr)){
							ArrayList<PDBResidue> nn = new ArrayList<PDBResidue>(res_neighbor.get(rr));
							if(nn.size() != 2){
								cadists.put(rr,1000.0);
							}else if(nn.get(0).getCA() == null || nn.get(1).getCA() == null){
								cadists.put(rr,1000.0);
							}else{
								cadists.put(rr,nn.get(0).getCA().distance(nn.get(1).getCA()));
								double[] angles = calcangles(rr,res_prev.get(rr),res_next.get(rr));
								if(angles == null){
									phiangles.put(rr,1000.0);
									psiangles.put(rr,1000.0);
									omegaangles.put(rr,1000.0);
									omega2angles.put(rr,1000.0);
								}else{
									phiangles.put(rr,angles[0]);
									psiangles.put(rr,angles[1]);
									omegaangles.put(rr,angles[2]);
									omega2angles.put(rr,angles[3]);
								}
							}
						}else{
							cadists.put(rr,1000.0);
							
						}
					}
					
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}

							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.equals("phi")){
								if(phiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(phiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("psi")){
								if(psiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(psiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega1")){
								if(omegaangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omegaangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega2")){
								if(omega2angles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omega2angles.get(rr)+"\t");
								}
							
							}else if(h.indexOf("neighborCA") > -1){
								if(cadists.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(cadists.get(rr)+"\t");
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}

		
	public static void calcAtomDistance_TopFiveNearestdist_ThreeAtoms_distdiff_phipsi(){
		//0.567 が最高
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			/*
			double tcbdist = 5.0;
			double tcbdist2 = 2.5;
			double tcbdist3 = 2.5;
			*/
			
			double tcbdist = 2.5;
			double tcbdist2 = 3.0;
			double tcbdist3 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_neighbor = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_next = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_prev = new HashMap<>();
			
			HashMap<PDBResidue,Double> phiangles = new HashMap<>();
			HashMap<PDBResidue,Double> psiangles = new HashMap<>();
			HashMap<PDBResidue,Double> omegaangles = new HashMap<>();
			HashMap<PDBResidue,Double> omega2angles = new HashMap<>(); 
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
					headertext.add(a+"_rankc"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankc"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}

			//headertext.add("neighborCA");
			headertext.add("phi");
			headertext.add("psi");
			headertext.add("omega1");
			headertext.add("omega2");
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_neighbor.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					for(int ii = 1;ii < c.residues.size();ii++){
						PDBResidue prev = null;
						PDBResidue curr = c.residues.get(ii);
						if(ii > 0){
							prev = c.residues.get(ii-1);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_a){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}


					for(int ii = 0;ii < c.residues.size();ii++){
						for(int jj = 0;jj < c.residues.size();jj++){
							if(ii == jj){//無くてもいい気がする
								continue;
							}
							PDBResidue curr = c.residues.get(ii);
							PDBResidue prev = c.residues.get(jj);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_b){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}
					
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs3 = new HashMap<>();
					//前後残基の CA 距離
					HashMap<PDBResidue,Double> cadists = new HashMap<>();
					
					
					for(PDBResidue rr:c.residues){
						if(res_neighbor.containsKey(rr)){
							ArrayList<PDBResidue> nn = new ArrayList<PDBResidue>(res_neighbor.get(rr));
							if(nn.size() != 2){
								cadists.put(rr,1000.0);
							}else if(nn.get(0).getCA() == null || nn.get(1).getCA() == null){
								cadists.put(rr,1000.0);
							}else{
								cadists.put(rr,nn.get(0).getCA().distance(nn.get(1).getCA()));
								double[] angles = calcangles(rr,res_prev.get(rr),res_next.get(rr));
								if(angles == null){
									phiangles.put(rr,1000.0);
									psiangles.put(rr,1000.0);
									omegaangles.put(rr,1000.0);
									omega2angles.put(rr,1000.0);
								}else{
									phiangles.put(rr,angles[0]);
									psiangles.put(rr,angles[1]);
									omegaangles.put(rr,angles[2]);
									omega2angles.put(rr,angles[3]);
								}
							}
						}else{
							cadists.put(rr,1000.0);
							
						}
					}
					
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
							
							if(tcbdist3 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist3;
								yy /= dist/tcbdist3;
								zz /= dist/tcbdist3;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs3.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}

							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								AtomDistance a3  = new AtomDistance(aa,badist,badist2-badist);
								a3.distance3 = tcbs3.get(rr).distance(aa);
								distrank.get(code).add(a3);
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rank[a-z]?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.equals("phi")){
								if(phiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(phiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("psi")){
								if(psiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(psiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega1")){
								if(omegaangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omegaangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega2")){
								if(omega2angles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omega2angles.get(rr)+"\t");
								}
							
							}else if(h.indexOf("neighborCA") > -1){
								if(cadists.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(cadists.get(rr)+"\t");
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else if(h.indexOf("_rankc") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance3))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	

	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_phipsi_centercb(){
		//0.567 が最高
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_neighbor = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_next = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_prev = new HashMap<>();
			
			HashMap<PDBResidue,Double> phiangles = new HashMap<>();
			HashMap<PDBResidue,Double> psiangles = new HashMap<>();
			HashMap<PDBResidue,Double> omegaangles = new HashMap<>();
			HashMap<PDBResidue,Double> omega2angles = new HashMap<>(); 
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}

			//headertext.add("neighborCA");
			headertext.add("phi");
			headertext.add("psi");
			headertext.add("omega1");
			headertext.add("omega2");
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_neighbor.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					for(int ii = 1;ii < c.residues.size();ii++){
						PDBResidue prev = null;
						PDBResidue curr = c.residues.get(ii);
						if(ii > 0){
							prev = c.residues.get(ii-1);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_a){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}


					for(int ii = 0;ii < c.residues.size();ii++){
						for(int jj = 0;jj < c.residues.size();jj++){
							if(ii == jj){//無くてもいい気がする
								continue;
							}
							PDBResidue curr = c.residues.get(ii);
							PDBResidue prev = c.residues.get(jj);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_b){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}
					
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					//前後残基の CA 距離
					HashMap<PDBResidue,Double> cadists = new HashMap<>();
					
					
					for(PDBResidue rr:c.residues){
						if(res_neighbor.containsKey(rr)){
							ArrayList<PDBResidue> nn = new ArrayList<PDBResidue>(res_neighbor.get(rr));
							if(nn.size() != 2){
								cadists.put(rr,1000.0);
							}else if(nn.get(0).getCA() == null || nn.get(1).getCA() == null){
								cadists.put(rr,1000.0);
							}else{
								cadists.put(rr,nn.get(0).getCA().distance(nn.get(1).getCA()));
								double[] angles = calcangles(rr,res_prev.get(rr),res_next.get(rr));
								if(angles == null){
									phiangles.put(rr,1000.0);
									psiangles.put(rr,1000.0);
									omegaangles.put(rr,1000.0);
									omega2angles.put(rr,1000.0);
								}else{
									phiangles.put(rr,angles[0]);
									psiangles.put(rr,angles[1]);
									omegaangles.put(rr,angles[2]);
									omega2angles.put(rr,angles[3]);
								}
							}
						}else{
							cadists.put(rr,1000.0);
							
						}
					}
					
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								
								Point3D dloc = calcCenterCB(rr);
								if(dloc == null){
									dloc = cb.loc;
								}
								
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								
								double xx = dloc.x-ca.loc.x;
								double yy = dloc.y-ca.loc.y;
								double zz = dloc.z-ca.loc.z;
								double dist = dloc.distance(ca.loc);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								//Point3D dloc = calcCenterCB(rr);
								Point3D dloc = cb.loc;
								if(dloc == null){
									dloc = cb.loc;
								}
								double xx = dloc.x-ca.loc.x;
								double yy = dloc.y-ca.loc.y;
								double zz = dloc.z-ca.loc.z;
								double dist = dloc.distance(ca.loc);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								Point3D dloc = calcCenterCB(rr);
								if(dloc == null){
									dloc = cb.loc;
								}
								double xx = dloc.x-ca.loc.x;
								double yy = dloc.y-ca.loc.y;
								double zz = dloc.z-ca.loc.z;
								double dist = dloc.distance(ca.loc);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}

							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.equals("phi")){
								if(phiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(phiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("psi")){
								if(psiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(psiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega1")){
								if(omegaangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omegaangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega2")){
								if(omega2angles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omega2angles.get(rr)+"\t");
								}
							
							}else if(h.indexOf("neighborCA") > -1){
								if(cadists.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(cadists.get(rr)+"\t");
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_phipsi(){
		//0.567 が最高
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_neighbor = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_next = new HashMap<>();
			HashMap<PDBResidue,PDBResidue> res_prev = new HashMap<>();
			
			HashMap<PDBResidue,Double> phiangles = new HashMap<>();
			HashMap<PDBResidue,Double> psiangles = new HashMap<>();
			HashMap<PDBResidue,Double> omegaangles = new HashMap<>();
			HashMap<PDBResidue,Double> omega2angles = new HashMap<>(); 
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}

			//headertext.add("neighborCA");
			headertext.add("phi");
			headertext.add("psi");
			headertext.add("omega1");
			headertext.add("omega2");
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_neighbor.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					for(int ii = 1;ii < c.residues.size();ii++){
						PDBResidue prev = null;
						PDBResidue curr = c.residues.get(ii);
						if(ii > 0){
							prev = c.residues.get(ii-1);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_a){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}


					for(int ii = 0;ii < c.residues.size();ii++){
						for(int jj = 0;jj < c.residues.size();jj++){
							if(ii == jj){//無くてもいい気がする
								continue;
							}
							PDBResidue curr = c.residues.get(ii);
							PDBResidue prev = c.residues.get(jj);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_b){
										res_neighbor.get(prev).add(curr);
										res_neighbor.get(curr).add(prev);
										res_prev.put(curr,prev);
										res_next.put(prev,curr);
									}
								}
							}
						}
					}
					
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					//前後残基の CA 距離
					HashMap<PDBResidue,Double> cadists = new HashMap<>();
					
					
					for(PDBResidue rr:c.residues){
						if(res_neighbor.containsKey(rr)){
							ArrayList<PDBResidue> nn = new ArrayList<PDBResidue>(res_neighbor.get(rr));
							if(nn.size() != 2){
								cadists.put(rr,1000.0);
							}else if(nn.get(0).getCA() == null || nn.get(1).getCA() == null){
								cadists.put(rr,1000.0);
							}else{
								cadists.put(rr,nn.get(0).getCA().distance(nn.get(1).getCA()));
								double[] angles = calcangles(rr,res_prev.get(rr),res_next.get(rr));
								if(angles == null){
									phiangles.put(rr,1000.0);
									psiangles.put(rr,1000.0);
									omegaangles.put(rr,1000.0);
									omega2angles.put(rr,1000.0);
								}else{
									phiangles.put(rr,angles[0]);
									psiangles.put(rr,angles[1]);
									omegaangles.put(rr,angles[2]);
									omega2angles.put(rr,angles[3]);
								}
							}
						}else{
							cadists.put(rr,1000.0);
							
						}
					}
					
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}

							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.equals("phi")){
								if(phiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(phiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("psi")){
								if(psiangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(psiangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega1")){
								if(omegaangles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omegaangles.get(rr)+"\t");
								}
							
							}else if(h.equals("omega2")){
								if(omega2angles.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(omega2angles.get(rr)+"\t");
								}
							
							}else if(h.indexOf("neighborCA") > -1){
								if(cadists.get(rr) == null){
									pw.write("1000\t");
								}else{
									pw.write(cadists.get(rr)+"\t");
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}

	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_centercb(){
		//あまり良くなかった
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.0;
			double tcbdist = 2.0;
			double tcbdist2 = 4.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}
					if(ignoreneighbor){
					for(int ii = 1;ii < c.residues.size();ii++){
						PDBResidue prev = null;
						PDBResidue curr = c.residues.get(ii);
						if(ii > 0){
							prev = c.residues.get(ii-1);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_a){
										res_ignore.get(prev).add(curr);
										res_ignore.get(curr).add(prev);
									}
								}
							}
						}
					}
					
					
					for(int ii = 0;ii < c.residues.size();ii++){
						for(int jj = 0;jj < c.residues.size();jj++){
							if(ii == jj){//無くてもいい気がする
								continue;
							}
							PDBResidue curr = c.residues.get(ii);
							PDBResidue prev = c.residues.get(jj);
							if(prev.getC() != null){
								if(curr.getN() != null){
									if(prev.getC().distance(curr.getN()) < covthreshold_b){
										res_ignore.get(prev).add(curr);
										res_ignore.get(curr).add(prev);
									}
								}
							}
						}
					}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							Point3D ccb = calcCenterCB(rr);
							if(ccb == null){
								ccb = cb.loc;
							}
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								
								double xx = ccb.x-ca.loc.x;
								double yy = ccb.y-ca.loc.y;
								double zz = ccb.z-ca.loc.z;
								double dist = ccb.distance(ca.loc);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = ccb.x-ca.loc.x;
								double yy = ccb.y-ca.loc.y;
								double zz = ccb.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = ccb.x-ca.loc.x;
								double yy = ccb.y-ca.loc.y;
								double zz = ccb.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}else if(cb != null){
							PDBAtom jcb = new PDBAtom();
							jcb.atom_code = "C";
							jcb.pdb_atom_code = jcbdist+"CB";

							jcb.loc.set(cb.loc);
							jcb.parent = rr;
							envatoms.add(jcb);
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	

	public static void calcAtomDistance_TopFiveNearestdist_oppositCB(){
		//下がるか同じくらい
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
						
					//if(ii < 2){ 自分の O がバックボーンとコンタクトするかを見ようとしたがあまり関係なさそうだった
					//	headertext.add(a+"_O_rank"+String.valueOf(ii+1));
					//	headertext.add(a+"_O_rankb"+String.valueOf(ii+1));
					//}
					if(a.equals("O")){
						headertext.add(jcbdist+a+"_rank"+String.valueOf(ii+1));
						headertext.add(jcbdist+a+"_rankb"+String.valueOf(ii+1));
					}
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
					
					headertext.add(a+"_"+jcbdist+"CB_X_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_X_rankb"+String.valueOf(ii+1));
					
				}
			}
			
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> txcb = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> txcb2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					HashSet<PDBResidue> skipper = new HashSet<>();
					
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						
						if(cc ==  null || oo == null || nn == null || ca == null|| cb == null){
							skipper.add(rr);
						}
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						
						if(jcbdist > 0  && ca != null && cb != null){
							PDBAtom jcb = new PDBAtom();
							jcb.atom_code = "C";
							jcb.pdb_atom_code = jcbdist+"CB";
							double xx = cb.loc.x-ca.loc.x;
							double yy = cb.loc.y-ca.loc.y;
							double zz = cb.loc.z-ca.loc.z;
							double dist = cb.distance(ca);
							xx /= dist/jcbdist;
							yy /= dist/jcbdist;
							zz /= dist/jcbdist;
							xx += ca.loc.x;
							yy += ca.loc.y;
							zz += ca.loc.z;
							jcb.loc.set(xx,yy,zz);
							jcb.parent = rr;
							envatoms.add(jcb);
						}else if(cb != null){
							if(jcbdist > 0){
								
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								jcb.loc.set(cb.loc);//CA が無いと計算できないので
								jcb.parent = rr;
								envatoms.add(jcb);
							}else{
								envatoms.add(cb);
							}
						}
						
						
						
						if(jcbdist > 0  && cc != null && oo != null){
							PDBAtom jcb = new PDBAtom();
							jcb.atom_code = "O";
							jcb.pdb_atom_code = jcbdist+"O";
							double xx = oo.loc.x-cc.loc.x;
							double yy = oo.loc.y-cc.loc.y;
							double zz = oo.loc.z-cc.loc.z;
							double dist = oo.distance(cc);
							xx /= dist/jcbdist;
							yy /= dist/jcbdist;
							zz /= dist/jcbdist;
							xx += cc.loc.x;
							yy += cc.loc.y;
							zz += cc.loc.z;
							jcb.loc.set(xx,yy,zz);
							jcb.parent = rr;
							envatoms.add(jcb);
						}else if(oo != null){
							if(jcbdist > 0){
								
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "O";
								jcb.pdb_atom_code = jcbdist+"O";
								jcb.loc.set(oo.loc);//C が無いと計算できないので
								jcb.parent = rr;
								envatoms.add(jcb);
							}
						}
						
						
						
						
						if(!skipper.contains(rr)){
							
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
							
							
							if(tcbdist > 0){
								//-----------PCBに関するもの
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"XCB";
								Point3D xcb = calcReverseCB(rr);
								double xx = xcb.x-ca.loc.x;
								double yy = xcb.y-ca.loc.y;
								double zz = xcb.z-ca.loc.z;
								double dist = xcb.distance(ca.loc);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								
								tcb.parent = rr;
								txcb.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								//-----------PCBに関するもの
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"XCB";
								Point3D xcb = calcReverseCB(rr);
								double xx = xcb.x-ca.loc.x;
								double yy = xcb.y-ca.loc.y;
								double zz = xcb.z-ca.loc.z;
								double dist = xcb.distance(ca.loc);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								
								tcb.parent = rr;
								txcb2.put(rr, tcb);
							}
							
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank_cbx = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						
						PDBAtom xcb = txcb.get(rr);
						PDBAtom xcb2 = txcb2.get(rr);
						
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1 || ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
						}
						
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else if(aa.pdb_atom_code.equals(jcbdist+"O")){
									code = jcbdist+"O_";
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								if(code.equals( jcbdist+"O_")){
									distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								}else{
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
									
								}
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}
										}
									}
								}
							}
							if(tcbdist > 0){
								double badist = xcb.distance(aa);
								double badist2 = xcb2.distance(aa);
								if(badist > mdist){
									continue;
								}
								String code = null;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_X_";
								}else{
									code = aa.pdb_atom_code+"_CB_X_";
								}
								if(!distrank_cbx.containsKey(code)){
									distrank_cbx.put(code,new ArrayList<AtomDistance>());
								}
								distrank_cbx.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								
							}
						}
						
						
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						for(String ss:distrank_cbx.keySet()){
							Collections.sort(distrank_cbx.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("CB_X_rank") > -1){
								if(h.indexOf("CB_X_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank_cbx.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank_cbx.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank_cbx.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank_cbx.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank_cbx.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank_cbx.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_withO(){
		//下がるか同じくらい
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 3;ii++){
					
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
						
					//if(ii < 2){ 自分の O がバックボーンとコンタクトするかを見ようとしたがあまり関係なさそうだった
					//	headertext.add(a+"_O_rank"+String.valueOf(ii+1));
					//	headertext.add(a+"_O_rankb"+String.valueOf(ii+1));
					//}
					if(a.equals("O")){
						headertext.add(jcbdist+a+"_rank"+String.valueOf(ii+1));
						headertext.add(jcbdist+a+"_rankb"+String.valueOf(ii+1));
					}
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 3;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
					
					//if(ii < 2){ 自分の O が他の CB とコンタクトするかを見ようとしたがあまり性能に関係なさそうだった
					
					//	headertext.add(a+"_"+jcbdist+"CB_O_rank"+String.valueOf(ii+1));
					//	headertext.add(a+"_"+jcbdist+"CB_O_rankb"+String.valueOf(ii+1));
					//}
				}
			}
			
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tos = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tos2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					HashSet<PDBResidue> skipper = new HashSet<>();
					
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						
						if(cc ==  null || oo == null || nn == null || ca == null|| cb == null){
							skipper.add(rr);
						}
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						
						if(jcbdist > 0  && ca != null && cb != null){
							PDBAtom jcb = new PDBAtom();
							jcb.atom_code = "C";
							jcb.pdb_atom_code = jcbdist+"CB";
							double xx = cb.loc.x-ca.loc.x;
							double yy = cb.loc.y-ca.loc.y;
							double zz = cb.loc.z-ca.loc.z;
							double dist = cb.distance(ca);
							xx /= dist/jcbdist;
							yy /= dist/jcbdist;
							zz /= dist/jcbdist;
							xx += ca.loc.x;
							yy += ca.loc.y;
							zz += ca.loc.z;
							jcb.loc.set(xx,yy,zz);
							jcb.parent = rr;
							envatoms.add(jcb);
						}else if(cb != null){
							if(jcbdist > 0){
								
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								jcb.loc.set(cb.loc);//CA が無いと計算できないので
								jcb.parent = rr;
								envatoms.add(jcb);
							}else{
								envatoms.add(cb);
							}
						}
						
						
						
						if(jcbdist > 0  && cc != null && oo != null){
							PDBAtom jcb = new PDBAtom();
							jcb.atom_code = "O";
							jcb.pdb_atom_code = jcbdist+"O";
							double xx = oo.loc.x-cc.loc.x;
							double yy = oo.loc.y-cc.loc.y;
							double zz = oo.loc.z-cc.loc.z;
							double dist = oo.distance(cc);
							xx /= dist/jcbdist;
							yy /= dist/jcbdist;
							zz /= dist/jcbdist;
							xx += cc.loc.x;
							yy += cc.loc.y;
							zz += cc.loc.z;
							jcb.loc.set(xx,yy,zz);
							jcb.parent = rr;
							envatoms.add(jcb);
						}else if(oo != null){
							if(jcbdist > 0){
								
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "O";
								jcb.pdb_atom_code = jcbdist+"O";
								jcb.loc.set(oo.loc);//C が無いと計算できないので
								jcb.parent = rr;
								envatoms.add(jcb);
							}
						}
						
						
						
						
						if(!skipper.contains(rr)){
							
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
							
							
							if(tcbdist > 0){
								//-----------Oに関するもの
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "O";
								tcb.pdb_atom_code = tcbdist+"O";
								
								double xx = oo.loc.x-cc.loc.x;
								double yy = oo.loc.y-cc.loc.y;
								double zz = oo.loc.z-cc.loc.z;
								double dist = oo.distance(cc);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += cc.loc.x;
								yy += cc.loc.y;
								zz += cc.loc.z;
								tcb.loc.set(xx,yy,zz);
								
								tcb.loc.set(oo.loc);
								
								tcb.parent = rr;
								tos.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								//-----------Oに関するもの
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "O";
								tcb.pdb_atom_code = tcbdist+"O";

								double xx = oo.loc.x-cc.loc.x;
								double yy = oo.loc.y-cc.loc.y;
								double zz = oo.loc.z-cc.loc.z;
								double dist = oo.distance(cc);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += cc.loc.x;
								yy += cc.loc.y;
								zz += cc.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;
								tos2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank_o = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						
						PDBAtom o = tos.get(rr);
						PDBAtom o2 = tos2.get(rr);
						
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1 || ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
						}
						
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else if(aa.pdb_atom_code.equals(jcbdist+"O")){
									code = jcbdist+"O_";
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								if(code.equals( jcbdist+"O_")){
									distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								}else{
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
									
								}
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}
										}
									}
								}
							}
							if(tcbdist > 0){
								double badist = o.distance(aa);
								double badist2 = o2.distance(aa);
								if(badist > mdist){
									continue;
								}
								String code = null;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_O_";
								}else{
									code = aa.pdb_atom_code+"_O_";
								}
								if(!distrank_o.containsKey(code)){
									distrank_o.put(code,new ArrayList<AtomDistance>());
								}
								distrank_o.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								
							}
						}
						
						
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						for(String ss:distrank_o.keySet()){
							Collections.sort(distrank_o.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_O_rank") > -1){
								if(h.indexOf("_O_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank_o.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank_o.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank_o.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank_o.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank_o.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank_o.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void calcAtomDistance_TopFiveNearestdist_FourAtoms_distdiff(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			//double jcbdist2 = 3.0;
			double tcbdist = 2.5;
			double tcbdist2 = 3.75;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								if(jcbflag){
									double badist2 = b2.distance(tcbs2.get(aa.parent));
									distrank.get(code).add(new AtomDistance(aa,badist,badist-badist2));
								}else{
									double badist2 = b2.distance(aa);
									distrank.get(code).add(new AtomDistance(aa,badist,badist-badist2));
								}
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<AtomDistance>> distrank = new HashMap<>();
						//HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								double badist2 = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<AtomDistance>());
								}
								distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2-badist));一番いい
								//distrank.get(code).add(new AtomDistance(aa,badist2,badist2-badist));かなり下がる
								//distrank.get(code).add(new AtomDistance(aa,badist,badist2));ちょっと下がる
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss),new AtomDistanceComparator());
						}
						
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance2))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank).distance))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	public static void calcAtomDistance_TopFiveNearestdist_TwoAtoms(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,5,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			double tcbdist2 = 5.0;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			//for(String a:atomnames){
			//	for(double d:distthreshold){
			//		headertext.add(a+"_"+d);
			//	}
			//}
			for(String a:atomnames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_rank"+String.valueOf(ii+1));
					headertext.add(a+"_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
					headertext.add(a+"_"+jcbdist+"CB_rankb"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					HashMap<PDBResidue,PDBAtom> tcbs2 = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
							if(tcbdist2 > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist2+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist2;
								yy /= dist/tcbdist2;
								zz /= dist/tcbdist2;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs2.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<Double>> distrank = new HashMap<>();
						HashMap<String,ArrayList<Double>> distrank2 = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						PDBAtom b2 = tcbs2.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							
							if(tcbdist > 0){
								double badist = b.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank.containsKey(code)){
									distrank.put(code,new ArrayList<Double>());
								}
								distrank.get(code).add(badist);
								double prevdist = 0;
								for(double dd:distthreshold){
									if(badist < dd){
									//if(badist >= prevdist && badist < dd){ //性能下がる
										prevdist = dd;
										String dcode = code+dd;
										if(contactcount.containsKey(dcode)){
											contactcount.put(dcode,1+contactcount.get(dcode));
										}
										if(jcbflag){

											String zcode = jcbdist+"CB_ALL_"+dd;
											if(contactcount.containsKey(zcode)){
												contactcount.put(zcode,1+contactcount.get(zcode));
											}
											if(groupmap.containsKey(aa.parent.getName())){
												String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
												if(contactcount.containsKey(gcode)){
													contactcount.put(gcode,1+contactcount.get(gcode));
												}
											}

										}
									}
								}
							}
							if(tcbdist2 > 0){
								double badist = b2.distance(aa);
								if(badist > mdist){
									continue;
								}


								String code = null;
								boolean jcbflag = false;
								if(aa.pdb_atom_code.equals(jcbdist+"CB")){
									code = aa.parent.getName()+"_"+jcbdist+"CB_";
									jcbflag = true;
								}else{
									code = aa.pdb_atom_code+"_";
								}
								if(!distrank2.containsKey(code)){
									distrank2.put(code,new ArrayList<Double>());
								}
								distrank2.get(code).add(badist);
							}
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss));
						}
						
						for(String ss:distrank2.keySet()){
							Collections.sort(distrank2.get(ss));
						}
						
						Pattern rpat = Pattern.compile("^(.+_)rankb?([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_rank") > -1){
								if(h.indexOf("_rankb") > -1){
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank2.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank2.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank2.get(code).get(rank)))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}else{
									Matcher mmat = rpat.matcher(h);
									if(mmat.find()){
										String code = mmat.group(1);
										int rank = Integer.parseInt(mmat.group(2))-1;
										if(!distrank.containsKey(code)){
											pw.write("1000\t");
										}else{
											if(distrank.get(code).size() > rank){
												pw.write(String.valueOf((float)((double)distrank.get(code).get(rank)))+"\t");
											}else{
												pw.write("1000\t");
											}
										}
									}else{
										System.err.println(h+" cannot be parsed.");
										throw new RuntimeException("???");
									}
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	public static void calcAtomDistance_TopFiveNearestdist(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,4,5,6,9,12};//最終Threshold 未満の Atom について、5 個までキョリの短い順に入る
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			
			
			
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			for(String a:atomnames){
				for(double d:distthreshold){
					headertext.add(a+"_"+d);
				}
			}
			
			for(String a:aanames){
				for(int ii = 0;ii < 5;ii++){
					headertext.add(a+"_"+jcbdist+"CB_rank"+String.valueOf(ii+1));
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				if(d.indexOf("_rank") == -1){
					contactcount_all.put(d,0);
				}
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						HashMap<String,ArrayList<Double>> distrank = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						for(String ss:headertext){
							if(ss.indexOf(jcbdist+"CB_rank") > -1){
							}else{
								contactcount.put(ss,0);
							}
									
						}
						double mdist = distthreshold[distthreshold.length-1];
						for(PDBAtom aa:envatoms){
							if(aa.isAlternative() && !aa.getAltCode().equals("A")){
								continue;
							}
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							double badist = b.distance(aa);
							if(badist > mdist){
								continue;
							}
							
							
							String code = null;
							boolean jcbflag = false;
							if(aa.pdb_atom_code.equals(jcbdist+"CB")){
								code = aa.parent.getName()+"_"+jcbdist+"CB_";
								jcbflag = true;
							}else{
								code = aa.pdb_atom_code+"_";
							}
							if(!distrank.containsKey(code)){
								distrank.put(code,new ArrayList<Double>());
							}
							distrank.get(code).add(badist);
							double prevdist = 0;
							for(double dd:distthreshold){
								if(badist < dd){
								//if(badist >= prevdist && badist < dd){ //性能下がる
									prevdist = dd;
									String dcode = code+dd;
									if(contactcount.containsKey(dcode)){
										contactcount.put(dcode,1+contactcount.get(dcode));
									}
									if(jcbflag){
										
										String zcode = jcbdist+"CB_ALL_"+dd;
										if(contactcount.containsKey(zcode)){
											contactcount.put(zcode,1+contactcount.get(zcode));
										}
										if(groupmap.containsKey(aa.parent.getName())){
											String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
											if(contactcount.containsKey(gcode)){
												contactcount.put(gcode,1+contactcount.get(gcode));
											}
										}
										
									}
								}
							}
							
						}
						for(String ss:distrank.keySet()){
							Collections.sort(distrank.get(ss));
						}
						
						Pattern rpat = Pattern.compile("^(.+_)rank([0-9]+)");
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else if(h.indexOf("_rank") > -1){
								Matcher mmat = rpat.matcher(h);
								if(mmat.find()){
									String code = mmat.group(1);
									int rank = Integer.parseInt(mmat.group(2))-1;
									if(!distrank.containsKey(code)){
										pw.write("1000\t");
									}else{
										if(distrank.get(code).size() > rank){
											pw.write(String.valueOf(distrank.get(code).get(rank))+"\t");
										}else{
											pw.write("1000\t");
										}
									}
								}else{
									System.err.println(h+" cannot be parsed.");
									throw new RuntimeException("???");
								}
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							if(contactcount_all.containsKey(h) && contactcount.containsKey(h)){
								contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
							}
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void calcAtomDistance_nearestdist(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {2,3,4,6,9,10};//最終項目については、最短距離が入っている
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			for(String a:atomnames){
				for(double d:distthreshold){
					headertext.add(a+"_"+d);
				}
			}
			
			for(String a:aanames){
				for(double d:distthreshold){
					headertext.add(a+"_"+jcbdist+"CB_"+d);
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				if(Pattern.compile(jcbdist+"CB_"+distthreshold[distthreshold.length-1]+"$")
						.matcher(d).find()){
					headertext.add(d+"_dist");
				}else{
				pw.write(d+"\t");
				}
				contactcount_all.put(d,0);
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						HashMap<String,Double> contactdist = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						for(String ss:headertext){
							contactcount.put(ss,0);
							if(ss.indexOf(distthreshold[distthreshold.length-1]+"CB") > -1 && ss.indexOf("CB_ALL") < 0){
								contactdist.put(ss+"_dist",1000.0);
							}
						}
						double mdist = distthreshold[distthreshold.length-1]+1;
						for(PDBAtom aa:envatoms){
							
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							double badist = b.distance(aa);
							if(badist > mdist){
								continue;
							}
							String code = null;
							boolean jcbflag = false;
							if(aa.pdb_atom_code.equals(jcbdist+"CB")){
								code = aa.parent.getName()+"_"+jcbdist+"CB_";
								jcbflag = true;
							}else{
								code = aa.pdb_atom_code+"_";
							}
							double prevdist = 0;
							for(double dd:distthreshold){
								if(badist < dd){
								//if(badist >= prevdist && badist < dd){ //性能下がる
									prevdist = dd;
									String dcode = code+dd;
									contactcount.put(dcode,1+contactcount.get(dcode));
									if(contactdist.containsKey(dcode+"_dist") && contactdist.get(dcode+"_dist") > badist){
										contactdist.put(dcode+"_dist", badist);
									}
									if(jcbflag){
										
										String zcode = jcbdist+"CB_ALL_"+dd;
										if(contactcount.containsKey(zcode)){
											contactcount.put(zcode,1+contactcount.get(zcode));
										}
										if(groupmap.containsKey(aa.parent.getName())){
											String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
											if(contactcount.containsKey(gcode)){
												contactcount.put(gcode,1+contactcount.get(gcode));
											}
										}
										
									}
								}
							}
							
						}
						for(String h:headertext){
							if(contactdist.containsKey(h+"_dist")){
								pw.write(contactdist.get(h+"_dist")+"\t");
							}else{
								pw.write(contactcount.get(h)+"\t");
							}
							contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void calcAtomDistance(){
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
				new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 

			double[] distthreshold = {1,2,3,4,5,6,7,8,9,10};
			double jcbdist = 2.5;
			double tcbdist = 2.5;
			
			boolean ignoreneighbor = false;//一つ隣を使わないようにすると性能が下がった
			double covthreshold_a = 2.2;//共有結合とみなすキョリ リスト上で隣にある場合
			double covthreshold_b = 1.6;//共有結合とみなすキョリ　リスト上で隣にない場合
			Pattern pat = Pattern.compile("\\.(pdb|ent)$");
			ArrayList<String> headertext = new ArrayList<>();
			String[] atomnames = {
				"C","CA","N","O"
			};
			String[] aanames = {
				"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
				,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
			};
			
			String[] groupnames = {
				"FYW"
				,"AG"
				,"RHK"
				,"DE"
				,"STNQ"
			};
			HashMap<String,String> groupmap = new HashMap<>();
			for(String s:groupnames){
				char[] cc = s.toCharArray();
				for(char c:cc){
					groupmap.put(PepProcess.one_to_three.get(c), s);
				}
			}
			
			HashMap<PDBResidue,HashSet<PDBResidue>> res_ignore = new HashMap<>();
			
			
			
			
			
			HashMap<String,Integer> contactcount_all = new HashMap<>();
			HashSet<String> validaas = new HashSet<>();
			for(String s:aanames){
				validaas.add(s);
			}
			for(String a:atomnames){
				for(double d:distthreshold){
					headertext.add(a+"_"+d);
				}
			}
			
			for(String a:aanames){
				for(double d:distthreshold){
					headertext.add(a+"_"+jcbdist+"CB_"+d);
				}
			}
			
			for(double d:distthreshold){
				headertext.add(jcbdist+"CB_ALL_"+d);
			}
			/*
			for(String s:groupnames){
				for(double d:distthreshold){
					headertext.add(s+"_"+jcbdist+"CB_"+d);
				}
			}
			*/
			
			pw.write("name\t");
			for(String d:headertext){
				pw.write(d+"\t");
				contactcount_all.put(d,0);
			}
			pw.write("target\n");
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				boolean invalidaa = false;
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					for(PDBResidue rr:c.residues){
						if(rr.isLigand() || rr.isMissing()){
							continue;
						}
						if(!validaas.contains(rr.getName())){
							invalidaa = true;
						}
					}
				}
				if(invalidaa){
					System.out.println(f.getName()+" has irregular aa.");
					continue;
				}
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					for(int ii = 0;ii < c.residues.size();ii++){
						res_ignore.put(c.residues.get(ii),new HashSet<PDBResidue>());
					}

					if(ignoreneighbor){
						for(int ii = 1;ii < c.residues.size();ii++){
							PDBResidue prev = null;
							PDBResidue curr = c.residues.get(ii);
							if(ii > 0){
								prev = c.residues.get(ii-1);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_a){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
						
						
						for(int ii = 0;ii < c.residues.size();ii++){
							for(int jj = 0;jj < c.residues.size();jj++){
								if(ii == jj){//無くてもいい気がする
									continue;
								}
								PDBResidue curr = c.residues.get(ii);
								PDBResidue prev = c.residues.get(jj);
								if(prev.getC() != null){
									if(curr.getN() != null){
										if(prev.getC().distance(curr.getN()) < covthreshold_b){
											res_ignore.get(prev).add(curr);
											res_ignore.get(curr).add(prev);
										}
									}
								}
							}
						}
					}
				}
				
				
				
				for(String cname:p.chains.keySet()){
					PDBChain c = p.chains.get(cname);
					ArrayList<PDBAtom> allatoms = new ArrayList<>();
					HashMap<PDBResidue,PDBAtom> tcbs = new HashMap<>();
					ArrayList<PDBAtom> backbones = new ArrayList<>();
					ArrayList<PDBAtom> envatoms = new ArrayList<>();
					for(PDBResidue rr:c.residues){
						if(rr.isLigand()){
							continue;
						}
						if(rr.isMissing()){
							continue;
						}
						PDBAtom oo = rr.getO();
						PDBAtom cc = rr.getC();
						PDBAtom nn = rr.getN();
						PDBAtom ca = rr.getCA();
						PDBAtom cb = rr.getCB();
						if(cc != null){
							envatoms.add(cc);
						}
						
						if(oo != null){
							envatoms.add(oo);
						}
						
						if(nn != null){
							envatoms.add(nn);
						}
						if(ca != null){
							envatoms.add(ca);
						}
						if(cb != null && ca != null){
							if(jcbdist > 0){
								PDBAtom jcb = new PDBAtom();
								jcb.atom_code = "C";
								jcb.pdb_atom_code = jcbdist+"CB";
								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/jcbdist;
								yy /= dist/jcbdist;
								zz /= dist/jcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								jcb.loc.set(xx,yy,zz);
								jcb.parent = rr;
								envatoms.add(jcb);

							}
							
							if(tcbdist > 0){
								PDBAtom tcb = new PDBAtom();
								tcb.atom_code = "C";
								tcb.pdb_atom_code = tcbdist+"CB";

								double xx = cb.loc.x-ca.loc.x;
								double yy = cb.loc.y-ca.loc.y;
								double zz = cb.loc.z-ca.loc.z;
								double dist = cb.distance(ca);
								xx /= dist/tcbdist;
								yy /= dist/tcbdist;
								zz /= dist/tcbdist;
								xx += ca.loc.x;
								yy += ca.loc.y;
								zz += ca.loc.z;
								tcb.loc.set(xx,yy,zz);
								tcb.parent = rr;

								tcbs.put(rr, tcb);
							}
						}
					}
					for(PDBResidue rr:tcbs.keySet()){
						
						HashMap<String,Integer> contactcount = new HashMap<>();
						pw.write(p.id+"_"+cname+"_"+rr.getRepresentativeCode()+"_"+tcbdist+"CB\t");
						PDBAtom b = tcbs.get(rr);
						for(String ss:headertext){
							contactcount.put(ss,0);
						}
						double mdist = distthreshold[distthreshold.length-1]+1;
						for(PDBAtom aa:envatoms){
							
							if(aa.parent == rr){
								continue;
							}
							if(res_ignore.get(rr).contains(aa.parent)){
								continue;
							}
							double badist = b.distance(aa);
							if(badist > mdist){
								continue;
							}
							String code = null;
							boolean jcbflag = false;
							if(aa.pdb_atom_code.equals(jcbdist+"CB")){
								code = aa.parent.getName()+"_"+jcbdist+"CB_";
								jcbflag = true;
							}else{
								code = aa.pdb_atom_code+"_";
							}
							double prevdist = 0;
							for(double dd:distthreshold){
								if(badist < dd){
								//if(badist >= prevdist && badist < dd){ //性能下がる
									prevdist = dd;
									String dcode = code+dd;
									contactcount.put(dcode,1+contactcount.get(dcode));
									if(jcbflag){
										
										String zcode = jcbdist+"CB_ALL_"+dd;
										if(contactcount.containsKey(zcode)){
											contactcount.put(zcode,1+contactcount.get(zcode));
										}
										if(groupmap.containsKey(aa.parent.getName())){
											String gcode = groupmap.get(aa.parent.getName())+"_"+jcbdist+"CB_"+dd;
											if(contactcount.containsKey(gcode)){
												contactcount.put(gcode,1+contactcount.get(gcode));
											}
										}
										
									}
								}
							}
							
						}
						for(String h:headertext){
							pw.write(contactcount.get(h)+"\t");
							contactcount_all.put(h,contactcount_all.get(h)+contactcount.get(h));
						}
						pw.write(rr.getName()+"_"+tcbdist+"CB\n");
					}
				}
				
			}
			for(String h:headertext){
				System.out.println(h+"\t"+contactcount_all.get(h));
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}

	public static void main(String[] args){
		//calcAtomDistance_nearestdist();
		//calcAtomDistance_TopFiveNearestdist();
		//calcAtomDistance_TopFiveNearestdist_TwoAtoms();
		//calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff();
		calcAtomDistance_TopFiveNearestdist_TwoAtoms_distdiff_phipsi_removing();
		//calcAtomDistance_TopFiveNearestdist_FourAtoms_distdiff();
	}
}

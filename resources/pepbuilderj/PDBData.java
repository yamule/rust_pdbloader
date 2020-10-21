package pepbuilderj;


import java.util.*;
import java.io.*;
import java.util.regex.*;


public class PDBData{
	public HashMap<String,PDBChain> chains = new HashMap<String,PDBChain>();
	public ArrayList<AtomConnection> connections = new ArrayList<>();
	public String title = "";
	public String id = "";
	public ArrayList<TransformUnit> transUnit = new ArrayList<>();
	public boolean scaleSection = false;
	public double scaleMatrix[][] = new double[3][4];
	public double crystVector[] = new double[3];
	PDBData(){
	
	}
	PDBData(PDBData p){
		scaleSection = p.scaleSection;
		for(int xx = 0;xx < scaleMatrix.length;xx++){
			for(int yy = 0;yy < scaleMatrix[0].length;yy++){
				scaleMatrix[xx][yy] = p.scaleMatrix[xx][yy];
			}	
		}
		crystVector[0] = p.crystVector[0];
		crystVector[1] = p.crystVector[1];
		crystVector[2] = p.crystVector[2];
		Iterator<String> ite = p.chains.keySet().iterator();
		while(ite.hasNext()){
			PDBChain c = p.chains.get(ite.next());
			PDBChain newc = c.copy();
			chains.put(newc.name,newc);
		}
		if(p.connections.size() > 0){
			HashMap<String,PDBAtom> nummap = new HashMap<String,PDBAtom>();
			for(String s:chains.keySet()){
				PDBChain c = chains.get(s);
				for(PDBResidue r:c.residues){
					for(PDBAtom a:r.atoms){
						if(nummap.containsKey(a.atom_number)){
							System.err.println("Duplicated atom id found.");
						}else{
							nummap.put(a.atom_number, a);
						}
					}
				}
			}
			for(AtomConnection ac:p.connections){
				if(!nummap.containsKey(ac.atomA.atom_number) || !nummap.containsKey(ac.atomB.atom_number)){
					System.err.println("Atom connection discrepancy found.");
				}else{
					connections.add(new AtomConnection(nummap.get(ac.atomA.atom_number),nummap.get(ac.atomB.atom_number)));
				}
			}
		}
	}
	public void setID(String i){
		id = i;
		
	}
	public String getID(){
		return id;
	}
	public void setTitle(String t){
		title = t;
	}
	public String getTitle(){
		return title;
	}
	
	/**
	 * Chain の DESCRIPTION を設定して表示する
	 */
	public static void mapChainName(String filename,PDBData pdb){
		try{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = null;
			Pattern pat = Pattern.compile("^COMPND[\\s]+[0-9]*[\\s]*(.+)[\\s]*$");
			Pattern chainpat = Pattern.compile("^COMPND[\\s]+[0-9]*[\\s]*CHAIN\\:[\\s]*([^\\s;]+);?");
			String descs = "";
			while((line = br.readLine()) != null){
				if(line.indexOf("COMPND") == 0){
					while(line.indexOf("COMPND") == 0){
						Matcher cmat = chainpat.matcher(line);
						if(cmat.find()){
							String c = cmat.group(1);
							if(pdb.chains.containsKey(c)){
								pdb.chains.get(c).setDesc(descs.replaceAll("[\\s]+"," "));
							}
							descs = "";
							while(line.indexOf("MOL_ID") == -1){
								line = br.readLine();
								if(line.indexOf("COMPND") != 0){
									break;
								}
							}
							
						}else{
							Matcher pmat = pat.matcher(line);
							 if(pmat.find()){
							 	descs += pmat.group(1)+" ";
							 }
							 
							line = br.readLine();
						}
					}
				}
				if(line.indexOf("ENDMDL") == 0){
					break;
				}
				
			}
			
			
			Iterator<String> ite = pdb.chains.keySet().iterator();
			while(ite.hasNext()){
				PDBChain c = pdb.chains.get(ite.next());
				System.out.println(c.desc+";;"+c.name);
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	
	public static PDBData loadPDBFile(String filename){
		PDBData ret = new PDBData();
		try{
			System.out.println(filename);
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = null;
			Pattern mispat = Pattern.compile("REMARK 465[\\s]+M[\\s]+RES[\\s]+C[\\s]+SSSEQI");
			ArrayList<String> transline = new ArrayList<>();
			ArrayList<String> scaleline = new ArrayList<>();
			ArrayList<String> crystline = new ArrayList<>();
			
			//REMARK 465     ALA B     1
			Pattern mispat2 = Pattern.compile("REMARK 465[\\s]+([A-Z]+)[\\s]+([A-Za-z0-9]+)[\\s]+([\\-0-9A-Z]+)");
			Pattern inspat = Pattern.compile("([\\-0-9]+)([A-Z]+)");
			// 1 -  6 Record name "ATOM "
			// 7 - 11 Integer serial Atom serial number.
			//13 - 16 Atom name Atom name.
			//17      Character altLoc Alternate location indicator.
			//18 - 20 Residue name resName Residue name.
			//22      Character chainID Chain identifier.
			//23 - 26 Integer resSeq Residue sequence number.
			//27      AChar iCode Code for insertion of residues.
			//31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.
			//39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.
			//47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.
			//55 - 60 Real(6.2) occupancy Occupancy.
			//61 - 66 Real(6.2) tempFactor Temperature factor.
			//77 - 78 LString(2) element Element symbol, right-justified.
			//79 - 80 LString(2) charge Charge on the atom.
			String tit = "";
			HashMap<Integer,PDBAtom> atommap = new HashMap<>();
			ArrayList<ArrayList<Integer>> co = new ArrayList<>();
			boolean connectionerrorflag = false;
			HashSet<String> terflag = new HashSet<>();
			while((line = br.readLine()) != null){
				
				if(line.length() < 6){
					continue;
				}
				
				String linelabel = line.substring(0,6);
				if(linelabel.equals("HEADER")){
					ret.setID(line.substring(62,66));
					continue;
				}
				if(linelabel.equals("ENDMDL")){
					break;
				}
				if(linelabel.equals("TITLE ")){
					tit += line.substring(10);
				}
				if(linelabel.equals("CONECT")){
					char[] ss = line.toCharArray();
					ArrayList<Integer> cal = new ArrayList<>();
					for(int ii = 6;ii < ss.length;){
						int start = ii;
						StringBuffer sb = new StringBuffer();
						for(;ii < start+5;ii++){
							if(ii < ss.length){
								sb.append(ss[ii]);
							}
						}
						String ls = sb.toString().replaceAll("[\\s]","");
						if(ls.length() > 0){
							cal.add(Integer.parseInt(ls));
						}
					}
					if(cal.size() > 0){
						co.add(cal);
					}
				}
				
				if(linelabel.equals("CRYST1")){
					crystline.add(line);
				}
				
				if(linelabel.equals("SCALE1")){
					scaleline.add(line);
				}
				if(linelabel.equals("SCALE2")){
					scaleline.add(line);
				}
				
				if(linelabel.equals("SCALE3")){
					scaleline.add(line);
				}
				if(linelabel.equals("REMARK")){
					if(line.indexOf("REMARK 350") == 0){//missing 部分の処理
						transline.add(line);
					}
					if(line.indexOf("REMARK 465") == 0){//missing 部分の処理
					
						boolean remflag = false;
						while((line = br.readLine()) != null){


							if(!remflag && mispat.matcher(line).find()){

								remflag = true;
							}else if(remflag){
								Matcher mat = mispat2.matcher(line);
								if(mat.find()){
									String rescode = mat.group(1);
									String chain = mat.group(2);
									String pos = mat.group(3);
									String inscode = "";

									Matcher pmat = inspat.matcher(pos);
									if(pmat.find()){
										pos = pmat.group(1);
										inscode = pmat.group(2);
									}



									if(ret.chains.get(chain) == null){
										ret.chains.put(chain,new PDBChain(ret,chain));
									}
									//PDBResidue(String n,int nn,String alt,String ins,boolean miss){
									ret.chains.get(chain).addMissingResidue(new PDBResidue(rescode,Integer.parseInt(pos),"",inscode,true));
								}else{
									//System.err.println(line+" was not processed.");
								}


							}

							if(line.indexOf("REMARK 465") != 0){
								linelabel = line.substring(0,6);
								break;
							}
						}
					}
					
				}
				//TER    1518      PHE B
				if(linelabel.equals("TER   ")){
					String chain_id = line.substring(21,22);
					terflag.add(chain_id);
				}
				if(linelabel.equals("ATOM  ")||linelabel.equals("HETATM")){
					String record_name = line.substring(0,6);
					String atom_serial_num = line.substring(6,11);
					String atom_name = line.substring(12,16);
					String alternate_location = line.substring(16,17);
					String residue_name = line.substring(17,20).toUpperCase();
					String chain_id = line.substring(21,22);
					String residue_sequence_number = line.substring(22,26);
					String insertion_code = line.substring(26,27);
					String x_coord = line.substring(30,38);
					String y_coord = line.substring(38,46);
					String z_coord = line.substring(46,54);
					String occupancy = "";
					boolean possibleligand = false;
					if(linelabel.equals("HETATM") && terflag.contains(chain_id)){
						possibleligand = true;
					}
					
					
					String temperature_factor = "";
					if(line.length() >= 60){
						occupancy = line.substring(54,60);
					}
					if(line.length() >= 61){
						temperature_factor = line.substring(60,Math.min(line.length(),66));
					}
					String element_symbol = "";
					if(line.length() >= 78){
						 element_symbol = line.substring(76,78).replaceAll(" ","");;
					}
					//String charge = line.substring(78,80);
					if(element_symbol.length() > 0){
					}else{
						element_symbol = atom_name.substring(0,2).replaceAll(" ","");
					}
					PDBAtom newp= new PDBAtom(atom_name
							,element_symbol
							,atom_serial_num
							,new Point3D(Float.parseFloat(x_coord)
							,Float.parseFloat(y_coord)
							,Float.parseFloat(z_coord)));
					newp.ligandFlag(possibleligand);
					newp.setLine(line.replaceAll("[\\r\\n]","")+"\n");
					if(temperature_factor.length() > 0){
						newp.setBFactor(Double.parseDouble(temperature_factor.replaceAll("[\\s]","")));
					}
					if(occupancy.length() > 0){
						newp.setOccupancy(Double.parseDouble(occupancy.replaceAll("[\\s]","")));
					}
					String si = atom_serial_num.replaceAll("[\\s]","");
					int acode = 0;
					if(si.length() > 0){
						acode = Integer.parseInt(si);
					}
					if(atommap.containsKey(acode)){
						//System.err.println("Atom code duplication\n"+line);
						connectionerrorflag = true;
					}else{
						atommap.put(acode,newp);
					}
					
					if( line.indexOf("HETATM") == 0){
						newp.setHETATM(true);
					}
					int residue_number = -1000;
					String altflag = "";
					try{
						residue_number  = Integer.parseInt(residue_sequence_number.replaceAll(" ",""));
					}catch(Exception exx){
						residue_number  = Integer.parseInt(residue_sequence_number.replaceAll("[^0-9\\-]",""));
						altflag = residue_sequence_number.replaceAll(" ","");
					}
					if(!alternate_location.equals(" ")){
						newp.setAltCode(alternate_location);
					}
					if(alternate_location.equals("A") || alternate_location.equals("1") || alternate_location.equals(" ")){
					}else{
						altflag = alternate_location;
					}
					
					if(ret.chains.get(chain_id) == null){
						ret.chains.put(chain_id,new PDBChain(ret,chain_id));
					}
					if(altflag.length() == 0){
						if(insertion_code.replaceAll("[\\s]","").length() == 0){
							ret.chains.get(chain_id).addAtom(newp,residue_name,residue_number);
						}else{
							ret.chains.get(chain_id).addSpecialAtom(newp,residue_name,residue_number,altflag,insertion_code);
						}
					}else{
						ret.chains.get(chain_id).addSpecialAtom(newp,residue_name,residue_number,altflag,insertion_code);
					}
				}
			}
			if(!connectionerrorflag && co.size() > 0 ){
				outer:for(ArrayList<Integer> i:co){
					for(int ii = 1;ii < i.size();ii++){
						try{
							ret.connections.add(new AtomConnection(atommap.get(i.get(0)),atommap.get(i.get(ii))));
						}catch(Exception e){
							ret.connections.clear();
							System.err.println("Error in connection field.");
							break outer;
						}
					}
				}
			}
			ret.setTitle(tit.replaceAll("^[\\s]+","").replaceAll("[\\s]+"," ").replaceAll("[\\s]+$",""));
			
		
			if(scaleline.size() > 0){
				if(scaleline.size() != 3){
					System.err.println("More than three SCALEn lines are found. The first 3 lines are used. ");
					
				}
				//SCALE1      0.000000 -0.012994  0.000000        0.23250                         
				//SCALE2      0.007258  0.000000  0.009856        0.00000                         
				//SCALE3     -0.008693  0.000000  0.006402        0.00000                         

				// 1 -  6         Record name   "SCALEn" n=1,  2, or 3
				//11 - 20         Real(10.6)    s[n][1]            Sn1
				//21 - 30         Real(10.6)    s[n][2]            Sn2
				//31 - 40         Real(10.6)    s[n][3]            Sn3
				//46 - 55         Real(10.5)    u[n]               Un
				ret.scaleSection = true;
				Pattern scalepat = Pattern.compile("SCALE[0-9].{4}(.{10})(.{10})(.{10}).{5}(.{10})");
				for(int ii = 0;ii < 3;ii++){
					String l = scaleline.get(ii);
					if(l.indexOf("SCALE"+String.valueOf(ii+1)) == -1){
						throw new RuntimeException("Confusion in SCALEn lines??? Please check.");
						
					}else{
						Matcher mmat = scalepat.matcher(l);
						if(mmat.find()){
							double d1 = Double.parseDouble(mmat.group(1).replaceAll(" ",""));
							double d2 = Double.parseDouble(mmat.group(2).replaceAll(" ",""));
							double d3 = Double.parseDouble(mmat.group(3).replaceAll(" ",""));
							double d4 = Double.parseDouble(mmat.group(4).replaceAll(" ",""));
							ret.scaleMatrix[ii][0] = d1;
							ret.scaleMatrix[ii][1] = d2;
							ret.scaleMatrix[ii][2] = d3;
							ret.scaleMatrix[ii][3] = d4;
						}else{
							System.err.println("Line "+l+" is not match the pattern "+"SCALE[0-9].{4}(.{10})(.{10})(.{10}).{5}(.{10})");
							System.err.println("Please check.");
							System.exit(-1);
						}
					}
				}
			}
			//1 -  6      Record name          "CRYST1"
			//7 - 15      Real(9.3)            a            a (Angstroms).
			//16 - 24      Real(9.3)            b            b (Angstroms).
			//25 - 33      Real(9.3)            c            c (Angstroms).
			//34 - 40      Real(7.2)            alpha        alpha (degrees).
			//41 - 47      Real(7.2)            beta         beta (degrees).
			//48 - 54      Real(7.2)            gamma        gamma (degrees).
			//56 - 66      LString              sGroup       Space group.
			//67 - 70      Integer              z            Z value.
			
			if(crystline.size() > 0){
				if(crystline.size() > 1){
					System.err.println("multiple cryst1 entries were found. The first was used.");
				}
				//Pattern spat = Pattern.compile("CRYST1(.{9})(.{9})(.{9})(.{7})(.{7})(.{7})(.{11})(.{4})");
				Pattern spat = Pattern.compile("CRYST1(.{9})(.{9})(.{9})(.{7})(.{7})(.{7})");
				Matcher smat = spat.matcher(crystline.get(0));
				if(smat.find()){
					ret.crystVector[0] = Double.parseDouble(smat.group(1).replaceAll(" ",""));
					ret.crystVector[1] = Double.parseDouble(smat.group(2).replaceAll(" ",""));
					ret.crystVector[2] = Double.parseDouble(smat.group(3).replaceAll(" ",""));
				}else{
					System.err.println("Cannot parse cryst1 line. The line was ignored.");
				}
				
			}
			
			if(transline.size() > 0){
				Pattern startpat = Pattern.compile("REMARK 350 BIOMOLECULE[\\s]*\\:[\\s]+([^\\s]+)");				TransformUnit currentunit = null;
				TransformSet currentset = null;
				for(int ii = 0;ii < transline.size();ii++){
					String str = transline.get(ii);
					Matcher smat = startpat.matcher(str);
					if(smat.find()){
						currentunit = new TransformUnit(smat.group(1));
						ret.transUnit.add(currentunit);
					}else if(str.indexOf("TO CHAINS:") > -1){
						StringBuffer cbuff = new StringBuffer();
						int st = ii;
						for(;ii < Math.min(transline.size(), st+5);ii++){//5行までしか対応していないが、複数行にわたる場合に対処しようとしている。
							String sstr = transline.get(ii);
							if(sstr.indexOf("BIOMT") > -1){
								ii--;
								break;
							}
							sstr = sstr.replaceAll(".+\\:","");
							cbuff.append(sstr+" ");
							
						}
						if(cbuff.length() == 0){
							System.err.println("FAILED TO LOAD BIOUNIT (REMARK 350 LINES).");
						}else{
							currentset = new TransformSet();
							currentunit.addTransformSet(currentset);
						}
						String[] pt = cbuff.toString().split(" ");
						Pattern npat = Pattern.compile("([a-zA-Z0-9]+)");
						for(String pp:pt){
							Matcher m = npat.matcher(pp);
							if(m.find()){
								PDBChain pc = ret.chains.get(m.group(1));
								if(pc != null){
									currentset.addChain(pc);
								}else{
									System.err.println("Chain "+m.group(1)+" was not found.");
								}
							}
						}
					}else{
						Pattern matpat = Pattern.compile("BIOMT([0-9]+)[\\s]+([0-9]+)[\\s]+([^\\s].+)");
						Pattern chkpat = Pattern.compile("REMARK +350 +.*[0-9][0-9]");
						int currentrow = 0;
						ArrayList<String> dline = new ArrayList<>();
						for(;ii < transline.size();ii++){
							String sst = transline.get(ii);
							
							if(sst.indexOf("TO CHAINS:") > -1 || startpat.matcher(sst).find()){
								ii--;
								
								break;
							}
							Matcher mmat = matpat.matcher(sst);
							if(mmat.find()){
								String code = mmat.group(1);
								String val = mmat.group(3);
								if(Integer.parseInt(code) != currentrow+1){
									System.err.println((currentrow+1)+" was expected.");
									System.err.print(transline.get(ii));
									continue;
								}
								dline.add(val);
								currentrow++;
								if(currentrow == 3){
									TransformMatrix tu = new TransformMatrix();
									for(int jj = 0;jj < dline.size();jj++){
										String[] s = dline.get(jj).split("[\\s]+");
										int cnum = 0;
										for(String ss:s){
											String sss = ss.replaceAll("[\\s]","");
											if(sss.length() > 0){
												tu.set(jj,cnum, Double.parseDouble(sss));
												cnum++;
												if(cnum == 4){
													break;
												}
											}
										}
									}
									currentset.addUnit(tu);
									currentrow = 0;
									dline.clear();
								}
								
							}else{
								//if(chkpat.matcher(sst).find()){
								//	System.err.println("Can not parse :");
								//	System.err.println(transline.get(ii));
								//}
							}	
						}
					}
				}
			}
			
			br.close();
		}catch(Exception e){
			System.err.println(filename);
			e.printStackTrace();
		}
		/*
		Iterator<String> ite = ret.chains.keySet().iterator();
		while(ite.hasNext()){
			PDBChain c = ret.chains.get(ite.next());
			//System.out.println(c.name);
			Collections.sort(c.residues, new PDBResidueComparator());
		}
		*/
		
		
		
		//mapChainName(filename,ret);//chain に Description を設定する場合コメントアウトを取る
		
		return ret;
	}
	public static void writeToFile(String filename,StringBuffer sb,boolean append){
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename,append),"UTF-8"))); 
			pw.print(sb);
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	

	
	public PDBData getTransformedEntry(String biomolecule){
		PDBData dum = new PDBData();
		int cou = 0;
		ArrayList<PDBChain> cc = null;
		for(int ii = 0;ii < transUnit.size();ii++){
			if(transUnit.get(ii).id.equals(biomolecule)){
				cc = transUnit.get(ii).getTransformedChains();
			}

		}
		if(cc != null){
			for(PDBChain c:cc){
				dum.chains.put(c.name+cou, c);
				c.parent = dum;
				cou++;
			}
		}else{
			System.err.println("Biomolecule "+biomolecule+" was not found.");
			System.exit(1);
		}
		return dum;
	}
	
	
	
	
	
	/**
	 * Calculate covalent bounds automatically.
	 */
	public void autoConnect(double boundthreshold){
		connections.clear();
		for(String cl:chains.keySet()){
			PDBChain c = chains.get(cl);
			c.removeWater();
			for(int ii = 0;ii < c.residues.size();ii++){
				PDBResidue r1 = c.residues.get(ii);
				for(int jj = 0;jj < r1.atoms.size();jj++){//bounds between atoms in same residue
					PDBAtom a = r1.atoms.get(jj);
					for(int kk = jj+1;kk < r1.atoms.size();kk++){
						PDBAtom b = r1.atoms.get(kk);
						if(a.getAltCode().equals(b.getAltCode()) || a.getAltCode().length() == 0 || b.getAltCode().length() == 0){
							if(a.pdb_atom_code.indexOf("H") == 0 && b.pdb_atom_code.indexOf("H") == 0){
								continue;
							}
							if(a.distance(b) < boundthreshold ){
								connections.add(new AtomConnection(a,b));
							}
						}
					}	
				}
				for(int jj = ii+1;jj < c.residues.size();jj++){//bounds between different residues
					PDBResidue r2 = c.residues.get(jj);
					outer:for(int ll = 0;ll < r1.atoms.size();ll++){
						PDBAtom a = r1.atoms.get(ll);
						for(int kk = 0;kk < r2.atoms.size();kk++){
							PDBAtom b = r2.atoms.get(kk);
							
							if(a.pdb_atom_code.indexOf("H") == 0 && b.pdb_atom_code.indexOf("H") == 0){
								continue;
							}
							if(a.distance(b) - r1.maxDist -r2.maxDist > boundthreshold){
								break outer;
							}
							if(a.distance(b) < boundthreshold){
								connections.add(new AtomConnection(a,b));
							}
						}	
					}
				}
			}
		}
	}
	
	public static void main(String[] args){
	}

}
class PDBResidueComparator implements Comparator<PDBResidue>{
	@SuppressWarnings("unchecked")
	public int compare(PDBResidue arg1, PDBResidue arg2){
		PDBResidue i1 = (PDBResidue)arg1;
		PDBResidue i2 = (PDBResidue)arg2;
		
		if(i1.getResidueNumber() > i2.getResidueNumber()){
			return 1;
		}
		if(i1.getResidueNumber() < i2.getResidueNumber()){
			return -1;
		}
		i1.getInsertionCode().compareTo(i2.getInsertionCode());
		return 0;
	}
	
}



class PDBChain{
	String name = "none";
	String desc = "none";
	PDBData parent = null;
	ArrayList<PDBResidue> missing = new ArrayList<PDBResidue>();//missing な残基が入っている
	ArrayList<PDBResidue> residues = new ArrayList<PDBResidue>();
//residuesには、missing, inserted と mapped の Residue が入っている。順番はファイルに出てきた順
	private HashMap<Integer,PDBResidue> mapped = new HashMap<Integer,PDBResidue>();//リファレンスと一致する残基が入っている。
	private HashMap<Integer,ArrayList<PDBResidue>> insertion = new HashMap<>();//リファレンスに無い残基が入っている
	String entityId = "";
	PDBChain(String n){
		name = n;
	}
	PDBChain(PDBData p,String n){
		parent = p;
		name = n;
	}
	
	public void setEntityId(String e){
		entityId = e;
	}
	public void changeName(String n){
		name = n;
	}
	public void changeParent(PDBData p){
		parent = p;
	}
	
	public PDBChain copy(){//insertion とかはまだコピーできてない。Fixit
		PDBChain ret = new PDBChain(name);
		ret.desc = this.desc;
		ret.parent = this.parent;
		for(PDBResidue r:residues){
			PDBResidue rr = r.getCopy();
			ret.residues.add(rr);
			rr.parent = ret;
		}
		
		return ret;
	}
	
	public void removeHETATM(){
		Iterator<PDBResidue> ite = residues.iterator();
		while(ite.hasNext()){
			PDBResidue r = ite.next();
			Iterator<PDBAtom> aite = r.atoms.iterator();
			while(aite.hasNext()){
				PDBAtom a = aite.next();
				if(a.isHETATM()){
					aite.remove();
				}
			}
			if(r.atoms.size() == 0){
				ite.remove();
			}
		}
	}
	public void removeWater(){
		Iterator<PDBResidue> ite = residues.iterator();
		while(ite.hasNext()){
			PDBResidue r = ite.next();
			if(r.getName().equals("HOH") || r.getName().equals("DOD")){
				ite.remove();
			}
		}
	}
	public String getName(){
		return name;
	}
	
	public void constructMap(){
		mapped.clear();
		for(PDBResidue r:residues){
			mapped.put(r.getResidueNumber(),r);
		}
	}
	public void addMissingResidue(PDBResidue r){
		missing.add(r);
		residues.add(r);
	}
	public PDBResidue getResidueAt(int i){
		return mapped.get(i);
	}
	public void addSpecialAtom(PDBAtom a,String rescode,int resnum,String altcode,String inscode){
		
		inscode = inscode.replaceAll("[\\s]","");
		altcode = altcode.replaceAll("[\\s]","");
		
		if(inscode != null && inscode.length() > 0){
		
			if(!insertion.containsKey(resnum)){
				ArrayList<PDBResidue> al = new ArrayList<PDBResidue>();
				insertion.put(resnum,al);
			}
			
			ArrayList<PDBResidue> pal = insertion.get(resnum);
			boolean flag = true;
			for(PDBResidue r:pal){
				if(r.getInsertionCode().equals(inscode)){
					r.addAtom(a,altcode);
					flag = false;
					break;
				}
			}
			if(flag){
				PDBResidue r = new PDBResidue(rescode,resnum,altcode,inscode);
				r.setParent(this);
				r.addAtom(a,altcode);
				residues.add(r);
				pal.add(r);
			}
		}else if(altcode != null && altcode.length() > 0){
			
			if(!mapped.containsKey(resnum)){
				PDBResidue r = new PDBResidue(rescode,resnum);
				r.setParent(this);
				mapped.put(resnum,r);
				residues.add(r);
				r.addAtom(a,altcode);
			}else{
				PDBResidue r = mapped.get(resnum);
				r.addAtom(a,altcode);
			}
		}else{
			System.err.println("unknown format.");
		}
		
	}
	
	public String getSequence(){
		return getSequence(true);
	}
	public String getSequence(boolean missing_lowercase){
		Collections.sort(residues, new PDBResidueComparator());
		StringBuffer sb = new StringBuffer();
		for(PDBResidue re:residues){
			String a = re.getName();
			if(PDBResidue.aaMap.get(a) == null){
				/*
				if(re.missing){
					sb.append("x");
				}else{
					sb.append("X");
				}
				*/
				//System.err.println(a+" was not considered as aa.");
			}else{
				if(re.missing && missing_lowercase){
					sb.append(PDBResidue.aaMap.get(a).toLowerCase());
				}else{
					sb.append(PDBResidue.aaMap.get(a));
				}
			}
		}
		return sb.toString();
	}
	
	
	public void addAtom(PDBAtom a,String rescode,int resnum){
		if(!mapped.containsKey(resnum)){
			PDBResidue r = new PDBResidue(rescode,resnum);
			r.setParent(this);
			mapped.put(resnum,r);
			residues.add(r);
		}
		
		mapped.get(resnum).addAtom(a);
	}
	
	
	
	public String getDesc(){
		return desc;
	}
	public void setDesc(String d){
		desc = d;
	}
}
class PDBWater extends PDBAtom{
	
	PDBWater(String pdb_c,String n,Point3D p){
		super(pdb_c,"HOH",n,p);
	}
	
}
class PDBResidue{
	public static final int COMMON_O = 0;
	public static final int COMMON_N = 1;
	public static final int COMMON_C = 2;
	public static final int COMMON_NH = 3;
	public static final int COMMON_CA = 4;
	public static final int COMMON_CB = 5;
	
	public PDBAtom[] commonAtom = new PDBAtom[6];
	private int resNum = 0;
	private String resName ="ALA";
	private String insertionCode = "";
	public static HashMap<String,String> aaMap = new HashMap<String,String>();
	boolean hetatm = false;//変な残基である場合
	boolean missing = false;
	String alternativeCode = "";
	PDBChain parent=null;
	double maxDist = -1000000;//最も遠い原子同士の距離。
	ArrayList<PDBAtom> atoms = new ArrayList<PDBAtom>();
	boolean isLigand = false;
	public static HashSet<String> aaSet = new HashSet<>();
	public static HashSet<String> validAASet = new HashSet<>();
	
	
	PDBResidue(){
		for(int ii = 0;ii < commonAtom.length;ii++){
			commonAtom[ii] = null;
		}
	}
	
	
	static{
		aaMap.put("ALA","A");
		aaMap.put("ARG","R");
		aaMap.put("ASN","N");
		aaMap.put("ASP","D");
		aaMap.put("CYS","C");
		aaMap.put("GLN","Q");
		aaMap.put("GLU","E");
		aaMap.put("GLY","G");
		aaMap.put("HIS","H");
		aaMap.put("ILE","I");
		aaMap.put("LEU","L");
		aaMap.put("LYS","K");
		aaMap.put("MET","M");
		aaMap.put("PHE","F");
		aaMap.put("PRO","P");
		aaMap.put("SER","S");
		aaMap.put("THR","T");
		aaMap.put("TRP","W");
		aaMap.put("TYR","Y");
		aaMap.put("VAL","V");
		validAASet.addAll(aaMap.keySet());
		aaMap.put("SEC","U");
		aaMap.put("UNK","X");
		aaSet.addAll(aaMap.keySet());
		
		aaMap.put("T","T");
		aaMap.put("A","A");
		aaMap.put("G","G");
		aaMap.put("C","C");
		aaMap.put("U","U");
		aaMap.put("N","N");
		
		aaMap.put("DT","T");
		aaMap.put("DA","A");
		aaMap.put("DG","G");
		aaMap.put("DC","C");
		aaMap.put("DU","U");
		aaMap.put("DN","N");
	}
	PDBResidue(String n,int nn){
		resName = n.replaceAll("[\\s]","");
		resNum = nn;
	}
	PDBResidue(String n,int nn,String alt){
		resName = n.replaceAll("[\\s]","");
		resNum = nn;
		alternativeCode = alt;
	}
	PDBResidue(String n,int nn,String alt,String ins){
		resName = n.replaceAll("[\\s]","");
		resNum = nn;
		alternativeCode = alt;
		insertionCode = ins;
	}
	PDBResidue(String n,int nn,String alt,String ins,boolean miss){
		resName = n.replaceAll("[\\s]","");
		resNum = nn;
		alternativeCode = alt;
		insertionCode = ins;
		missing = miss;
	}
	
	public boolean isMissing(){
		return missing;
	}
	public PDBResidue getCopy(){
		PDBResidue ret = new PDBResidue();
		ret.resNum = resNum;
		ret.resName = resName;
		ret.insertionCode = insertionCode;
		ret.hetatm = hetatm;//変な残基である場合
		ret.missing = missing;
		ret.alternativeCode = alternativeCode;
		ret.parent=parent;
		for(PDBAtom a:atoms){
			PDBAtom pp = a.getCopy();
			pp.parent = ret;
			ret.atoms.add(pp);
		}
		return ret;
	}
	
	
	
	public void setHETATM(boolean b){
		hetatm = b;
	}
	public boolean isHETATM(){
		return hetatm;
	}
	public boolean isAA(){
		return aaSet.contains(resName.toUpperCase());
	}
	public boolean isValidAA(){
		return validAASet.contains(resName.toUpperCase());
	}
	public void setName(String xt){
		resName = xt.replaceAll("[\\s]","");
	}
	public void setParent(PDBChain p){
		parent = p;
	}
	public String getName(){
		return resName;
	}
	public String getOneLetterName(){
		if(resName.length() == 1){
			return resName;
		}
		if(aaMap.containsKey(resName)){
			return aaMap.get(resName);
		}
		return "X";
	}
	public void setResidueNumber(int res){
		resNum = res;
	}
	public int getResidueNumber(){
		return resNum;
	}
	public void removeAlt(){
		for(int ii = 0;ii < commonAtom.length;ii++){
			commonAtom[ii] = null;
		}
		Iterator<PDBAtom> aite = atoms.iterator();
		while(aite.hasNext()){
			PDBAtom a = aite.next();
			if(a.getAltCode().length() == 0
			||a.getAltCode().equals("A")
			||a.getAltCode().equals(" ")
			){
			}else{
				aite.remove();
			}
		}
	}
	
	
	public void addAtom(PDBAtom a,String altcode){
		a.setAltCode(altcode);
		addAtom(a);
	}
	public void addAtom(PDBAtom a){
		for(PDBAtom aa:atoms){
			maxDist = Math.max(maxDist, aa.distance(a));
		}
		if(a.isLigand()){
			ligandFlag(true);
		}
		atoms.add(a);
		if(a.isHETATM()){
			setHETATM(true);
		}
		a.setParent(this);
	}
	public boolean isLigand(){
		return isLigand;
	}
	public void ligandFlag(boolean b){
		isLigand = b;
	}
	
	/**
	 * pdb_atom_code が s と等しいAtomを返す
	 * @param s
	 * @return 
	 */
	public PDBAtom getAtomByName(String s){
		PDBAtom candidate = null;
		for(PDBAtom a:atoms){
			if(a.pdb_atom_code.equals(s)){
				if(a.isAlternative()){
					candidate = a;
				}else{
					return a;
				}
			}
			
		}
		return candidate;
	}
	
	
	/**
	 * 挿入の記号が入っている場合に対応するため、String で残基番号の情報を返す
	 */
	public String getPositionCode(){
		return resNum+insertionCode;
	}
	/**
	 * 残基を特定する情報を表すコードを返す。
	 */
	public String getRepresentativeCode(){
		if(alternativeCode.length() == 0){
			return resName+resNum+insertionCode;
		}else{
			return resName+resNum+insertionCode+"_"+alternativeCode;
			
		}
	}
	public void addAtom(String pdb_c,String rcode,String n,Point3D p){
		PDBAtom a = new PDBAtom(pdb_c,rcode,n,p);
		a.setParent(this);
		atoms.add(a);
		if(a.isHETATM()){
			setHETATM(true);
		}
	}
	
	
	public double debugGCode(){
		Point3D ca = null;
		Point3D oc = null;
		Point3D n = null;
		Point3D cb = null;
		
		
		for(PDBAtom a:atoms){
			if(a.pdb_atom_code.equals("C")){
				oc = a.loc;
			}
			if(a.pdb_atom_code.equals("N")){
				n = a.loc;
			}
			if(a.pdb_atom_code.indexOf("CA") == 0){
				ca = a.loc;
			}
			if(a.pdb_atom_code.indexOf("CB") == 0){
				cb = a.loc;
			}
		}
		Point3D p3 = new Point3D();
		p3.x = n.x/2+oc.x/2;
		p3.y = n.y/2+oc.y/2;
		p3.z = n.z/2+oc.z/2;
		
		return p3.distance(cb);
		
	}
	public Point3D calcPseudoCB(){
		Point3D ca = null;
		Point3D oc = null;
		Point3D n = null;
		
		
		for(PDBAtom a:atoms){
			if(a.pdb_atom_code.equals("C")){
				oc = a.loc;
			}
			if(a.pdb_atom_code.equals("N")){
				n = a.loc;
			}
			if(a.pdb_atom_code.indexOf("CA") == 0){
				ca = a.loc;
			}
		}
		
		if(ca == null || oc == null || n == null){
			System.err.println(getRepresentativeCode()+" has some missing atom and cannnot calculate CB.");
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
		
		Point3D ret = Point3D.rotate(pb,vec,-2.186025);
		//Point3D ret = Point3D.rotate(pb,vec,-2.22);
		
		double rlen = Math.sqrt(ret.x*ret.x+ret.y*ret.y+ret.z*ret.z);
		if(rlen == 0){
			System.err.println("CB was not able to be calculated.");
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
	public String getAlternativeCode(){
		return alternativeCode;
		
	}
	public String getInsertionCode(){
		return insertionCode;
		
	}
	public void setAlternativeCode(String s){
		alternativeCode = s;
	}
	
	public void setInsertionCode(String s){
		insertionCode = s;
	}
	/**
	 * CB を返す。複数ある場合は最初のCB。 Glysin の場合は Pseudo CB を返す。
	 */
	public PDBAtom getCB(){
		if(commonAtom[COMMON_CB] != null){
			return commonAtom[COMMON_CB];
		}
		if(!resName.equalsIgnoreCase("GLY")){
			
			for(PDBAtom a:atoms){
				if(a.pdb_atom_code.indexOf("CB") == 0){
					commonAtom[COMMON_CB] = a;
					return a;
				}
			}
		}else{
			Point3D p3 = calcPseudoCB();
			if(p3 == null){
				return null;
			}
			PDBAtom ret =  new PDBAtom("CB","C","-1",p3);
			ret.setParent(this);
			commonAtom[COMMON_CB] = ret;
			return ret;
		}
		
		return null;
	}
	public PDBAtom getCA(){
		if(commonAtom[COMMON_CA] != null){
			return commonAtom[COMMON_CA];
		}
		PDBAtom ret = this.getAtomByName("CA");
		commonAtom[COMMON_CA] = ret;
		return ret;
	}
	
	public PDBAtom getC(){
		if(commonAtom[COMMON_C] != null){
			return commonAtom[COMMON_C];
		}
		PDBAtom ret = this.getAtomByName("C");
		commonAtom[COMMON_C] = ret;
		return ret;
	}
	
	
	public PDBAtom getO(){
		if(commonAtom[COMMON_O] != null){
			return commonAtom[COMMON_O];
		}
		PDBAtom ret = this.getAtomByName("O");
		commonAtom[COMMON_O] = ret;
		return ret;
	}
	
	public PDBAtom getN(){
		if(commonAtom[COMMON_N] != null){
			return commonAtom[COMMON_N];
		}
		PDBAtom ret = this.getAtomByName("N");
		commonAtom[COMMON_N] = ret;
		return ret;
	}
	public PDBAtom getN_H(){
		if(commonAtom[COMMON_NH] != null){
			return commonAtom[COMMON_NH];
		}
		//PDBAtom ret = this.getAtomByName("H");
		PDBAtom ret = calcNH(this);
		commonAtom[COMMON_NH] = ret;
		return ret;
		
	}
	
	public static PDBAtom calcNH(PDBResidue res){
		
		PDBAtom h =  new PDBAtom("H","H","-1",res.getN().loc);
		h.setParent(res);
		if(res.getC() == null 
				|| res.getO() == null
				|| res.getN() == null){
			return null;
		}
		
		Point3D c = res.getC().loc;
		Point3D o = res.getO().loc;
		Point3D n = res.getN().loc;

		double codist = c.distance(o);

		///ちょっと違う気がするが、DSSP に従う
		double xx =  (c.x - o.x)/codist + n.x;
		double yy =  (c.y - o.y)/codist + n.y;
		double zz =  (c.z - o.z)/codist + n.z;
		h.loc.set(xx,yy,zz);
		return h;
	}
	
	public PDBResidue copy_obsolete(){
		PDBResidue ret = new PDBResidue(resName,resNum);
		ret.alternativeCode = alternativeCode;
		for(PDBAtom a:atoms){
			ret.addAtom(a.getCopy());
			a.parent = ret;
		}
		
		return ret;
	}
}



class PDBAtom{
	public static int BASIC_SPHERE_POLYGON = 12;
	public static float BASIC_SPHERE_RADIUS = 0.4f;
	String atom_code ="H";//元素記号
	String pdb_atom_code = "H0";//場所等のデータが入った記号
	
	String place_code_1 = "";
	String place_code_2 = "";
	String atom_number = "";
	
	double bfactor = 0.0;
	double occupancy = 1.0;
	double vdwRadius = 1.8;
	PDBResidue parent = null;
	Point3D loc = new Point3D();
	boolean hetatm  = false;
	boolean isLigand = false;
	private String altcode = "";
	String atomLine ="";
	PDBAtom(){
	}
	
	PDBAtom(String pdb_c,String rcode,String n,Point3D p){
		
		atom_code = rcode;
		pdb_atom_code = pdb_c.replaceAll(" ","");
		atom_number = n.replaceAll(" ","");
		if(pdb_atom_code.length() >= 2){
			char[] cc =  pdb_atom_code.toCharArray();
			if(cc.length > 2){
				place_code_1 = String.valueOf(cc[1]);
				place_code_2 = String.valueOf(cc[2]);
			}
			if(cc.length > 1){
				place_code_1 = String.valueOf(cc[1]);
			}
		}
		
		loc.set(p);
	}
	public void setBFactor(double b){
		bfactor = b;
	}
	public void ligandFlag(boolean b){
		isLigand = b;
	}
	public boolean isLigand(){
		return isLigand;
	}
	public void setOccupancy(double b){
		occupancy = b;
	}
		
	public void setAltCode(String a){
		altcode = a.replaceAll("[\\s]","");
	}
	public void setLine(String s){
		atomLine = s;
	}
	
	public boolean isBackBoneAtom(){
		return pdb_atom_code.equals("CA") || pdb_atom_code.equals("N") || pdb_atom_code.equals("O") || pdb_atom_code.equals("C");
	}
	
	public String getAltCode(){
		return altcode;
	}
	public void setParent(PDBResidue p){
		parent = p;
	}
	public boolean isAlternative(){
		return altcode != null && altcode.length() != 0;
	}
	public PDBAtom getCopy(){
		PDBAtom ret = new PDBAtom();
		ret.loc = new Point3D(loc);
		ret.hetatm = hetatm;
		ret.parent = parent;
		ret.pdb_atom_code = pdb_atom_code;
		ret.atom_number = atom_number;
		ret.atom_code = atom_code;
		ret.place_code_1 = place_code_1;
		ret.place_code_2 = place_code_2;
		
		ret.altcode = altcode;
		ret.bfactor = bfactor;
		return ret;
	}
	public void setHETATM(boolean b){
		hetatm = b;
	}
	public boolean isHETATM(){
		return hetatm;
	}
	public void setPDBATOMCode(String c){
		pdb_atom_code = c;
	}
	public Point3D getLocation(){
		return loc;
	}
	public double distance(Point3D a){
		return loc.distance(a);
	}
	public double distance(PDBAtom a){
		return loc.distance(a.loc);
	}
	
	
	public static String getStringWithLength(double d,int len,int fnum){
		String s = String.format("%."+fnum+"f", d);
		return getStringWithLength(s,len);
	}
	public static String getStringWithLength(int i,int len){
		return getStringWithLength(String.valueOf(i),len);
	}
	public static String getStringWithLength(String str,int len){
		String s = "            "+str;
		return s.substring(s.length()-len,s.length());
	}
	public String makePDBString(){
		StringBuffer sb = new StringBuffer();
		if(hetatm){
			sb.append("HETATM");
		}else{
			sb.append("ATOM  ");
		}
		sb.append(getStringWithLength(this.atom_number,5));
		sb.append(" ");
		sb.append(getStringWithLength(this.pdb_atom_code,4));
		if(altcode.length() > 0){
			sb.append(altcode);
		}else{
			sb.append(" ");
		}
		sb.append(getStringWithLength(this.parent.getName(),3));
		sb.append(" ");
		sb.append(getStringWithLength(this.parent.parent.name,1));
		sb.append(getStringWithLength(this.parent.getResidueNumber(),4));
		String ins = this.parent.getInsertionCode();
		if(ins.length() > 0){
			sb.append(ins);
		}else{
			sb.append(" ");
		}
		sb.append("   ");
		sb.append(getStringWithLength(loc.x,8,3));
		sb.append(getStringWithLength(loc.y,8,3));
		sb.append(getStringWithLength(loc.z,8,3));
		sb.append(getStringWithLength(this.occupancy,6,2));
		sb.append(getStringWithLength(this.bfactor,6,2));
		sb.append("          ");
		sb.append(getStringWithLength(this.atom_code,2));
		sb.append("  ");
		
		return sb.toString();
		
	}
	
}



class Point3D implements Location3D{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	
	public void set(double xx,double yy,double zz){
		x = xx;
		y = yy;
		z = zz;
		
	}
	public void set(Point3D pp){
		x = pp.x;
		y = pp.y;
		z = pp.z;
		
	}
	public void set(double[] dd){
		set(dd[0],dd[1],dd[2]);
	}
	public double distance(Point3D pp){
		return Math.sqrt((x-pp.x)*(x-pp.x)+(y-pp.y)*(y-pp.y)+(z-pp.z)*(z-pp.z));
		
	}
	Point3D(double xx,double yy,double zz){
		set(xx,yy,zz);
	}
	Point3D(){
		
	}
	Point3D(Point3D pp){
		this(pp.x,pp.y,pp.z);
	}
	public String toString(){
	
		return x+","+y+","+z;
	}
	/**
	 * vec を軸として target を radian 度回転する。
	 * 
	 */
	public static Point3D rotate(Point3D target,Point3D vec,double radian){
		double len = vec.x*vec.x+vec.y*vec.y+vec.z*vec.z;
		if(len == 0){
			System.err.println("can not calculate rotated point with zero vector.");
			return target;
		}
		len = Math.sqrt(len);
		double vx = vec.x;
		double vy = vec.y;
		double vz = vec.z;
		if(len != 1){
			vx /= len;
			vy /= len;
			vz /= len;
		}
		Point3D ret = new Point3D();
		double cos = Math.cos(radian);
		double sin = Math.sin(radian);
		
		ret.x = (vx*vx*(1-cos)+cos)*target.x +(vx*vy*(1-cos)-vz*sin)*target.y +(vz*vx*(1-cos)+vy*sin)*target.z;
		ret.y = (vx*vy*(1-cos)+vz*sin)*target.x +(vy*vy*(1-cos)+cos)*target.y +(vz*vy*(1-cos)-vx*sin)*target.z;
		ret.z = (vx*vz*(1-cos)-vy*sin)*target.x +(vz*vy*(1-cos)+vx*sin)*target.y +(vz*vz*(1-cos)+cos)*target.z;
		return ret;
	}
	
	public double getX(){
		return x;
	}
	public double getY(){
		return y;
	}
	public double getZ(){
		return z;
	}
	
	/**
	 * 単位ベクトルにする
	 */
	public void standarize(){
		double len = x*x+y*y+z*z;
		if(len == 0){
		}else{
			len = Math.sqrt(len);
			x /= len;
			y /= len;
			z /= len; 
		}
	}
}


class AtomConnection{
	PDBAtom atomA = null;
	PDBAtom atomB = null;
	AtomConnection(PDBAtom a,PDBAtom b){
		atomA = a;
		atomB = b;
	}
}

/**
 * Biological unit を作成するためのクラス
 Load するときに TransformMatrix を作成しておき、loadtargetChains に Chain を つっこんでおく。
 * 適用する Chain がファイル内で指定されているのでこのような構成。
 getTransformedChains で変換された Chain が得られる。
 普通変換前の TransformMatrix も入っているので元の Chain と合わせると重複になる。
 * @author yamule
 */
class TransformUnit{
	String id = "";
	TransformUnit(String i){
		id = i;
	}
	ArrayList<TransformSet> transformSet = new ArrayList<>();
	public void addTransformSet(TransformSet t){
		transformSet.add(t);
	}
	public ArrayList<PDBChain> getTransformedChains(){
		ArrayList<PDBChain> ret = new ArrayList<>();
		for(TransformSet t:transformSet){
			ret.addAll(t.getTransformedChains());
		}
		return ret;
	}
	
	public ArrayList<PDBChain> getTransformedChains(Collection<PDBChain> chains){
		ArrayList<PDBChain> ret = new ArrayList<>();
		for(TransformSet t:transformSet){
			ret.addAll(t.getTransformedChains());
		}
		return ret;
	}
}

class TransformSet{
	ArrayList<PDBChain> targetChains = new ArrayList<>();
	ArrayList<TransformMatrix> units = new ArrayList<>();
	public void addChain(PDBChain c){
		targetChains.add(c);
	}
	public void addUnit(TransformMatrix u){
		units.add(u);
	}
	
	public ArrayList<PDBChain> getTransformedChains(){
		return getTransformedChains(targetChains);
	}
	public ArrayList<PDBChain> getTransformedChains(Collection<PDBChain> tchains){
		ArrayList<PDBChain> ret = new ArrayList<>();
		for(TransformMatrix t:units){
			for(PDBChain c:tchains){
				PDBChain cc = c.copy();
				for(PDBResidue r:cc.residues){
					for(PDBAtom a:r.atoms){
						Point3D p3 = t.transform(a.loc);
						a.loc.set(p3);
					}
				}
				ret.add(cc);
			}
		}
		return ret;
	}
}

class TransformMatrix{
	double mat[][] = new double[3][4];
	
	public void set(int r,int c,double d){
		mat[r][c] = d;
	}
	
	
	
	public Point3D transform(Point3D v){
		Point3D ret = new Point3D(0,0,0);
		double res[] = new double[3];
		double pt[] = new double[3];
		pt[0] = v.x;
		pt[1] = v.y;
		pt[2] = v.z;
		
		for(int ii = 0;ii < 3;ii++){
			res[ii] = pt[0]*mat[ii][0]
					+pt[1]*mat[ii][1]
					+pt[2]*mat[ii][2]
					+mat[ii][3];
		}
		ret.x = res[0];
		ret.y = res[1];
		ret.z = res[2];
		return ret;
	}
}
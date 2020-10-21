/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
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
public class DefaultAtomCoordinate {
	//20171116 に 1.5Å以下？未満？ Monomer でフィルタリングして最初の単一チェーンのエントリ 5omt からとった
	//PDBData loadPDB からコピペしてきたやつ。面倒だったので
	public static PDBData asPDBData(){
		PDBData ret = new PDBData();
		try{
			Pattern mispat = Pattern.compile("REMARK 465[\\s]+M[\\s]+RES[\\s]+C[\\s]+SSSEQI");
			ArrayList<String> transline = new ArrayList<>();
			ArrayList<String> scaleline = new ArrayList<>();
			ArrayList<String> crystline = new ArrayList<>();
			
			Pattern mispat2 = Pattern.compile("REMARK 465[\\s]+([A-Z]+)[\\s]+([A-Za-z0-9]+)[\\s]+([\\-0-9A-Z]+)");
			Pattern inspat = Pattern.compile("([\\-0-9]+)([A-Z]+)");
			String tit = "";
			HashMap<Integer,PDBAtom> atommap = new HashMap<>();
			ArrayList<ArrayList<Integer>> co = new ArrayList<>();
			boolean connectionerrorflag = false;
			for(int ll = 0;ll < atomlines.length;ll++){
				String line = atomlines[ll];
				if(line.length() < 6){
					continue;
				}
				
				String linelabel = line.substring(0,6);
				
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
					newp.setLine(line.replaceAll("[\\r\\n]","")+"\n");
					if(temperature_factor.length() > 0){
						newp.setBFactor(Double.parseDouble(temperature_factor.replaceAll("[\\s]","")));
					}
					if(occupancy.length() > 0){
						newp.setOccupancy(Double.parseDouble(occupancy.replaceAll("[\\s]","")));
					}
					int acode = Integer.parseInt(atom_serial_num.replaceAll("[\\s]",""));
					if(atommap.containsKey(acode)){
						System.err.println("Atom code duplication\n"+line);
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
				Pattern startpat = Pattern.compile("REMARK 350 BIOMOLECULE[\\s]*\\:[\\s]([^\\s]+)");
				TransformUnit currentunit = null;
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
							
							if(sst.indexOf("CHAIN") > -1 || sst.indexOf("BIOMOLECULE") > -1){
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
			
		}catch(Exception e){
			e.printStackTrace();
		}
		Iterator<String> ite = ret.chains.keySet().iterator();
		while(ite.hasNext()){
			PDBChain c = ret.chains.get(ite.next());
			//System.out.println(c.name);
			Collections.sort(c.residues, new PDBResidueComparator());
		}
		
		
		return ret;
	}
	
	static String[] atomlines ={ 
		
"ATOM      6  N   ARG A  35      12.506  36.565   7.716  1.00 18.00           N  ",
"ATOM      7  CA  ARG A  35      12.654  35.467   8.680  1.00 16.94           C  ",
"ATOM      8  C   ARG A  35      12.701  34.068   8.058  1.00 15.64           C  ",
"ATOM      9  O   ARG A  35      13.272  33.117   8.613  1.00 16.57           O  ",
"ATOM     10  CB  ARG A  35      11.497  35.573   9.689  1.00 22.02           C  ",
"ATOM     11  CG  ARG A  35      11.719  34.920  10.993  1.00 21.39           C  ",
"ATOM     12  CD  ARG A  35      10.696  35.325  12.058  1.00 20.47           C  ",
"ATOM     13  NE  ARG A  35       9.381  34.700  11.909  1.00 18.84           N  ",
"ATOM     14  CZ  ARG A  35       9.044  33.472  12.332  1.00 19.49           C  ",
"ATOM     15  NH1 ARG A  35       9.925  32.690  12.910  1.00 20.87           N  ",
"ATOM     16  NH2 ARG A  35       7.804  33.038  12.140  1.00 24.42           N  ",
"ATOM     17  N   TYR A  36      12.071  33.937   6.914  1.00 12.70           N  ",
"ATOM     18  CA  TYR A  36      11.982  32.660   6.183  1.00 12.66           C  ",
"ATOM     19  C   TYR A  36      11.837  32.934   4.715  1.00 13.35           C  ",
"ATOM     20  O   TYR A  36      11.511  34.064   4.303  1.00 15.85           O  ",
"ATOM     21  CB  TYR A  36      10.830  31.796   6.721  1.00 14.81           C  ",
"ATOM     22  CG  TYR A  36       9.509  32.488   6.731  1.00 14.38           C  ",
"ATOM     23  CD1 TYR A  36       8.685  32.490   5.630  1.00 17.53           C  ",
"ATOM     24  CD2 TYR A  36       9.081  33.149   7.878  1.00 14.85           C  ",
"ATOM     25  CE1 TYR A  36       7.463  33.160   5.661  1.00 19.34           C  ",
"ATOM     26  CE2 TYR A  36       7.896  33.820   7.917  1.00 18.93           C  ",
"ATOM     27  CZ  TYR A  36       7.094  33.845   6.803  1.00 21.32           C  ",
"ATOM     28  OH  TYR A  36       5.866  34.508   6.901  1.00 27.25           O  ",
"ATOM     29  N   ASP A  37      12.113  31.908   3.927  1.00 11.65           N  ",
"ATOM     30  CA  ASP A  37      12.075  31.991   2.490  1.00 12.55           C  ",
"ATOM     31  C   ASP A  37      10.817  31.403   1.857  1.00 13.91           C  ",
"ATOM     32  O   ASP A  37      10.380  31.858   0.801  1.00 17.76           O  ",
"ATOM     33  CB  ASP A  37      13.320  31.285   1.905  1.00 12.52           C  ",
"ATOM     34  CG  ASP A  37      14.614  31.912   2.379  1.00 13.26           C  ",
"ATOM     35  OD1 ASP A  37      14.776  33.141   2.167  1.00 15.66           O  ",
"ATOM     36  OD2 ASP A  37      15.460  31.192   2.989  1.00 14.12           O  ",
"ATOM     45  N   VAL A  39       7.064  28.776   2.988  1.00 12.94           N  ",
"ATOM     46  CA  VAL A  39       6.205  28.170   4.001  1.00 13.18           C  ",
"ATOM     47  C   VAL A  39       5.850  26.760   3.579  1.00 12.90           C  ",
"ATOM     48  O   VAL A  39       5.751  26.465   2.388  1.00 17.73           O  ",
"ATOM     49  CB  VAL A  39       4.942  29.010   4.256  1.00 17.00           C  ",
"ATOM     50  CG1 VAL A  39       5.284  30.440   4.537  1.00 18.59           C  ",
"ATOM     51  CG2 VAL A  39       3.899  28.876   3.171  1.00 20.54           C  ",
"ATOM     52  N   LEU A  40       5.644  25.906   4.583  1.00 11.67           N  ",
"ATOM     53  CA  LEU A  40       5.170  24.533   4.404  1.00 12.74           C  ",
"ATOM     54  C   LEU A  40       4.005  24.329   5.383  1.00 10.91           C  ",
"ATOM     55  O   LEU A  40       4.148  24.504   6.585  1.00 12.22           O  ",
"ATOM     56  CB  LEU A  40       6.299  23.578   4.709  1.00 15.21           C  ",
"ATOM     57  CG  LEU A  40       6.042  22.094   4.494  1.00 17.02           C  ",
"ATOM     58  CD1 LEU A  40       7.328  21.386   4.125  1.00 19.92           C  ",
"ATOM     59  CD2 LEU A  40       5.406  21.403   5.700  1.00 19.32           C  ",
"ATOM     72  N   PHE A  42       2.081  21.697   7.179  1.00 12.09           N  ",
"ATOM     73  CA  PHE A  42       2.342  20.304   7.482  1.00 11.25           C  ",
"ATOM     74  C   PHE A  42       1.065  19.456   7.281  1.00 12.19           C  ",
"ATOM     75  O   PHE A  42      -0.004  19.879   7.679  1.00 12.37           O  ",
"ATOM     76  CB  PHE A  42       2.828  20.171   8.916  1.00 11.73           C  ",
"ATOM     77  CG  PHE A  42       3.246  18.755   9.251  1.00 11.75           C  ",
"ATOM     78  CD1 PHE A  42       4.510  18.308   8.930  1.00 13.02           C  ",
"ATOM     79  CD2 PHE A  42       2.341  17.865   9.767  1.00 11.29           C  ",
"ATOM     80  CE1 PHE A  42       4.870  16.994   9.181  1.00 14.43           C  ",
"ATOM     81  CE2 PHE A  42       2.697  16.555   9.991  1.00 12.96           C  ",
"ATOM     82  CZ  PHE A  42       3.971  16.138   9.700  1.00 14.58           C  ",
"ATOM     83  N   PRO A  43       1.175  18.237   6.699  1.00 11.92           N  ",
"ATOM     84  CA  PRO A  43       0.035  17.374   6.405  1.00 12.93           C  ",
"ATOM     85  C   PRO A  43      -0.477  16.654   7.633  1.00 12.74           C  ",
"ATOM     86  O   PRO A  43      -0.318  15.424   7.791  1.00 12.88           O  ",
"ATOM     87  CB  PRO A  43       0.595  16.419   5.336  1.00 14.05           C  ",
"ATOM     88  CG  PRO A  43       2.033  16.292   5.735  1.00 14.13           C  ",
"ATOM     89  CD  PRO A  43       2.404  17.685   6.116  1.00 13.89           C  ",
"ATOM     90  N   ALA A  44      -1.058  17.428   8.548  1.00 12.12           N  ",
"ATOM     91  CA  ALA A  44      -1.484  16.884   9.813  1.00 14.63           C  ",
"ATOM     92  C   ALA A  44      -2.581  15.833   9.675  1.00 15.11           C  ",
"ATOM     93  O   ALA A  44      -2.648  14.919  10.489  1.00 18.33           O  ",
"ATOM     94  CB  ALA A  44      -1.945  18.015  10.712  1.00 18.87           C  ",
"ATOM     95  N   SER A  45      -3.401  15.902   8.632  1.00 15.52           N  ",
"ATOM     96  CA  SER A  45      -4.445  14.866   8.500  1.00 18.43           C  ",
"ATOM     97  C   SER A  45      -3.841  13.498   8.191  1.00 18.39           C  ",
"ATOM     98  O   SER A  45      -4.376  12.456   8.589  1.00 23.63           O  ",
"ATOM     99  CB  SER A  45      -5.476  15.237   7.423  1.00 20.57           C  ",
"ATOM    100  OG  SER A  45      -4.934  15.192   6.083  1.00 27.19           O  ",
"ATOM    131  N   GLU A  49      -0.012  12.281  14.444  1.00 12.58           N  ",
"ATOM    132  CA  GLU A  49       0.803  12.238  15.665  1.00 12.03           C  ",
"ATOM    133  C   GLU A  49       2.056  13.078  15.507  1.00 11.40           C  ",
"ATOM    134  O   GLU A  49       2.442  13.818  16.436  1.00 12.12           O  ",
"ATOM    135  CB  GLU A  49       1.196  10.806  16.029  1.00 14.74           C  ",
"ATOM    136  CG  GLU A  49       0.115   9.950  16.550  1.00 17.40           C  ",
"ATOM    137  CD  GLU A  49       0.706   8.607  17.015  1.00 20.70           C  ",
"ATOM    138  OE1 GLU A  49       1.732   8.490  17.785  1.00 20.18           O  ",
"ATOM    139  OE2 GLU A  49       0.149   7.640  16.556  1.00 25.98           O  ",
"ATOM    140  N   THR A  50       2.709  12.987  14.352  1.00 10.98           N  ",
"ATOM    141  CA  THR A  50       3.863  13.862  14.078  1.00 10.30           C  ",
"ATOM    142  C   THR A  50       3.436  15.324  14.094  1.00  9.95           C  ",
"ATOM    143  O   THR A  50       4.073  16.177  14.693  1.00 10.96           O  ",
"ATOM    144  CB  THR A  50       4.548  13.441  12.792  1.00 10.51           C  ",
"ATOM    145  OG1 THR A  50       4.958  12.066  12.896  1.00 11.72           O  ",
"ATOM    146  CG2 THR A  50       5.755  14.255  12.518  1.00 11.56           C  ",
"ATOM    147  N   GLY A  51       2.310  15.611  13.431  1.00 10.02           N  ",
"ATOM    148  CA  GLY A  51       1.837  16.973  13.396  1.00  9.90           C  ",
"ATOM    149  C   GLY A  51       1.524  17.539  14.752  1.00  9.93           C  ",
"ATOM    150  O   GLY A  51       1.802  18.708  15.023  1.00 10.97           O  ",
"ATOM    156  N   HIS A  53       2.903  16.721  17.566  1.00 10.29           N  ",
"ATOM    157  CA  HIS A  53       4.177  17.033  18.189  1.00  9.80           C  ",
"ATOM    158  C   HIS A  53       4.722  18.368  17.714  1.00 10.19           C  ",
"ATOM    159  O   HIS A  53       5.118  19.215  18.538  1.00 10.33           O  ",
"ATOM    160  CB  HIS A  53       5.176  15.915  17.990  1.00 10.64           C  ",
"ATOM    161  CG  HIS A  53       6.469  16.165  18.686  1.00 10.42           C  ",
"ATOM    162  ND1 HIS A  53       7.681  15.705  18.204  1.00 11.81           N  ",
"ATOM    163  CD2 HIS A  53       6.742  16.842  19.832  1.00 11.10           C  ",
"ATOM    164  CE1 HIS A  53       8.642  16.113  19.022  1.00 11.59           C  ",
"ATOM    165  NE2 HIS A  53       8.096  16.797  20.011  1.00 11.06           N  ",
"ATOM    166  N   ILE A  54       4.769  18.555  16.393  1.00  9.71           N  ",
"ATOM    167  CA  ILE A  54       5.320  19.790  15.859  1.00  9.20           C  ",
"ATOM    168  C   ILE A  54       4.518  20.989  16.404  1.00  9.88           C  ",
"ATOM    169  O   ILE A  54       5.096  21.987  16.840  1.00 10.12           O  ",
"ATOM    170  CB  ILE A  54       5.320  19.775  14.331  1.00  9.72           C  ",
"ATOM    171  CG1 ILE A  54       6.195  18.630  13.821  1.00 11.25           C  ",
"ATOM    172  CG2 ILE A  54       5.745  21.116  13.832  1.00 11.39           C  ",
"ATOM    173  CD1 ILE A  54       6.080  18.414  12.336  1.00 13.38           C  ",
"ATOM    201  N   LYS A  59       4.538  25.672  20.117  1.00 10.95           N  ",
"ATOM    202  CA  LYS A  59       3.902  26.422  21.224  1.00 11.52           C  ",
"ATOM    203  C   LYS A  59       4.717  26.379  22.540  1.00 11.95           C  ",
"ATOM    204  O   LYS A  59       4.455  27.166  23.416  1.00 17.26           O  ",
"ATOM    205  CB  LYS A  59       2.476  25.925  21.442  1.00 13.42           C  ",
"ATOM    206  CG  LYS A  59       1.540  26.275  20.292  1.00 16.54           C  ",
"ATOM    207  CD  LYS A  59       0.136  25.849  20.670  1.00 25.52           C  ",
"ATOM    208  CE  LYS A  59      -0.884  26.259  19.640  1.00 35.12           C  ",
"ATOM    209  NZ  LYS A  59      -0.896  25.251  18.579  1.00 39.77           N  ",
"ATOM    249  N   CYS A  66      12.759  25.969   9.802  1.00 11.30           N  ",
"ATOM    250  CA  CYS A  66      13.823  25.027   9.566  1.00 11.27           C  ",
"ATOM    251  C   CYS A  66      14.787  25.666   8.556  1.00 10.91           C  ",
"ATOM    252  O   CYS A  66      14.408  25.837   7.391  1.00 11.68           O  ",
"ATOM    253  CB  CYS A  66      13.284  23.728   8.966  1.00 13.28           C  ",
"ATOM    254  SG  CYS A  66      14.616  22.561   8.414  1.00 16.12           S  ",
"ATOM    344  N   GLN A  78      22.242  12.435  17.465  1.00 28.13           N  ",
"ATOM    345  CA  GLN A  78      23.023  11.925  18.636  1.00 33.97           C  ",
"ATOM    346  C   GLN A  78      22.686  12.800  19.822  1.00 33.47           C  ",
"ATOM    347  O   GLN A  78      22.396  12.314  20.901  1.00 42.10           O  ",
"ATOM    348  CB  GLN A  78      24.557  12.020  18.451  1.00 39.30           C  ",
"ATOM    349  CG  GLN A  78      25.307  10.885  17.763  1.00 49.80           C  ",
"ATOM    350  CD  GLN A  78      25.063   9.518  18.374  1.00 50.32           C  ",
"ATOM    351  OE1 GLN A  78      24.083   8.837  18.049  1.00 62.39           O  ",
"ATOM    352  NE2 GLN A  78      25.966   9.098  19.249  1.00 58.57           N  ",
"ATOM    478  N   TRP A  95      11.629  17.644  17.143  1.00  9.84           N  ",
"ATOM    479  CA  TRP A  95      11.333  18.957  16.588  1.00  9.54           C  ",
"ATOM    480  C   TRP A  95      12.117  19.968  17.428  1.00  9.80           C  ",
"ATOM    481  O   TRP A  95      11.882  20.028  18.647  1.00 10.80           O  ",
"ATOM    482  CB  TRP A  95       9.857  19.299  16.555  1.00  9.22           C  ",
"ATOM    483  CG  TRP A  95       9.661  20.600  15.879  1.00  9.81           C  ",
"ATOM    484  CD1 TRP A  95       9.702  21.861  16.448  1.00 10.50           C  ",
"ATOM    485  CD2 TRP A  95       9.546  20.792  14.476  1.00  9.53           C  ",
"ATOM    486  NE1 TRP A  95       9.585  22.813  15.461  1.00 10.90           N  ",
"ATOM    487  CE2 TRP A  95       9.526  22.187  14.239  1.00  9.99           C  ",
"ATOM    488  CE3 TRP A  95       9.511  19.910  13.377  1.00 10.03           C  ",
"ATOM    489  CZ2 TRP A  95       9.403  22.711  12.962  1.00 11.57           C  ",
"ATOM    490  CZ3 TRP A  95       9.409  20.444  12.113  1.00 12.14           C  ",
"ATOM    491  CH2 TRP A  95       9.324  21.828  11.920  1.00 13.19           C  ",
"ATOM    499  N   MET A  97      13.740  19.413  13.598  1.00  9.08           N  ",
"ATOM    500  CA  MET A  97      14.249  18.213  12.948  1.00  9.36           C  ",
"ATOM    501  C   MET A  97      15.768  18.220  12.860  1.00  9.59           C  ",
"ATOM    502  O   MET A  97      16.431  19.237  12.795  1.00 10.61           O  ",
"ATOM    503  CB  MET A  97      13.706  18.115  11.511  1.00 10.06           C  ",
"ATOM    504  CG  MET A  97      12.191  18.017  11.480  1.00 11.12           C  ",
"ATOM    505  SD  MET A  97      11.524  17.886   9.820  1.00 13.65           S  ",
"ATOM    506  CE  MET A  97      12.018  19.454   9.127  1.00 15.04           C  ",
"ATOM    645  N   ASN A 117       5.332   8.338  11.981  1.00 13.93           N  ",
"ATOM    646  CA  ASN A 117       6.576   9.054  11.916  1.00 13.56           C  ",
"ATOM    647  C   ASN A 117       7.392   8.775  10.657  1.00 12.94           C  ",
"ATOM    648  O   ASN A 117       8.003   9.691  10.077  1.00 13.77           O  ",
"ATOM    649  CB  ASN A 117       7.414   8.815  13.144  1.00 14.47           C  ",
"ATOM    650  CG  ASN A 117       8.365   9.951  13.408  1.00 15.92           C  ",
"ATOM    651  OD1 ASN A 117       9.536   9.720  13.557  1.00 17.13           O  ",
"ATOM    652  ND2 ASN A 117       7.860  11.188  13.460  1.00 15.47           N  ",	
	};
	static String [] sourcelines = {
"ATOM      6  N   ARG A  35      12.506  36.565   7.716  1.00 18.00           N  ",
"ATOM      7  CA  ARG A  35      12.654  35.467   8.680  1.00 16.94           C  ",
"ATOM      8  C   ARG A  35      12.701  34.068   8.058  1.00 15.64           C  ",
"ATOM      9  O   ARG A  35      13.272  33.117   8.613  1.00 16.57           O  ",
"ATOM     10  CB  ARG A  35      11.497  35.573   9.689  1.00 22.02           C  ",
"ATOM     11  CG  ARG A  35      11.719  34.920  10.993  1.00 21.39           C  ",
"ATOM     12  CD  ARG A  35      10.696  35.325  12.058  1.00 20.47           C  ",
"ATOM     13  NE  ARG A  35       9.381  34.700  11.909  1.00 18.84           N  ",
"ATOM     14  CZ  ARG A  35       9.044  33.472  12.332  1.00 19.49           C  ",
"ATOM     15  NH1 ARG A  35       9.925  32.690  12.910  1.00 20.87           N  ",
"ATOM     16  NH2 ARG A  35       7.804  33.038  12.140  1.00 24.42           N  ",
"ATOM     17  N   TYR A  36      12.071  33.937   6.914  1.00 12.70           N  ",
"ATOM     18  CA  TYR A  36      11.982  32.660   6.183  1.00 12.66           C  ",
"ATOM     19  C   TYR A  36      11.837  32.934   4.715  1.00 13.35           C  ",
"ATOM     20  O   TYR A  36      11.511  34.064   4.303  1.00 15.85           O  ",
"ATOM     21  CB  TYR A  36      10.830  31.796   6.721  1.00 14.81           C  ",
"ATOM     22  CG  TYR A  36       9.509  32.488   6.731  1.00 14.38           C  ",
"ATOM     23  CD1 TYR A  36       8.685  32.490   5.630  1.00 17.53           C  ",
"ATOM     24  CD2 TYR A  36       9.081  33.149   7.878  1.00 14.85           C  ",
"ATOM     25  CE1 TYR A  36       7.463  33.160   5.661  1.00 19.34           C  ",
"ATOM     26  CE2 TYR A  36       7.896  33.820   7.917  1.00 18.93           C  ",
"ATOM     27  CZ  TYR A  36       7.094  33.845   6.803  1.00 21.32           C  ",
"ATOM     28  OH  TYR A  36       5.866  34.508   6.901  1.00 27.25           O  ",
"ATOM     29  N   ASP A  37      12.113  31.908   3.927  1.00 11.65           N  ",
"ATOM     30  CA  ASP A  37      12.075  31.991   2.490  1.00 12.55           C  ",
"ATOM     31  C   ASP A  37      10.817  31.403   1.857  1.00 13.91           C  ",
"ATOM     32  O   ASP A  37      10.380  31.858   0.801  1.00 17.76           O  ",
"ATOM     33  CB  ASP A  37      13.320  31.285   1.905  1.00 12.52           C  ",
"ATOM     34  CG  ASP A  37      14.614  31.912   2.379  1.00 13.26           C  ",
"ATOM     35  OD1 ASP A  37      14.776  33.141   2.167  1.00 15.66           O  ",
"ATOM     36  OD2 ASP A  37      15.460  31.192   2.989  1.00 14.12           O  ",
"ATOM     37  N   ASP A  38      10.249  30.389   2.481  1.00 13.41           N  ",
"ATOM     38  CA  ASP A  38       9.029  29.737   2.038  1.00 13.80           C  ",
"ATOM     39  C   ASP A  38       8.311  29.177   3.223  1.00 11.72           C  ",
"ATOM     40  O   ASP A  38       8.865  29.074   4.318  1.00 12.45           O  ",
"ATOM     41  CB  ASP A  38       9.235  28.708   0.912  1.00 17.86           C  ",
"ATOM     42  CG  ASP A  38       7.978  28.564  -0.064  1.00 22.48           C  ",
"ATOM     43  OD1 ASP A  38       6.839  29.147   0.148  1.00 22.75           O  ",
"ATOM     44  OD2 ASP A  38       8.121  27.808  -1.077  1.00 26.00           O  ",
"ATOM     45  N   VAL A  39       7.064  28.776   2.988  1.00 12.94           N  ",
"ATOM     46  CA  VAL A  39       6.205  28.170   4.001  1.00 13.18           C  ",
"ATOM     47  C   VAL A  39       5.850  26.760   3.579  1.00 12.90           C  ",
"ATOM     48  O   VAL A  39       5.751  26.465   2.388  1.00 17.73           O  ",
"ATOM     49  CB  VAL A  39       4.942  29.010   4.256  1.00 17.00           C  ",
"ATOM     50  CG1 VAL A  39       5.284  30.440   4.537  1.00 18.59           C  ",
"ATOM     51  CG2 VAL A  39       3.899  28.876   3.171  1.00 20.54           C  ",
"ATOM     52  N   LEU A  40       5.644  25.906   4.583  1.00 11.67           N  ",
"ATOM     53  CA  LEU A  40       5.170  24.533   4.404  1.00 12.74           C  ",
"ATOM     54  C   LEU A  40       4.005  24.329   5.383  1.00 10.91           C  ",
"ATOM     55  O   LEU A  40       4.148  24.504   6.585  1.00 12.22           O  ",
"ATOM     56  CB  LEU A  40       6.299  23.578   4.709  1.00 15.21           C  ",
"ATOM     57  CG  LEU A  40       6.042  22.094   4.494  1.00 17.02           C  ",
"ATOM     58  CD1 LEU A  40       7.328  21.386   4.125  1.00 19.92           C  ",
"ATOM     59  CD2 LEU A  40       5.406  21.403   5.700  1.00 19.32           C  ",
"ATOM     60  N   TYR A  41       2.854  23.930   4.852  1.00 11.26           N  ",
"ATOM     61  CA  TYR A  41       1.697  23.554   5.682  1.00 11.67           C  ",
"ATOM     62  C   TYR A  41       1.803  22.046   5.928  1.00 12.04           C  ",
"ATOM     63  O   TYR A  41       1.675  21.213   5.022  1.00 13.27           O  ",
"ATOM     64  CB  TYR A  41       0.412  23.917   4.963  1.00 13.18           C  ",
"ATOM     65  CG  TYR A  41       0.282  25.391   4.678  1.00 12.81           C  ",
"ATOM     66  CD1 TYR A  41       0.026  26.279   5.709  1.00 15.19           C  ",
"ATOM     67  CD2 TYR A  41       0.412  25.894   3.404  1.00 14.60           C  ",
"ATOM     68  CE1 TYR A  41      -0.075  27.641   5.488  1.00 16.92           C  ",
"ATOM     69  CE2 TYR A  41       0.330  27.246   3.205  1.00 16.10           C  ",
"ATOM     70  CZ  TYR A  41       0.082  28.108   4.249  1.00 16.84           C  ",
"ATOM     71  OH  TYR A  41      -0.028  29.475   3.995  1.00 21.92           O  ",
"ATOM     72  N   PHE A  42       2.081  21.697   7.179  1.00 12.09           N  ",
"ATOM     73  CA  PHE A  42       2.342  20.304   7.482  1.00 11.25           C  ",
"ATOM     74  C   PHE A  42       1.065  19.456   7.281  1.00 12.19           C  ",
"ATOM     75  O   PHE A  42      -0.004  19.879   7.679  1.00 12.37           O  ",
"ATOM     76  CB  PHE A  42       2.828  20.171   8.916  1.00 11.73           C  ",
"ATOM     77  CG  PHE A  42       3.246  18.755   9.251  1.00 11.75           C  ",
"ATOM     78  CD1 PHE A  42       4.510  18.308   8.930  1.00 13.02           C  ",
"ATOM     79  CD2 PHE A  42       2.341  17.865   9.767  1.00 11.29           C  ",
"ATOM     80  CE1 PHE A  42       4.870  16.994   9.181  1.00 14.43           C  ",
"ATOM     81  CE2 PHE A  42       2.697  16.555   9.991  1.00 12.96           C  ",
"ATOM     82  CZ  PHE A  42       3.971  16.138   9.700  1.00 14.58           C  ",
"ATOM     83  N   PRO A  43       1.175  18.237   6.699  1.00 11.92           N  ",
"ATOM     84  CA  PRO A  43       0.035  17.374   6.405  1.00 12.93           C  ",
"ATOM     85  C   PRO A  43      -0.477  16.654   7.633  1.00 12.74           C  ",
"ATOM     86  O   PRO A  43      -0.318  15.424   7.791  1.00 12.88           O  ",
"ATOM     87  CB  PRO A  43       0.595  16.419   5.336  1.00 14.05           C  ",
"ATOM     88  CG  PRO A  43       2.033  16.292   5.735  1.00 14.13           C  ",
"ATOM     89  CD  PRO A  43       2.404  17.685   6.116  1.00 13.89           C  ",
"ATOM     90  N   ALA A  44      -1.058  17.428   8.548  1.00 12.12           N  ",
"ATOM     91  CA  ALA A  44      -1.484  16.884   9.813  1.00 14.63           C  ",
"ATOM     92  C   ALA A  44      -2.581  15.833   9.675  1.00 15.11           C  ",
"ATOM     93  O   ALA A  44      -2.648  14.919  10.489  1.00 18.33           O  ",
"ATOM     94  CB  ALA A  44      -1.945  18.015  10.712  1.00 18.87           C  ",
"ATOM     95  N   SER A  45      -3.401  15.902   8.632  1.00 15.52           N  ",
"ATOM     96  CA  SER A  45      -4.445  14.866   8.500  1.00 18.43           C  ",
"ATOM     97  C   SER A  45      -3.841  13.498   8.191  1.00 18.39           C  ",
"ATOM     98  O   SER A  45      -4.376  12.456   8.589  1.00 23.63           O  ",
"ATOM     99  CB  SER A  45      -5.476  15.237   7.423  1.00 20.57           C  ",
"ATOM    100  OG  SER A  45      -4.934  15.192   6.083  1.00 27.19           O  ",
"ATOM    101  N   ARG A  46      -2.744  13.504   7.448  1.00 16.63           N  ",
"ATOM    102  CA  ARG A  46      -2.094  12.269   7.022  1.00 16.41           C  ",
"ATOM    103  C   ARG A  46      -1.213  11.696   8.159  1.00 17.24           C  ",
"ATOM    104  O   ARG A  46      -1.156  10.470   8.345  1.00 18.21           O  ",
"ATOM    105  CB  ARG A  46      -1.219  12.477   5.777  1.00 23.50           C  ",
"ATOM    106  CG  ARG A  46      -1.929  12.619   4.434  1.00 28.14           C  ",
"ATOM    107  CD  ARG A  46      -0.971  13.013   3.313  1.00 38.97           C  ",
"ATOM    108  NE  ARG A  46       0.052  12.007   2.967  1.00 49.36           N  ",
"ATOM    109  CZ  ARG A  46      -0.164  10.902   2.246  1.00 51.28           C  ",
"ATOM    110  NH1 ARG A  46      -1.381  10.591   1.798  1.00 62.38           N  ",
"ATOM    111  NH2 ARG A  46       0.840  10.085   1.964  1.00 50.85           N  ",
"ATOM    112  N   TYR A  47      -0.545  12.578   8.912  1.00 12.74           N  ",
"ATOM    113  CA  TYR A  47       0.444  12.162   9.919  1.00 12.44           C  ",
"ATOM    114  C   TYR A  47       0.107  12.923  11.186  1.00 11.86           C  ",
"ATOM    115  O   TYR A  47       0.852  13.810  11.616  1.00 11.57           O  ",
"ATOM    116  CB  TYR A  47       1.880  12.407   9.459  1.00 13.30           C  ",
"ATOM    117  CG  TYR A  47       2.138  11.885   8.087  1.00 14.23           C  ",
"ATOM    118  CD1 TYR A  47       2.112  10.535   7.841  1.00 16.56           C  ",
"ATOM    119  CD2 TYR A  47       2.401  12.743   7.032  1.00 15.43           C  ",
"ATOM    120  CE1 TYR A  47       2.257  10.081   6.552  1.00 19.33           C  ",
"ATOM    121  CE2 TYR A  47       2.550  12.286   5.747  1.00 20.35           C  ",
"ATOM    122  CZ  TYR A  47       2.512  10.955   5.504  1.00 22.07           C  ",
"ATOM    123  OH  TYR A  47       2.689  10.482   4.187  1.00 30.14           O  ",
"ATOM    124  N   PRO A  48      -1.035  12.594  11.806  1.00 11.98           N  ",
"ATOM    125  CA  PRO A  48      -1.496  13.404  12.900  1.00 13.00           C  ",
"ATOM    126  C   PRO A  48      -0.659  13.389  14.143  1.00 11.52           C  ",
"ATOM    127  O   PRO A  48      -0.632  14.389  14.868  1.00 13.19           O  ",
"ATOM    128  CB  PRO A  48      -2.908  12.851  13.184  1.00 15.18           C  ",
"ATOM    129  CG  PRO A  48      -2.880  11.477  12.676  1.00 15.97           C  ",
"ATOM    130  CD  PRO A  48      -2.006  11.539  11.460  1.00 14.77           C  ",
"ATOM    131  N   GLU A  49      -0.012  12.281  14.444  1.00 12.58           N  ",
"ATOM    132  CA  GLU A  49       0.803  12.238  15.665  1.00 12.03           C  ",
"ATOM    133  C   GLU A  49       2.056  13.078  15.507  1.00 11.40           C  ",
"ATOM    134  O   GLU A  49       2.442  13.818  16.436  1.00 12.12           O  ",
"ATOM    135  CB  GLU A  49       1.196  10.806  16.029  1.00 14.74           C  ",
"ATOM    136  CG  GLU A  49       0.115   9.950  16.550  1.00 17.40           C  ",
"ATOM    137  CD  GLU A  49       0.706   8.607  17.015  1.00 20.70           C  ",
"ATOM    138  OE1 GLU A  49       1.732   8.490  17.785  1.00 20.18           O  ",
"ATOM    139  OE2 GLU A  49       0.149   7.640  16.556  1.00 25.98           O  ",
"ATOM    140  N   THR A  50       2.709  12.987  14.352  1.00 10.98           N  ",
"ATOM    141  CA  THR A  50       3.863  13.862  14.078  1.00 10.30           C  ",
"ATOM    142  C   THR A  50       3.436  15.324  14.094  1.00  9.95           C  ",
"ATOM    143  O   THR A  50       4.073  16.177  14.693  1.00 10.96           O  ",
"ATOM    144  CB  THR A  50       4.548  13.441  12.792  1.00 10.51           C  ",
"ATOM    145  OG1 THR A  50       4.958  12.066  12.896  1.00 11.72           O  ",
"ATOM    146  CG2 THR A  50       5.755  14.255  12.518  1.00 11.56           C  ",
"ATOM    147  N   GLY A  51       2.310  15.611  13.431  1.00 10.02           N  ",
"ATOM    148  CA  GLY A  51       1.837  16.973  13.396  1.00  9.90           C  ",
"ATOM    149  C   GLY A  51       1.524  17.539  14.752  1.00  9.93           C  ",
"ATOM    150  O   GLY A  51       1.802  18.708  15.023  1.00 10.97           O  ",
"ATOM    151  N   ALA A  52       0.898  16.735  15.617  1.00 10.43           N  ",
"ATOM    152  CA  ALA A  52       0.595  17.195  16.961  1.00 11.37           C  ",
"ATOM    153  C   ALA A  52       1.864  17.533  17.738  1.00 10.84           C  ",
"ATOM    154  O   ALA A  52       1.898  18.497  18.516  1.00 12.56           O  ",
"ATOM    155  CB  ALA A  52      -0.227  16.163  17.724  1.00 12.94           C  ",
"ATOM    156  N   HIS A  53       2.903  16.721  17.566  1.00 10.29           N  ",
"ATOM    157  CA  HIS A  53       4.177  17.033  18.189  1.00  9.80           C  ",
"ATOM    158  C   HIS A  53       4.722  18.368  17.714  1.00 10.19           C  ",
"ATOM    159  O   HIS A  53       5.118  19.215  18.538  1.00 10.33           O  ",
"ATOM    160  CB  HIS A  53       5.176  15.915  17.990  1.00 10.64           C  ",
"ATOM    161  CG  HIS A  53       6.469  16.165  18.686  1.00 10.42           C  ",
"ATOM    162  ND1 HIS A  53       7.681  15.705  18.204  1.00 11.81           N  ",
"ATOM    163  CD2 HIS A  53       6.742  16.842  19.832  1.00 11.10           C  ",
"ATOM    164  CE1 HIS A  53       8.642  16.113  19.022  1.00 11.59           C  ",
"ATOM    165  NE2 HIS A  53       8.096  16.797  20.011  1.00 11.06           N  ",
"ATOM    166  N   ILE A  54       4.769  18.555  16.393  1.00  9.71           N  ",
"ATOM    167  CA  ILE A  54       5.320  19.790  15.859  1.00  9.20           C  ",
"ATOM    168  C   ILE A  54       4.518  20.989  16.404  1.00  9.88           C  ",
"ATOM    169  O   ILE A  54       5.096  21.987  16.840  1.00 10.12           O  ",
"ATOM    170  CB  ILE A  54       5.320  19.775  14.331  1.00  9.72           C  ",
"ATOM    171  CG1 ILE A  54       6.195  18.630  13.821  1.00 11.25           C  ",
"ATOM    172  CG2 ILE A  54       5.745  21.116  13.832  1.00 11.39           C  ",
"ATOM    173  CD1 ILE A  54       6.080  18.414  12.336  1.00 13.38           C  ",
"ATOM    174  N   SER A  55       3.200  20.890  16.375  1.00 10.18           N  ",
"ATOM    175  CA  SER A  55       2.353  21.984  16.849  1.00 10.89           C  ",
"ATOM    176  C   SER A  55       2.662  22.331  18.311  1.00 10.70           C  ",
"ATOM    177  O   SER A  55       2.850  23.491  18.650  1.00 12.40           O  ",
"ATOM    178  CB  SER A  55       0.881  21.600  16.715  1.00 12.88           C  ",
"ATOM    179  OG  SER A  55       0.025  22.603  17.231  1.00 17.85           O  ",
"ATOM    180  N   ASP A  56       2.733  21.316  19.166  1.00 10.26           N  ",
"ATOM    181  CA  ASP A  56       2.988  21.516  20.587  1.00 11.35           C  ",
"ATOM    182  C   ASP A  56       4.407  22.067  20.802  1.00 11.77           C  ",
"ATOM    183  O   ASP A  56       4.624  22.914  21.654  1.00 12.21           O  ",
"ATOM    184  CB  ASP A  56       2.805  20.208  21.368  1.00 13.08           C  ",
"ATOM    185  CG  ASP A  56       1.362  19.807  21.569  1.00 17.02           C  ",
"ATOM    186  OD1 ASP A  56       0.462  20.617  21.320  1.00 21.40           O  ",
"ATOM    187  OD2 ASP A  56       1.167  18.673  22.018  1.00 19.66           O  ",
"ATOM    188  N   ALA A  57       5.371  21.538  20.059  1.00 10.48           N  ",
"ATOM    189  CA  ALA A  57       6.753  21.961  20.233  1.00 11.12           C  ",
"ATOM    190  C   ALA A  57       6.873  23.441  19.868  1.00 10.47           C  ",
"ATOM    191  O   ALA A  57       7.565  24.187  20.570  1.00 10.98           O  ",
"ATOM    192  CB  ALA A  57       7.700  21.102  19.403  1.00 10.34           C  ",
"ATOM    193  N   ILE A  58       6.241  23.859  18.780  1.00  9.42           N  ",
"ATOM    194  CA  ILE A  58       6.291  25.295  18.426  1.00  9.54           C  ",
"ATOM    195  C   ILE A  58       5.635  26.136  19.540  1.00  9.46           C  ",
"ATOM    196  O   ILE A  58       6.162  27.215  19.880  1.00 11.33           O  ",
"ATOM    197  CB  ILE A  58       5.634  25.485  17.025  1.00  9.92           C  ",
"ATOM    198  CG1 ILE A  58       6.548  24.896  15.929  1.00 10.99           C  ",
"ATOM    199  CG2 ILE A  58       5.302  26.954  16.779  1.00 10.16           C  ",
"ATOM    200  CD1 ILE A  58       5.874  24.802  14.571  1.00 11.68           C  ",
"ATOM    201  N   LYS A  59       4.538  25.672  20.117  1.00 10.95           N  ",
"ATOM    202  CA  LYS A  59       3.902  26.422  21.224  1.00 11.52           C  ",
"ATOM    203  C   LYS A  59       4.717  26.379  22.540  1.00 11.95           C  ",
"ATOM    204  O   LYS A  59       4.455  27.166  23.416  1.00 17.26           O  ",
"ATOM    205  CB  LYS A  59       2.476  25.925  21.442  1.00 13.42           C  ",
"ATOM    206  CG  LYS A  59       1.540  26.275  20.292  1.00 16.54           C  ",
"ATOM    207  CD  LYS A  59       0.136  25.849  20.670  1.00 25.52           C  ",
"ATOM    208  CE  LYS A  59      -0.884  26.259  19.640  1.00 35.12           C  ",
"ATOM    209  NZ  LYS A  59      -0.896  25.251  18.579  1.00 39.77           N  ",
"ATOM    210  N   ALA A  60       5.728  25.515  22.605  1.00 11.33           N  ",
"ATOM    211  CA  ALA A  60       6.710  25.459  23.672  1.00 11.93           C  ",
"ATOM    212  C   ALA A  60       7.990  26.240  23.375  1.00 11.71           C  ",
"ATOM    213  O   ALA A  60       9.001  26.178  24.112  1.00 14.41           O  ",
"ATOM    214  CB  ALA A  60       6.981  24.022  24.088  1.00 13.04           C  ",
"ATOM    215  N   GLY A  61       7.986  26.982  22.258  1.00 11.50           N  ",
"ATOM    216  CA  GLY A  61       9.079  27.865  21.931  1.00 12.16           C  ",
"ATOM    217  C   GLY A  61      10.151  27.293  21.061  1.00 11.43           C  ",
"ATOM    218  O   GLY A  61      11.162  27.951  20.816  1.00 13.18           O  ",
"ATOM    219  N   HIS A  62       9.969  26.070  20.571  1.00 10.94           N  ",
"ATOM    220  CA  HIS A  62      10.929  25.525  19.610  1.00 11.03           C  ",
"ATOM    221  C   HIS A  62      10.756  26.205  18.262  1.00 10.70           C  ",
"ATOM    222  O   HIS A  62       9.640  26.566  17.858  1.00 12.02           O  ",
"ATOM    223  CB  HIS A  62      10.759  24.010  19.506  1.00 11.18           C  ",
"ATOM    224  CG  HIS A  62      11.055  23.313  20.781  1.00 11.13           C  ",
"ATOM    225  ND1 HIS A  62      12.334  22.928  21.111  1.00 13.27           N  ",
"ATOM    226  CD2 HIS A  62      10.265  22.982  21.834  1.00 12.34           C  ",
"ATOM    227  CE1 HIS A  62      12.312  22.375  22.319  1.00 12.93           C  ",
"ATOM    228  NE2 HIS A  62      11.068  22.392  22.774  1.00 13.17           N  ",
"ATOM    229  N   ALA A  63      11.858  26.320  17.534  1.00 11.39           N  ",
"ATOM    230  CA  ALA A  63      11.856  27.070  16.256  1.00 10.78           C  ",
"ATOM    231  C   ALA A  63      10.781  26.564  15.309  1.00 10.72           C  ",
"ATOM    232  O   ALA A  63      10.693  25.358  15.054  1.00 11.29           O  ",
"ATOM    233  CB  ALA A  63      13.216  26.984  15.623  1.00 12.62           C  ",
"ATOM    234  N   ASP A  64      10.034  27.471  14.725  1.00 11.28           N  ",
"ATOM    235  CA  ASP A  64       9.091  27.162  13.677  1.00 10.83           C  ",
"ATOM    236  C   ASP A  64       9.692  27.226  12.279  1.00 10.45           C  ",
"ATOM    237  O   ASP A  64       9.081  26.710  11.328  1.00 12.17           O  ",
"ATOM    238  CB  ASP A  64       7.817  28.036  13.756  1.00 13.95           C  ",
"ATOM    239  CG  ASP A  64       8.071  29.586  13.722  1.00 17.09           C  ",
"ATOM    240  OD1 ASP A  64       9.236  30.027  13.750  1.00 19.44           O  ",
"ATOM    241  OD2 ASP A  64       7.086  30.352  13.671  1.00 23.24           O  ",
"ATOM    242  N   VAL A  65      10.881  27.809  12.155  1.00 10.43           N  ",
"ATOM    243  CA  VAL A  65      11.608  27.824  10.889  1.00 10.68           C  ",
"ATOM    244  C   VAL A  65      12.696  26.734  10.896  1.00 11.28           C  ",
"ATOM    245  O   VAL A  65      13.458  26.641  11.865  1.00 13.03           O  ",
"ATOM    246  CB  VAL A  65      12.224  29.210  10.637  1.00 12.85           C  ",
"ATOM    247  CG1 VAL A  65      13.120  29.179   9.404  1.00 14.99           C  ",
"ATOM    248  CG2 VAL A  65      11.121  30.236  10.490  1.00 13.64           C  ",
"ATOM    249  N   CYS A  66      12.759  25.969   9.802  1.00 11.30           N  ",
"ATOM    250  CA  CYS A  66      13.823  25.027   9.566  1.00 11.27           C  ",
"ATOM    251  C   CYS A  66      14.787  25.666   8.556  1.00 10.91           C  ",
"ATOM    252  O   CYS A  66      14.408  25.837   7.391  1.00 11.68           O  ",
"ATOM    253  CB  CYS A  66      13.284  23.728   8.966  1.00 13.28           C  ",
"ATOM    254  SG  CYS A  66      14.616  22.561   8.414  1.00 16.12           S  ",
"ATOM    255  N   THR A  67      16.022  25.920   8.999  1.00 10.44           N  ",
"ATOM    256  CA  THR A  67      17.094  26.336   8.109  1.00 10.70           C  ",
"ATOM    257  C   THR A  67      17.772  25.025   7.686  1.00 10.34           C  ",
"ATOM    258  O   THR A  67      18.448  24.385   8.485  1.00 10.56           O  ",
"ATOM    259  CB  THR A  67      18.089  27.286   8.785  1.00 12.01           C  ",
"ATOM    260  OG1 THR A  67      17.355  28.389   9.330  1.00 13.39           O  ",
"ATOM    261  CG2 THR A  67      19.158  27.766   7.816  1.00 14.19           C  ",
"ATOM    262  N   ILE A  68      17.525  24.614   6.455  1.00  9.48           N  ",
"ATOM    263  CA  ILE A  68      17.918  23.275   6.034  1.00  9.61           C  ",
"ATOM    264  C   ILE A  68      19.439  23.133   6.090  1.00 10.82           C  ",
"ATOM    265  O   ILE A  68      20.166  23.902   5.482  1.00 11.97           O  ",
"ATOM    266  CB  ILE A  68      17.399  22.962   4.626  1.00  9.85           C  ",
"ATOM    267  CG1 ILE A  68      15.867  23.009   4.575  1.00 10.87           C  ",
"ATOM    268  CG2 ILE A  68      17.881  21.590   4.174  1.00 10.80           C  ",
"ATOM    269  CD1 ILE A  68      15.309  22.891   3.195  1.00 11.84           C  ",
"ATOM    270  N   GLU A  69      19.876  22.067   6.769  1.00 10.89           N  ",
"ATOM    271  CA  GLU A  69      21.298  21.754   6.832  1.00 11.62           C  ",
"ATOM    272  C   GLU A  69      21.379  20.225   6.948  1.00 12.16           C  ",
"ATOM    273  O   GLU A  69      21.375  19.675   8.032  1.00 13.82           O  ",
"ATOM    274  CB  GLU A  69      21.990  22.458   7.991  1.00 12.76           C  ",
"ATOM    275  CG  GLU A  69      23.522  22.525   7.867  1.00 14.60           C  ",
"ATOM    276  CD  GLU A  69      24.149  21.136   7.750  1.00 16.46           C  ",
"ATOM    277  OE1 GLU A  69      24.231  20.492   8.809  1.00 20.20           O  ",
"ATOM    278  OE2 GLU A  69      24.563  20.714   6.621  1.00 17.90           O  ",
"ATOM    279  N   ARG A  70      21.459  19.572   5.829  1.00 12.05           N  ",
"ATOM    280  CA  ARG A  70      21.260  18.139   5.780  1.00 12.03           C  ",
"ATOM    281  C   ARG A  70      22.436  17.335   6.298  1.00 12.75           C  ",
"ATOM    282  O   ARG A  70      22.257  16.230   6.813  1.00 15.49           O  ",
"ATOM    283  CB  ARG A  70      20.972  17.738   4.339  1.00 11.97           C  ",
"ATOM    284  CG  ARG A  70      19.609  18.189   3.813  1.00 11.77           C  ",
"ATOM    285  CD  ARG A  70      19.436  17.996   2.317  1.00 12.61           C  ",
"ATOM    286  NE  ARG A  70      18.079  18.271   1.909  1.00 12.33           N  ",
"ATOM    287  CZ  ARG A  70      17.645  19.364   1.289  1.00 12.42           C  ",
"ATOM    288  NH1 ARG A  70      18.457  20.357   0.948  1.00 13.23           N  ",
"ATOM    289  NH2 ARG A  70      16.370  19.471   0.940  1.00 13.99           N  ",
"ATOM    290  N   SER A  71      23.639  17.868   6.169  1.00 14.67           N  ",
"ATOM    291  CA  SER A  71      24.867  17.095   6.447  1.00 17.25           C  ",
"ATOM    292  C   SER A  71      25.043  16.782   7.938  1.00 17.94           C  ",
"ATOM    293  O   SER A  71      25.630  15.782   8.332  1.00 24.25           O  ",
"ATOM    294  CB  SER A  71      26.050  17.836   5.847  1.00 19.82           C  ",
"ATOM    295  OG  SER A  71      26.422  18.898   6.675  1.00 28.70           O  ",
"ATOM    296  N   GLY A  72      24.525  17.649   8.769  1.00 17.08           N  ",
"ATOM    297  CA  GLY A  72      24.703  17.555  10.214  1.00 17.60           C  ",
"ATOM    298  C   GLY A  72      23.704  16.677  10.947  1.00 16.43           C  ",
"ATOM    299  O   GLY A  72      23.699  16.659  12.169  1.00 17.97           O  ",
"ATOM    300  N   ALA A  73      22.801  16.020  10.213  1.00 15.51           N  ",
"ATOM    301  CA  ALA A  73      21.659  15.387  10.839  1.00 16.62           C  ",
"ATOM    302  C   ALA A  73      22.030  14.254  11.795  1.00 17.19           C  ",
"ATOM    303  O   ALA A  73      21.434  14.120  12.859  1.00 19.38           O  ",
"ATOM    304  CB  ALA A  73      20.683  14.932   9.790  1.00 17.21           C  ",
"ATOM    305  N   ASP A  74      23.043  13.452  11.458  1.00 22.28           N  ",
"ATOM    306  CA  ASP A  74      23.475  12.400  12.402  1.00 25.73           C  ",
"ATOM    307  C   ASP A  74      23.965  12.923  13.713  1.00 24.12           C  ",
"ATOM    308  O   ASP A  74      23.589  12.368  14.760  1.00 30.09           O  ",
"ATOM    309  CB  ASP A  74      24.572  11.553  11.796  1.00 29.13           C  ",
"ATOM    310  CG  ASP A  74      24.055  10.665  10.711  1.00 34.40           C  ",
"ATOM    311  OD1 ASP A  74      22.837  10.404  10.641  1.00 41.43           O  ",
"ATOM    312  OD2 ASP A  74      24.874  10.179   9.929  1.00 46.88           O  ",
"ATOM    313  N   LYS A  75      24.827  13.943  13.680  1.00 23.81           N  ",
"ATOM    314  CA  LYS A  75      25.326  14.552  14.916  1.00 22.46           C  ",
"ATOM    315  C   LYS A  75      24.122  15.105  15.679  1.00 24.80           C  ",
"ATOM    316  O   LYS A  75      24.027  14.902  16.889  1.00 29.50           O  ",
"ATOM    317  CB  LYS A  75      26.287  15.737  14.646  1.00 30.01           C  ",
"ATOM    318  CG  LYS A  75      27.731  15.371  14.317  1.00 40.91           C  ",
"ATOM    319  CD  LYS A  75      28.624  16.616  14.208  1.00 47.87           C  ",
"ATOM    320  CE  LYS A  75      28.179  17.642  13.140  1.00 45.85           C  ",
"ATOM    321  NZ  LYS A  75      28.111  17.176  11.704  1.00 34.24           N  ",
"ATOM    322  N   ARG A  76      23.179  15.757  14.988  1.00 17.99           N  ",
"ATOM    323  CA  ARG A  76      22.072  16.329  15.695  1.00 15.56           C  ",
"ATOM    324  C   ARG A  76      21.136  15.295  16.333  1.00 15.13           C  ",
"ATOM    325  O   ARG A  76      20.656  15.517  17.436  1.00 18.16           O  ",
"ATOM    326  CB  ARG A  76      21.266  17.267  14.790  1.00 16.28           C  ",
"ATOM    327  CG  ARG A  76      22.025  18.527  14.412  1.00 15.16           C  ",
"ATOM    328  CD  ARG A  76      21.168  19.564  13.760  1.00 14.67           C  ",
"ATOM    329  NE  ARG A  76      22.026  20.488  13.035  1.00 15.24           N  ",
"ATOM    330  CZ  ARG A  76      22.323  20.448  11.753  1.00 14.54           C  ",
"ATOM    331  NH1 ARG A  76      21.763  19.576  10.912  1.00 14.55           N  ",
"ATOM    332  NH2 ARG A  76      23.199  21.338  11.272  1.00 17.41           N  ",
"ATOM    333  N   ARG A  77      20.884  14.172  15.676  1.00 18.60           N  ",
"ATOM    334  CA  ARG A  77      20.182  13.077  16.339  1.00 22.77           C  ",
"ATOM    335  C   ARG A  77      20.908  12.549  17.550  1.00 25.58           C  ",
"ATOM    336  O   ARG A  77      20.246  12.195  18.552  1.00 27.36           O  ",
"ATOM    337  CB  ARG A  77      20.015  11.913  15.391  1.00 25.93           C  ",
"ATOM    338  CG  ARG A  77      18.879  12.150  14.484  1.00 24.45           C  ",
"ATOM    339  CD  ARG A  77      18.265  10.889  13.850  1.00 30.79           C  ",
"ATOM    340  NE  ARG A  77      17.816  11.266  12.504  1.00 33.79           N  ",
"ATOM    341  CZ  ARG A  77      18.654  11.407  11.479  1.00 30.20           C  ",
"ATOM    342  NH1 ARG A  77      19.939  11.132  11.620  1.00 34.91           N  ",
"ATOM    343  NH2 ARG A  77      18.212  11.803  10.324  1.00 36.26           N  ",
"ATOM    344  N   GLN A  78      22.242  12.435  17.465  1.00 28.13           N  ",
"ATOM    345  CA  GLN A  78      23.023  11.925  18.636  1.00 33.97           C  ",
"ATOM    346  C   GLN A  78      22.686  12.800  19.822  1.00 33.47           C  ",
"ATOM    347  O   GLN A  78      22.396  12.314  20.901  1.00 42.10           O  ",
"ATOM    348  CB  GLN A  78      24.557  12.020  18.451  1.00 39.30           C  ",
"ATOM    349  CG  GLN A  78      25.307  10.885  17.763  1.00 49.80           C  ",
"ATOM    350  CD  GLN A  78      25.063   9.518  18.374  1.00 50.32           C  ",
"ATOM    351  OE1 GLN A  78      24.083   8.837  18.049  1.00 62.39           O  ",
"ATOM    352  NE2 GLN A  78      25.966   9.098  19.249  1.00 58.57           N  ",
"ATOM    353  N   GLU A  79      22.690  14.109  19.587  1.00 27.21           N  ",
"ATOM    354  CA  GLU A  79      22.455  15.081  20.621  1.00 30.44           C  ",
"ATOM    355  C   GLU A  79      21.007  15.037  21.105  1.00 23.43           C  ",
"ATOM    356  O   GLU A  79      20.763  14.951  22.296  1.00 25.22           O  ",
"ATOM    357  CB  GLU A  79      22.806  16.487  20.121  1.00 33.74           C  ",
"ATOM    358  CG  GLU A  79      24.287  16.687  19.740  1.00 43.00           C  ",
"ATOM    359  CD  GLU A  79      25.284  16.603  20.899  1.00 56.18           C  ",
"ATOM    360  OE1 GLU A  79      24.888  16.357  22.063  1.00 71.67           O  ",
"ATOM    361  OE2 GLU A  79      26.497  16.788  20.641  1.00 74.16           O  ",
"ATOM    362  N   SER A  80      20.040  15.059  20.205  1.00 18.47           N  ",
"ATOM    363  CA  SER A  80      18.649  15.213  20.665  1.00 18.60           C  ",
"ATOM    364  C   SER A  80      18.162  13.992  21.383  1.00 16.37           C  ",
"ATOM    365  O   SER A  80      17.353  14.107  22.307  1.00 17.66           O  ",
"ATOM    366  CB  SER A  80      17.705  15.476  19.522  1.00 19.06           C  ",
"ATOM    367  OG  SER A  80      17.866  14.451  18.633  1.00 21.31           O  ",
"ATOM    368  N   LEU A  81      18.668  12.817  20.971  1.00 16.04           N  ",
"ATOM    369  CA  LEU A  81      18.158  11.532  21.481  1.00 17.24           C  ",
"ATOM    370  C   LEU A  81      18.954  10.962  22.641  1.00 19.79           C  ",
"ATOM    371  O   LEU A  81      18.628   9.881  23.118  1.00 20.53           O  ",
"ATOM    372  CB  LEU A  81      18.056  10.522  20.341  1.00 19.18           C  ",
"ATOM    373  CG  LEU A  81      17.161  10.944  19.178  1.00 20.40           C  ",
"ATOM    374  CD1 LEU A  81      17.027   9.781  18.206  1.00 25.01           C  ",
"ATOM    375  CD2 LEU A  81      15.775  11.364  19.615  1.00 21.17           C  ",
"ATOM    376  N   LYS A  82      20.006  11.647  23.061  1.00 21.79           N  ",
"ATOM    377  CA  LYS A  82      20.885  11.130  24.114  1.00 26.48           C  ",
"ATOM    378  C   LYS A  82      20.052  10.746  25.331  1.00 24.96           C  ",
"ATOM    379  O   LYS A  82      19.281  11.545  25.825  1.00 26.26           O  ",
"ATOM    380  CB  LYS A  82      21.930  12.200  24.505  1.00 29.81           C  ",
"ATOM    381  CG  LYS A  82      23.039  11.712  25.439  1.00 32.71           C  ",
"ATOM    382  CD  LYS A  82      23.763  12.886  26.086  1.00 41.51           C  ",
"ATOM    383  CE  LYS A  82      24.834  12.397  27.051  1.00 44.72           C  ",
"ATOM    384  NZ  LYS A  82      25.988  11.788  26.340  1.00 48.96           N  ",
"ATOM    385  N   GLY A  83      20.213   9.518  25.818  1.00 25.97           N  ",
"ATOM    386  CA  GLY A  83      19.551   9.094  27.050  1.00 26.22           C  ",
"ATOM    387  C   GLY A  83      18.055   8.791  26.967  1.00 24.32           C  ",
"ATOM    388  O   GLY A  83      17.444   8.441  27.968  1.00 34.77           O  ",
"ATOM    389  N   ILE A  84      17.470   8.888  25.771  1.00 23.61           N  ",
"ATOM    390  CA  ILE A  84      16.050   8.622  25.575  1.00 21.65           C  ",
"ATOM    391  C   ILE A  84      15.903   7.170  25.095  1.00 21.50           C  ",
"ATOM    392  O   ILE A  84      16.314   6.833  23.990  1.00 22.78           O  ",
"ATOM    393  CB  ILE A  84      15.422   9.595  24.579  1.00 21.87           C  ",
"ATOM    394  CG1 ILE A  84      15.500  11.008  25.147  1.00 23.17           C  ",
"ATOM    395  CG2 ILE A  84      13.962   9.221  24.328  1.00 20.45           C  ",
"ATOM    396  CD1 ILE A  84      14.953  12.075  24.247  1.00 22.91           C  ",
"ATOM    397  N   PRO A  85      15.299   6.311  25.930  1.00 24.84           N  ",
"ATOM    398  CA  PRO A  85      15.221   4.913  25.530  1.00 27.60           C  ",
"ATOM    399  C   PRO A  85      14.257   4.697  24.383  1.00 25.12           C  ",
"ATOM    400  O   PRO A  85      13.374   5.516  24.139  1.00 27.39           O  ",
"ATOM    401  CB  PRO A  85      14.690   4.214  26.775  1.00 29.13           C  ",
"ATOM    402  CG  PRO A  85      13.952   5.267  27.529  1.00 30.25           C  ",
"ATOM    403  CD  PRO A  85      14.646   6.563  27.230  1.00 27.51           C  ",
"ATOM    404  N   THR A  86      14.452   3.596  23.683  1.00 29.56           N  ",
"ATOM    405  CA  THR A  86      13.554   3.221  22.611  1.00 33.08           C  ",
"ATOM    406  C   THR A  86      12.330   2.602  23.210  1.00 33.21           C  ",
"ATOM    407  O   THR A  86      12.324   2.194  24.373  1.00 36.34           O  ",
"ATOM    408  CB  THR A  86      14.197   2.295  21.592  1.00 34.29           C  ",
"ATOM    409  OG1 THR A  86      14.518   1.079  22.229  1.00 42.69           O  ",
"ATOM    410  CG2 THR A  86      15.476   2.948  21.018  1.00 35.37           C  ",
"ATOM    411  N   LYS A  87      11.254   2.632  22.432  1.00 29.53           N  ",
"ATOM    412  CA  LYS A  87       9.984   2.043  22.881  1.00 32.34           C  ",
"ATOM    413  C   LYS A  87       9.364   1.296  21.710  1.00 34.07           C  ",
"ATOM    414  O   LYS A  87       9.349   1.793  20.575  1.00 32.47           O  ",
"ATOM    415  CB  LYS A  87       9.106   3.168  23.395  1.00 33.39           C  ",
"ATOM    416  CG  LYS A  87       7.703   2.860  23.883  1.00 43.91           C  ",
"ATOM    417  CD  LYS A  87       7.144   4.140  24.515  1.00 50.08           C  ",
"ATOM    418  CE  LYS A  87       5.636   4.251  24.439  1.00 52.22           C  ",
"ATOM    419  NZ  LYS A  87       4.988   3.311  25.387  1.00 57.09           N  ",
"ATOM    420  N   PRO A  88       8.880   0.068  21.960  1.00 32.37           N  ",
"ATOM    421  CA  PRO A  88       8.492  -0.727  20.801  1.00 26.65           C  ",
"ATOM    422  C   PRO A  88       7.238  -0.181  20.137  1.00 26.81           C  ",
"ATOM    423  O   PRO A  88       6.312   0.240  20.821  1.00 27.34           O  ",
"ATOM    424  CB  PRO A  88       8.295  -2.145  21.398  1.00 30.26           C  ",
"ATOM    425  CG  PRO A  88       8.031  -1.896  22.860  1.00 30.89           C  ",
"ATOM    426  CD  PRO A  88       8.934  -0.763  23.185  1.00 30.92           C  ",
"ATOM    427  N   GLY A  89       7.235  -0.162  18.810  1.00 23.11           N  ",
"ATOM    428  CA  GLY A  89       6.103   0.383  18.050  1.00 19.81           C  ",
"ATOM    429  C   GLY A  89       6.183   1.898  17.797  1.00 17.43           C  ",
"ATOM    430  O   GLY A  89       5.324   2.442  17.097  1.00 17.38           O  ",
"ATOM    431  N   PHE A  90       7.162   2.575  18.429  1.00 17.35           N  ",
"ATOM    432  CA  PHE A  90       7.278   4.025  18.383  1.00 15.80           C  ",
"ATOM    433  C   PHE A  90       8.653   4.417  17.929  1.00 16.98           C  ",
"ATOM    434  O   PHE A  90       9.650   3.746  18.233  1.00 22.11           O  ",
"ATOM    435  CB  PHE A  90       7.014   4.639  19.757  1.00 18.53           C  ",
"ATOM    436  CG  PHE A  90       5.623   4.480  20.219  1.00 19.08           C  ",
"ATOM    437  CD1 PHE A  90       5.213   3.306  20.792  1.00 27.96           C  ",
"ATOM    438  CD2 PHE A  90       4.706   5.475  20.030  1.00 19.74           C  ",
"ATOM    439  CE1 PHE A  90       3.898   3.122  21.192  1.00 37.34           C  ",
"ATOM    440  CE2 PHE A  90       3.399   5.330  20.413  1.00 22.15           C  ",
"ATOM    441  CZ  PHE A  90       2.992   4.150  20.998  1.00 31.01           C  ",
"ATOM    442  N   ASP A  91       8.702   5.520  17.184  1.00 13.09           N  ",
"ATOM    443  CA  ASP A  91       9.954   6.185  16.916  1.00 13.60           C  ",
"ATOM    444  C   ASP A  91      10.054   7.442  17.785  1.00 12.49           C  ",
"ATOM    445  O   ASP A  91       9.068   7.932  18.265  1.00 14.49           O  ",
"ATOM    446  CB  ASP A  91      10.059   6.532  15.435  1.00 16.51           C  ",
"ATOM    447  CG  ASP A  91      10.183   5.284  14.537  1.00 17.56           C  ",
"ATOM    448  OD1 ASP A  91      10.466   4.151  15.022  1.00 19.63           O  ",
"ATOM    449  OD2 ASP A  91       9.968   5.428  13.349  1.00 18.72           O  ",
"ATOM    450  N   ARG A  92      11.274   7.920  17.992  1.00 13.92           N  ",
"ATOM    451  CA  ARG A  92      11.517   9.134  18.796  1.00 13.52           C  ",
"ATOM    452  C   ARG A  92      11.585  10.319  17.864  1.00 14.72           C  ",
"ATOM    453  O   ARG A  92      12.614  10.550  17.239  1.00 20.23           O  ",
"ATOM    454  CB  ARG A  92      12.823   8.985  19.568  1.00 14.16           C  ",
"ATOM    455  CG  ARG A  92      12.767   7.879  20.618  1.00 16.67           C  ",
"ATOM    456  CD  ARG A  92      14.092   7.444  21.159  1.00 18.43           C  ",
"ATOM    457  NE  ARG A  92      14.848   6.812  20.092  1.00 19.44           N  ",
"ATOM    458  CZ  ARG A  92      16.137   6.533  20.163  1.00 22.11           C  ",
"ATOM    459  NH1 ARG A  92      16.834   6.788  21.269  1.00 25.51           N  ",
"ATOM    460  NH2 ARG A  92      16.741   6.035  19.122  1.00 21.57           N  ",
"ATOM    461  N   ASP A  93      10.470  11.029  17.732  1.00 11.71           N  ",
"ATOM    462  CA  ASP A  93      10.360  12.157  16.827  1.00 12.53           C  ",
"ATOM    463  C   ASP A  93      10.990  13.377  17.446  1.00 10.96           C  ",
"ATOM    464  O   ASP A  93      10.770  13.661  18.635  1.00 13.28           O  ",
"ATOM    465  CB  ASP A  93       8.883  12.492  16.564  1.00 12.78           C  ",
"ATOM    466  CG  ASP A  93       8.702  13.577  15.474  1.00 14.23           C  ",
"ATOM    467  OD1 ASP A  93       9.212  13.366  14.329  1.00 16.83           O  ",
"ATOM    468  OD2 ASP A  93       8.098  14.632  15.788  1.00 15.24           O  ",
"ATOM    469  N   GLU A  94      11.690  14.148  16.615  1.00 11.49           N  ",
"ATOM    470  CA  GLU A  94      12.389  15.342  17.062  1.00 11.18           C  ",
"ATOM    471  C   GLU A  94      11.832  16.580  16.369  1.00 10.09           C  ",
"ATOM    472  O   GLU A  94      11.668  16.591  15.141  1.00 11.69           O  ",
"ATOM    473  CB  GLU A  94      13.862  15.270  16.627  1.00 14.85           C  ",
"ATOM    474  CG  GLU A  94      14.665  14.085  17.018  1.00 16.79           C  ",
"ATOM    475  CD  GLU A  94      15.868  13.874  16.102  1.00 21.15           C  ",
"ATOM    476  OE1 GLU A  94      16.840  14.690  16.143  1.00 20.28           O  ",
"ATOM    477  OE2 GLU A  94      15.787  12.952  15.261  1.00 28.95           O  ",
"ATOM    478  N   TRP A  95      11.629  17.644  17.143  1.00  9.84           N  ",
"ATOM    479  CA  TRP A  95      11.333  18.957  16.588  1.00  9.54           C  ",
"ATOM    480  C   TRP A  95      12.117  19.968  17.428  1.00  9.80           C  ",
"ATOM    481  O   TRP A  95      11.882  20.028  18.647  1.00 10.80           O  ",
"ATOM    482  CB  TRP A  95       9.857  19.299  16.555  1.00  9.22           C  ",
"ATOM    483  CG  TRP A  95       9.661  20.600  15.879  1.00  9.81           C  ",
"ATOM    484  CD1 TRP A  95       9.702  21.861  16.448  1.00 10.50           C  ",
"ATOM    485  CD2 TRP A  95       9.546  20.792  14.476  1.00  9.53           C  ",
"ATOM    486  NE1 TRP A  95       9.585  22.813  15.461  1.00 10.90           N  ",
"ATOM    487  CE2 TRP A  95       9.526  22.187  14.239  1.00  9.99           C  ",
"ATOM    488  CE3 TRP A  95       9.511  19.910  13.377  1.00 10.03           C  ",
"ATOM    489  CZ2 TRP A  95       9.403  22.711  12.962  1.00 11.57           C  ",
"ATOM    490  CZ3 TRP A  95       9.409  20.444  12.113  1.00 12.14           C  ",
"ATOM    491  CH2 TRP A  95       9.324  21.828  11.920  1.00 13.19           C  ",
"ATOM    492  N   PRO A  96      12.965  20.806  16.840  1.00  9.33           N  ",
"ATOM    493  CA  PRO A  96      13.261  20.882  15.415  1.00  9.35           C  ",
"ATOM    494  C   PRO A  96      13.885  19.588  14.891  1.00  9.49           C  ",
"ATOM    495  O   PRO A  96      14.501  18.835  15.647  1.00 10.41           O  ",
"ATOM    496  CB  PRO A  96      14.182  22.101  15.297  1.00  9.90           C  ",
"ATOM    497  CG  PRO A  96      13.901  22.890  16.566  1.00  9.96           C  ",
"ATOM    498  CD  PRO A  96      13.648  21.869  17.602  1.00  9.77           C  ",
"ATOM    499  N   MET A  97      13.740  19.413  13.598  1.00  9.08           N  ",
"ATOM    500  CA  MET A  97      14.249  18.213  12.948  1.00  9.36           C  ",
"ATOM    501  C   MET A  97      15.768  18.220  12.860  1.00  9.59           C  ",
"ATOM    502  O   MET A  97      16.431  19.237  12.795  1.00 10.61           O  ",
"ATOM    503  CB  MET A  97      13.706  18.115  11.511  1.00 10.06           C  ",
"ATOM    504  CG  MET A  97      12.191  18.017  11.480  1.00 11.12           C  ",
"ATOM    505  SD  MET A  97      11.524  17.886   9.820  1.00 13.65           S  ",
"ATOM    506  CE  MET A  97      12.018  19.454   9.127  1.00 15.04           C  ",
"ATOM    507  N   ALA A  98      16.330  17.018  12.784  1.00 10.05           N  ",
"ATOM    508  CA  ALA A  98      17.770  16.857  12.714  1.00 10.21           C  ",
"ATOM    509  C   ALA A  98      18.353  17.422  11.410  1.00  9.81           C  ",
"ATOM    510  O   ALA A  98      19.523  17.783  11.398  1.00 11.53           O  ",
"ATOM    511  CB  ALA A  98      18.128  15.397  12.839  1.00 11.62           C  ",
"ATOM    512  N   MET A  99      17.557  17.538  10.330  1.00 10.35           N  ",
"ATOM    513  CA  MET A  99      17.998  18.125   9.072  1.00 10.92           C  ",
"ATOM    514  C   MET A  99      17.973  19.635   9.087  1.00 10.93           C  ",
"ATOM    515  O   MET A  99      18.169  20.230   8.023  1.00 12.02           O  ",
"ATOM    516  CB  MET A  99      17.193  17.579   7.868  1.00 10.94           C  ",
"ATOM    517  CG  MET A  99      15.744  17.952   7.903  1.00 11.76           C  ",
"ATOM    518  SD  MET A  99      14.835  17.519   6.433  1.00 13.31           S  ",
"ATOM    519  CE  MET A  99      15.564  18.628   5.220  1.00 14.59           C  ",
"ATOM    520  N   CYS A 100      17.750  20.268  10.237  1.00 11.15           N  ",
"ATOM    521  CA  CYS A 100      17.708  21.723  10.354  1.00 11.84           C  ",
"ATOM    522  C   CYS A 100      18.800  22.202  11.289  1.00 11.06           C  ",
"ATOM    523  O   CYS A 100      19.101  21.547  12.281  1.00 12.39           O  ",
"ATOM    524  CB  CYS A 100      16.372  22.102  10.905  1.00 13.43           C  ",
"ATOM    525  SG  CYS A 100      14.948  21.413  10.078  1.00 15.57           S  ",
"ATOM    526  N   GLU A 101      19.315  23.389  11.007  1.00 11.14           N  ",
"ATOM    527  CA  GLU A 101      20.298  23.987  11.913  1.00 12.77           C  ",
"ATOM    528  C   GLU A 101      19.807  24.073  13.355  1.00 12.30           C  ",
"ATOM    529  O   GLU A 101      20.575  24.014  14.352  1.00 15.64           O  ",
"ATOM    530  CB  GLU A 101      20.616  25.416  11.417  1.00 15.76           C  ",
"ATOM    531  CG  GLU A 101      21.426  25.445  10.162  1.00 20.52           C  ",
"ATOM    532  CD  GLU A 101      22.057  26.781   9.782  1.00 27.66           C  ",
"ATOM    533  OE1 GLU A 101      21.862  27.803  10.500  1.00 46.89           O  ",
"ATOM    534  OE2 GLU A 101      22.746  26.761   8.732  1.00 34.25           O  ",
"ATOM    535  N   GLU A 102      18.497  24.279  13.464  1.00 11.51           N  ",
"ATOM    536  CA  GLU A 102      17.830  24.488  14.764  1.00 12.78           C  ",
"ATOM    537  C   GLU A 102      17.673  23.220  15.617  1.00 11.73           C  ",
"ATOM    538  O   GLU A 102      17.292  23.300  16.766  1.00 14.10           O  ",
"ATOM    539  CB  GLU A 102      16.435  25.079  14.507  1.00 13.57           C  ",
"ATOM    540  CG  GLU A 102      16.421  26.498  13.880  1.00 14.02           C  ",
"ATOM    541  CD  GLU A 102      16.634  26.599  12.361  1.00 14.19           C  ",
"ATOM    542  OE1 GLU A 102      16.714  25.545  11.670  1.00 13.05           O  ",
"ATOM    543  OE2 GLU A 102      16.624  27.756  11.830  1.00 17.79           O  ",
"ATOM    544  N   GLY A 103      17.928  22.063  15.007  1.00 11.37           N  ",
"ATOM    545  CA  GLY A 103      17.852  20.788  15.671  1.00 12.56           C  ",
"ATOM    546  C   GLY A 103      19.005  20.468  16.581  1.00 12.51           C  ",
"ATOM    547  O   GLY A 103      19.884  21.317  16.857  1.00 15.27           O  ",
"ATOM    548  N   GLY A 104      19.005  19.221  17.060  1.00 12.32           N  ",
"ATOM    549  CA  GLY A 104      20.057  18.792  17.975  1.00 13.54           C  ",
"ATOM    550  C   GLY A 104      19.748  19.066  19.439  1.00 14.38           C  ",
"ATOM    551  O   GLY A 104      18.636  18.802  19.922  1.00 16.59           O  ",
"ATOM    552  N   LYS A 105      20.747  19.546  20.165  1.00 19.07           N  ",
"ATOM    553  CA  LYS A 105      20.615  19.783  21.607  1.00 20.74           C  ",
"ATOM    554  C   LYS A 105      19.403  20.637  21.872  1.00 19.01           C  ",
"ATOM    555  O   LYS A 105      19.158  21.644  21.197  1.00 20.03           O  ",
"ATOM    556  CB  LYS A 105      21.880  20.442  22.189  1.00 28.64           C  ",
"ATOM    557  CG  LYS A 105      23.023  19.476  22.379  1.00 38.36           C  ",
"ATOM    558  CD  LYS A 105      24.016  20.013  23.401  1.00 44.93           C  ",
"ATOM    559  CE  LYS A 105      25.173  19.042  23.620  1.00 55.13           C  ",
"ATOM    560  NZ  LYS A 105      26.180  19.555  24.595  1.00 64.22           N  ",
"ATOM    561  N   GLY A 106      18.606  20.193  22.823  1.00 17.68           N  ",
"ATOM    562  CA  GLY A 106      17.441  20.945  23.237  1.00 16.32           C  ",
"ATOM    563  C   GLY A 106      16.151  20.634  22.505  1.00 13.39           C  ",
"ATOM    564  O   GLY A 106      15.095  21.078  22.946  1.00 15.24           O  ",
"ATOM    565  N   ALA A 107      16.226  19.855  21.430  1.00 12.51           N  ",
"ATOM    566  CA  ALA A 107      15.023  19.606  20.688  1.00 11.60           C  ",
"ATOM    567  C   ALA A 107      13.980  18.862  21.500  1.00 11.17           C  ",
"ATOM    568  O   ALA A 107      14.312  18.034  22.346  1.00 13.33           O  ",
"ATOM    569  CB  ALA A 107      15.316  18.831  19.433  1.00 12.40           C  ",
"ATOM    570  N   SER A 108      12.719  19.123  21.168  1.00 10.43           N  ",
"ATOM    571  CA  SER A 108      11.624  18.342  21.744  1.00 11.65           C  ",
"ATOM    572  C   SER A 108      11.615  16.947  21.152  1.00 10.93           C  ",
"ATOM    573  O   SER A 108      11.651  16.813  19.927  1.00 11.52           O  ",
"ATOM    574  CB  SER A 108      10.287  18.999  21.427  1.00 11.72           C  ",
"ATOM    575  OG  SER A 108       9.144  18.335  21.961  1.00 12.22           O  ",
"ATOM    576  N   VAL A 109      11.500  15.956  22.016  1.00 11.59           N  ",
"ATOM    577  CA  VAL A 109      11.419  14.562  21.577  1.00 12.55           C  ",
"ATOM    578  C   VAL A 109      10.132  13.942  22.089  1.00 12.52           C  ",
"ATOM    579  O   VAL A 109       9.781  14.114  23.264  1.00 15.42           O  ",
"ATOM    580  CB  VAL A 109      12.631  13.766  22.017  1.00 13.71           C  ",
"ATOM    581  CG1 VAL A 109      12.519  12.320  21.585  1.00 13.98           C  ",
"ATOM    582  CG2 VAL A 109      13.915  14.389  21.483  1.00 15.14           C  ",
"ATOM    583  N   ARG A 110       9.417  13.251  21.189  1.00 12.54           N  ",
"ATOM    584  CA  ARG A 110       8.157  12.611  21.576  1.00 12.56           C  ",
"ATOM    585  C   ARG A 110       8.108  11.245  20.918  1.00 11.92           C  ",
"ATOM    586  O   ARG A 110       8.374  11.157  19.709  1.00 12.99           O  ",
"ATOM    587  CB  ARG A 110       6.950  13.449  21.178  1.00 13.28           C  ",
"ATOM    588  CG  ARG A 110       5.618  12.850  21.653  1.00 14.13           C  ",
"ATOM    589  CD  ARG A 110       4.404  13.600  21.168  1.00 14.43           C  ",
"ATOM    590  NE  ARG A 110       4.321  14.920  21.751  1.00 15.45           N  ",
"ATOM    591  CZ  ARG A 110       3.265  15.710  21.626  1.00 14.78           C  ",
"ATOM    592  NH1 ARG A 110       2.249  15.373  20.826  1.00 15.41           N  ",
"ATOM    593  NH2 ARG A 110       3.258  16.882  22.256  1.00 16.93           N  ",
"ATOM    594  N   TYR A 111       7.615  10.227  21.647  1.00 12.42           N  ",
"ATOM    595  CA  TYR A 111       7.342   8.940  21.011  1.00 13.19           C  ",
"ATOM    596  C   TYR A 111       6.146   9.094  20.068  1.00 12.82           C  ",
"ATOM    597  O   TYR A 111       5.082   9.562  20.510  1.00 15.26           O  ",
"ATOM    598  CB  TYR A 111       7.044   7.866  22.029  1.00 14.69           C  ",
"ATOM    599  CG  TYR A 111       8.216   7.534  22.921  1.00 15.70           C  ",
"ATOM    600  CD1 TYR A 111       9.375   6.952  22.385  1.00 16.75           C  ",
"ATOM    601  CD2 TYR A 111       8.193   7.813  24.306  1.00 19.33           C  ",
"ATOM    602  CE1 TYR A 111      10.479   6.638  23.212  1.00 18.51           C  ",
"ATOM    603  CE2 TYR A 111       9.284   7.506  25.128  1.00 21.39           C  ",
"ATOM    604  CZ  TYR A 111      10.416   6.924  24.587  1.00 20.46           C  ",
"ATOM    605  OH  TYR A 111      11.499   6.630  25.402  1.00 28.21           O  ",
"ATOM    606  N   VAL A 112       6.328   8.702  18.806  1.00 12.29           N  ",
"ATOM    607  CA  VAL A 112       5.273   8.769  17.785  1.00 12.46           C  ",
"ATOM    608  C   VAL A 112       5.200   7.400  17.106  1.00 12.74           C  ",
"ATOM    609  O   VAL A 112       6.230   6.796  16.848  1.00 12.62           O  ",
"ATOM    610  CB  VAL A 112       5.630   9.852  16.729  1.00 13.12           C  ",
"ATOM    611  CG1 VAL A 112       4.683   9.796  15.525  1.00 15.41           C  ",
"ATOM    612  CG2 VAL A 112       5.573  11.251  17.357  1.00 13.35           C  ",
"ATOM    613  N   SER A 113       3.996   6.881  16.880  1.00 12.52           N  ",
"ATOM    614  CA ASER A 113       3.838   5.595  16.196  0.50 12.63           C  ",
"ATOM    615  CA BSER A 113       3.898   5.566  16.257  0.50 13.28           C  ",
"ATOM    616  C   SER A 113       4.715   5.512  14.975  1.00 13.17           C  ",
"ATOM    617  O   SER A 113       4.735   6.445  14.176  1.00 14.04           O  ",
"ATOM    618  CB ASER A 113       2.388   5.427  15.775  0.50 14.80           C  ",
"ATOM    619  CB BSER A 113       2.446   5.140  16.044  0.50 16.70           C  ",
"ATOM    620  OG ASER A 113       1.549   5.465  16.919  0.50 16.38           O  ",
"ATOM    621  OG BSER A 113       1.796   5.910  15.062  0.50 20.79           O  ",
"ATOM    622  N   SER A 114       5.388   4.376  14.779  1.00 13.42           N  ",
"ATOM    623  CA  SER A 114       6.397   4.283  13.745  1.00 15.23           C  ",
"ATOM    624  C   SER A 114       5.846   4.539  12.359  1.00 15.66           C  ",
"ATOM    625  O   SER A 114       6.513   5.228  11.569  1.00 17.22           O  ",
"ATOM    626  CB  SER A 114       7.150   2.952  13.802  1.00 16.60           C  ",
"ATOM    627  OG  SER A 114       7.796   2.815  15.040  1.00 21.79           O  ",
"ATOM    628  N   SER A 115       4.650   4.048  12.030  1.00 16.37           N  ",
"ATOM    629  CA ASER A 115       4.123   4.267  10.683  0.50 17.27           C  ",
"ATOM    630  CA BSER A 115       4.118   4.259  10.684  0.50 17.51           C  ",
"ATOM    631  C   SER A 115       3.849   5.742  10.414  1.00 17.90           C  ",
"ATOM    632  O   SER A 115       4.121   6.232   9.327  1.00 19.11           O  ",
"ATOM    633  CB ASER A 115       2.850   3.447  10.460  0.50 19.42           C  ",
"ATOM    634  CB BSER A 115       2.832   3.449  10.483  0.50 22.15           C  ",
"ATOM    635  OG ASER A 115       3.164   2.075  10.497  0.50 28.14           O  ",
"ATOM    636  OG BSER A 115       1.872   3.794  11.463  0.50 24.69           O  ",
"ATOM    637  N   ASP A 116       3.323   6.441  11.420  1.00 16.13           N  ",
"ATOM    638  CA  ASP A 116       3.056   7.876  11.302  1.00 14.28           C  ",
"ATOM    639  C   ASP A 116       4.344   8.625  11.130  1.00 13.58           C  ",
"ATOM    640  O   ASP A 116       4.460   9.465  10.207  1.00 14.20           O  ",
"ATOM    641  CB  ASP A 116       2.289   8.319  12.561  1.00 13.60           C  ",
"ATOM    642  CG  ASP A 116       1.829   9.763  12.562  1.00 14.01           C  ",
"ATOM    643  OD1 ASP A 116       2.648  10.689  12.560  1.00 13.86           O  ",
"ATOM    644  OD2 ASP A 116       0.600   9.992  12.763  1.00 15.25           O  ",
"ATOM    645  N   ASN A 117       5.332   8.338  11.981  1.00 13.93           N  ",
"ATOM    646  CA  ASN A 117       6.576   9.054  11.916  1.00 13.56           C  ",
"ATOM    647  C   ASN A 117       7.392   8.775  10.657  1.00 12.94           C  ",
"ATOM    648  O   ASN A 117       8.003   9.691  10.077  1.00 13.77           O  ",
"ATOM    649  CB  ASN A 117       7.414   8.815  13.144  1.00 14.47           C  ",
"ATOM    650  CG  ASN A 117       8.365   9.951  13.408  1.00 15.92           C  ",
"ATOM    651  OD1 ASN A 117       9.536   9.720  13.557  1.00 17.13           O  ",
"ATOM    652  ND2 ASN A 117       7.860  11.188  13.460  1.00 15.47           N  ",
"ATOM    653  N   ARG A 118       7.423   7.517  10.216  1.00 14.72           N  ",
"ATOM    654  CA  ARG A 118       8.194   7.162   9.041  1.00 14.46           C  ",
"ATOM    655  C   ARG A 118       7.550   7.731   7.774  1.00 14.85           C  ",
"ATOM    656  O   ARG A 118       8.239   8.213   6.876  1.00 17.06           O  ",
"ATOM    657  CB  ARG A 118       8.391   5.640   8.981  1.00 14.32           C  ",
"ATOM    658  CG  ARG A 118       9.311   5.151  10.088  1.00 15.03           C  ",
"ATOM    659  CD  ARG A 118       9.355   3.651  10.253  1.00 15.54           C  ",
"ATOM    660  NE  ARG A 118       9.969   3.290  11.508  1.00 15.78           N  ",
"ATOM    661  CZ  ARG A 118      10.217   2.048  11.892  1.00 14.04           C  ",
"ATOM    662  NH1 ARG A 118       9.990   1.049  11.071  1.00 15.68           N  ",
"ATOM    663  NH2 ARG A 118      10.779   1.833  13.065  1.00 14.66           N  ",
"ATOM    664  N   GLY A 119       6.232   7.744   7.769  1.00 15.60           N  ",
"ATOM    665  CA  GLY A 119       5.533   8.352   6.661  1.00 14.46           C  ",
"ATOM    666  C   GLY A 119       5.800   9.855   6.603  1.00 14.21           C  ",
"ATOM    667  O   GLY A 119       6.071  10.390   5.546  1.00 14.70           O  ",
"ATOM    668  N   ALA A 120       5.738  10.517   7.763  1.00 13.18           N  ",
"ATOM    669  CA  ALA A 120       5.980  11.952   7.812  1.00 11.46           C  ",
"ATOM    670  C   ALA A 120       7.395  12.218   7.348  1.00 11.46           C  ",
"ATOM    671  O   ALA A 120       7.671  13.185   6.629  1.00 12.50           O  ",
"ATOM    672  CB  ALA A 120       5.766  12.452   9.234  1.00 12.03           C  ",
"ATOM    673  N   GLY A 121       8.343  11.404   7.823  1.00 12.44           N  ",
"ATOM    674  CA  GLY A 121       9.743  11.624   7.438  1.00 13.23           C  ",
"ATOM    675  C   GLY A 121       9.974  11.563   5.929  1.00 14.48           C  ",
"ATOM    676  O   GLY A 121      10.723  12.373   5.362  1.00 14.88           O  ",
"ATOM    677  N   SER A 122       9.347  10.575   5.277  1.00 15.19           N  ",
"ATOM    678  CA ASER A 122       9.466  10.430   3.829  0.50 16.59           C  ",
"ATOM    679  CA BSER A 122       9.474  10.417   3.837  0.50 16.64           C  ",
"ATOM    680  C   SER A 122       8.781  11.590   3.104  1.00 13.83           C  ",
"ATOM    681  O   SER A 122       9.252  12.100   2.079  1.00 15.49           O  ",
"ATOM    682  CB ASER A 122       8.889   9.118   3.318  0.50 20.40           C  ",
"ATOM    683  CB BSER A 122       8.879   9.080   3.400  0.50 21.01           C  ",
"ATOM    684  OG ASER A 122       9.071   9.069   1.896  0.50 23.62           O  ",
"ATOM    685  OG BSER A 122       9.579   7.991   3.983  0.50 25.92           O  ",
"ATOM    686  N   TRP A 123       7.620  11.983   3.607  1.00 12.35           N  ",
"ATOM    687  CA  TRP A 123       6.910  13.108   3.039  1.00 11.58           C  ",
"ATOM    688  C   TRP A 123       7.756  14.373   3.098  1.00 10.90           C  ",
"ATOM    689  O   TRP A 123       7.862  15.092   2.107  1.00 11.38           O  ",
"ATOM    690  CB  TRP A 123       5.572  13.297   3.764  1.00 13.41           C  ",
"ATOM    691  CG  TRP A 123       4.726  14.329   3.129  1.00 14.14           C  ",
"ATOM    692  CD1 TRP A 123       3.768  14.114   2.209  1.00 17.46           C  ",
"ATOM    693  CD2 TRP A 123       4.788  15.752   3.322  1.00 13.80           C  ",
"ATOM    694  NE1 TRP A 123       3.192  15.298   1.849  1.00 20.80           N  ",
"ATOM    695  CE2 TRP A 123       3.847  16.329   2.459  1.00 16.83           C  ",
"ATOM    696  CE3 TRP A 123       5.561  16.595   4.116  1.00 15.01           C  ",
"ATOM    697  CZ2 TRP A 123       3.603  17.700   2.433  1.00 19.10           C  ",
"ATOM    698  CZ3 TRP A 123       5.356  17.950   4.029  1.00 17.44           C  ",
"ATOM    699  CH2 TRP A 123       4.414  18.488   3.176  1.00 17.16           C  ",
"ATOM    700  N   VAL A 124       8.349  14.643   4.278  1.00 11.08           N  ",
"ATOM    701  CA  VAL A 124       9.151  15.843   4.434  1.00 12.47           C  ",
"ATOM    702  C   VAL A 124      10.385  15.790   3.546  1.00 12.55           C  ",
"ATOM    703  O   VAL A 124      10.730  16.775   2.891  1.00 12.99           O  ",
"ATOM    704  CB  VAL A 124       9.522  16.062   5.912  1.00 14.46           C  ",
"ATOM    705  CG1 VAL A 124      10.623  17.111   6.068  1.00 17.84           C  ",
"ATOM    706  CG2 VAL A 124       8.292  16.409   6.706  1.00 14.62           C  ",
"ATOM    707  N   GLY A 125      11.062  14.648   3.506  1.00 13.90           N  ",
"ATOM    708  CA  GLY A 125      12.275  14.592   2.669  1.00 15.88           C  ",
"ATOM    709  C   GLY A 125      11.985  14.842   1.214  1.00 18.17           C  ",
"ATOM    710  O   GLY A 125      12.692  15.590   0.544  1.00 20.86           O  ",
"ATOM    711  N   ASN A 126      10.889  14.268   0.718  1.00 14.70           N  ",
"ATOM    712  CA  ASN A 126      10.505  14.529  -0.645  1.00 14.53           C  ",
"ATOM    713  C   ASN A 126      10.055  15.943  -0.915  1.00 13.73           C  ",
"ATOM    714  O   ASN A 126      10.437  16.536  -1.911  1.00 16.96           O  ",
"ATOM    715  CB  ASN A 126       9.502  13.467  -1.095  1.00 16.67           C  ",
"ATOM    716  CG  ASN A 126      10.237  12.311  -1.649  1.00 18.18           C  ",
"ATOM    717  OD1 ASN A 126      10.878  12.435  -2.721  1.00 23.91           O  ",
"ATOM    718  ND2 ASN A 126      10.263  11.239  -0.894  1.00 20.90           N  ",
"ATOM    719  N   ARG A 127       9.230  16.489  -0.015  1.00 11.73           N  ",
"ATOM    720  CA  ARG A 127       8.751  17.830  -0.217  1.00 12.22           C  ",
"ATOM    721  C   ARG A 127       9.887  18.840  -0.180  1.00 15.17           C  ",
"ATOM    722  O   ARG A 127       9.835  19.828  -0.906  1.00 22.87           O  ",
"ATOM    723  CB  ARG A 127       7.672  18.164   0.810  1.00 13.32           C  ",
"ATOM    724  CG  ARG A 127       6.970  19.508   0.616  1.00 16.13           C  ",
"ATOM    725  CD  ARG A 127       6.023  19.442  -0.629  1.00 19.86           C  ",
"ATOM    726  NE  ARG A 127       5.297  20.697  -0.880  1.00 25.19           N  ",
"ATOM    727  CZ  ARG A 127       5.831  21.752  -1.483  1.00 32.19           C  ",
"ATOM    728  NH1 ARG A 127       7.111  21.739  -1.874  1.00 30.62           N  ",
"ATOM    729  NH2 ARG A 127       5.106  22.858  -1.667  1.00 36.83           N  ",
"ATOM    730  N   LEU A 128      10.836  18.654   0.731  1.00 14.49           N  ",
"ATOM    731  CA  LEU A 128      11.936  19.616   0.870  1.00 15.93           C  ",
"ATOM    732  C   LEU A 128      13.030  19.428  -0.129  1.00 16.71           C  ",
"ATOM    733  O   LEU A 128      13.902  20.269  -0.252  1.00 16.56           O  ",
"ATOM    734  CB  LEU A 128      12.477  19.602   2.303  1.00 16.46           C  ",
"ATOM    735  CG  LEU A 128      11.522  20.105   3.415  1.00 17.79           C  ",
"ATOM    736  CD1 LEU A 128      12.235  20.170   4.779  1.00 20.21           C  ",
"ATOM    737  CD2 LEU A 128      10.906  21.443   3.023  1.00 20.05           C  ",
"ATOM    738  N   ASN A 129      12.944  18.373  -0.922  1.00 20.86           N  ",
"ATOM    739  CA  ASN A 129      13.954  18.156  -1.944  1.00 18.59           C  ",
"ATOM    740  C   ASN A 129      13.965  19.204  -3.045  1.00 18.29           C  ",
"ATOM    741  O   ASN A 129      14.925  19.279  -3.836  1.00 17.87           O  ",
"ATOM    742  CB  ASN A 129      13.785  16.777  -2.584  1.00 19.64           C  ",
"ATOM    743  CG  ASN A 129      14.965  16.419  -3.482  1.00 17.97           C  ",
"ATOM    744  OD1 ASN A 129      16.090  16.467  -3.061  1.00 19.31           O  ",
"ATOM    745  ND2 ASN A 129      14.687  16.026  -4.691  1.00 20.02           N  ",
"ATOM    746  N   GLY A 130      12.942  20.058  -3.104  1.00 21.23           N  ",
"ATOM    747  CA  GLY A 130      12.963  21.192  -3.975  1.00 22.76           C  ",
"ATOM    748  C   GLY A 130      13.728  22.404  -3.455  1.00 17.72           C  ",
"ATOM    749  O   GLY A 130      13.856  23.387  -4.143  1.00 24.02           O  ",
"ATOM    750  N   TYR A 131      14.301  22.302  -2.263  1.00 14.58           N  ",
"ATOM    751  CA  TYR A 131      14.966  23.447  -1.649  1.00 12.97           C  ",
"ATOM    752  C   TYR A 131      16.414  23.085  -1.303  1.00 11.82           C  ",
"ATOM    753  O   TYR A 131      16.689  22.061  -0.680  1.00 12.81           O  ",
"ATOM    754  CB  TYR A 131      14.248  23.825  -0.357  1.00 12.95           C  ",
"ATOM    755  CG  TYR A 131      12.811  24.217  -0.520  1.00 13.31           C  ",
"ATOM    756  CD1 TYR A 131      12.449  25.442  -1.081  1.00 14.15           C  ",
"ATOM    757  CD2 TYR A 131      11.833  23.408  -0.063  1.00 15.35           C  ",
"ATOM    758  CE1 TYR A 131      11.109  25.811  -1.194  1.00 17.22           C  ",
"ATOM    759  CE2 TYR A 131      10.490  23.761  -0.143  1.00 18.54           C  ",
"ATOM    760  CZ  TYR A 131      10.145  24.956  -0.716  1.00 16.85           C  ",
"ATOM    761  OH  TYR A 131       8.839  25.280  -0.770  1.00 21.97           O  ",
"ATOM    762  N   ALA A 132      17.315  23.993  -1.657  1.00 12.61           N  ",
"ATOM    763  CA  ALA A 132      18.717  23.839  -1.338  1.00 12.72           C  ",
"ATOM    764  C   ALA A 132      18.978  23.979   0.149  1.00 11.82           C  ",
"ATOM    765  O   ALA A 132      18.273  24.722   0.865  1.00 12.64           O  ",
"ATOM    766  CB  ALA A 132      19.531  24.882  -2.058  1.00 13.38           C  ",
"ATOM    767  N   ASP A 133      20.065  23.376   0.613  1.00 12.30           N  ",
"ATOM    768  CA  ASP A 133      20.536  23.667   1.949  1.00 11.87           C  ",
"ATOM    769  C   ASP A 133      20.681  25.186   2.140  1.00 12.02           C  ",
"ATOM    770  O   ASP A 133      21.067  25.895   1.216  1.00 13.85           O  ",
"ATOM    771  CB  ASP A 133      21.854  23.001   2.254  1.00 12.70           C  ",
"ATOM    772  CG  ASP A 133      21.742  21.485   2.417  1.00 12.76           C  ",
"ATOM    773  OD1 ASP A 133      21.290  20.827   1.454  1.00 14.08           O  ",
"ATOM    774  OD2 ASP A 133      22.141  21.003   3.506  1.00 15.29           O  ",
"ATOM    775  N   GLY A 134      20.332  25.642   3.340  1.00 10.96           N  ",
"ATOM    776  CA  GLY A 134      20.339  27.057   3.637  1.00 10.82           C  ",
"ATOM    777  C   GLY A 134      19.004  27.738   3.519  1.00  9.88           C  ",
"ATOM    778  O   GLY A 134      18.798  28.853   4.032  1.00 11.29           O  ",
"ATOM    779  N   THR A 135      18.082  27.125   2.795  1.00 11.09           N  ",
"ATOM    780  CA  THR A 135      16.739  27.674   2.656  1.00 10.56           C  ",
"ATOM    781  C   THR A 135      16.050  27.604   3.998  1.00 10.53           C  ",
"ATOM    782  O   THR A 135      16.170  26.577   4.695  1.00 11.21           O  ",
"ATOM    783  CB  THR A 135      15.945  26.902   1.590  1.00 12.03           C  ",
"ATOM    784  OG1 THR A 135      16.666  26.850   0.378  1.00 14.91           O  ",
"ATOM    785  CG2 THR A 135      14.611  27.569   1.326  1.00 13.10           C  ",
"ATOM    786  N   ARG A 136      15.302  28.662   4.340  1.00  9.86           N  ",
"ATOM    787  CA  ARG A 136      14.577  28.745   5.593  1.00  9.69           C  ",
"ATOM    788  C   ARG A 136      13.087  28.505   5.333  1.00 10.23           C  ",
"ATOM    789  O   ARG A 136      12.468  29.264   4.610  1.00 10.95           O  ",
"ATOM    790  CB  ARG A 136      14.816  30.131   6.213  1.00 10.70           C  ",
"ATOM    791  CG  ARG A 136      16.298  30.375   6.506  1.00 14.43           C  ",
"ATOM    792  CD  ARG A 136      16.600  31.799   6.739  1.00 15.29           C  ",
"ATOM    793  NE  ARG A 136      16.349  32.581   5.539  1.00 14.36           N  ",
"ATOM    794  CZ  ARG A 136      16.289  33.907   5.517  1.00 18.41           C  ",
"ATOM    795  NH1 ARG A 136      16.557  34.594   6.614  1.00 21.09           N  ",
"ATOM    796  NH2 ARG A 136      15.990  34.510   4.377  1.00 23.25           N  ",
"ATOM    797  N   ILE A 137      12.556  27.411   5.875  1.00  9.33           N  ",
"ATOM    798  CA  ILE A 137      11.172  27.020   5.676  1.00  9.16           C  ",
"ATOM    799  C   ILE A 137      10.399  27.182   6.973  1.00  9.81           C  ",
"ATOM    800  O   ILE A 137      10.699  26.528   7.953  1.00 10.75           O  ",
"ATOM    801  CB  ILE A 137      11.098  25.575   5.168  1.00 10.64           C  ",
"ATOM    802  CG1 ILE A 137      11.918  25.390   3.855  1.00 13.68           C  ",
"ATOM    803  CG2 ILE A 137       9.660  25.113   4.960  1.00 12.91           C  ",
"ATOM    804  CD1 ILE A 137      11.447  26.173   2.648  1.00 17.58           C  ",
"ATOM    805  N   LEU A 138       9.344  28.017   6.923  1.00 10.00           N  ",
"ATOM    806  CA  LEU A 138       8.411  28.131   8.041  1.00 10.19           C  ",
"ATOM    807  C   LEU A 138       7.403  26.987   7.967  1.00  9.61           C  ",
"ATOM    808  O   LEU A 138       6.697  26.854   6.957  1.00 12.39           O  ",
"ATOM    809  CB  LEU A 138       7.676  29.466   7.957  1.00 12.22           C  ",
"ATOM    810  CG  LEU A 138       6.581  29.648   9.000  1.00 13.38           C  ",
"ATOM    811  CD1 LEU A 138       7.101  29.722  10.386  1.00 14.98           C  ",
"ATOM    812  CD2 LEU A 138       5.796  30.926   8.681  1.00 15.40           C  ",
"ATOM    813  N   PHE A 139       7.340  26.201   9.047  1.00  9.98           N  ",
"ATOM    814  CA  PHE A 139       6.334  25.151   9.148  1.00 10.43           C  ",
"ATOM    815  C   PHE A 139       5.121  25.707   9.865  1.00 11.82           C  ",
"ATOM    816  O   PHE A 139       5.267  26.324  10.915  1.00 12.80           O  ",
"ATOM    817  CB  PHE A 139       6.859  23.953   9.938  1.00 11.48           C  ",
"ATOM    818  CG  PHE A 139       7.660  22.970   9.115  1.00 12.52           C  ",
"ATOM    819  CD1 PHE A 139       8.838  23.329   8.475  1.00 13.68           C  ",
"ATOM    820  CD2 PHE A 139       7.255  21.621   9.088  1.00 15.98           C  ",
"ATOM    821  CE1 PHE A 139       9.552  22.396   7.715  1.00 14.85           C  ",
"ATOM    822  CE2 PHE A 139       7.977  20.671   8.379  1.00 18.52           C  ",
"ATOM    823  CZ  PHE A 139       9.117  21.072   7.676  1.00 16.51           C  ",
"ATOM    824  N   ILE A 140       3.946  25.457   9.308  1.00 11.60           N  ",
"ATOM    825  CA  ILE A 140       2.682  25.887   9.872  1.00 13.47           C  ",
"ATOM    826  C   ILE A 140       1.814  24.660  10.106  1.00 13.09           C  ",
"ATOM    827  O   ILE A 140       1.565  23.920   9.163  1.00 13.54           O  ",
"ATOM    828  CB  ILE A 140       1.987  26.868   8.917  1.00 16.68           C  ",
"ATOM    829  CG1 ILE A 140       2.860  28.093   8.725  1.00 19.80           C  ",
"ATOM    830  CG2 ILE A 140       0.604  27.198   9.414  1.00 18.09           C  ",
"ATOM    831  CD1 ILE A 140       2.270  29.181   7.878  1.00 24.51           C  ",
"ATOM    832  N   VAL A 141       1.362  24.441  11.339  1.00 14.51           N  ",
"ATOM    833  CA  VAL A 141       0.448  23.346  11.711  1.00 16.70           C  ",
"ATOM    834  C   VAL A 141      -0.652  24.007  12.523  1.00 22.92           C  ",
"ATOM    835  O   VAL A 141      -1.796  23.784  12.263  1.00 37.89           O  ",
"ATOM    836  CB  VAL A 141       1.086  22.283  12.643  1.00 23.19           C  ",
"ATOM    837  CG1 VAL A 141       0.316  20.973  12.582  1.00 25.66           C  ",
"ATOM    838  CG2 VAL A 141       2.551  22.087  12.399  1.00 26.05           C  ",
"ATOM    839  N   GLN A 142      -0.202  24.729  13.537  1.00 30.17           N  ",
"ATOM    840  CA  GLN A 142      -0.849  25.863  14.249  1.00 43.44           C  ",
"ATOM    841  C   GLN A 142      -1.266  25.367  15.610  1.00 53.97           C  ",
"ATOM    842  O   GLN A 142      -2.348  25.603  16.175  1.00 73.17           O  ",
"ATOM    843  CB  GLN A 142      -1.971  26.585  13.489  1.00 51.01           C  ",
"ATOM    844  CG  GLN A 142      -1.444  27.425  12.340  1.00 51.77           C  ",
"ATOM    845  CD  GLN A 142      -2.436  27.561  11.197  1.00 56.41           C  ",
"ATOM    846  OE1 GLN A 142      -3.152  26.612  10.853  1.00 68.90           O  ",
"ATOM    847  NE2 GLN A 142      -2.473  28.741  10.588  1.00 44.85           N  ",
"ATOM    848  OXT GLN A 142      -0.406  24.689  16.162  1.00 41.66           O  "

	};
	public static void main(String[] args){
		HashMap<String,String> firstcode = new HashMap<>();
		for(int ii = 0;ii < sourcelines.length;ii++){
			String[] pt = sourcelines[ii].split("");
			String aacode = pt[17]+pt[18]+pt[19];
			String rnum = pt[23]+pt[24]+pt[25];
			if(!firstcode.containsKey(aacode)){
				
				System.out.println("\""+sourcelines[ii]+"\",");
				firstcode.put(aacode,rnum);	
			}else{
				if(rnum.equals(firstcode.get(aacode))){
					System.out.println("\""+sourcelines[ii]+"\",");
				}
			}
		}
		if(firstcode.size() != 20){
			throw new RuntimeException("residue ga tari nai");
		}
	}
	
	
}

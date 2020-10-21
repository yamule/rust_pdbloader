
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Pattern;
import static pepbuilderj.TemplateBaseModeller.makeAtomLines;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class Preparesubmit {
	public static void main_forgromacs(String[] args){
	//public static void main(String[] args){
		//String infilename = "C:\\dummy\\vbox_share\\casp13\\queries\\T1004\\try2\\T1004_5m9f_m1.pdb.1.c.pdb.0.cc.pdb";
		String infilename = 
				"C:\\dummy\\vbox_share\\casp13\\queries\\T1009\\model1.pfas.10565_minim.pdb.pr.pdb.mul.63.pdb.0.b.pdb";
		

		String outfilename = infilename+".z.pdb";
		PDBData p = PDBData.loadPDBFile(infilename);
		ArrayList<PDBResidue> allres = new ArrayList<>();
		
		int offset = 0;
		int resnum = 1;
		for(String c:p.chains.keySet()){
			PDBChain cc = p.chains.get(c);
			cc.name = "A";
			if(c.equals("B")){
				offset += 1000;
				resnum = 1;
			}
			if(c.equals("C")){
				offset += 1000;
				resnum = 1;
			}
			
			if(c.equals("D")){
				offset += 1000;
				resnum = 1;
			}
			
			
			if(c.equals("E")){
				offset += 1000;
				resnum = 1;
			}
			
			if(c.equals("F")){
				offset += 1000;
				resnum = 1;
			}
			
			if(c.equals("G")){
				offset += 1000;
				resnum = 1;
			}
			if(c.equals("H")){
				offset += 1000;
				resnum = 1;
			}
			
			for(PDBResidue r:cc.residues){
				//r.setResidueNumber(resnum++);
				r.setResidueNumber(r.getResidueNumber()+offset);
				if(PDBResidue.validAASet.contains(r.getName())){
					Iterator<PDBAtom> a = r.atoms.iterator();
					while(a.hasNext()){
						PDBAtom at = a.next();
						if(at.pdb_atom_code.equals("CD")&&r.getName().equals("ILE")){
							at.pdb_atom_code = "CD1";
						}
						if(at.pdb_atom_code.indexOf("H") == 0
								|| at.pdb_atom_code.indexOf("1H") == 0
								|| at.pdb_atom_code.indexOf("2H") == 0
								|| at.pdb_atom_code.indexOf("3H") == 0
								|| at.pdb_atom_code.indexOf("OW") == 0
								
								){
							a.remove();
						}
					}
					allres.add(r);
				}
			}
		}
		writeToFile(makeAtomLines(allres,allres.get(0).getResidueNumber(),1,"A",true),
					outfilename);

	}
	public static void main_fromgromacsinputmultimer(String[] args){
	//public static void main(String[] args){
		
		String infilename = "C:\\dummy\\vbox_share\\casp13\\queries\\T1009\\step42c.pdb";
		String outfilename = infilename+".pr.pdb";
		PDBData p = PDBData.loadPDBFile(infilename);
		ArrayList<PDBResidue> allres = new ArrayList<>();
		
		PDBChain ac = new PDBChain("A");
		PDBChain bc = new PDBChain("B");
		PDBChain cc = new PDBChain("C");
		PDBChain dc = new PDBChain("D");
		PDBChain ec = new PDBChain("E");
		PDBChain fc = new PDBChain("F");
		PDBChain gc = new PDBChain("G");
		PDBChain hc = new PDBChain("H");
		for(String c:p.chains.keySet()){
			PDBChain ccc = p.chains.get(c);
					
			for(PDBResidue r:ccc.residues){
				if(r.getResidueNumber() < 1000){
					r.setParent(ac);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				
				if(r.getResidueNumber() > 7000){
					r.setParent(hc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(r.getResidueNumber() > 6000){
					r.setParent(gc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(r.getResidueNumber() > 5000){
					r.setParent(fc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(r.getResidueNumber() > 4000){
					r.setParent(ec);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(r.getResidueNumber() > 3000){
					r.setParent(dc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				
				if(r.getResidueNumber() > 2000){
					r.setParent(cc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(r.getResidueNumber() > 1000){
					r.setParent(bc);
					r.setResidueNumber(r.getResidueNumber()%1000);
				}
				if(PDBResidue.validAASet.contains(r.getName())){
					Iterator<PDBAtom> a = r.atoms.iterator();
					while(a.hasNext()){
						PDBAtom at = a.next();
						if(at.pdb_atom_code.equals("O1")){
							at.pdb_atom_code = "O";
						}
						if(at.pdb_atom_code.equals("CD")&&r.getName().equals("ILE")){
							at.pdb_atom_code = "CD1";
						}
						if(at.pdb_atom_code.indexOf("H") == 0
								|| at.pdb_atom_code.indexOf("1H") == 0
								|| at.pdb_atom_code.indexOf("2H") == 0
								|| at.pdb_atom_code.indexOf("3H") == 0
								|| at.pdb_atom_code.indexOf("OW") == 0
								|| at.pdb_atom_code.indexOf("O2") == 0
								){
							a.remove();
						}
					}
					allres.add(r);
				}
				//if(r.getResidueNumber() > 300){
				//	r.atoms.clear();
				//}
			}
		}
		AddScoreOnPDB as = new AddScoreOnPDB();
		as.addScoreOnBFactor(allres);
		System.out.println(allres.size());
		ArrayList<String> head = new ArrayList<>();
//PFRMAT     TS
//TARGET     T0949
//AUTHOR     3711-3160-9340
//METHOD     Template based modelling
//MODEL      1
//PARENT     1CC3_B//チェーンも必要
		writeToFile(makeAtomLines(allres,allres.get(0).getResidueNumber(),1,"A",true),
					outfilename);
//TER
//END
	}
	public static void main(String[] args){
	//public static void main_normal(String[] args){
	
	
		String dirname = "C:\\dummy\\vbox_share\\casp13\\queries\\T1015s1\\singles\\submitter\\fin";
		File dir = new File(dirname);
		if(dir.isDirectory()){
			File[] flist = dir.listFiles();
			//flist = new File[1];
			//flist[0] = new File("C:\\dummy\\vbox_share\\casp13\\queries\\T1013\\test1.pdb");
			for(File f:flist){
				String infilename = f.getAbsolutePath();

						if(Pattern.compile("pdb$").matcher(infilename).find()){
				String outfilename = infilename+".pr.pdb";
				PDBData p = PDBData.loadPDBFile(infilename);
				ArrayList<PDBResidue> allres = new ArrayList<>();

				int offset = 0;
				for(String c:p.chains.keySet()){
					PDBChain cc = p.chains.get(c);
					if(cc.name.equals(" ") || cc.name.equals("")){ 
						cc.name = "A";
					}
					for(PDBResidue r:cc.residues){
						r.setResidueNumber(r.getResidueNumber()+offset);
						if(PDBResidue.validAASet.contains(r.getName())){
							Iterator<PDBAtom> a = r.atoms.iterator();
							while(a.hasNext()){
								PDBAtom at = a.next();
								if(at.pdb_atom_code.equals("O1")){
									at.pdb_atom_code = "O";
								}

								if(at.pdb_atom_code.equals("CD")&&r.getName().equals("ILE")){
									at.pdb_atom_code = "CD1";
								}
								if(at.pdb_atom_code.indexOf("H") == 0
										|| at.pdb_atom_code.indexOf("1H") == 0
										|| at.pdb_atom_code.indexOf("2H") == 0
										|| at.pdb_atom_code.indexOf("3H") == 0
										|| at.pdb_atom_code.indexOf("OW") == 0
										|| at.pdb_atom_code.indexOf("O2") == 0
										){
									a.remove();
								}
							}
							allres.add(r);
						}
					}
				}
				AddScoreOnPDB as = new AddScoreOnPDB();
				as.addScoreOnBFactor(allres);
				ArrayList<String> head = new ArrayList<>();
		//PFRMAT     TS
		//TARGET     T0949
		//AUTHOR     3711-3160-9340
		//METHOD     Template based modelling
		//MODEL      1
		//PARENT     1CC3_B//チェーンも必要
				StringBuffer sb = makeAtomLines(allres,allres.get(0).getResidueNumber(),1,"A",true);
				sb.append("TER\nEND\n");
				writeToFile(sb,
							outfilename);
		//TER
		//END
					}
				}
		}
	}
}

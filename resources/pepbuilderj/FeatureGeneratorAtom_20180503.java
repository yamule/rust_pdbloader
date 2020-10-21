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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class FeatureGeneratorAtom_20180503 implements FeatureGenerator{
	
	//バックボーンからそこに来るべき残基を予測する
	double jcbdist = 2.5;
	double tcbdist = 2.5;
	double tcb2dist = 5.0;
	
	double maxdist = 12.0;
	
	//CA から jcbdist A 先にある偽の CB
	HashMap<PDBResidue,PDBAtom> jCB = new HashMap<>();
	
	//CA から tcbdist A 先にある注目している残基の偽の CB
	HashMap<PDBResidue,PDBAtom> tCB = new HashMap<>();
	
	//CA から tcbdist2 A 先にある注目している残基の偽の CB 
	HashMap<PDBResidue,PDBAtom> t2CB = new HashMap<>();
	
	
	
	//Feature に使われる Atom だけ入っている
	public static final int BACKATOM_C = 0;
	public static final int BACKATOM_O = 1;
	public static final int BACKATOM_N = 2;
	public static final int BACKATOM_CA = 3;
	//public static final int BACKATOM_JCB = 4;
	ArrayList<ArrayList<PDBAtom>> backboneAtoms = new ArrayList<>();
	static HashSet<String> validResidue = new HashSet<>();
	
	static String[] aaname = {
		"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
		,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
	};
	ArrayList<ArrayList<PDBAtom>> selectedAtoms = new ArrayList<>();
	static HashMap<String,Integer> selectedAtoms_index = new HashMap<>();
	static ArrayList<String> selectedAtoms_label = new ArrayList<>();
	static{
		String[] att = {"ASN_OD1","ASN_CB","ASP_OD1","ASP_CB","CYS_SG","CYS_CB","PRO_CD"
				,"PRO_CB","PRO_CG","HIS_CE1","HIS_CB","HIS_CG","TRP_NE1","TRP_CB","TRP_CG"
				,"THR_OG1","THR_CB","LEU_CD2","LEU_CB","GLU_OE1","GLU_CB","ILE_CG2","ILE_CD1"
				,"ILE_CB","ARG_CB","ARG_NH1","ARG_CG","PHE_CE1","PHE_CB","PHE_CG","VAL_CG2","VAL_CB","LYS_CB","LYS_CG","GLN_OE1","GLN_CB","GLN_CG","ALA_CB","SER_OG","SER_CB","MET_CE","MET_CG","MET_CB","TYR_OH","TYR_CB","TYR_CG","GLY_CA"
			};
		
		for(int ii = 0;ii < att.length;ii++){
			selectedAtoms_label.add(att[ii]);
			selectedAtoms_index.put(att[ii],ii);
		}
		
		for(String aa:aaname){
			validResidue.add(aa);
		}
	}
	
	public ArrayList<String> getResourceFilePath(){
		String treefile="resources/tree_simple.0.dat";
		String scorefile="resources/scores.0.dat";
		ArrayList<String> ret = new ArrayList<>();
		ret.add(treefile);
		ret.add(scorefile);
		return ret;
	}
	FeatureGeneratorAtom_20180503(){
		
		
	}
	
	public String transCode(String s){
		String[] az = s.split("[^A-Za-z0-9]");
		return az[0]+"_"+az[1];
	}
	
	public void prepare(Collection<PDBResidue> res){
		backboneAtoms.clear();
		tCB.clear();
		
		for(int ii = 0;ii < 4;ii++){
			backboneAtoms.add(new ArrayList<PDBAtom>());
		}
		selectedAtoms.clear();
		for(int ii = 0;ii < selectedAtoms_label.size();ii++){
			selectedAtoms.add(new ArrayList<PDBAtom>());
		}
		
		for(PDBResidue r:res){
			if(r.isLigand || r.isMissing()){
				continue;
			}
			PDBAtom ca = r.getCA();
			PDBAtom c = r.getC();
			PDBAtom n = r.getN();
			PDBAtom o = r.getO();
			if(ca !=  null){
				backboneAtoms.get(BACKATOM_CA).add(ca);
			}
			
			if(n !=  null){
				backboneAtoms.get(BACKATOM_N).add(n);
			}
			
			if(o !=  null){
				backboneAtoms.get(BACKATOM_O).add(o);
			}
			if(c !=  null){
				backboneAtoms.get(BACKATOM_C).add(c);
			}
			for(PDBAtom a:r.atoms){
				if(a.isAlternative()){
				}else{
					String code = r.getName()+"_"+a.pdb_atom_code;
					if(selectedAtoms_index.containsKey(code)){
						selectedAtoms.get(selectedAtoms_index.get(code)).add(a);
					}
				}
			}
		}
	}
	
	public static ArrayList<String> getHeader(){
		
		ArrayList<String> headers = new ArrayList<>();
		//headers.add("CA_count3");
		//headers.add("CA_count6");
		//headers.add("CA_count9");
		//headers.add("CA_count3diff");
		//headers.add("CA_count6diff");
		//headers.add("CA_count9diff");
		String[] b = {"","","",""};
		b[BACKATOM_CA] = "CA";
		b[BACKATOM_C] = "C";
		b[BACKATOM_N] = "N";
		b[BACKATOM_O] = "O";
		for(String bb:b){
			for(int ii = 0;ii < 3;ii++){
				headers.add(bb+"_dist"+String.valueOf(ii));
			}
		}
		for(String r:selectedAtoms_label){
			for(int ii = 0;ii < 3;ii++){
				headers.add(r+"_dist"+String.valueOf(ii));
			}
		}
		return headers;
	}
	
	
	public ArrayList<PDBAtom> getTargetAtoms(PDBResidue target){
		ArrayList<PDBAtom> ret = new ArrayList<>();
		
		for(PDBAtom aa:target.atoms){
			if(aa.isAlternative()){
				continue;
			}else{
				if(selectedAtoms_label.contains(target.getName()+"_"+aa.pdb_atom_code)){
					ret.add(aa);
				}
			}
		}
		return ret;
	}
	
	//二つは別々の generateFeatures から来たものでないと Parent が同一だったりするかも
	public ArrayList<Double> mergeFeatures(ArrayList<Double> f1,ArrayList<Double> f2){
		ArrayList<Double> ret = new ArrayList<>();
		//backbone について近い奴三つ
		for(int ii = 0;ii < 4;ii++){
			ArrayList<Double> ss = new ArrayList<>();
			for(int jj = 0;jj < 3;jj++){
				ss.add(f1.get(ii*3+jj));
				ss.add(f2.get(ii*3+jj));
			}
			Collections.sort(ss);
			for(int jj = 0;jj < 3;jj++){
				ret.add(ss.get(jj));
			}
			
		}
		int offset = ret.size();
		for(int ii = 0;ii < selectedAtoms_label.size();ii++){
			ArrayList<Double> ss = new ArrayList<>();
			for(int jj = 0;jj < 3;jj++){
				ss.add(f1.get(offset+ii*3+jj));
				ss.add(f2.get(offset+ii*3+jj));
			}
			Collections.sort(ss);
			for(int jj = 0;jj < 3;jj++){
				ret.add(ss.get(jj));
			}		
		}
		return ret;
		
	}
	
	
	
	public ArrayList<Double> generateFeatures_merge(PDBAtom target,PDBAtom target2){ArrayList<Double> ret = new ArrayList<>();
		if(true){	
			throw new RuntimeException("使わない予定");
		}
	//使わない予定。
		//backbone について近い奴三つ
		for(int ii = 0;ii < 4;ii++){
			ArrayList<PDBAtom> al = backboneAtoms.get(ii);
			ArrayList<AtomDistance> dlis = new ArrayList<>();
			for(PDBAtom a:al){
				if(a.parent == target.parent){
					continue;
				}
				double ddist = a.distance(target);
				if(ddist < 12.0){
					dlis.add(new AtomDistance(a,ddist,ddist));
				}
			}
			for(PDBAtom a:al){
				if(a.parent == target2.parent){
					continue;
				}
				double ddist = a.distance(target2);
				if(ddist < 12.0){
					dlis.add(new AtomDistance(a,ddist,ddist));
				}
			}
			
			Collections.sort(dlis,new AtomDistanceComparator());
			for(int jj = 0;jj < 3;jj++){
				if(dlis.size() > jj){
					ret.add(dlis.get(jj).distance);
				}else{
					ret.add(1000.0);
				}
			}
		}
		for(int ii = 0;ii < selectedAtoms_label.size();ii++){
			ArrayList<PDBAtom> aal = selectedAtoms.get(ii);
			ArrayList<AtomDistance> dlis = new ArrayList<>();
			for(PDBAtom aa:aal){
				if(aa.parent == target.parent){
					continue;
				}
				double tdis = aa.distance(target);
				dlis.add(new AtomDistance(aa,tdis,tdis));
			}
			
			
			for(PDBAtom aa:aal){
				if(aa.parent == target2.parent){
					continue;
				}
				double tdis = aa.distance(target2);
				dlis.add(new AtomDistance(aa,tdis,tdis));
			}
			
			Collections.sort(dlis,new AtomDistanceComparator());
			for(int jj = 0;jj < 3;jj++){
				if(dlis.size() > jj){
					ret.add(dlis.get(jj).distance);
				}else{
					ret.add(1000.0);
				}
			}
		}
		return ret;
	}
	
	
	public double[] generateFeaturesA_merge(PDBAtom target,PDBAtom target2){
		return listToArray(generateFeatures_merge(target,target2));
	}
	
	
	public double[] generateFeaturesA(PDBAtom target){
		return listToArray(generateFeatures(target));
	}
	public ArrayList<Double> generateFeatures(PDBAtom target){
		ArrayList<Double> ret = new ArrayList<>();
		
		//backbone について近い奴三つ
		for(int ii = 0;ii < 4;ii++){
			ArrayList<PDBAtom> al = backboneAtoms.get(ii);
			ArrayList<AtomDistance> dlis = new ArrayList<>();
			for(PDBAtom a:al){
				if(a.parent == target.parent){
					continue;
				}
				double ddist = a.distance(target);
				if(ddist < 12.0){
					dlis.add(new AtomDistance(a,ddist,ddist));
				}
			}
			Collections.sort(dlis,new AtomDistanceComparator());
			for(int jj = 0;jj < 3;jj++){
				if(dlis.size() > jj){
					ret.add(dlis.get(jj).distance);
				}else{
					ret.add(1000.0);
				}
			}
		}
		for(int ii = 0;ii < selectedAtoms_label.size();ii++){
			ArrayList<PDBAtom> aal = selectedAtoms.get(ii);
			ArrayList<AtomDistance> dlis = new ArrayList<>();
			for(PDBAtom aa:aal){
				if(aa.parent == target.parent){
					continue;
				}
				double tdis = aa.distance(target);
				dlis.add(new AtomDistance(aa,tdis,tdis));
			}
			Collections.sort(dlis,new AtomDistanceComparator());
			for(int jj = 0;jj < 3;jj++){
				if(dlis.size() > jj){
					ret.add(dlis.get(jj).distance);
				}else{
					ret.add(1000.0);
				}
			}
		}
		return ret;
	}
	public static double[] listToArray(ArrayList<Double> d){
		if(d == null){
			return null;
		}
		double[] ret = new double[d.size()];
		for(int ii = 0;ii < d.size();ii++){
			ret[ii] = d.get(ii);
		}
		return ret;
	}
	
	public static ArrayList<String> countAtoms_PDB(String treefile,File[] lis){
		RandomForestProcess rg = new RandomForestProcess();
		rg.loadForestFromFile_SimpleFormat(treefile);
		HashMap<String,Integer> allanswers = new HashMap<>();//回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		for(File f:lis){
			if(!Pattern.compile("\\.(pdb|ent)$").matcher(f.getPath()).find()){
				continue;
			}
			System.out.println(f.getPath());
			PDBData p = PDBData.loadPDBFile(f.getPath());
			FeatureGeneratorAtom_20180503 gen = new FeatureGeneratorAtom_20180503();
			for(String cc:p.chains.keySet()){
				PDBChain pc = p.chains.get(cc);
				gen.prepare(pc.residues);
				for(PDBResidue r:pc.residues){
					if(r.isLigand || r.isMissing()){
						continue;
					}
					for(PDBAtom aa:r.atoms){
						if(aa.isAlternative()){
							continue;
						}
						if(selectedAtoms_label.contains(r.getName()+"_"+aa.pdb_atom_code)){
							
					
							ArrayList<Double> res = gen.generateFeatures(aa);
							if(res == null){
								continue;
							}
							String ans = r.getName()+"_"+aa.pdb_atom_code;
							String pred = rg.predict_Majority(listToArray(res));

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
					}
				}
			}
		}
//#count all: 195323
//#count LEU_CB: 18664
//trueatom_predictedas ALA_CB@5:	ALA_CB	342	GLY_PCB	197	SER_CB	108	LEU_CB	63	PHE_CB	48	VAL_CB	44	THR_CB	43	PRO_CB	40	TYR_CB	35	CYS_CB	28	ILE_CB	24	MET_CB	17	HIS_CB	17	ASN_CB	16	ASP_CB	16	LYS_CB	15	GLN_CB	12	TRP_CB	12	GLU_CB	10	ARG_CB	3

		ArrayList<String> lines = new ArrayList<>();
		int allcount = 0;
		for(String s:allanswers.keySet()){
			lines.add("#count "+s+":"+allanswers.get(s));
			allcount += allanswers.get(s);
		}
		lines.add("#count all:"+allcount);
		int posicount = 0;
		int negacount = 0;
		int zerocount = 0;
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
			
			double predcount = 0;
			for(VSorter v:al){
				sb.append("\t"+labels.get(v.index)+"\t"+v.val);
				predcount+=v.val;
			}
			
			for(String hh:hs.keySet()){
				double backd = allanswers.get(hh)/(double)allcount;
				double predd = hs.get(hh)/(double)predcount;
				if(predd/backd > Math.sqrt(2)){
					posicount+=hs.get(hh);
				}else if(predd/backd < 1.0/Math.sqrt(2)){
					negacount+=hs.get(hh);
				}else{
					zerocount+=hs.get(hh);
				}
			}
			lines.add(sb.toString());
		}
		System.out.println("posi\t"+posicount+"\t"+"zero\t" + zerocount + "\t"+"nega\t"+negacount);
		return lines;
	}
	
	public static ArrayList<String> countAtoms(String treefile,String tablefile){
		RandomForestProcess rg = new RandomForestProcess();
		rg.loadForestFromFile_SimpleFormat(treefile);
		FuzzyDecisionTree co = FuzzyDecisionTree.loadTable(tablefile,true,true);
		HashMap<String,Integer> allanswers = new HashMap<>();//回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		for(DirtySample d:co.samples){
			String ans = d.classLabel;
			String pred = rg.predict_Majority(listToArray(d.values));
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
//#count all: 195323
//#count LEU_CB: 18664
//trueatom_predictedas ALA_CB@5:	ALA_CB	342	GLY_PCB	197	SER_CB	108	LEU_CB	63	PHE_CB	48	VAL_CB	44	THR_CB	43	PRO_CB	40	TYR_CB	35	CYS_CB	28	ILE_CB	24	MET_CB	17	HIS_CB	17	ASN_CB	16	ASP_CB	16	LYS_CB	15	GLN_CB	12	TRP_CB	12	GLU_CB	10	ARG_CB	3

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
	
	public static void makeScoreTable(String[] args){
		
		ArrayList<String> res = countAtoms(
		"C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\0506\\tree_simple.0.dat_0"
		,"C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\0506\\table_merged.dat");
		for(String s:res){
			System.out.println(s);
		}
		printStrings("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\0506\\scores.0.dat",res);
		
		
		String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\";
		File dir = new File(sampledirname);
			File[] lis = dir.listFiles();
		ArrayList<String> dres = countAtoms_PDB(
		"C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\0506\\tree_simple.0.dat_0"
		,lis);
		
		printStrings("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\0506\\scores.0.dat_c",dres);
		
	}
	public static void makeInputTable(String[] args){
		String outfilename = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_maketree_sidechain\\table_merged.dat";
		String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered\\";
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
			new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 
			
			
			
			File dir = new File(sampledirname);
			File[] lis = dir.listFiles();
			
			ArrayList<String> header = FeatureGeneratorAtom_20180503.getHeader();
			pw.write("name\t");
			for(String h:header){
				pw.write(h+"\t");
			}
			pw.write("target\n");
			
			for(File f:lis){
				if(!Pattern.compile("\\.(pdb|ent)$").matcher(f.getPath()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				FeatureGeneratorAtom_20180503 gen = new FeatureGeneratorAtom_20180503();
				for(String cc:p.chains.keySet()){
					PDBChain pc = p.chains.get(cc);
					gen.prepare(pc.residues);
					for(PDBResidue r:pc.residues){
						if(r.isLigand || r.isMissing()){
							continue;
						}
						
						ArrayList<PDBAtom> vals = new ArrayList<>();
						for(PDBAtom aa:r.atoms){
							if(aa.isAlternative()){
								continue;
							}
							if(selectedAtoms_label.contains(r.getName()+"_"+aa.pdb_atom_code)){
								vals.add(aa);
							}
						}
						for(PDBAtom aa:vals){
							ArrayList<Double> res = gen.generateFeatures(aa);
							if(res == null){
								continue;
							}
							pw.write(p.getID()+"_"+cc+"_"+r.getRepresentativeCode()+"_"+aa.pdb_atom_code+"\t");
							for(Double d:res){
								pw.write(String.valueOf((float)((double)d))+"\t");
							}
							pw.write(r.getName()+"_"+aa.pdb_atom_code+"\n");
						}
					}
				}
			}
			pw.close();
			
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void main(String[] args){
		//makeInputTable(args);
		makeScoreTable(args);
	}
}

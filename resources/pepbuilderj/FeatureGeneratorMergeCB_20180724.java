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
import java.util.Iterator;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class FeatureGeneratorMergeCB_20180724 implements FeatureGenerator{
	
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
	public static final int ENVATOM_C = 0;
	public static final int ENVATOM_O = 1;
	public static final int ENVATOM_N = 2;
	public static final int ENVATOM_CA = 3;
	public static final int ENVATOM_JCB = 4;
	ArrayList<ArrayList<PDBAtom>> envAtoms = new ArrayList<>();
	HashSet<String> validResidue = new HashSet<>();
	
	static String[] aaname = {
		"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
		,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
	};
	
	public ArrayList<String> getResourceFilePath(){
		ArrayList<String> ret = new ArrayList<>();
		String treefile="resources/mergecb_0724/reftree.0.0.fin.dat";
		String scorefile="resources/mergecb_0724/reftree.0.0.fin.dat.scores";
		ret.add(treefile);
		ret.add(scorefile);
		
		return ret;
	}
	public String transCode(String s){
		return s.split("_")[0];
	}
	
	FeatureGeneratorMergeCB_20180724(){
		
		
		for(String aa:aaname){
			validResidue.add(aa);
		}
	}
	
	
	
	public ArrayList<Double> mergeFeatures(ArrayList<Double> a,ArrayList<Double> b){
		return null;
	}
	
	public ArrayList<Double> generateFeatures_merge(PDBAtom target,PDBAtom target2){
		//not supported currently;
		return null;
	}
	public double[] generateFeaturesA_merge(PDBAtom target,PDBAtom target2){
		//not supported currently;
		return null;
	}
	
	
	public PDBAtom calcJumpCB(PDBResidue r,String code,double distance){
		PDBAtom ca = r.getCA();
		PDBAtom cb = r.getCB();
		
		if(cb == null){
			return null;
		}
		if(ca == null && cb != null){
			PDBAtom tcb = new PDBAtom();
			tcb.atom_code = "C";
			tcb.pdb_atom_code = code;
			tcb.loc.set(cb.loc);
			tcb.parent = r;
			return tcb;
		}
		
		PDBAtom tcb = new PDBAtom();
		tcb.atom_code = "C";
		tcb.pdb_atom_code = code;

		double xx = cb.loc.x-ca.loc.x;
		double yy = cb.loc.y-ca.loc.y;
		double zz = cb.loc.z-ca.loc.z;
		double dist = cb.distance(ca);
		xx /= dist/distance;
		yy /= dist/distance;
		zz /= dist/distance;
		xx += ca.loc.x;
		yy += ca.loc.y;
		zz += ca.loc.z;
		tcb.loc.set(xx,yy,zz);
		tcb.parent = r;
		return tcb;
	}
	
	
	
	public void prepare(Collection<PDBResidue> res){
		envAtoms.clear();
		jCB.clear();
		tCB.clear();
		t2CB.clear();
		
		for(int ii = 0;ii < 5;ii++){
			envAtoms.add(new ArrayList<PDBAtom>());
		}
		
		for(PDBResidue r:res){
			if(r.isLigand || r.isMissing()){
				continue;
			}
			r.removeAlt();
			
			
			PDBAtom ca = r.getCA();
			PDBAtom c = r.getC();
			PDBAtom n = r.getN();
			PDBAtom o = r.getO();
			if(ca !=  null){
				envAtoms.get(ENVATOM_CA).add(ca);
			}
			
			if(n !=  null){
				envAtoms.get(ENVATOM_N).add(n);
			}
			
			if(o !=  null){
				envAtoms.get(ENVATOM_O).add(o);
			}
			if(c !=  null){
				envAtoms.get(ENVATOM_C).add(c);
			}
			
			PDBAtom jcb = calcJumpCB(r,"jCB",jcbdist);
			if(jcb != null){
				jCB.put(r, jcb);
				envAtoms.get(ENVATOM_JCB).add(jcb);
			}
			PDBAtom tcb = calcJumpCB(r,"tCB",tcbdist);
			if(tcb != null){
				tCB.put(r, tcb);
			}
			PDBAtom tcb2 = calcJumpCB(r,"t2CB",tcb2dist);
			if(tcb2 != null){
				t2CB.put(r, tcb2);
			}
		}
	}
	
	/**
	 * ここはマニュアルで修正の必要がある。
	 * @return 
	 */
	public static ArrayList<String> getHeader(){
		
		ArrayList<String> headers = new ArrayList<>();
		headers.add("CA_count3");
		headers.add("CA_count6");
		headers.add("CA_count9");
		headers.add("JCB_count3");
		headers.add("JCB_count6");
		headers.add("JCB_count9");
		headers.add("JCB2_count3");
		headers.add("JCB2_count6");
		headers.add("JCB2_count9");
		//headers.add("CA_count3diff");
		//headers.add("CA_count6diff");
		//headers.add("CA_count9diff");
		String[] b = {"","","","",""};
		b[ENVATOM_CA] = "CA";
		b[ENVATOM_C] = "C";
		b[ENVATOM_N] = "N";
		b[ENVATOM_O] = "O";
		b[ENVATOM_JCB] = "JCB";
		for(String bb:b){
			for(int ii = 0;ii < 3;ii++){
				headers.add(bb+"_dist"+String.valueOf(ii));
				headers.add(bb+"_distdiff"+String.valueOf(ii));
			}
		}
		
		
		return headers;
	}
	
	
	public void updateResidue(PDBResidue r){
		
		PDBAtom jcb = calcJumpCB(r,"jCB",jcbdist);
		if(jcb != null){
			jCB.get(r).loc.set(jcb.loc);
		}
		PDBAtom tcb = calcJumpCB(r,"tCB",tcbdist);
		if(tcb != null){
			tCB.get(r).loc.set(tcb.loc);
		}
		PDBAtom tcb2 = calcJumpCB(r,"t2CB",tcb2dist);
		if(tcb2 != null){
			t2CB.get(r).loc.set(tcb2.loc);
		}
	}
	
	
	
	
	public ArrayList<PDBAtom> getTargetAtoms(PDBResidue target){
		ArrayList<PDBAtom> ret = new ArrayList<>();
		PDBAtom tcb1 = null; 
		if(tCB.containsKey(target)){
			tcb1 = tCB.get(target);
		}else{
			tcb1 = calcJumpCB(target,"tCB",tcbdist);
		}
		if(tcb1 != null){
			ret.add(tcb1);
		}
		return ret;
	}
	
	public double[] generateFeaturesA(PDBAtom target){
		return listToArray(generateFeatures(target));
	}
	/**
	 * ここはマニュアルで修正の必要がある。
	 * @return 
	 */
	public ArrayList<Double> generateFeatures(PDBAtom tcbb){
		ArrayList<Double> ret = new ArrayList<>();
		
		PDBResidue target = tcbb.parent;
		/*
		Fixme
		TCB を渡して親である残基オブジェクトをとるようにしているがあまりにも汚い。。。
		
		*/
		
		
		//3, 6, 12 A 内にある CA の数
		PDBAtom tcb1 = null; 
		PDBAtom tcb2 = null; 
		if(tCB.containsKey(target)){
			tcb1 = tCB.get(target);
		}else{
			tcb1 = calcJumpCB(target,"tCB",tcbdist);
		}
		if(t2CB.containsKey(target)){
			tcb2 = t2CB.get(target);
		}else{
			tcb2 = calcJumpCB(target,"t2CB",tcbdist);
		}
		if(tcb1 == null){
			return null;
		}
		int cacount[] = {0,0,0};
		
		for(PDBAtom a:envAtoms.get(ENVATOM_CA)){
			if(a.parent == target){
				continue;
			}
			double ddist = tcb1.distance(a);
			if(ddist < 9.0){
				cacount[2]++;
			}
			if(ddist < 6.0){
				cacount[1]++;
			}
			if(ddist < 3.0){
				cacount[0]++;
			}
		}
		ret.add((double)cacount[0]);
		ret.add((double)cacount[1]);
		ret.add((double)cacount[2]);
		
		
		
		int jacbcount[] = {0,0,0};
		
		for(PDBAtom a:envAtoms.get(ENVATOM_JCB)){
			if(a.parent == target){
				continue;
			}
			double ddist = tcb1.distance(a);
			if(ddist < 9.0){
				jacbcount[2]++;
			}
			if(ddist < 6.0){
				jacbcount[1]++;
			}
			if(ddist < 3.0){
				jacbcount[0]++;
			}
		}
		ret.add((double)jacbcount[0]);
		ret.add((double)jacbcount[1]);
		ret.add((double)jacbcount[2]);
		
		int[] jacbcount2 = {0,0,0};
		
		for(PDBAtom a:envAtoms.get(ENVATOM_JCB)){
			if(a.parent == target){
				continue;
			}
			double ddist = tcb2.distance(a);
			if(ddist < 9.0){
				jacbcount2[2]++;
			}
			if(ddist < 6.0){
				jacbcount2[1]++;
			}
			if(ddist < 3.0){
				jacbcount2[0]++;
			}
		}
		ret.add((double)jacbcount2[0]);
		ret.add((double)jacbcount2[1]);
		ret.add((double)jacbcount2[2]);
		
		
		
		
		//backbone について近い奴三つ
		for(int ii = 0;ii < 5;ii++){
			ArrayList<PDBAtom> al = envAtoms.get(ii);
			ArrayList<AtomDistance> dlis = new ArrayList<>();
			for(PDBAtom a:al){
				if(a.parent == target){
					continue;
				}
				double ddist = a.distance(tcb1);
				if(ddist < 12.0){
					double dd = a.distance(tcb2);
					dlis.add(new AtomDistance(a,ddist,dd-ddist));
				}
			}
			Collections.sort(dlis,new AtomDistanceComparator());
			for(int jj = 0;jj < 3;jj++){
				if(dlis.size() > jj){
					ret.add(dlis.get(jj).distance);
					ret.add(dlis.get(jj).distance2);
					
					
				}else{
					ret.add(1000.0);
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
	
	public static ArrayList<String> countAtoms_PDB(String treefile
			,File[] lis,File[] lis_test){
		RandomForestProcess rg = new RandomForestProcess();
		rg.loadForestFromFile_SimpleFormat(treefile);
		HashMap<String,Integer> allanswers = new HashMap<>();//回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		HashMap<String,HashMap<String,Integer>> pred_answers_test = new HashMap<>();//予測ラベルと予測された回答ラベルの数
		for(File f:lis){
			if(!Pattern.compile("\\.(pdb|ent)$").matcher(f.getPath()).find()){
				continue;
			}
			System.out.println(f.getPath());
			PDBData p = PDBData.loadPDBFile(f.getPath());
			FeatureGeneratorMergeCB_20180724 gen = new FeatureGeneratorMergeCB_20180724();
			for(String cc:p.chains.keySet()){
				PDBChain pc = p.chains.get(cc);
				gen.prepare(pc.residues);
				for(PDBResidue r:pc.residues){
					if(r.isLigand || r.isMissing() || r.getCB() == null){
						continue;
					}
					ArrayList<Double> res = gen.generateFeatures(r.getCB());
					if(res == null){
						continue;
					}
					String ans = r.getName()+"_"+gen.tcbdist+"CB";
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
		
		if(lis_test != null){
			for(File f:lis_test){
				if(!Pattern.compile("\\.(pdb|ent)$").matcher(f.getPath()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				FeatureGeneratorMergeCB_20180724 gen = new FeatureGeneratorMergeCB_20180724();
				for(String cc:p.chains.keySet()){
					PDBChain pc = p.chains.get(cc);
					gen.prepare(pc.residues);
					for(PDBResidue r:pc.residues){
						if(r.isLigand || r.isMissing() || r.getCB() == null){
							continue;
						}
						ArrayList<Double> res = gen.generateFeatures(r.getCB());
						if(res == null){
							continue;
						}
						String ans = r.getName()+"_"+gen.tcbdist+"CB";
						String pred = rg.predict_Majority(listToArray(res));

						//if(!allanswers.containsKey(ans)){
						//	allanswers.put(ans,0);
						//}
						if(!pred_answers_test.containsKey(pred)){
							pred_answers_test.put(pred,new HashMap<String,Integer>());
						}
						if(!pred_answers_test.get(pred).containsKey(ans)){
							pred_answers_test.get(pred).put(ans,0);
						}
						//allanswers.put(ans,1+allanswers.get(ans));
						pred_answers_test.get(pred).put(ans,1+pred_answers_test.get(pred).get(ans));
					}
				}
			}
			
		}else{
			pred_answers_test = pred_answers;
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
			HashMap<String,Integer> hs_test = pred_answers_test.get(s);
			
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
			if(hs_test != null){
				for(String hh:hs.keySet()){
					double backd = allanswers.get(hh)/(double)allcount;
					double predd = hs.get(hh)/(double)predcount;
					int hcou = 0;
					if(hs_test.get(hh) == null){
					}else{
						hcou = hs_test.get(hh);
					}
					if(predd/backd > Math.sqrt(2)){
						posicount+=hcou;
					}else if(predd/backd < 1.0/Math.sqrt(2)){
						negacount+=hcou;
					}else{
						zerocount+=hcou;
					}
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
	
	public static void makeScoreTable(String treefile,String outscore){
			//ArrayList<String> res = countAtoms(
			//"C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\reftree_noresinfo."+ii+".fin.dat"
			//,"C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\table_mergecb.dat");
			//for(String s:res){
			//	System.out.println(s);
			//}
			//printStrings("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_jumpcb_noresinfo\\0511\\scores.0.dat",res);


		String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\pdb_train";
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();

		String testdirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\pdb_test_1";
		File testdir = new File(testdirname);
		File[] lis_test= testdir.listFiles();

		ArrayList<String> dres = countAtoms_PDB(
		treefile
		,lis,lis_test);

		printStrings(outscore,dres);
		
	}
	public static void makeInputTable(String[] args){
		String outfilename = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\table_mergecb.dat";
		String sampledirname = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered_plusligand\\";
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
			new OutputStreamWriter(new FileOutputStream(outfilename,false),"UTF-8"))); 
			
			
			
			File dir = new File(sampledirname);
			File[] lis = dir.listFiles();
			
			ArrayList<String> header = FeatureGeneratorMergeCB_20180724.getHeader();
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
				FeatureGeneratorMergeCB_20180724 gen = new FeatureGeneratorMergeCB_20180724();
				for(String cc:p.chains.keySet()){
					PDBChain pc = p.chains.get(cc);
					gen.prepare(pc.residues);
					for(PDBResidue r:pc.residues){
						if(r.isLigand || r.isMissing() || r.getCB() == null){
							continue;
						}
						r.removeAlt();
						ArrayList<Double> res = gen.generateFeatures(r.getCB());
						if(res == null){
							continue;
						}
						pw.write(p.getID()+"_"+cc+"_"+r.getRepresentativeCode()+"_JCB\t");
						for(Double d:res){
							pw.write(String.valueOf((float)((double)d))+"\t");
						}
						pw.write(r.getName()+"_"+gen.tcbdist+"CB\n");
					}
				}
			}
			pw.close();
			
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public static void createScoringLine(String pdbdir, String treefile,String scorefile,String outfile){
		FeatureGeneratorMergeCB_20180724 gcb = new FeatureGeneratorMergeCB_20180724();
		FuzzyDecisionTreeScoring_generator fgen = new FuzzyDecisionTreeScoring_generator(gcb,treefile,scorefile);
		String sampledirname = pdbdir;
		File dir = new File(sampledirname);
		File[] lis = dir.listFiles();
		ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
		scoring.add(fgen);
		ArrayList<String> rstring = new ArrayList<>();
		for(File f:lis){
			PDBData pdb = PDBData.loadPDBFile(f.getAbsolutePath());
			ArrayList<PDBResidue> res = new ArrayList<>();
			
			for(String s:pdb.chains.keySet()){
				for(PDBResidue r:pdb.chains.get(s).residues){
					if(r.isValidAA()){
						r.removeAlt();
						res.add(r);
					}
				}
			}
			fgen.prepare(res);
			double[] d = FuzzyDecisionTreeScoring_generator.calcResidueScores(res,scoring);
			rstring.add("file:"+f.getAbsolutePath());
			for(int ii = 0;ii < d.length;ii++){
				rstring.add("r:"+res.get(ii).getRepresentativeCode()+"\t"+String.valueOf((float)d[ii]));
			}
		}
		printStrings(outfile,rstring);
	}
	
	
	public static void main(String[] args){
		//makeInputTable(args);
		for(int ii = 0;ii < 10;ii++){
			String treedir = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\fintree\\";
			String basedir = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_mergecb\\";
			File f = new File(treedir+"reftree.0."+ii+".fin.dat");
			if(!f.exists()){
				continue;
			}
			
			createScoringLine(
					basedir+"pdb_train",
					f.getAbsolutePath()
					,treedir+"reftree.0."+ii+".fin.dat.scores"
			,treedir+"train_pred_res."+ii+".dat");
		}
	}
}

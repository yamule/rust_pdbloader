/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static pepbuilderj.BackBoneSample.parseBlock;

/**
 *
 * @author kimidori
 */
public class AtomDistancePenalty {
	public static final int TYPE_ALL = 0;
	public static final int TYPE_BACKBONE = 1;
	public static final int TYPE_SIDECHAIN = 2;
	
	public static final int BACKBONE_N_CA = 0;
	public static final int BACKBONE_CA_C = 1;
	public static final int BACKBONE_C_O = 2;
	public static final int BACKBONE_N_C = 3;
	public static final int BACKBONE_PREVC_CA = 4;
	public static final int BACKBONE_PREVCA_N = 5;
	public static final int BACKBONE_PREVC_N = 6;
	
	public static final int MIN95 = 0;
	public static final int MAX95 = 1;
	public static final int MEAN5_95 = 2;
	
	
	
	double backbone_average[] = new double[7];
	double[] dist_prevc_nextn = new double[3];
	
	HashMap<String,Double> prev_current_backbone_min = new HashMap<>();
	HashMap<String,Double> prev_current_backbone_max = new HashMap<>();
	HashMap<String,Double> pairsMinDist = new HashMap<>();
	double trelance_ratio = 0.1;
	double trelance_ratio_peptideBonds = 0.1;
	
	double crashPenalty = -5;
	double peptideBondBrokenPenalty = -15000;
	
	
	
	HashMap<PDBAtom,HashMap<PDBAtom,Double>> backboneNextInteraction_min = new HashMap<>();
	HashMap<PDBAtom,HashMap<PDBAtom,Double>> backboneNextInteraction_max = new HashMap<>();
	HashMap<PDBAtom,ArrayList<PDBAtom>> peptideBonds = new HashMap<>();
	
	AtomDistancePenalty(){
		try{
			InputStream iss = AtomDistancePenalty.class.getResourceAsStream("resources/atom_mindist.dat");
			BufferedReader br = new BufferedReader(new InputStreamReader(iss));
			String ln = null;
			ArrayList<String> buff = new ArrayList<>();
			Pattern lpat = Pattern.compile("([A-Z]+)\\:([A-Z]+)");
			while((ln = br.readLine()) != null){
				String[] pt = ln.toUpperCase().replaceAll("[\r\n]","").split("\t");
				Matcher mat1 = lpat.matcher(pt[0]);
				Matcher mat2 = lpat.matcher(pt[1]);
				HashMap<String,String> mapp = MyWorld.lineToHash(ln);
				if(pt[0].equals("!DIST_PREVC_NEXTN")){
					dist_prevc_nextn[MIN95] = Double.parseDouble(mapp.get("min95"));
					dist_prevc_nextn[MAX95] = Double.parseDouble(mapp.get("max95"));
					dist_prevc_nextn[MEAN5_95] = Double.parseDouble(mapp.get("mean5_95"));
				}
				if(mat1.find() && mat2.find()){
					String m1 = mat1.group(1);
					String m1b = mat1.group(2);
					String m2 = mat2.group(1);
					String m2b = mat2.group(2);
					
					String m1c = m1+"_"+m1b;
					String m2c = m2+"_"+m2b;
					if((m1c.equals("CURRENT_N") && m2c.equals("CURRENT_CA"))
						||	(m1c.equals("CURRENT_CA") && m2c.equals("CURRENT_N"))){
						backbone_average[BACKBONE_N_CA] = Double.parseDouble(mapp.get("mean5_95"));
					}
					if((m1c.equals("CURRENT_CA") && m2c.equals("CURRENT_C"))
						||	(m1c.equals("CURRENT_C") && m2c.equals("CURRENT_CA"))){
						backbone_average[BACKBONE_CA_C] = Double.parseDouble(mapp.get("mean5_95"));
					}
					if((m1c.equals("CURRENT_C") && m2c.equals("CURRENT_O"))
						||	(m1c.equals("CURRENT_O") && m2c.equals("CURRENT_C"))){
						backbone_average[BACKBONE_C_O] = Double.parseDouble(mapp.get("mean5_95"));
					}
					if((m1c.equals("CURRENT_N") && m2c.equals("CURRENT_C"))
						||	(m1c.equals("CURRENT_C") && m2c.equals("CURRENT_N"))){
						backbone_average[BACKBONE_N_C] = Double.parseDouble(mapp.get("mean5_95"));
					}
					
					if((m1c.equals("PREV_C") && m2c.equals("CURRENT_CA"))
						||	(m1c.equals("CURRENT_CA") && m2c.equals("PREV_C"))){
						backbone_average[BACKBONE_PREVC_CA] = Double.parseDouble(mapp.get("mean5_95"));
					}
					if((m1c.equals("PREV_CA") && m2c.equals("CURRENT_N"))
						||	(m1c.equals("CURRENT_N") && m2c.equals("PREV_CA"))){
						backbone_average[BACKBONE_PREVCA_N] = Double.parseDouble(mapp.get("mean5_95"));
					}
					if((m1c.equals("PREV_C") && m2c.equals("CURRENT_N"))
						||	(m1c.equals("CURRENT_N") && m2c.equals("PREV_C"))){
						backbone_average[BACKBONE_PREVC_N] = Double.parseDouble(mapp.get("mean5_95"));
					}


					
					if(m1.equals("CURRENT")){
						if(m2.equals("PREV")){
							prev_current_backbone_min.put(mat2.group(2)+"_"+mat1.group(2),Double.parseDouble(mapp.get("min95")));
							prev_current_backbone_max.put(mat2.group(2)+"_"+mat1.group(2),Double.parseDouble(mapp.get("max95")));
						}else if(!m2.equals("CURRENT")){
							System.out.println("??? unexpected format "+ln);
						}
					}else if(m1.equals("PREV")){
						if(m2.equals("CURRENT")){
							prev_current_backbone_min.put(m1b+"_"+m2b,Double.parseDouble(mapp.get("min95")));
							prev_current_backbone_max.put(m1b+"_"+m2b,Double.parseDouble(mapp.get("max95")));
						}else{
							System.out.println("??? unexpected format "+ln);
						}
					}else{
						pairsMinDist.put(pt[0]+"_"+pt[1],Double.parseDouble(mapp.get("min95")));
						pairsMinDist.put(pt[1]+"_"+pt[0],Double.parseDouble(mapp.get("min95")));
					}
					
				}
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	public double dockingCrashPenalty(ArrayList<PDBResidue> allres
	,ArrayList<PDBResidue> allres2
	,double stopThreshold){
		double ret = 0;
		if(this.peptideBondBrokenPenalty > 0){
			this.peptideBondBrokenPenalty *= -1;
		}
		if(this.crashPenalty > 0){
			this.crashPenalty *= -1;
		}
		for(int jj = 0;jj < allres.size();jj++){
		
			PDBResidue p = allres.get(jj);
			jouter:for(int ii = 0;ii < allres2.size();ii++){
				PDBResidue aa = allres2.get(ii);
				if(aa.getCA().distance(p.getCA()) > 20){
					continue;
				}
				for(PDBAtom a:aa.atoms){
					String acode = a.parent.getName()+":"+a.pdb_atom_code;
					if(a.parent == p){
						continue;
					}
					for(PDBAtom r:p.atoms){
						double dist = r.distance(a);
						if(dist > 6.0){
							continue;
						}
						if(backboneNextInteraction_min.containsKey(r) && backboneNextInteraction_min.get(r).containsKey(a)){
								if(backboneNextInteraction_min.get(r).get(a) > (trelance_ratio+1.0)*dist){
								ret += crashPenalty;
							}
						}else{
							if(dist < 1.8){// CYS_CYS が 2.0 くらいなので
							//	ret += crashPenalty;
							//	continue;
							//}
							//String rcode = r.parent.getName()+":"+r.pdb_atom_code;
							//if(pairsMinDist.get(rcode+"_"+acode) > (trelance_ratio+1.0)*dist){
								ret += crashPenalty;
								if(ret < stopThreshold){
									return ret;
								}
								break jouter;
							}
						}
					}
				}
			}
		}
		return ret;
	}
	
	
	
	/**
	 * 今から置こうとしている Residue のペナルティを計算するため、allatoms はいちいち渡すようにする
	 * @param allatoms
	 * @param p
	 * @return 
	 */
	public double calcCrashPenalty(ArrayList<PDBResidue> allres,PDBResidue targetresidue,int type){
		double ret = 0;
		if(this.peptideBondBrokenPenalty > 0){
			this.peptideBondBrokenPenalty *= -1;
		}
		if(this.crashPenalty > 0){
			this.crashPenalty *= -1;
		}
		for(int ii = 0;ii < allres.size();ii++){
			PDBResidue aa = allres.get(ii);
			if(aa == targetresidue){
				continue;
			}
			if(aa.getCA().distance(targetresidue.getCA()) > 20){
				continue;
			}
			for(PDBAtom a:aa.atoms){
				String acode = a.parent.getName()+":"+a.pdb_atom_code;
				if(a.parent == targetresidue){
					continue;
				}
				for(PDBAtom r:targetresidue.atoms){
					if(type == TYPE_BACKBONE && !backboneNextInteraction_min.containsKey(r)){
						continue;
					}else if(type == TYPE_SIDECHAIN && backboneNextInteraction_min.containsKey(r)){
						continue;
					}
					double dist = r.distance(a);
					if(dist > 4.0){
						continue;
					}
					if(peptideBonds.containsKey(r)){
						ArrayList<PDBAtom> pz = peptideBonds.get(r);
						for(PDBAtom z:pz){
							double dz = r.distance(z);
							double sdist = prev_current_backbone_max.get("C_N");
							if(dz> sdist*(trelance_ratio+1.0)){
								ret += Math.max(peptideBondBrokenPenalty
										,peptideBondBrokenPenalty*(dz-sdist)/sdist);
							}
						}
					}
					if(backboneNextInteraction_min.containsKey(r) && backboneNextInteraction_min.get(r).containsKey(a)){

						if(backboneNextInteraction_min.get(r).get(a) > (trelance_ratio+1.0)*dist){
							ret += crashPenalty;
						}
					}else{
						
						if(dist < 1.8){// CYS_CYS が 2.0 くらいなので
							ret += crashPenalty;
							continue;
						}
						String rcode = r.parent.getName()+":"+r.pdb_atom_code;
						if(pairsMinDist.get(rcode+"_"+acode) > (trelance_ratio+1.0)*dist){
							ret += crashPenalty;
						}
					}
				}
			}
		}
		return ret;
	}
	
	
	
	public void prepare(ArrayList<PDBResidue> al){
		backboneNextInteraction_min.clear();
		for(int ii = 0;ii < al.size();ii++){
			PDBResidue rr = al.get(ii);
			if(!backboneNextInteraction_min.containsKey(rr.getC())){
				backboneNextInteraction_min.put(rr.getC(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getCA())){
				backboneNextInteraction_min.put(rr.getCA(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getN())){
				backboneNextInteraction_min.put(rr.getN(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getO())){
				backboneNextInteraction_min.put(rr.getO(),new HashMap<PDBAtom,Double>());
			}
			if(!peptideBonds.containsKey(rr.getN())){//HashMap である必要はないが。。。
				peptideBonds.put(rr.getN(),new ArrayList<PDBAtom>());
			}
			
			if(!peptideBonds.containsKey(rr.getC())){
				peptideBonds.put(rr.getC(),new ArrayList<PDBAtom>());
			}
			if(ii > 0){
				PDBResidue r1 = al.get(ii-1);
				peptideBonds.get(rr.getN()).add(r1.getC());
				
				backboneNextInteraction_min.get(rr.getC()).put(r1.getC(),prev_current_backbone_min.get("C_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getN(),prev_current_backbone_min.get("N_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getCA(),prev_current_backbone_min.get("CA_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getO(),prev_current_backbone_min.get("O_C"));
				
				
				
				
				backboneNextInteraction_min.get(rr.getN()).put(r1.getC(),prev_current_backbone_min.get("C_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getN(),prev_current_backbone_min.get("N_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getCA(),prev_current_backbone_min.get("CA_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getO(),prev_current_backbone_min.get("O_N"));
				
				
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getC(),prev_current_backbone_min.get("C_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getN(),prev_current_backbone_min.get("N_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getCA(),prev_current_backbone_min.get("CA_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getO(),prev_current_backbone_min.get("O_CA"));
				
				backboneNextInteraction_min.get(rr.getO()).put(r1.getC(),prev_current_backbone_min.get("C_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getN(),prev_current_backbone_min.get("N_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getCA(),prev_current_backbone_min.get("CA_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getO(),prev_current_backbone_min.get("O_O"));
				
				
				
				
				
				
			}
			
			if(ii < al.size()-1){
				PDBResidue r1 = al.get(ii+1);
				
				peptideBonds.get(rr.getC()).add(r1.getN());
				
				backboneNextInteraction_min.get(rr.getC()).put(r1.getC(),prev_current_backbone_min.get("C_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getN(),prev_current_backbone_min.get("C_N"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getCA(),prev_current_backbone_min.get("C_CA"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getO(),prev_current_backbone_min.get("C_O"));
				
				backboneNextInteraction_min.get(rr.getN()).put(r1.getC(),prev_current_backbone_min.get("N_C"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getN(),prev_current_backbone_min.get("N_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getCA(),prev_current_backbone_min.get("N_CA"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getO(),prev_current_backbone_min.get("N_O"));
				
				
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getC(),prev_current_backbone_min.get("CA_C"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getN(),prev_current_backbone_min.get("CA_N"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getCA(),prev_current_backbone_min.get("CA_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getO(),prev_current_backbone_min.get("CA_O"));
				
			}
		}
		
	}
	public void prepare_f(ArrayList<FloatingResidue> al){
		backboneNextInteraction_min.clear();
		for(int ii = 0;ii < al.size();ii++){
			FloatingResidue rr = al.get(ii);
			if(!backboneNextInteraction_min.containsKey(rr.getC())){
				backboneNextInteraction_min.put(rr.getC(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getCA())){
				backboneNextInteraction_min.put(rr.getCA(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getN())){
				backboneNextInteraction_min.put(rr.getN(),new HashMap<PDBAtom,Double>());
			}
			if(!backboneNextInteraction_min.containsKey(rr.getO())){
				backboneNextInteraction_min.put(rr.getO(),new HashMap<PDBAtom,Double>());
			}
			if(!peptideBonds.containsKey(rr.getN())){//HashMap である必要はないが。。。
				peptideBonds.put(rr.getN(),new ArrayList<PDBAtom>());
			}
			
			if(!peptideBonds.containsKey(rr.getC())){
				peptideBonds.put(rr.getC(),new ArrayList<PDBAtom>());
			}
			if(rr.prev != null){
				PDBResidue r1 =rr.prev;
				peptideBonds.get(rr.getN()).add(r1.getC());
				
				backboneNextInteraction_min.get(rr.getC()).put(r1.getC(),prev_current_backbone_min.get("C_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getN(),prev_current_backbone_min.get("N_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getCA(),prev_current_backbone_min.get("CA_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getO(),prev_current_backbone_min.get("O_C"));
				
				
				
				
				backboneNextInteraction_min.get(rr.getN()).put(r1.getC(),prev_current_backbone_min.get("C_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getN(),prev_current_backbone_min.get("N_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getCA(),prev_current_backbone_min.get("CA_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getO(),prev_current_backbone_min.get("O_N"));
				
				
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getC(),prev_current_backbone_min.get("C_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getN(),prev_current_backbone_min.get("N_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getCA(),prev_current_backbone_min.get("CA_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getO(),prev_current_backbone_min.get("O_CA"));
				
				backboneNextInteraction_min.get(rr.getO()).put(r1.getC(),prev_current_backbone_min.get("C_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getN(),prev_current_backbone_min.get("N_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getCA(),prev_current_backbone_min.get("CA_O"));
				backboneNextInteraction_min.get(rr.getO()).put(r1.getO(),prev_current_backbone_min.get("O_O"));
				
				
				
				
				
				
			}
			
			if(rr.next != null){
				PDBResidue r1 = rr.next;
				
				peptideBonds.get(rr.getC()).add(r1.getN());
				
				backboneNextInteraction_min.get(rr.getC()).put(r1.getC(),prev_current_backbone_min.get("C_C"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getN(),prev_current_backbone_min.get("C_N"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getCA(),prev_current_backbone_min.get("C_CA"));
				backboneNextInteraction_min.get(rr.getC()).put(r1.getO(),prev_current_backbone_min.get("C_O"));
				
				backboneNextInteraction_min.get(rr.getN()).put(r1.getC(),prev_current_backbone_min.get("N_C"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getN(),prev_current_backbone_min.get("N_N"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getCA(),prev_current_backbone_min.get("N_CA"));
				backboneNextInteraction_min.get(rr.getN()).put(r1.getO(),prev_current_backbone_min.get("N_O"));
				
				
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getC(),prev_current_backbone_min.get("CA_C"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getN(),prev_current_backbone_min.get("CA_N"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getCA(),prev_current_backbone_min.get("CA_CA"));
				backboneNextInteraction_min.get(rr.getCA()).put(r1.getO(),prev_current_backbone_min.get("CA_O"));
				
			}
		}
		
	}
	public static void main(String[] args){
		AtomDistancePenalty ap = new AtomDistancePenalty();
		System.out.println("test");
	}
}

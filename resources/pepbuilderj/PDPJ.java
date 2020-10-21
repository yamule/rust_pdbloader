/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author kimidori
 */
public class PDPJ {
	public static final double CONTACT_THRESHOLD = 8.0;
	public static ArrayList<ArrayList<PDBResidue>> performAllPDP(ArrayList<PDBResidue> al){
		 ArrayList<ArrayList<PDBResidue>> ret = new ArrayList<>();
		 HashMap<PDBResidue,HashMap<PDBResidue,ClosePair>> distpair = new HashMap<>();
		 
		for(int ii  = 0; ii < al.size();ii++){
			PDBResidue p1 = al.get(ii);
			if(p1.getCA() == null){
				continue;
			}
			distpair.put(p1,new HashMap<PDBResidue,ClosePair>());
			
			for(int jj = ii+1;jj < al.size();jj++){
				PDBResidue p2 = al.get(jj);
				if(p2.getCA() != null){
					double d = p1.getCA().distance(p2.getCA());
					if(d < CONTACT_THRESHOLD){
						distpair.get(p1).put(p2,new ClosePair(p1,p2,d));
					}
				}
			}
		}
		ArrayList<ArrayList<PDBResidue>> updated = new ArrayList<>();
		updated.add(al);
		ArrayList<ArrayList<PDBResidue>> candidate = new ArrayList<>();
		while(updated.size() > 0){
			ArrayList<PDBResidue> dom = updated.remove(0);
			ArrayList<ArrayList<PDBResidue>> rs = cutDomain(dom,distpair);
			if(rs == null){
				candidate.add(dom);
			}else{
				updated.addAll(rs);
			}
		}
		boolean chk = true;
		while(chk){
			double maxcontact = 2;
			int c1 = -1;
			int c2 = -1;
			chk = false;
			for(int ii = 0;ii < candidate.size();ii++){
				for(int jj = ii+1;jj < candidate.size();jj++){
					double d = calcInterDomainContact(candidate.get(ii),candidate.get(jj),distpair);
					if(d > maxcontact){
						c1 = ii;
						c2 = jj;
					}
				}
			}
			if(c1 > -1){
				ArrayList<PDBResidue> r = candidate.get(c1);
				ArrayList<PDBResidue> r2 = candidate.remove(c2);
				r.addAll(r2);
				chk = true;
			}
		}
		return candidate;
	}
	
	public static double calcInterDomainContact(ArrayList<PDBResidue> al
			,ArrayList<PDBResidue> al2
	,HashMap<PDBResidue,HashMap<PDBResidue,ClosePair>> distpair){
		int count2 = 0;
		for(PDBResidue p1:al){
			HashMap<PDBResidue,ClosePair> hs = null;
			if(!distpair.containsKey(p1)){
				continue;
			}
			hs = distpair.get(p1);
			for(PDBResidue p2:al2){
				if(!hs.containsKey(p2)){
					continue;
				}
				if(hs.get(p2).dist < CONTACT_THRESHOLD){
					count2++;
				}
			}
		}
		
		return count2/(Math.pow(al.size(),0.43)*Math.pow(al2.size(),0.43));
	}
	public static double calcInterDomainContact(ArrayList<PDBResidue> al,int start,int end
			,HashMap<PDBResidue,HashMap<PDBResidue,ClosePair>> distpair){
		int count1 = 0;
		int count2 = 0;
		
		int all1 = 0;
		int all2 = 0;
		start = Math.max(0,start);
		end = Math.min(end,al.size()-1);
		for(int ii  = 0; ii < al.size();ii++){
			PDBResidue p1 = al.get(ii);
			HashMap<PDBResidue,ClosePair> hs = null;
			if(!distpair.containsKey(p1)){
				continue;
			}
			hs = distpair.get(p1);
			for(int jj = ii+1;jj < al.size();jj++){
				PDBResidue p2 = al.get(jj);
				if(!hs.containsKey(p2)){
					continue;
				}
				if(hs.get(p2).dist < CONTACT_THRESHOLD){
					if(((ii-start)*(ii-end) <= 0 && (jj-start)*(jj-end) <= 0)
						|| ((ii-start)*(ii-end) > 0 && (jj-start)*(jj-end) > 0)
							){//intra domain
						count1++;
					}else{//inter domain
						count2++;
					}
				}
			}
		}
		return count2/(Math.pow(al.size()-end+start-1,0.43)*Math.pow(end-start+1,0.43));
	}
	
	
	public static ArrayList<ArrayList<PDBResidue>> cutDomain(ArrayList<PDBResidue> al
			,HashMap<PDBResidue,HashMap<PDBResidue,ClosePair>> distpair){
		 ArrayList<ArrayList<PDBResidue>> ret = new ArrayList<>();
		 int splitstart = -1;
		 int splitend = -1;
		 
		 int cutlen_min = 20;//BioJava が 12 より大きくないと切らないようにしていたが PDB は 20 っぽいのでそれに従う
		 
		 double normalizedcontact = 10000;
		 double ncsum = 0;
		 int count = 0;
		 for(int ii = cutlen_min+1;ii < al.size()-cutlen_min-1;ii++){
			 double nc = calcInterDomainContact(al,0,ii,distpair);
			 if(nc < normalizedcontact){
				splitstart = 0;
				splitend = ii;
				normalizedcontact = nc;
			 }
			ncsum+=nc;
			count++;
		 }
		 double dthreshold = ncsum/count*0.5;
		 for(int ii = cutlen_min+1;ii < al.size()-cutlen_min-1;ii++){
			 PDBResidue p1 = al.get(ii);
			 
			HashMap<PDBResidue,ClosePair> hs = null;
			if(!distpair.containsKey(p1)){
				continue;
			}
			hs = distpair.get(p1);
			 for(int jj = ii+1;jj < al.size()-cutlen_min-1;jj++){
				PDBResidue p2 = al.get(jj);
				if(!hs.containsKey(p2)){
					continue;
				}
				if(hs.get(p2).dist < CONTACT_THRESHOLD){
					if(p2.getResidueNumber()-p1.getResidueNumber()
					> 34 && jj-ii+1
							> 34){
						double nc = calcInterDomainContact(al,ii,jj,distpair);
						
						if(nc < normalizedcontact){
						   splitstart = ii;
						   splitend = jj;
						   normalizedcontact = nc;
						}
						if(ii > cutlen_min+2 && jj <  al.size()-cutlen_min-2){
							nc = calcInterDomainContact(al,ii-1,jj+1,distpair);

							if(nc < normalizedcontact){
							   splitstart = ii-1;
							   splitend = jj+1;
							   normalizedcontact = nc;
							}
						}
					}
				}
				
			 }
		 }
		// if(normalizedcontact < Math.pow(al.size()+splitstart-splitend-1,2/3.0) 
		//	 && normalizedcontact < Math.pow(splitend-splitstart+1,2/3.0)){
		if(normalizedcontact < dthreshold){
		//System.out.println(splitstart+";"+splitend);
			ArrayList<PDBResidue> a1 = new ArrayList<>();
			ArrayList<PDBResidue> a2 = new ArrayList<>();
			for(int ii = 0;ii < al.size();ii++){
				if((ii-splitstart)*(ii-splitend) <= 0){
					a1.add(al.get(ii));
				}else{
					a2.add(al.get(ii));
				}
			}
			ret.add(a1);
			ret.add(a2);
			return ret;
		 }
		 return null;
	}
	public static void main(String[] args){
		//途中
		PDBData p = PDBData.loadPDBFile("C:\\Users\\kimidori\\Downloads\\4up3.pdb");
		for(String c:p.chains.keySet()){
			PDBChain cc = p.chains.get(c);
			ArrayList<PDBResidue> rr = new ArrayList<>();
			for(PDBResidue r:cc.residues){
				if(!r.isLigand() && !r.isMissing()){
					rr.add(r);
				}
			}
			ArrayList<ArrayList<PDBResidue>> res = performAllPDP(rr);
			for(ArrayList<PDBResidue> r:res){
				System.out.println(r.get(0).getResidueNumber()+"-"+r.get(r.size()-1).getResidueNumber()+"-");
			}
			System.out.println(c+";"+res.size());
		}
	}
	
	
}

class ClosePair{
	PDBResidue r1;
	PDBResidue r2;
	double dist;
	ClosePair(PDBResidue rr1,PDBResidue rr2,double d){
		r1 = rr1;
		r2 = rr2;
		dist = d;
	}
}
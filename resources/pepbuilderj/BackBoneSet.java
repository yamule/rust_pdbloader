/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashMap;
import static pepbuilderj.PepProcess.phi;
import static pepbuilderj.PepProcess.psi;

/**
 *
 * @author kimidori
 */
public class BackBoneSet{
	HashMap<String,ArrayList<BackBoneSample>> backbones = new HashMap<>();
	String backboneDir = "resources/sampledresidues/";
	BackBoneSet(){
		String[] aaname = {
			"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
			,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
		};
		for(String a:aaname){
			backbones.put(a,BackBoneSample.load(
					BackBoneSample.class.getResourceAsStream(backboneDir+a+".backbones.dat")));
		}
	}
	
	public BackBoneSample getRandomly(String r){
		ArrayList<BackBoneSample> sal = backbones.get(r);
		int code = (int)(Math.random()*sal.size());
		return sal.get(code);
	}
	public BackBoneSample getRandomly(FloatingResidue r){
		
		ArrayList<BackBoneSample> sal = backbones.get(r.getName());
		int code = (int)(Math.random()*sal.size());
		r.backboneIndex = code;
		return sal.get(code);
		
	}
	public BackBoneSample getNext(FloatingResidue f){
		ArrayList<BackBoneSample> sal = backbones.get(f.getName());
		if(f.backboneIndex >= sal.size()-1){
			f.backboneIndex = -1;
		}
		f.backboneIndex++;
		return sal.get(f.backboneIndex);
	}
	
	public BackBoneSample getOneAt(FloatingResidue f,int index){
		f.backboneIndex = index-1;
		return getNext(f);
	}
	
	
	
	public double getProb(ArrayList<PDBResidue> al,int ii){
		ArrayList<BackBoneSample> sal = backbones.get(al.get(ii).getName());
		
		if(ii == 0 || ii == al.size()-1){
			return 0;
		}
		
		
		double mdist =Double.MAX_VALUE;
		BackBoneSample minsample = null;
		double x = phi(al.get(ii),al.get(ii-1))+180;
		double y = psi(al.get(ii),al.get(ii+1))+180;
		
		double distx = 360;
		double disty = 360;
		for(BackBoneSample b:sal){
			double xd = Math.abs(x-b.phi);
			if(xd <= distx && xd < 20 ){
				distx = xd;
			}else{
				continue;
			}
			
			double yd = Math.abs(y-b.psi);
			if(yd <= disty && yd < 20 ){
				disty = yd;
			}else{
				continue;
			}
			
			
			minsample = b;
		}
		if(minsample == null){
			return 0;
		}
		return minsample.prob;
	}
	
	public static void main(String[] args){
		BackBoneSet bb = new BackBoneSet();
	}
}

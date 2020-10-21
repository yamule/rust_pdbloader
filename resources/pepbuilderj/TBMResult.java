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
//TBM の結果
//テンプレート側がギャップの場合 mapped に false が入っている。

public class TBMResult{
	ArrayList<PDBResidue> residues = new ArrayList<>();
	ArrayList<Boolean> mapped = new ArrayList<>();
	ArrayList<PDBResidue> template_residues = new ArrayList<>();
	HashMap<PDBResidue,PDBResidue> target_template_map = new HashMap<>();
	HashMap<PDBResidue,PDBResidue> template_target_map = new HashMap<>();
	double score = 0;
	TBMResult(){
	}
	public void setTemplateResidues(ArrayList<PDBResidue> c){
		for(PDBResidue p:c){
			if(p == null){
				continue;
			}
			if(p.isMissing() || p.isLigand()){
			}else{
				template_residues.add(p);
			}
		}
	}
	public void addResidue(PDBResidue target,PDBResidue template,boolean flag){
		residues.add(target);
		if(template != null){
			target_template_map.put(target,template);
			template_target_map.put(template, target);
			
		}
		mapped.add(flag);
	}

}

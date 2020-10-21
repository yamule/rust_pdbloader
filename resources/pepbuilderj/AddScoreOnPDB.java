/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import static pepbuilderj.TemplateBaseModeller.makeAtomLines;
import static pepbuilderj.TemplateBaseModeller.writeToFile;

/**
 *
 * @author kimidori
 */
public class AddScoreOnPDB{
	ArrayList<FuzzyDecisionTreeScoring_generator> scoring = new ArrayList<>();
	AddScoreOnPDB(){
	FuzzyDecisionTreeScoring_generator j 
				= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorAtom_20180503());
		scoring.add(j);
	}
	public void addScoreOnBFactor(ArrayList<PDBResidue> res){
		for(FuzzyDecisionTreeScoring_generator g:scoring){
			g.prepare(res);
		}
		double[] scores = FuzzyDecisionTreeScoring_generator.calcResidueScores(res, scoring);
		
		double min = scores[0];
		double max = scores[0];
		
		double scoresum = 0;
		for(int ii = 0;ii < res.size();ii++){
			if(min == FuzzyDecisionTreeScoring_generator.SCORE_INVALIDRESIDUE){
				min = scores[ii];
			}
			if(max == FuzzyDecisionTreeScoring_generator.SCORE_INVALIDRESIDUE){
				max = scores[ii];
			}
			scoresum += scores[ii];
			min = Math.min(scores[ii],min);
			max = Math.max(scores[ii],max);
		}
		for(int ii = 0;ii < res.size();ii++){
			PDBResidue r = res.get(ii);
			for(PDBAtom a:r.atoms){
				if(FuzzyDecisionTreeScoring_generator.SCORE_INVALIDRESIDUE
						== scores[ii]){
					scores[ii] = min;
				}
				//(8pi^2*d^2)/3
				double bscore = 1.0-scores[ii]/max;
				a.setBFactor(8*Math.pow(Math.PI,2.0)*Math.pow(bscore,2));
			}
		}
		System.out.println("score\t"+scoresum+"\t"+res.size());
	}
	public static void main(String[] args){
		PDBData d = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\bioo\\docking_try\\3KB5.pdb");
		ArrayList<PDBResidue> ress = new ArrayList<>();
		for(String c:d.chains.keySet()){
			ArrayList<PDBResidue> r = d.chains.get(c).residues;
			for(PDBResidue rr:r){
				if(rr.isLigand() || rr.isMissing()){
					continue;
				}
				ress.add(rr);
			}
		}
		AddScoreOnPDB aop = new AddScoreOnPDB();
		aop.addScoreOnBFactor(ress);
	}
}

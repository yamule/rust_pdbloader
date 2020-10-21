/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.Collection;

/**
 *
 * @author kimidori
 */
public interface FeatureGenerator {
	public void prepare(Collection<PDBResidue> res);
	public ArrayList<Double> generateFeatures(PDBAtom target);
	public double[] generateFeaturesA(PDBAtom target);
	
	public ArrayList<Double> generateFeatures_merge(PDBAtom target,PDBAtom target2);
	public ArrayList<Double> mergeFeatures(ArrayList<Double> target,ArrayList<Double> target2);
	public double[] generateFeaturesA_merge(PDBAtom target,PDBAtom target2);
	
	public ArrayList<PDBAtom> getTargetAtoms(PDBResidue r);
	public ArrayList<String> getResourceFilePath();
	/**
	 * 正解ラベル、予測ラベルを与えると突き合わせの文字列を返す
	 * @param b
	 * @return 
	 */
	public String transCode(String b);
}

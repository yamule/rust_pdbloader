/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.HashMap;

/**
 *
 * @author yamule
 */
public class PepBuilderJ {
	public static HashMap<String,String> parseArg(String[] args){
		HashMap<String,String> ret = new HashMap<>();
		for(int ii = 0;ii < args.length;ii++){
			if(args[ii].indexOf("-") == 0){
				if(args.length > ii+1){
					ret.put(args[ii],args[ii+1]);
				}else{
					ret.put(args[ii],"");
				}
			}
		}
		return ret;
	}
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		//	なんかすごいループができるのでチェック
		//String[] test = {"-model","-in"
		//		,"C:\\dummy\\vbox_share\\casp13\\queries\\T0957\\s1\\test11pfas.fas"
		//,"-out","C:\\dummy\\vbox_share\\casp13\\queries\\test100.fas.pdb"};
		
		
		//String[] test = {"-model","-in"
		//		,"C:\\dummy\\vbox_share\\casp13\\queries\\T0958\\model1.pfas.fas"
		//,"-out","C:\\dummy\\vbox_share\\casp13\\queries\\T0958\\model1.pdb"};
		//String[] test = "-model -in C:\\dummy\\vbox_share\\casp13\\queries\\realigntest\\yoke\\mul.8.2ihy.model_1.chain_B.pdb.align.fas -out testout.pdb".split("[\\s]+");
		//args = test;
		
		HashMap<String,String> arghash = parseArg(args);
		
		if(arghash.containsKey("-model")){
			int realignnum = 0;
			int realignwindow = 1;
			if(arghash.containsKey("-realign")){
				realignnum = Integer.parseInt(arghash.get("-realign"));
			}
			if(arghash.containsKey("-realign_window")){
				realignwindow = Integer.parseInt(arghash.get("-realign_window"));
			}
			int rotemerrefine = 0;
			if(arghash.containsKey("-check_rotemer")){
				rotemerrefine = Integer.parseInt(arghash.get("-check_rotemer"));
			}
			TemplateBaseModeller.model(arghash.get("-in"),arghash.get("-out"),arghash.containsKey("-nocorrupt")
					,realignnum,realignwindow,rotemerrefine);
			
			
			
			
		}else if(arghash.containsKey("-score")){
			FuzzyDecisionTreeScoring_generator fs = new FuzzyDecisionTreeScoring_generator(
			new FeatureGeneratorAtom_20180503());
			double[] res = fs.calcPDBScore(arghash.get("-in"));
			for(int ii = 0;ii < res.length;ii++){
				System.out.print(res[ii]+"\t");
			}
		}else if(arghash.containsKey("-scoresum_cb")){
			FuzzyDecisionTreeScoring_generator fs = new FuzzyDecisionTreeScoring_generator(
			new FeatureGeneratorCB_20180501());
			double[] res = fs.calcPDBScore(arghash.get("-in"));
			double ss = 0;
			for(int ii = 0;ii < res.length;ii++){
				ss += res[ii];
			}
			System.out.println(arghash.get("-in")+"\t"+ss);
		}else if(arghash.containsKey("-scoresum")){
			FuzzyDecisionTreeScoring_generator fs = new FuzzyDecisionTreeScoring_generator(
			new FeatureGeneratorAtom_20180503());
			double[] res = fs.calcPDBScore(arghash.get("-in"));
			double ss = 0;
			for(int ii = 0;ii < res.length;ii++){
				ss += res[ii];
			}
			System.out.println(arghash.get("-in")+"\t"+ss);
		}else if(arghash.containsKey("-refinement")){
			ResidueRefine.main(args);
		}
		
		//else{
		//	LatticeWorld.main(args);
		//}
	}
	
}

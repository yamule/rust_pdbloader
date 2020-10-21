/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.File;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class CalcContacts {
	public static void main(String[] args){
		//モノマーの internal コンタクト数を数える
		String dirname = "C:\\dummy\\vbox_share\\casp13\\queries\\T1015s1\\threadingresult";
		File dir = new File(dirname);
		if(dir.isDirectory()){
			File[] flist = dir.listFiles();
			for(File f:flist){
				String filename = f.getAbsolutePath();
				if(Pattern.compile("pdb$").matcher(filename).find()){
					PDBData d = PDBData.loadPDBFile(filename);
					int count =0;
					int crash = 0;
					for(String s:d.chains.keySet()){
						PDBChain c = d.chains.get(s);
						for(PDBResidue r:c.residues){
							if(r.isLigand() || r.getCA() == null){
								continue;
							}
							for(PDBResidue r2:c.residues){
								if(r2.isLigand() || r2.getC() == null){
									continue;
								}
								if(r.getResidueNumber() >= r2.getResidueNumber()){
									continue;
								}
								if(r == r2){
									continue;
								}
								if(Math.abs(r.getResidueNumber()-r2.getResidueNumber()) > 10 && r.getCA().distance(r2.getCA()) < 8.0){
									count++;
								}
								if(Math.abs(r.getResidueNumber()-r2.getResidueNumber()) > 5 && r.getCA().distance(r2.getCA()) < 3.0){
									crash++;
								}
								
							}
						}
					}
					System.out.println(filename+"\t"+count+"\t"+crash);
				}
			}
		}
	}
}

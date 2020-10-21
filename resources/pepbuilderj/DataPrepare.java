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
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 *
 * @author kimidori
 */
public class DataPrepare {
	/**
	 * 
	 */
	public static void MakeResidueFragments(){
		String outdir_ = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\residue_fragments";
		String indir_ = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\onechain_plain_filtered";
		String resfilename = "C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\distout.dat";
		File outdir = new File(outdir_);
		if(!outdir.exists()){
			outdir.mkdir();
		}
		File indir =new File(indir_);
		File[] lis = indir.listFiles();
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(
					new OutputStreamWriter(new FileOutputStream(resfilename,false),"UTF-8"))); 

			//td <- read.table("distout.dat",header=F);
			//pca = prcomp(td, scale=T)
			//plot(pca$x[,"PC1"],pca$x[,"PC2"])
			Pattern pat = Pattern.compile("\\.pdb$");
			HashMap<String,Integer> code_count = new HashMap<>();
			for(File f:lis){
				if(!pat.matcher(f.getName()).find()){
					continue;
				}
				System.out.println(f.getPath());
				PDBData p = PDBData.loadPDBFile(f.getPath());
				for(String s:p.chains.keySet()){
					PDBChain c = p.chains.get(s);
					ArrayList<PDBResidue> rr = new ArrayList<>(c.residues);
					boolean[] cbreak = PepProcess.checkChainBreak(rr);
					for(int ii = 1;ii < rr.size()-1;ii++){
						if(!cbreak[ii] && !cbreak[ii-1]){
							PDBResidue r = rr.get(ii);
							if(!r.getName().equals("ALA")){
								continue;
							}
							PDBAtom prevc = rr.get(ii-1).getC();
							PDBAtom nexn = rr.get(ii+1).getN();
							ArrayList<PDBAtom> atoms = new ArrayList<>();
							atoms.add(prevc);
							atoms.add(r.getN());
							atoms.add(r.getCA());
							atoms.add(r.getC());
							atoms.add(nexn);
							for(int kk = 0;kk < atoms.size();kk++){
								if(atoms.get(kk) == null){
									System.out.println("");
								}
							}
							double[] angles = new double[2];
							for(int kk = 0;kk < 2;kk++){
								angles[kk] = PepProcess.dihedralAngle(
										atoms.get(kk).loc
										, atoms.get(kk+1).loc
										, atoms.get(kk+2).loc
										, atoms.get(kk+3).loc
								);
							}
							
							
							/*
							pw.write(angles[0]+"\t");
							pw.write(angles[1]+"\t");
							for(int jj = 0;jj < atoms.size();jj++){
								for(int kk = jj+1;kk < atoms.size();kk++){
									double dd = atoms.get(jj).distance(atoms.get(kk));
									pw.write(dd+"\t");
								}
							}
							pw.write("\n");
							*/
							
							StringBuffer rcode = new StringBuffer();
							rcode.append((int)((angles[0]+180)/20)+"\t"+(int)((angles[1]+180)/20));
							//for(int jj = 0;jj < atoms.size();jj++){
							//	for(int kk = jj+1;kk < atoms.size();kk++){
							//		double dd = atoms.get(jj).distance(atoms.get(kk));
							//		rcode.append("\t"+(int)(dd*10));
							//	}
							//}
							rcode.append("\t"+String.valueOf((int)(prevc.distance(nexn)*4)));
							String rc = rcode.toString();
							if(!code_count.containsKey(rc)){
								code_count.put(rc,0);
							}
							code_count.put(rc,code_count.get(rc)+1);
						}
					}
				}
			}
			for(String s:code_count.keySet()){
				if(code_count.get(s) > 1){
					pw.write(s+"\t"+code_count.get(s)+"\n");
				}
			}
			
			
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
		
		
	}
	public static void main(String[] args){
		MakeResidueFragments();
	}
	
	
	
}

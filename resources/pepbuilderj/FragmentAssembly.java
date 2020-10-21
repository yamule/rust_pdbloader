/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 *
 * @author kimidori
 */
public class FragmentAssembly {
	
	public static ArrayList<ArrayList<PDBResidue>> getUniqueResidue(ArrayList<String> filelist){
		HashSet<Integer> resused = new HashSet<>();
		ArrayList<ArrayList<PDBResidue>> ret = new ArrayList<>();
		for(String ss:filelist){
			PDBData d = PDBData.loadPDBFile(ss);
			for(String c:d.chains.keySet()){
				ArrayList<PDBResidue> res = new ArrayList<>();
				for(PDBResidue r:d.chains.get(c).residues){
					if(resused.contains(r.getResidueNumber())){
						if(res.size() > 5){
							for(PDBResidue rr:res){
								resused.add(rr.getResidueNumber());
							}
							ret.add(res);
						}
						res = new ArrayList<>();
					}else{
						
						res.add(r);
					}
				}
				if(res.size() > 5){
					ret.add(res);
					for(PDBResidue r:res){
						resused.add(r.getResidueNumber());
					}
				}
			}
		}
		return ret;
	}
	public static void main(String[] args){
		ChainBuilder cb = new ChainBuilder();
		ArrayList<String> st = new ArrayList<>();
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.0.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.1.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.2.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.3.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.4.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.5.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.6.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.7.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.8.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.9.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.10.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.11.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.12.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.13.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.14.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.15.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.16.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.17.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.18.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.19.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.20.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.21.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.22.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.23.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.24.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.25.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.26.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.27.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.28.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.29.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.30.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.31.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.32.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.33.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.34.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.35.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.36.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.37.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.38.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.39.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.40.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.41.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.42.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.43.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.44.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.45.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.46.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.47.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.48.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.49.th.pdb");
st.add("C:\\dummy\\vbox_share\\casp13\\queries\\T0991\\threadingresult\\T0991.50.th.pdb");

		for(int kk = 0;kk < 200;kk++){
			if(st.size() == 0){
				break;
			}
			if(kk != 0){
				for(int ii = 0;ii < st.size();ii++){
					int p = (int)(Math.random()*st.size());
					int q = (int)(Math.random()*st.size());
					String pp = st.get(p);
					String qq = st.get(q);
					st.set(p,qq);
					st.set(q,pp);
				}
				
			}
			ArrayList<VSorter> sorter = new ArrayList<>();
			ArrayList<ArrayList<PDBResidue>> blis =  getUniqueResidue(st);
			for(int ii = 0;ii < blis.size();ii++){
				ArrayList<PDBResidue> p = blis.get(ii);
				sorter.add(new VSorter(p.get(0).getResidueNumber(),ii));
			}
			Collections.sort(sorter,new VComparator());
			ArrayList<ArrayList<FloatingResidue>> llis = new ArrayList<>();
			for(VSorter vv:sorter){
				ArrayList<FloatingResidue> fl = ChainBuilder.changeToFloating(blis.get(vv.index));
				cb.prepare(fl);
				for(FloatingResidue r:fl){
					cb.refiner.maxRotamer(r, -5,0.1, true, cb.sidechains, false);
				}
				llis.add(fl);
			}
			ArrayList<FloatingResidue> fr = cb.sequencialDock(llis);

			ChainBuilder.saveChains(fr,"teststr."+kk+".pdb");
		}
	}
}

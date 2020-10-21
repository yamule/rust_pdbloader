/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

/**
 *
 * @author kimidori
 */
import java.util.*;
import java.util.regex.*;
import java.io.*;
import static pepbuilderj.ChainBuilder.saveOneChain;
import static pepbuilderj.FuzzyTreeAlign.blosum;
import static pepbuilderj.Refine3.changeBackBone;
import static pepbuilderj.Refine3.getUnmappedRegions;



public class TemplateBaseModeller{
	
	public static final double BONDLENGTH_STATIC = 1.8;
	public static HashMap<String,String> aaCode = new HashMap<>();
		static{
		aaCode.put("ALA","A");
		aaCode.put("ARG","R");
		aaCode.put("ASN","N");
		aaCode.put("ASP","D");
		aaCode.put("CYS","C");
		aaCode.put("GLN","Q");
		aaCode.put("GLU","E");
		aaCode.put("GLY","G");
		aaCode.put("HIS","H");
		aaCode.put("ILE","I");
		aaCode.put("LEU","L");
		aaCode.put("LYS","K");
		aaCode.put("MET","M");
		aaCode.put("PHE","F");
		aaCode.put("PRO","P");
		aaCode.put("SER","S");
		aaCode.put("THR","T");
		aaCode.put("TRP","W");
		aaCode.put("TYR","Y");
		aaCode.put("VAL","V");
		
		//aaCode.put("SEC","U");
		//aaCode.put("UNK","X");
		
		
	}
	
	HashMap<String,ArrayList<PDBResidue>> sources = new HashMap<>();//マップの素になる残基の座標
	ResidueRefine refiner = new ResidueRefine();
	
	public TemplateBaseModeller(){
		PDBData source = DefaultAtomCoordinate.asPDBData();
		for(String s:source.chains.keySet()){
			for(PDBResidue res:source.chains.get(s).residues){
				if(res.isValidAA()){
					if(sources.get(res.getName()) == null){
						sources.put(res.getName(),new ArrayList<PDBResidue>());
						sources.put(PDBResidue.aaMap.get(res.getName()),new ArrayList<PDBResidue>());
					}
					sources.get(res.getName()).add(res);
					sources.get(PDBResidue.aaMap.get(res.getName())).add(res);
				}
			}
		}
	}
	
	/**
	 * TBM に認識される構造の Fasta 配列を返す。
	 * @return 
	 */
	
	
	public static ArrayList<String> getFastaSequenceForTBM(PDBData pdb){
		ArrayList<String> ret = new ArrayList<>();
		for(String c:pdb.chains.keySet()){
			PDBChain cc = pdb.chains.get(c);
			ArrayList<PDBResidue> r = prepareTemplateResidues(cc.residues);
			StringBuffer sb = new StringBuffer();
			sb.append(">"+pdb.getID()+":"+c+"\n");
			int cou = 0;
			for(PDBResidue rr:r){
				sb.append(aaCode.get(rr.getName()));
				cou++;
				if(cou%60 == 0){
					sb.append("\n");
				}
			}
			sb.append("\n");
			if(cou > 0){
				ret.add(sb.toString());
			}else{
				System.err.println("Chain "+c+" does not have valid aa.");
			}
		}
		return ret;
		
	}
	
	
	/**
	 * remove ligand and residue which do not have OCNCA
	 * @return 
	 */
	public static ArrayList<PDBResidue> prepareTemplateResidues(ArrayList<PDBResidue> al){
		ArrayList<PDBResidue> ret = new ArrayList<>();
		for(PDBResidue rr:al){
			if(rr.isMissing() || rr.isLigand()){
				continue;
			}
			if(rr.getC() == null || rr.getCA() == null || rr.getO() == null || rr.getN() == null){
				continue;
			}
			//分からない奴は全部 UNKNOWN にしてしまう
			if(!aaCode.containsKey(rr.getName())){
				Iterator<PDBAtom> aite = rr.atoms.iterator();
				boolean hascb = false;
				//その他の残基に変えたい場合ここ変更
				System.err.println(rr.getRepresentativeCode()+" was changed to ALA.");
				rr.setName("ALA");
				while(aite.hasNext()){
					PDBAtom a = aite.next();
					if(a.pdb_atom_code.equals("CA")
						|| a.pdb_atom_code.equals("C")
						|| a.pdb_atom_code.equals("N")
						|| a.pdb_atom_code.equals("O")
						|| a.pdb_atom_code.indexOf("CB") == 0){
						if(a.pdb_atom_code.indexOf("CB") == 0){
							hascb = true;
						}
					}else{
						aite.remove();
					}
				}
				if(!hascb){
					System.err.println("Pseudo CB was added. the atom number cannot used without sorting.");
					rr.addAtom(rr.getCB());
				}
				ret.add(rr);
			}else{
				ret.add(rr);
			}
		}
		return PepProcess.makeFilteredAA(ret,true);
	}
	
	
	/**
	 * Missing でもなく Valid な AA のみを抽出して返す
	 */
	public static ArrayList<PDBResidue> getValidAAResidues(ArrayList<PDBResidue> res){
		ArrayList<PDBResidue> ret = new ArrayList<>();
		
		int check = -1;
		for(PDBResidue r:res){
			if(r.isValidAA() && !r.isMissing()){
				ret.add(r);
			}else{
				if(!r.isValidAA()){
					if(ret.size() > 0){
						check = ret.size();//途中に HETATM があった場合等
					}
				}
			}
		}
		if(check > -1 && ret.size() != check){
			System.err.println("Template has invalid residues in its inside.");
			return null;
		}
		return ret;
	}
	
	
	
	
	
	public TBMResult model(String src,String temp,PDBChain template){
		String[] ss = src.replaceAll("[\\s]","").replaceAll("\\.","-").split("");
		String[] tt = temp.replaceAll("[\\s]","").replaceAll("\\.","-").split("");
		ArrayList<String> sa = new ArrayList<>();
		ArrayList<String> ta = new ArrayList<>();
		for(String s:ss){
			if(s.length() > 0){
				sa.add(s);
			}
		}
		for(String s:tt){
			if(s.length() > 0){
				ta.add(s);
			}
		}
		
		return model(sa,ta,template);
	}
	
	
	/**
	 * Alignment のデータ Temp を与え、Align されてない Residue については必要ないのでテンプレートから削除する
	 * 
	 * 
	 * @param temp
	 * @param template 
	 */
	public static void filtAlignedResidues(ArrayList<String> temp,ArrayList<PDBResidue> pos_temp){

		StringBuffer resinchain = new StringBuffer();
		StringBuffer resinalignment = new StringBuffer();
		HashMap<Integer,Integer> aapos_posinalignment = new HashMap<>();
		int cou = 0;
		for(int ii = 0;ii < temp.size();ii++){
			String t = temp.get(ii);
			if(t.equals("-") || t.equals(".")){
			}else{
				resinalignment.append(t.toUpperCase());
				aapos_posinalignment.put(cou,ii);
				cou++;
			}
		}
		Iterator<PDBResidue> rite = pos_temp.iterator();
		while(rite.hasNext()){
			PDBResidue r = rite.next();
			if(r.isLigand() || r.isMissing()){
				rite.remove();
			}
		}
		for(PDBResidue r:pos_temp){
			if(aaCode.containsKey(r.getName())){
				resinchain.append(aaCode.get(r.getName()));
			}else{
				resinchain.append("X");
			}
		}
		ArrayList<PDBResidue> ret = new ArrayList<>();
		if(!resinalignment.equals(resinchain)){
			//アラインメント中の文字と PDB ファイルの中の文字が違う場合 SW を試みる
			SmithWaterman sw = new SmithWaterman();
			sw.smat = new BLOSUM62_X0();
			SWResult res = sw.align(resinalignment.toString(), resinchain.toString());
			//System.out.println(SmithWaterman.listToString(res.qseq));
			//System.out.println(SmithWaterman.listToString(res.sseq));
			int rlen = res.qseq.size();
			int qcou = 0;
			int scou = 0;
			for(int ii = 0;ii < rlen;ii++){
				if(res.qseq.get(ii).equals('-')){//アラインメント側でギャップになっているとき PDB 側は null を入れてあとで削除
					pos_temp.set(scou,null);
				}
				if(!res.qseq.get(ii).equals('-')
						&& !res.sseq.get(ii).equals('-')){
					ret.add(pos_temp.get(scou));
				}
				
				if(res.sseq.get(ii).equals('-')){//PDB 側でギャップになっているときアラインメント側はギャップに
					temp.set(aapos_posinalignment.get(qcou),"-");
				}
				if(!res.qseq.get(ii).equals('-')){
					qcou++;
				}
				if(!res.sseq.get(ii).equals('-')){
					scou++;
				}
			}
		}
		Iterator<PDBResidue> ite = pos_temp.iterator();
		while(ite.hasNext()){//null が入っているものは削除する
			PDBResidue p = ite.next();
			if(p == null){
				ite.remove();
			}
		}
		//for(int ii = 0;ii < pos_temp.size();ii++){
		//	System.out.println(pos_temp.get(ii).getName()+";;"+ret.get(ii).getName());
		//}
		
	}
	
	
	/**
	 * アラインメントされている Backbone と CB だけ合わせる。
	 * @param src
	 * @param tt
	 * @param template
	 * @param fix
	 * @return 
	 */
	
	public TBMResult model(ArrayList<String> src,ArrayList<String> tt,PDBChain template){
		return model(src,tt,template.residues);
	}
	public TBMResult model(ArrayList<String> src,ArrayList<String> tt,ArrayList<PDBResidue> ress){
		ArrayList<ArrayList<PDBResidue>> ret = new ArrayList<>();
		ArrayList<PDBResidue> pos_temp = prepareTemplateResidues(ress);
		ArrayList<String> temp = new ArrayList<>(tt);
		
		filtAlignedResidues(temp,pos_temp);
		
		
		
		
		
				
		ArrayList<String> src_aa_str = new ArrayList<>();
		ArrayList<PDBResidue> pr = new ArrayList<>();
		int tcou = 0;
		for(int ii = 0;ii < temp.size();ii++){
			if(!src.get(ii).equals("-")){
				src_aa_str.add(src.get(ii));
				if(!temp.get(ii).equals("-")){
					pr.add(pos_temp.get(tcou));
					if(!temp.get(ii).equals(PDBResidue.aaMap.get(pos_temp.get(tcou).getName()))){
						System.err.println(ii+" alignment discrepancy \t"+temp.get(ii)
								+"\t"+pos_temp.get(tcou).getName());
					}
				}else{
					pr.add(null);//インサーションになる部分は null
				}
			}
			
			if(!temp.get(ii).equals("-")){
					tcou++;
			}
		}
		TBMResult res = new TBMResult();
		res.setTemplateResidues(pos_temp);
		for(int ii = 0;ii < src_aa_str.size();ii++){
			String code = src_aa_str.get(ii);
			if(sources.get(code) != null){
				//PDBResidue rr = sources.get(code).get((int)(sources.get(code).size()*Math.random()));
				PDBResidue rr = sources.get(code).get(0);
				PDBResidue p = pr.get(ii);
				if(p == null){//ターゲットに余分残基があるので適当に置く
					res.addResidue(rr.getCopy(),null,false);
				}else{
					PDBResidue c  = rr.getCopy();
					placeResidueTo(c,p);
					res.addResidue(c,p,true);
				}
			}else{
				System.err.println(code+" is not defined in source material file.");
			}
		}
		for(int ii = 0;ii < res.residues.size();ii++){
			res.residues.get(ii).setResidueNumber(ii);
		}
		return res;
	}
	
	
	/**
	 * あるベクトルの Y 軸中心の回転 sin, cos, Z 軸中心の回転 sin, cos を得る
	 * x 軸方向のベクトルを掛けると元の線分と同じベクトルになる
	 * 
	 */
	public static double[] getRotationAngle(double xx,double yy,double zz){
		double lenxz = xx*xx+zz*zz;
		double ysin = 0;
		double ycos = 0;
		
		
		double px = xx;
		double py = yy;
		double pz = zz;
		if(lenxz != 0){
			lenxz = Math.sqrt(lenxz);
			ycos = xx/lenxz;
			ysin = zz/lenxz;
			px = xx*ycos+zz*ysin;
			pz = 0;
		}
		
		 
		double zsin = 0;
		double zcos = 0;
		
		double lenyy = px*px+py*py;
		if(lenyy != 0){
			lenyy = Math.sqrt(lenyy);
			zsin = py/lenyy;
			zcos = px/lenyy;
		}
		double ret[] = {ycos,ysin,zcos,zsin};
		return ret;
	}
	
	public static double[] rotateWith(double xx,double yy,double zz,double[] angle){
		double px = xx*angle[2]-yy*angle[3];
		double py = xx*angle[3]+yy*angle[2];
		double rx = px*angle[0]-zz*angle[1];
		double rz = px*angle[1]+zz*angle[0];
		double ret[] = {rx,py,rz};
		return ret;
	}
	public static double[] rotateXZero(double xx,double yy,double zz,double[] angle){
		double rx = xx*angle[0]+zz*angle[1];
		double rz = -xx*angle[1]+zz*angle[0];
		
		double px = rx*angle[2]+yy*angle[3];
		double py = -rx*angle[3]+yy*angle[2];
		
		double ret[] = {px,py,rz};
		return ret;
	}
	
	
	/**
	 * pd1 の位置に pd2 を移動する
	 * 
	 */
	
	public static void placeResidueTo(PDBResidue pd1,PDBResidue pd2){
		
		
		if(pd1.getCA() == null 
		|| pd2.getCA() == null ){
			return;
		}
		
		
		Point3D p1 = pd1.getCA().loc;
		Point3D p2;
		Point3D p3 = pd2.getCA().loc;
		Point3D p4;
		if(pd1.getCB() != null){
			p2 = pd1.getCB().loc;
		}else{
			p2 = pd1.calcPseudoCB();
			
		}
		if(pd2.getCB() != null){
			p4 = pd2.getCB().loc;
		}else{
			p4 = pd2.calcPseudoCB();	
		}
		
		if( pd1.getCB() == null || pd2.getCB() == null){
			for(int ii = 0;ii < pd1.atoms.size();ii++){
				if(pd1.atoms.get(ii) != pd1.getCA()){
					PDBAtom a = pd1.atoms.get(ii);
					a.loc.x -= p1.x;
					a.loc.y -= p1.y;
					a.loc.z -= p1.z;
					
					a.loc.x += p3.x;
					a.loc.y += p3.y;
					a.loc.z += p3.z;
				}
			}
		}else{
			double angle[] = getRotationAngle(p4.x-p3.x,p4.y-p3.y,p4.z-p3.z);
			double angle2[] = getRotationAngle(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
			
			for(int ii = 0;ii < pd1.atoms.size();ii++){
				if(pd1.atoms.get(ii) != pd1.getCA()){
					PDBAtom a = pd1.atoms.get(ii);
					a.loc.x -= p1.x;
					a.loc.y -= p1.y;
					a.loc.z -= p1.z;
					double d1[] = rotateXZero(a.loc.x,a.loc.y,a.loc.z,angle2);
					double dd[] = rotateWith(d1[0],d1[1],d1[2],angle);
					a.loc.set(dd);
					
					a.loc.x += p3.x;
					a.loc.y += p3.y;
					a.loc.z += p3.z;
					
				}
			}
			
			pd1.getCA().loc.set(pd2.getCA().loc);
			if(pd1.getC() != null && pd2.getC() != null){
				pd1.getC().loc.set(pd2.getC().loc);
			}
			if(pd1.getO() != null && pd2.getO() != null){
				pd1.getO().loc.set(pd2.getO().loc);
			}
			if(pd1.getN() != null && pd2.getN() != null){
				pd1.getN().loc.set(pd2.getN().loc);
			}
			
			ArrayList<PDBAtom> sorted_tar = new ArrayList<>();
			ArrayList<PDBAtom> sorted_ref = new ArrayList<>(); 
			String[] gz = {"G","D","E","Z","H"};
			for(int ii = 0;ii < gz.length;ii++){
				for(PDBAtom a:pd1.atoms){
					if(a.pdb_atom_code.indexOf(gz[ii]) > -1){
						sorted_tar.add(a);
					}
				}
				for(PDBAtom a:pd2.atoms){
					if(a.pdb_atom_code.indexOf(gz[ii]) > -1){
						sorted_ref.add(a);
					}
				}
			}
		}
	}
	
	public static StringBuffer makeAtomLines(ArrayList<PDBResidue> al,int atomstart,int residuestart,String chainname
	){
		return  makeAtomLines(al,atomstart,residuestart,chainname,false);
	}
	public static StringBuffer makeAtomLines(ArrayList<PDBResidue> al,int atomstart,int residuestart,String chainname
	,boolean nocorrection){
		StringBuffer ret = new StringBuffer();
		String str[] = new String[16];
		int len[] = {6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6};//使用するカラム数
		int aa = 0;
		int rescount = 0;
		for(PDBResidue res:al){
			for(PDBAtom atom:res.atoms){
				str[0] = "ATOM  ";
				if(!nocorrection){
					str[1] = String.valueOf(atomstart+aa);
					str[2] = " ";
					str[3] = atom.pdb_atom_code;
					str[4] = "";//res.getAlternativeCode();
					str[5] = res.getName();
					str[6] = "";//res.getInsertionCode();
					str[7] = chainname;//res.parent.getName();
					str[8] = String.valueOf(rescount+residuestart);//String.valueOf(res.getResidueNumber());
				}else{
					
					str[1] = String.valueOf(atom.atom_number);
					str[2] = " ";
					str[3] = atom.pdb_atom_code;
					str[4] =  atom.getAltCode();
					str[5] = res.getName();
					str[6] = res.getInsertionCode();
					if(res.parent == null){
						str[7] = "X";
					}else{
						str[7] = res.parent.getName();
					}
					str[8] = String.valueOf(res.getResidueNumber());
				}
				
				
				str[9] = " ";
				str[10] = "   ";
				str[11] = String.format("%.3f",atom.loc.x);
				str[12] = String.format("%.3f",atom.loc.y);
				str[13] = String.format("%.3f",atom.loc.z);
				str[14] = String.format("%6.2f",atom.occupancy);
				str[15] = String.format("%6.2f",atom.bfactor);
				aa++;
				for(int kk = 0;kk < str.length;kk++){
					String ks = str[kk];
					if(ks.length() == len[kk]){
					}else if(ks.length() > len[kk]){
						ks = str[kk].substring(0,len[kk]);
					}else{
						ks = "                          "+str[kk];
						ks = ks.substring(ks.length()-len[kk]);
					}
					ret.append(ks);
				}
				ret.append("\n");
			}
			rescount++;
		}
		return ret;
	}
	
	
	public static ArrayList<PDBResidue> fillSpace(PDBResidue pd1,ArrayList<PDBResidue> crushcheck){
			
		ArrayList<PDBAtom> sorted_tar = new ArrayList<>();
		String[] gz = {"G","D","E","Z","H"};
		for(int ii = 0;ii < gz.length;ii++){
			for(PDBAtom a:pd1.atoms){
				if(a.pdb_atom_code.indexOf(gz[ii]) > -1){
					sorted_tar.add(a);
				}
			}
		}
		
		
		ArrayList<PDBAtom> checker = new ArrayList<>();
		
		for(PDBResidue res:crushcheck){
			if(res == pd1){
				checker.add(pd1.getN());
				checker.add(pd1.getO());
				checker.add(pd1.getC());
				continue;
			}
			if(res.getCB().distance(pd1.getCB()) < 16){//CB が 16 ANGSTROM 内にあるときだけチェックする
				checker.addAll(res.atoms);
			}
		}
		
		PDBAtom prev_tar = pd1.getCA();
		PDBAtom tar_current = pd1.getCB();
		ArrayList<PDBAtom> sorted_tar2 = new ArrayList<>(sorted_tar);
		
		double score = 0;
		while(true){
			if(sorted_tar.size() == 0){
				break;
			}
			score = Math.min(score,rotateAndFillSpace(prev_tar,tar_current,sorted_tar,checker));
			
			PDBAtom p = sorted_tar.remove(0);
			if(p.place_code_2.length() > 0){//二股に分かれている場合それ以上回転できないとする。これは実際は正しくないので暇のある時に変更のこと。
				break;
			}
			prev_tar = tar_current;
			tar_current = p;
		}
		
		
		if(score < 0){
			HashSet<PDBResidue> p = new HashSet<>();
			for(PDBAtom a:pd1.atoms){
				if("ONCA".indexOf(a.pdb_atom_code) < 0){
					for(PDBAtom c:checker){
						if(c.distance(a) < BONDLENGTH_STATIC){
							p.add(c.parent);
						}
					}
				}
				
			}
			
			return new ArrayList<PDBResidue>(p);
		}
		return null;
	}
	/**
	 * クラッシュする残基があった場合そのリストを返す
	 * 
	 */
	
	public static ArrayList<PDBResidue> refineCrush(PDBResidue pd1,ArrayList<PDBResidue> crushcheck){
		
	
		ArrayList<PDBAtom> sorted_tar = new ArrayList<>();
		String[] gz = {"G","D","E","Z","H"};
		for(int ii = 0;ii < gz.length;ii++){
			for(PDBAtom a:pd1.atoms){
				if(a.pdb_atom_code.indexOf(gz[ii]) > -1){
					sorted_tar.add(a);
				}
			}
		}
		
		
		ArrayList<PDBAtom> checker = new ArrayList<>();
		
		for(PDBResidue res:crushcheck){
			if(res == pd1){
				checker.add(pd1.getN());
				checker.add(pd1.getO());
				checker.add(pd1.getC());
				continue;
			}
			if(res.getCB().distance(pd1.getCB()) < 16){//CB が 16 ANGSTROM 内にあるときだけチェックする
				checker.addAll(res.atoms);
			}
		}
		
		PDBAtom prev_tar = pd1.getCA();
		PDBAtom tar_current = pd1.getCB();
		
		boolean flag = true;
		while(true){
			if(sorted_tar.size() == 0){
				break;
			}
			if(rotateAndCheckCrush(prev_tar,tar_current,sorted_tar,checker) == 0.0){
				
				break;
			}else{
				flag = false;
			}
			PDBAtom p = sorted_tar.remove(0);
			if(p.place_code_2.length() > 0){//二股に分かれている場合それ以上回転できないとする。これは実際は正しくないので暇のある時に変更のこと。
				break;
			}
			prev_tar = tar_current;
			tar_current = p;
		}
		
		if(!flag){
			HashSet<PDBResidue> p = new HashSet<>();
			for(PDBAtom a:pd1.atoms){
				if("ONCACB".indexOf(a.pdb_atom_code) < 0){
					for(PDBAtom c:checker){
						if(c.distance(a) < BONDLENGTH_STATIC){
							//System.out.println(a.parent.getName()+":"+a.pdb_atom_code+";"+c.parent.getName()+":"+c.pdb_atom_code+";"+a.parent.getCB().distance(c.parent.getCB()));
							p.add(c.parent);
						}
					}
				}
			}
			return new ArrayList<PDBResidue>(p);
		}
		
		return null;
	}
	
	/**
	 * 他残基の原子が BONDLENGTH_STATIC オングストロームに無いように原子を移動させる
	 * @param tar1 基準となるベクトルの開始点
	 * @param tar2 終了点
	 * @param sorted_tar 基準より後ろにあるため依存して回転する原子
	 * @param checker クラッシュをチェックする原子
	 */
	public static double rotateAndCheckCrush(PDBAtom tar1,PDBAtom tar2,ArrayList<PDBAtom> sorted_tar,ArrayList<PDBAtom> checker){
		ArrayList<HashMap<PDBAtom,Point3D>> rotated = new ArrayList<>();
		
		Point3D vec = new Point3D(tar2.loc);
		vec.x -= tar1.loc.x;
		vec.y -= tar1.loc.y;
		vec.z -= tar1.loc.z;
		
		double mindist = Double.MAX_VALUE;
		HashMap<PDBAtom,Point3D> minhash = null;
		pouter:for(double dd = 0.0;dd <= 3.2;dd+=0.1){
			for(int s = -1;s <= 1;s+=2){//回転の小さい方から調べる
				double currentdist = 0;
				HashMap<PDBAtom,Point3D> hm = new HashMap<>();
				for(PDBAtom t:sorted_tar){
					Point3D od = new Point3D(t.loc);
					od.x -= tar2.loc.x;
					od.y -= tar2.loc.y;
					od.z -= tar2.loc.z;
					Point3D pd = Point3D.rotate(od,vec,dd*s);
					pd.x += tar2.loc.x;
					pd.y += tar2.loc.y;
					pd.z += tar2.loc.z;
					hm.put(t,pd);
					for(PDBAtom rr:checker){
						currentdist += (pd.distance(rr.loc) < BONDLENGTH_STATIC)?(1):(0);
					}
					
				}
				if(mindist > currentdist){
					mindist = currentdist;
					minhash = hm;
				}
				
				if(currentdist == 0.0){
					break pouter;
				}
			}
		}
		for(PDBAtom a:minhash.keySet()){
			a.loc.set(minhash.get(a));
		}
		return mindist;
	}

	/**
	 * 他残基の原子がなるたけ 5.0 オングストロームにあるように原子を移動させる
	 * @param tar1 基準となるベクトルの開始点
	 * @param tar2 終了点
	 * @param sorted_tar 基準より後ろにあるため依存して回転する原子
	 * @param checker クラッシュをチェックする原子
	 */
	public static double rotateAndFillSpace(PDBAtom tar1,PDBAtom tar2,ArrayList<PDBAtom> sorted_tar,ArrayList<PDBAtom> checker){
		ArrayList<HashMap<PDBAtom,Point3D>> rotated = new ArrayList<>();
		
		Point3D vec = new Point3D(tar2.loc);
		vec.x -= tar1.loc.x;
		vec.y -= tar1.loc.y;
		vec.z -= tar1.loc.z;
		
		double maxscore = -1000;
		HashMap<PDBAtom,Point3D> minhash = null;
		for(double dd = 0.0;dd <= 3.2;dd+=0.2){
			for(int s = -1;s <= 1;s+=2){//回転の小さい方から調べる
				double currentscore = 0;
				HashMap<PDBAtom,Point3D> hm = new HashMap<>();
				for(PDBAtom t:sorted_tar){
					Point3D od = new Point3D(t.loc);
					od.x -= tar2.loc.x;
					od.y -= tar2.loc.y;
					od.z -= tar2.loc.z;
					Point3D pd = Point3D.rotate(od,vec,dd*s);
					pd.x += tar2.loc.x;
					pd.y += tar2.loc.y;
					pd.z += tar2.loc.z;
					hm.put(t,pd);
					for(PDBAtom rr:checker){
						double dist = pd.distance(rr.loc);
						
						currentscore += ( dist < 5.0)?((dist<3.0)?(-40):(1)):(0);
					}
				}
				if(maxscore < currentscore){
					maxscore = currentscore;
					minhash = hm;
				}
			}
		}
		if(minhash == null){
			return maxscore;
		}
		for(PDBAtom a:minhash.keySet()){
			a.loc.set(minhash.get(a));
		}
		return maxscore;
	}
	/**
	 * ある残基を回転してある残基と原子の位置が最も近いものに合わせる
	 * ラベルが同じ残基が何もなかった場合に False を返す
	 * tar1->tar2 ref1->ref2 は合っている
	 * 
	 */
	
	
	public static boolean rotateAndSetMin(PDBAtom tar1,PDBAtom tar2,PDBAtom ref1,PDBAtom ref2,ArrayList<PDBAtom> sorted_tar,ArrayList<PDBAtom> sorted_ref){
		HashMap<PDBAtom,ArrayList<PDBAtom>> corresp = new HashMap<>();
		
		for(PDBAtom t:sorted_tar){
			String tcode = t.place_code_1;
			ArrayList<PDBAtom> pa = new ArrayList<>();
			
			for(PDBAtom r:sorted_ref){
				if(t.pdb_atom_code.equals(r.pdb_atom_code)){//同じ名前の原子があった場合最優先する
					pa.add(r);
				}
			}
			
			
			if(pa.size() == 0){
				for(PDBAtom r:sorted_ref){
					String rcode = r.place_code_1;
					
					if(tcode.equals(rcode)){
						pa.add(r);
					}
				}
			}
			corresp.put(t,pa);
		}
		if(corresp.size() == 0){
			return false;
		}
		ArrayList<HashMap<PDBAtom,Point3D>> rotated = new ArrayList<>();
		
		Point3D vec = new Point3D(tar2.loc);
		vec.x -= tar1.loc.x;
		vec.y -= tar1.loc.y;
		vec.z -= tar1.loc.z;
		
		double mindist = Double.MAX_VALUE;
		HashMap<PDBAtom,Point3D> minhash = null;
		for(double dd = 0.0;dd <= 6.4;dd+=0.1){
			double currentdist = 0;
			HashMap<PDBAtom,Point3D> hm = new HashMap<>();
			for(PDBAtom t:sorted_tar){
				Point3D od = new Point3D(t.loc);
				od.x -= tar2.loc.x;
				od.y -= tar2.loc.y;
				od.z -= tar2.loc.z;
				Point3D pd = Point3D.rotate(od,vec,dd);
				pd.x += tar2.loc.x;
				pd.y += tar2.loc.y;
				pd.z += tar2.loc.z;
				hm.put(t,pd);
				if(corresp.containsKey(t)){
					for(PDBAtom rr:corresp.get(t)){
						currentdist += pd.distance(rr.loc);
					}
				}
			}
			
			if(mindist > currentdist){
				mindist = currentdist;
				minhash = hm;
			}
			rotated.add(hm);
		}
		
		for(PDBAtom a:minhash.keySet()){
			a.loc.set(minhash.get(a));
		}
		return true;
		
	}
	
	
	
	
	
	public static void writeToFile(StringBuffer sb,String filename){
		try{
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename,false),"UTF-8")));
			pw.write(sb.toString());
			pw.close();
		}catch(Exception ex){
			ex.printStackTrace();
		}
	}
	public static void writeToFile(ArrayList<String> al,String filename){
		try{
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename,false),"UTF-8")));
			for(String s:al){
				pw.write(s);
				pw.write("\n");
			}
			pw.close();
		}catch(Exception ex){
			ex.printStackTrace();
		}
	}
	
		
	
	public static void printArray(double[] d){
		for(int ii = 0;ii < d.length;ii++){
			System.out.print(d[ii]+",");
		}
		System.out.println("");
	}
	
	
	public static ArrayList<String> positionArrayWithOffset(ArrayList<PDBResidue> res,int offset){
		ArrayList<String> ret = new ArrayList<>();
		for(int ii = offset;ii < res.size();ii++){
			ret.add(res.get(ii).getPositionCode());
		}
		return ret;
	}
	
	
	public static ArrayList<String> aaArrayWithOffset(ArrayList<PDBResidue> res,int offset){
		ArrayList<String> ret = new ArrayList<>();
		for(int ii = offset;ii < res.size();ii++){
			ret.add(res.get(ii).getName());
		}
		return ret;
	}
	
	/**
	 * 
	 * @param target
	 * @param template
	 * @param pc
	 * @param nocorrupt マップされた残基は位置を変えない
	 * @return 
	 */
	public static SWResult alignmentRefineModel(String target,String template
			,PDBChain pc){
			TemplateBaseModeller tbm  = new TemplateBaseModeller();
			TBMResult res = tbm.model(target
					,template,pc);
			FuzzyTreeAlign fta = new FuzzyTreeAlign();
			//Alanin に変更するか
			//res.template_residues = PepProcess.makeFilteredAA(res.template_residues);
			
			PSSMScoring testpssm = fta.calcPSSM(res.template_residues,res.template_residues);
			//>sp|Q8I2J4.1|PROF_PLAF7 RecName: Full=Profilin
			for(int ii = 0;ii < res.template_residues.size();ii++){
				PDBResidue d = res.template_residues.get(ii);
				if(res.template_target_map.containsKey(d)){
					PDBResidue ta = res.template_target_map.get(d);
					char tcc = PDBResidue.aaMap.get(ta.getName()).charAt(0);
					//FuzzyTreeAlign.fillWithBLOSUM(testpssm,tcc,ii);
					int cpos = testpssm.pssm.char_index_map[tcc];
					testpssm.pssm.scores.get(ii).set(cpos
							,(double)blosum.scoreMat[tcc][tcc]);
					//System.out.println(blosum.scoreMat[tcc][tcc]);
				}
			}
			
			SmithWaterman sw = new SmithWaterman();
			sw.penalO =30;
			sw.penalE =0;
			SWResult swres = sw.align_with_PSSM(target.replaceAll("[^A-Za-z]",""),testpssm);
			return swres;
	}
	public static void model(String infile,String outfilename,boolean nocorrupt
			,int realignnum,int realignwindow,int rotemerrefine){
		String targetfas = infile;

		ArrayList<TBMEntry> ls = TBMEntry.loadFasta(targetfas);
		HashMap<Integer,TBMEntry> targets = new HashMap<>();
		ArrayList<TBMEntry> templates = new ArrayList<>();
		for(TBMEntry t:ls){
			if(t.name.equals("target")){
				targets.put(t.id, t);
			}else{
				templates.add(t);
			}
		}
		TemplateBaseModeller tbm = new TemplateBaseModeller();
		ArrayList<String> messages = new ArrayList<>();
		ArrayList<FloatingResidue> allres = new ArrayList<>();
		for(TBMEntry t:templates){
			PDBData pdb1 = PDBData.loadPDBFile(t.file);
			TBMEntry target = targets.get(t.targetId);
			
			//SWResult swres = alignmentRefineModel(target.sequence,t.sequence,pdb1.chains.get(t.chainName));
			//System.out.println("--------------------");
			//System.out.println(target.sequence);
			//System.out.println(t.sequence);
			//System.out.println(swres.score);
			//target.sequence = SmithWaterman.listToString(swres.qseq);
			//t.sequence = SmithWaterman.listToString(swres.sseq);
			
			//System.out.println(target.sequence);
			//System.out.println(t.sequence);
			//System.out.println("#--------------------");
			
			if(realignnum > 0){
				ArrayList<PDBResidue> pos_temp = prepareTemplateResidues(
						PepProcess.makeFilteredAA(pdb1.chains.get(t.chainName).residues));

				ArrayList<String> temp = new ArrayList<>();
				String[] cc = t.sequence.split("");
				for(String c:cc){
					if(c.length() > 0){
						temp.add(c);
					}
				}
				filtAlignedResidues(temp,pos_temp);
				ArrayList<Character> tempchar = new ArrayList<>();
				for(String c:temp){
					tempchar.add(c.toCharArray()[0]);
				}
				t.sequence = SmithWaterman.listToString(tempchar);
				for(int ii = 0;ii < realignnum;ii++){
					FuzzyTreeAlign fz = new FuzzyTreeAlign();
					SWResult sres = fz.refineGap(SmithWaterman.toCharList(target.sequence)
							, SmithWaterman.toCharList(t.sequence), pos_temp,realignwindow);
					//System.out.println(target.sequence);
					//System.out.println(t.sequence);
					//System.out.println(SmithWaterman.listToString(sres.qseq));
					//System.out.println(SmithWaterman.listToString(sres.sseq));
					target.sequence = SmithWaterman.listToString(sres.qseq);
					t.sequence = SmithWaterman.listToString(sres.sseq);
				}
			}
			
			TBMResult res = tbm.model(target.sequence
					,t.sequence,pdb1.chains.get(t.chainName));
			
			/*
			ChainBuilder.savePDBChains(res.residues,"testres1.pdb");
			
			
			TBMResult res2 = AliFine.refined(res);
			for(PDBResidue r:res2.residues){
				if(res2.target_template_map.containsKey(r)){
					TemplateBaseModeller.placeResidueTo(r, res2.target_template_map.get(r));
				}
			}
			ChainBuilder.savePDBChains(res2.residues,"testres2.pdb");
			*/
			ArrayList<FloatingResidue> modelled = new ArrayList<>();
			ArrayList<FloatingResidue> unmapped = new ArrayList<>();
			for(int ii = 0;ii < res.residues.size();ii++){
				PDBResidue p = null;
				PDBResidue n = null;
				if(ii > 0 && res.mapped.get(ii-1)){
					p = res.residues.get(ii-1);
				}
				if(ii < res.residues.size()-1 && res.mapped.get(ii+1)){
					n = res.residues.get(ii+1);
				}
				modelled.add(new FloatingResidue(res.residues.get(ii),p,n));
				
			}
			HashSet<FloatingResidue> mapped = new HashSet<>();
			HashSet<FloatingResidue> flag = new HashSet<>();
			for(int ii = 0;ii < res.residues.size();ii++){
				if(ii > 0){
					modelled.get(ii).prev = modelled.get(ii-1);
				}
				if(ii < res.residues.size()-1){
					modelled.get(ii).next = modelled.get(ii+1);
				}
				
				if(!res.mapped.get(ii)){
					//Map されてない残基の前後もマップされてないことにする
					if(ii > 0 && !flag.contains(modelled.get(ii-1))){
						unmapped.add(modelled.get(ii-1));
						flag.add(modelled.get(ii-1));
					}
					if(ii < modelled.size()-1 && !flag.contains(modelled.get(ii+1))){
						unmapped.add(modelled.get(ii+1));
						flag.add(modelled.get(ii+1));
					}
					if(!flag.contains(modelled.get(ii))){
						unmapped.add(modelled.get(ii));
					}
					flag.add(modelled.get(ii));
				}
			}
			mapped.addAll(modelled);
			mapped.removeAll(unmapped);
			//tbm.fillGaps(unmapped, mapped,100);
			
			PDBChain parentchain = new PDBChain(target.chainName);
			for(int ii = 0;ii < modelled.size();ii++){
				modelled.get(ii).parent = parentchain;
				modelled.get(ii).setResidueNumber(ii+target.start);
			}
			
			
			messages.add("#residue stats");
			for(int ii = 0;ii < modelled.size();ii++){
				PDBResidue r = modelled.get(ii);
				if(mapped.contains(r)){
					messages.add("chain_name:"+r.parent.getName()+"\t"+"residue_name:"+r.getName()
							+"\t"+"residue_number:"+r.getResidueNumber()+"\t"+"ok");
				}else{
					messages.add("chain_name:"+r.parent.getName()+"\t"+"residue_name:"+r.getName()
							+"\t"+"residue_number:"+r.getResidueNumber()+"\t"+"ng");
				}
			}
			
			
			ChainBuilder cb = new ChainBuilder();
			HashSet<Integer> processed = new HashSet<>();
			
			for(FloatingResidue fr:modelled){
				if(fr.backboneIndex == -1){
					fr.calcPseudoNeighbour(cb.backbones.getRandomly(fr.getName()));
				}
			}
			for(FloatingResidue fr:modelled){
				fr.saveLoc();
				if(fr.sidechainIndex == -1){
					fr.changeSideChain(cb.sidechains.getNext(fr));
				}
			}
			//Refine3.fixAllUnmapped(modelled, unmapped_withneighbour, 100);
			
			
			//ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(modelled,unmapped_withneighbour);
			ArrayList<PDBResidue> modelled_r = new ArrayList<>();
			modelled_r.addAll(modelled);
			ChainBuilder cb2 = new ChainBuilder();
			
			FuzzyDecisionTreeScoring_generator j 
					= new FuzzyDecisionTreeScoring_generator(new FeatureGeneratorCB_20180501());
			cb2.scoring.clear();
			cb2.scoring.add(j);
			cb2.prepare(modelled);
			
			//fixme unmapped の前後が mapped の場合それは含める
			if(!nocorrupt){
				ArrayList<ArrayList<FloatingResidue>> regions = getUnmappedRegions(modelled,new HashSet<FloatingResidue> (unmapped));
				for(ArrayList<FloatingResidue> pd:regions){
					ArrayList<PDBResidue> pd_r = new ArrayList<>();	
					pd_r.addAll(pd);
					//System.out.println(pd.size()+";;"+modelled.size());
					changeBackBone(cb2,modelled,modelled_r,pd,pd_r,true,10);
				} 
				
				
				HashSet<FloatingResidue> unmapped_withneighbour = new HashSet<>();
				unmapped_withneighbour.addAll(unmapped);

				//fixme 露出残基をとって、露出している方から削っていく
				//ChainBuilder.saveChains(modelled,outfilename+".p.pdb");
				for(int ii = 0;ii < modelled.size()-1;ii++){
					if(modelled.get(ii).getC().distance(modelled.get(ii+1).getN()) > 1.4){
						//System.out.println(";;"+modelled.get(ii).getRepresentativeCode());
						unmapped_withneighbour.addAll(cb.getGapFillingResidues(modelled.get(ii)));
					}

				}
				regions = getUnmappedRegions(modelled,new HashSet<FloatingResidue> (unmapped_withneighbour));

				for(ArrayList<FloatingResidue> pd:regions){
					ArrayList<PDBResidue> pd_r = new ArrayList<>();	
					pd_r.addAll(pd);
					//System.out.println(pd.get(0).getRepresentativeCode()+";"+pd.size()+";;"+modelled.size());
					changeBackBone(cb2,modelled,modelled_r,pd,pd_r,true,50*pd.size());
				}
			}
			allres.addAll(modelled);
		}
		if(rotemerrefine > 0){
			ChainBuilder cb = new ChainBuilder();
			cb.prepare(allres);
			for(int ii = 0;ii < rotemerrefine;ii++){
				for(FloatingResidue r:allres){
					cb.refiner.maxRotamer(r, -5,0.3, true, cb.sidechains, false);
				}
			}
		}
		ChainBuilder.saveChains(allres,outfilename);
		writeToFile(messages,outfilename+".res");
	}
	
	
}

class TBMEntry{
	String name = "";
	String file ="";
	String chainName="";
	String sequence = "";
	int id = -1;
	int targetId = -1;
	int start = 1;
	
	TBMEntry(){
		
	}
	
	public static ArrayList<TBMEntry> loadFasta(String filename){
		ArrayList<TBMEntry> ret = new ArrayList<>();
		try{
			System.out.println(filename);
			FileInputStream iss = new FileInputStream(new File(filename));
			BufferedReader br = new BufferedReader(new InputStreamReader(iss));
			String ln = null;
			StringBuffer buff = new StringBuffer();
			String currentname = "";
			String currentfilename = "";
			String currentchainname = "";
			int currentid =  -1;
			int currenttargetid = -1;
			int currentstartnum = 1;
			Pattern lpat = Pattern.compile(">[\\s]*([^\\s]+)");
			Pattern descpat = Pattern.compile(">[\\s]*[^\\s]+[\\s]+([^\\s].+)");
			Pattern labelpat = Pattern.compile("[\\s]*([^\\s]+)=(.+)");
			Pattern quotepat = Pattern.compile("^(['\"])");
			while((ln = br.readLine()) != null){
				Matcher mat = lpat.matcher(ln);
				if(mat.find()){
					if(buff.length() > 0){
						TBMEntry t = new TBMEntry();
						t.name = currentname;
						t.file = currentfilename;
						t.chainName = currentchainname;
						t.sequence = buff.toString();
						t.id = currentid;
						t.targetId = currenttargetid;
						t.start = currentstartnum;
						ret.add(t);
						
					}
					
					buff = new StringBuffer();
					currentname = "";
					currentfilename = "";
					currentchainname = "";
					currentid =  -1;
					currenttargetid = -1;
					currentstartnum = 1;
					currentname = mat.group(1).replaceAll("[\\s]","");
					String desc = "";
					Matcher mmat = descpat.matcher(ln);
					if(mmat.find()){
						desc = mmat.group(1);
					}
					Matcher pmat = labelpat.matcher(desc);
					while(pmat.find()){
						String code = pmat.group(1);
						String tail = pmat.group(2);
						String res = "";
						Matcher qmat = quotepat.matcher(tail);
						if(qmat.find()){
							String sep = qmat.group(1);
							Matcher smat = Pattern.compile(sep+"([^"+sep+"]+)"+sep).matcher(tail);
							if(smat.find()){
								res = smat.group(1);
								desc = tail.replaceFirst(sep+"([^"+sep+"]+)"+sep,"");
							}else{
								throw new RuntimeException("parsing failed:\n"+ln);
							}
						}else{
							Matcher smat = Pattern.compile("^([^\\s]+)").matcher(tail);
							if(smat.find()){
								res = smat.group(1);
								desc = tail.replaceFirst("^([^\\s]+)","");
							}
						
						}
						pmat = labelpat.matcher(desc);
						
						if(code.equals("file")){
							currentfilename = res;
						}else if(code.equals("chain")){
							currentchainname = res;
						}else if(code.equals("id")){
							currentid = Integer.parseInt(res);
						}else if(code.equals("targetid")){
							currenttargetid = Integer.parseInt(res);
						}else if(code.equals("start")){
							currentstartnum = Integer.parseInt(res);
						}else{
							System.err.println("unknown label "+code+" is not used.");
						}
					}
					
				}else{
					buff.append(ln.replaceAll("[\\s]",""));
				}
			}
			
			
			if(buff.length() > 0){
				TBMEntry t = new TBMEntry();
				t.name = currentname;
				t.file = currentfilename;
				t.chainName = currentchainname;
				t.sequence = buff.toString();
				t.id = currentid;
				t.targetId = currenttargetid;
				ret.add(t);
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}

}

/*
use strict;
use warnings;

my $sepdir_path = "E:/ubuntu4/MaHMMPDB/MaHMMPDB/sepdb/";
my $pepbuilder = "java -jar C:/Users/kimidori/Documents/NetBeansProjects/PepBuilderJ/dist/PepBuilderJ.jar";

my $infile = 'C:\dummy\vbox_share\casp13\queries\T0950_mfas\T0950.jackali.mfa_0.a3m.hhr_hhhit0.fas';
my $outfile = $infile.".t";
my ($name,$desc,$seq) = loadseq($infile);
my @nn = @{$name};

open(OUT,"> $outfile");
for(my $ii = 0;$ii <= $#nn;$ii++){
	if($nn[$ii] =~ /^T[Ss0R][0-9][0-9]/){
		print OUT ">target chain=\"A\" id=0\n";
		print OUT join("",@{${$seq}[$ii]})."\n";
	}else{
		$nn[$ii] =~ /^...._([^\.]+)/;
		my $chainname = $1;
		my $filename = getSepPDBPath($nn[$ii]).".pdb";
		if(!-f $filename){
			die;
		}
		print OUT ">".$nn[$ii]."  chain=\"$chainname\" id=$ii targetid=0 file=\"$filename\"\n";
		print OUT join("",@{${$seq}[$ii]})."\n";
		
	}
}
close(OUT);

system($pepbuilder." -model -in $outfile -out $outfile".".pdb ");



#ドメインごとに区切った奴
sub getSepPDBPath{
	my $fullname = $_[0];
	$fullname =~ /^.(..)/;
	my $sepdirname = $1;
	
	#chain id までが入っている
	$fullname =~ /^(...._[^\.]+)/;
	my $cdirname = $1;
	my $fullpath = $sepdir_path."/".$sepdirname."/".$cdirname."/".$fullname;
	
	return $fullpath;
}



sub loadseq{
	my $infile = $_[0];
	
	my @name;
	my @desc;
	my @seq;
	open(FFIN,$infile) or die;
	while(my $ss = <FFIN>){
		$ss =~ s/[\r\n]//g;
		if($ss =~ /^[\s]*>[\s]*([^\s]+)/){
			push(@name,$1);
			push(@desc,"");
			push(@seq,"");
			
			if($ss =~ /^[\s]*>[\s]*[^\s]+[\s]+([^\s].+)/){
				$desc[$#name] = $1;
			}
		}else{
			if($#name > -1){
				$seq[$#name] .= $ss;
			}
		}
	}
	close(FFIN);
	for(my $ii = 0;$ii <= $#seq;$ii++){
		$seq[$ii] =~ s/[\s]//g;
		my @a = split(//,$seq[$ii]);
		$seq[$ii] = \@a;
	}
	
	return \@name,\@desc,\@seq;
}

*/
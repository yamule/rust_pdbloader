/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author kimidori
 */
public class PepProcess {


	public static HashMap<Character,String> one_to_three = new HashMap<>();
	public static HashMap<String,Character> three_to_one = new HashMap<>();
	
	static{
		three_to_one.put("ALA",'A');
		three_to_one.put("ARG",'R');
		three_to_one.put("ASN",'N');
		three_to_one.put("ASP",'D');
		three_to_one.put("CYS",'C');
		three_to_one.put("GLN",'Q');
		three_to_one.put("GLU",'E');
		three_to_one.put("GLY",'G');
		three_to_one.put("HIS",'H');
		three_to_one.put("ILE",'I');
		three_to_one.put("LEU",'L');
		three_to_one.put("LYS",'K');
		three_to_one.put("MET",'M');
		three_to_one.put("PHE",'F');
		three_to_one.put("PRO",'P');
		three_to_one.put("SER",'S');
		three_to_one.put("THR",'T');
		three_to_one.put("TRP",'W');
		three_to_one.put("TYR",'Y');
		three_to_one.put("VAL",'V');
		three_to_one.put("UNK",'X');
		for(String s:three_to_one.keySet()){
			one_to_three.put(three_to_one.get(s), s);
		}
	}
	
	
	public static String getAACode(char c){
		return one_to_three.get(c);
	}
	public static char getAALetter(String s){
		return three_to_one.get(s);
	}
	public static double len(Point3D a){
		double len = a.x*a.x+a.y*a.y+a.z*a.z;
		if(len > 0){
			return Math.sqrt(len);
		}
		return 0;
	}
	
	public static double len(double x,double y){
		double len = x*x+y*y;
		if(len > 0){
			return Math.sqrt(len);
		}
		return len;
	}
	public static Point3D subtracted(Point3D a,Point3D b){
		return new Point3D(a.x-b.x,a.y-b.y,a.z-b.z);
	}
	public static Point3D added(Point3D a,Point3D b){
		return new Point3D(a.x+b.x,a.y+b.y,a.z+b.z);
	}
	
	public static Point3D norm(Point3D p1,Point3D p2){
		return PepProcess.norm(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
	}
	public static Point3D normalized(Point3D p){
		double len = p.x*p.x+p.y*p.y+p.z*p.z;
		Point3D ret =  new Point3D(p);
		
		if(len == 0){
			return ret;
		}
		len = Math.sqrt(len);
		ret.x /= len;
		ret.y /= len;
		ret.z /= len;
		return ret;
	}
	public static double dot(Point3D p1,Point3D p2){
		return p1.x*p2.x + p1.y*p2.y+p1.z*p2.z;
	}
	public static Point3D norm(double v1x,double v1y,double v1z,double v2x,double v2y,double v2z){
		
		Point3D  ret = new Point3D(
		(v1y*v2z-v1z*v2y),
		(v1z*v2x-v1x*v2z),
		(v1x*v2y-v1y*v2x));
		
		return normalized(ret);
		
	}
	
	public static void adjustVector3D_Dot(Point3D targetstart_,
			Point3D targetend_,
			Point3D targetend2_,
			Point3D refstart_,
			Point3D refend_,
			Point3D refend2_,
			ArrayList<Point3D> al){
		
		Point3D targetstart = new Point3D(targetstart_);
		Point3D targetend = new Point3D(targetend_);
		Point3D targetend2 = new Point3D(targetend2_);
		Point3D refstart = new Point3D(refstart_);
		Point3D refend = new Point3D(refend_);
		Point3D refend2 = new Point3D(refend2_);
		
		
		ArrayList<Point3D> tmpp = new ArrayList<>(al);
		Point3D t = new Point3D(targetend2);
		tmpp.add(t);
		adjustVector3D_Dot(targetstart,targetend,refstart,refend,tmpp);
		tmpp.remove(tmpp.size()-1);
		Point3D p = subtracted(t,refstart);
		Point3D q = subtracted(refend,refstart);
		Point3D r = subtracted(refend2,refstart);
		
		Point3D c1 = norm(q,p);
		Point3D c2 = norm(q,r);
		
		Point3D c3 = norm(q,c1);
		
		double d1 = dot(c1,c2);
		double d2 = dot(c2,c3);
		
		double rot = Math.atan2(d2,d1);
		
		
		for(Point3D pp:al){
			pp.x -= refstart.x;
			pp.y -= refstart.y;
			pp.z -= refstart.z;
			Point3D ccc = Point3D.rotate(pp,q,rot);
			
			ccc.x += refstart.x;
			ccc.y += refstart.y;
			ccc.z += refstart.z;
			pp.set(ccc);
		}
	}
	/**
	 * targetstart -> targetend 
	 * というベクトルが
	 * refstart->refend
	 * になるような回転を、al 内のPoint に適用する。
	 * 内積と、ベクトルを軸にする回転の公式を使っている。
	 * ちょっと遅い。
	 * @param targetstart
	 * @param targetend
	 * @param refstart
	 * @param refend
	 * @param al 
	 */
	public static void adjustVector3D_Dot(Point3D targetstart_,
			Point3D targetend_,
			Point3D refstart_,
			Point3D refend_,
			ArrayList<Point3D> al){
		
		Point3D targetstart = new Point3D(targetstart_);
		Point3D targetend = new Point3D(targetend_);
		Point3D refstart = new Point3D(refstart_);
		Point3D refend = new Point3D(refend_);
		
		Point3D a = subtracted(targetend,targetstart);
		Point3D b = subtracted(refend,refstart);
		double blen = b.x*b.x+b.y*b.y+b.z*b.z;
		double alen = a.x*a.x+a.y*a.y+a.z*a.z;
		if(alen > 0 && blen > 0){
			alen = Math.sqrt(alen);
			blen = Math.sqrt(alen);
		}else{
			for(Point3D pp:al){
				pp.x -= targetstart.x;
				pp.y -= targetstart.y;
				pp.z -= targetstart.z;
				
				pp.x += refstart.x;
				pp.y += refstart.y;
				pp.z += refstart.z;
			}
			return;
		}
		
		
		Point3D p = norm(a,b);
		Point3D c = norm(p,a);
		b.x /= blen;
		b.y /= blen;
		b.z /= blen;
		
		
		a.x /= alen;
		a.y /= alen;
		a.z /= alen;
		
		double s = dot(c,b);
		double t = dot(b,a);
		double rot = Math.atan2(s, t);
		
		
		for(Point3D pp:al){
			pp.x -= targetstart.x;
			pp.y -= targetstart.y;
			pp.z -= targetstart.z;
			Point3D ccc = Point3D.rotate(pp,p,rot);
			pp.set(ccc);
			pp.x += refstart.x;
			pp.y += refstart.y;
			pp.z += refstart.z;
		}
		
		
	}
	
	
	/**
	 * 
	 * targetstart -> targetend 
	 * というベクトルが
	 * refstart->refend
	 * になるような回転を、al 内のPoint に適用する。
	 * 一度 target のベクトルを X 要素以外ゼロのベクトルに直して、ref のベクトルに合わせる計算を
	 * べた書きしている。
	 * @param targetstart
	 * @param targetend
	 * @param refstart
	 * @param refend
	 * @param al 
	 */
	public static void adjustVector3D(Point3D targetstart,
			Point3D targetend,
			Point3D refstart,
			Point3D refend,
			ArrayList<Point3D> al){
		
		Point3D a = subtracted(targetend,targetstart);
		Point3D b = subtracted(refend,refstart);
		double tx = targetstart.x;
		double ty = targetstart.y;
		double tz = targetstart.z;
		
		double rx = refstart.x;
		double ry = refstart.y;
		double rz = refstart.z;
		
		
		
		
		
		double alen = len(a);
		double blen = len(b);
		if(alen <= 0 || blen <= 0){
			for(Point3D pp:al){
				pp.x -= tx;
				pp.y -= ty;
				pp.z -= tz;
				
				pp.x += rx;
				pp.y += ry;
				pp.z += rz;
			}
			return;
		}
		
		double a2len = len(a.x,a.y);
		double b2len = len(b.x,b.y);
		
		boolean a2zeroflag = false;
		
		double cost = 0;
		double sint = 0;
		if(a2len > 0){
			cost = a.x/a2len;
			sint = a.y/a2len;
		}else{
			a2zeroflag = true;
		}
		double costk = a2len/alen;
		double sintk = a.z/alen;
		
		
		boolean b2zeroflag = false;
		double coss = 0;
		double sins = 0;
		if(b2len > 0){
			sins = b.y/b2len;
			coss = b.x/b2len;
		}else{
			b2zeroflag = true;
		}
		
		double cossk = b2len/blen;
		double sinsk = b.z/blen;
		
		for(Point3D ppz:al){
			ppz.x -= tx;
			ppz.y -= ty;
			ppz.z -= tz;
			
			
			double dddx = ppz.x;
			double dddy = ppz.y;
			if(!a2zeroflag){
				dddx = cost*ppz.x+sint*ppz.y;
				dddy = -sint*ppz.x+cost*ppz.y;
			}
			Point3D ppp = new Point3D(
			costk*dddx+sintk*ppz.z,
					dddy,
			-sintk*dddx+costk*ppz.z
			);
			
			
			double zzx = ppp.x*cossk-ppp.z*sinsk;
			double zzz = ppp.z*cossk+ppp.x*sinsk;
			if(!b2zeroflag){
				ppp.x = coss*zzx-sins*ppp.y;
				ppp.y = sins*zzx+coss*ppp.y;
			}
			ppp.z = zzz;
			
			ppp.x += rx;
			ppp.y += ry;
			ppp.z += rz;
			ppz.set(ppp);
		}
		return;
	}
	
	/**
	 *  targetstart -> targetend 
	 *  targetstart -> targetend2
	 * という面が
	 * refstart->refend
	 * refstart-> refend2
	 * と同じ方向の法線ベクトルを持つようになるような回転を、al 内のPoint に適用する。
	 * 一度 target のベクトルを Z 要素ゼロの面に直して、ref のベクトルに合わせる計算を
	 * べた書きしている。
	 * @param targetstart
	 * @param targetend
	 * @param targetend2
	 * @param refstart
	 * @param refend
	 * @param refend2
	 * @param al 
	 */
	
	public static void adjustVector3D(Point3D targetstart,
			Point3D targetend,
			Point3D targetend2,
			Point3D refstart,
			Point3D refend,
			Point3D refend2,
			ArrayList<Point3D> al){
		
		Point3D a = subtracted(targetend,targetstart);
		//Point3D a2 = subtracted(targetend2,targetstart);
		Point3D b = subtracted(refend,refstart);
		//Point3D b2 = subtracted(refend2,refstart);
		
		double[] ttrot = getRotValues(targetstart,targetend,targetend2);
		double[] rrrot = getRotValues(refstart,refend,refend2);
		//double[] ttrot = {trot[0],trot[1],trot[2],trot[3],1,0};
		//double[] rrrot = {rrot[0],rrot[1],rrot[2],rrot[3],1,0};
		double tx = targetstart.x;
		double ty = targetstart.y;
		double tz = targetstart.z;
		
		double rx = refstart.x;
		double ry = refstart.y;
		double rz = refstart.z;
		
		
		for(Point3D ppz:al){
			ppz.x -= tx;
			ppz.y -= ty;
			ppz.z -= tz;
			
			Point3D ppp = new Point3D(ppz);
			apply3D(ppp,ttrot,true);
			apply3D(ppp,rrrot,false);
			
			
			ppp.x += rx;
			ppp.y += ry;
			ppp.z += rz;
			ppz.set(ppp);
		}
	}
	
	
	public static final int COST = 0;
	public static final int SINT = 1;
	public static final int COSTK = 2;
	public static final int SINTK = 3;
	public static final int COSTU = 4;
	public static final int SINTU = 5;
	/**
	 *  X  要素以外ゼロのベクトルを
	 * targetstart-> end というベクトル方向に合わせるための回転要素を計算し返す。
	 * @param targetstart
	 * @param targetend
	 * @return 
	 */
	public static double[] getRotValues(Point3D targetstart,
			Point3D targetend){
		Point3D a = subtracted(targetend,targetstart);
		double cost = 1;
		double sint = 0;
		double costk = 1;
		double sintk = 0;
		double alen = len(a);
		double a2len = len(a.x,a.y);
		if(a2len > 0){
			cost = a.x/a2len;
			sint = a.y/a2len;
		}else{
			//a2zeroflag = true;
		}
		if(alen > 0){
			costk = a2len/alen;
			sintk = a.z/alen;
		}
		double[] ret = new double[4];
		ret[COST] = cost;
		ret[SINT] = sint;
		ret[COSTK] = costk;
		ret[SINTK] = sintk;
		return ret;
	}
	/**
	 * X,Y 平面（Z 要素がゼロ）を
	 * tst->ten
	 * tst->ten2 という三角形の面と同じ法線を持つよう変換するための回転の要素を返す。
	 * ここで計算した値と変換したい点のリストを apply3D 関数に渡すと変換してくれる。
	 * @param tst
	 * @param ten
	 * @param ten2
	 * @return 
	 */
	
	public static double[] getRotValues(Point3D tst,
			Point3D ten,Point3D ten2){
		double rota[] =getRotValues(tst,ten); 
		
		double tx = ten2.x-tst.x;
		double ty = ten2.y-tst.y;
		double tz = ten2.z-tst.z;	
		
		double dddx = rota[COST]*tx+rota[SINT]*ty;
		double dddy = -rota[SINT]*tx+rota[COST]*ty;
		double ddx = rota[COSTK]*dddx+rota[SINTK]*tz;
		double ddz = -rota[SINTK]*dddx+rota[COSTK]*tz;
		
		double len3 = len(dddy,ddz);
		
		double costu = 1;
		double sintu = 0;
		if(len3 > 0){
			costu = dddy/len3;
			sintu = ddz/len3;
		}
		double ret[] = new double[6];
		for(int ii = 0;ii < 4;ii++){
			ret[ii] = rota[ii];
		}
		ret[COSTU] = costu;
		ret[SINTU] = sintu;
		return ret;
	}
	
	/**
	 * getRotValues 関数で得た要素を利用し、与えられた点を回転する。
	 * revflag を true にすると、rotvalue を計算するために使用した面の Z 要素ゼロになるような回転を適用する。
	 * @param p
	 * @param rot
	 * @param revflag 
	 */
	public static void apply3D(Point3D p,double[] rot,boolean revflag){
		if(revflag){
			double dx = rot[COST]*p.x+rot[SINT]*p.y;
			double dy = -rot[SINT]*p.x+rot[COST]*p.y;
			double ddx = rot[COSTK]*dx+rot[SINTK]*p.z;
			double dz = -rot[SINTK]*dx+rot[COSTK]*p.z;
			double ddy = rot[COSTU]*dy+rot[SINTU]*dz;
			double ddz = -rot[SINTU]*dy+rot[COSTU]*dz;
			p.set(ddx,ddy,ddz);
			return;
		}else{
			double dy = rot[COSTU]*p.y-rot[SINTU]*p.z;
			double dz = rot[SINTU]*p.y+rot[COSTU]*p.z;
			double dx = rot[COSTK]*p.x-rot[SINTK]*dz;
			double ddz = rot[SINTK]*p.x+rot[COSTK]*dz;

			double ddx = rot[COST]*dx-rot[SINT]*dy;
			double ddy = rot[SINT]*dx+rot[COST]*dy;
			p.set(ddx,ddy,ddz);
		}
	}
	
	public static double dihedralAngle(Point3D a,Point3D b,Point3D c,Point3D d){
		Point3D ab = subtracted(a,b);
		Point3D bc = subtracted(b,c);
		Point3D dc = subtracted(d,c);
		
		Point3D p = norm(bc,ab);
		
		Point3D q = norm(bc,dc);
		Point3D r = norm(bc,q);
		
		double s = dot(r,p);
		double t = dot(q,p);
		
		//if(s > 0 && t > 0){//DSSP ではこうしているが・・・？あとでチェック
			return Math.atan2(s, t)*180/Math.PI;
		//}
		//return 0;
		
	}
	
	
	public static double phi(PDBResidue a,PDBResidue prev){

		return dihedralAngle(prev.getC().loc, a.getN().loc, a.getCA().loc, a.getC().loc);
	}
	
	public static double omega(PDBResidue a,PDBResidue nex){

		return dihedralAngle(a.getCA().loc, a.getC().loc, nex.getN().loc, nex.getCA().loc);
	}

	public static double psi(PDBResidue a,PDBResidue nex){

		return dihedralAngle(a.getN().loc, a.getCA().loc, a.getC().loc, nex.getN().loc);
	}
	/**
	 * Missing とか削除した Residue の配列を返す。
	 * @param residue
	 * @return 
	 */
	public static ArrayList<PDBResidue> makeFilteredAA(ArrayList<PDBResidue> al){
		return makeFilteredAA(al,false);
	}
	public static ArrayList<PDBResidue> makeFilteredAA(ArrayList<PDBResidue> al,boolean removehoxt){
		ArrayList<PDBResidue> ret = new ArrayList<>();
		for(PDBResidue r:al){
			
			if(r.isLigand() || r.isMissing()){
				continue;
			}
			if(r.getC() == null || r.getCA() == null || r.getO() == null || r.getN() == null){
				continue;
			}
			if((!r.getName().equals("GLY")) && r.getCB() == null){
				continue;
			}
			//AltAtomについても削除する
			PDBResidue ff = r.getCopy();
			ret.add(ff);
			Iterator<PDBAtom> ite = ff.atoms.iterator();
			HashMap<String,ArrayList<PDBAtom>> atomcount = new HashMap<>();
			while(ite.hasNext()){
				PDBAtom a = ite.next();
				if(removehoxt){
					if(a.pdb_atom_code.indexOf("H") == 0){
						ite.remove();
						continue;
					}else if(a.pdb_atom_code.equals("OXT")){
						ite.remove();
						continue;
					}
				}
				
				if(!atomcount.containsKey(a.pdb_atom_code)){
					atomcount.put(a.pdb_atom_code,new ArrayList<PDBAtom>());
				}
				atomcount.get(a.pdb_atom_code).add(a);
			}
			//二つ以上同じ原子がある場合（AltPosition）、一つを残して削除する
			HashSet<PDBAtom> store = new HashSet<>();
			boolean chkflag = false;
			for(String cou:atomcount.keySet()){
				ArrayList<PDBAtom> aal = atomcount.get(cou);
				if(aal.size() == 1){
					store.add(aal.get(0));
				}else{
					Collections.sort(aal,new AtomListComparator());
					//for(PDBAtom ttt:aal){//適当なアルファベットを入れてチェックし確認した。
					//	System.out.println(ttt.getAltCode());
					//}
					store.add(aal.get(0));
					chkflag = true;
				}
			}
			if(chkflag){
				ite = ff.atoms.iterator();
				while(ite.hasNext()){
					PDBAtom a = ite.next();
					if(!store.contains(a)){
						ite.remove();
					}
				}
			}
		}
		
		return ret;
	}
	/**
	 * 
	 * @param residues
	 * @return 
	 */
	public static boolean[] checkChainBreak(ArrayList<PDBResidue> al){
		boolean[] ret = new boolean[al.size()];
		ret[al.size()-1] = true;
		for(int ii = 0;ii < al.size();ii++){
			ret[ii] = false;
			PDBResidue r = al.get(ii);
			
			if(r.isLigand() || r.isMissing() ||
				r.getC() == null || r.getCA() == null || r.getO() == null || r.getN() == null){
				if(ii != 0){
					ret[ii-1] = true;
				}
				ret[ii] = true;
				continue;
			}
			if(ii == al.size()-1){
				break;
			}
			
			
			PDBResidue nex = al.get(ii+1);
			int dff = nex.getResidueNumber() - r.getResidueNumber();
			if(dff != 1 && dff != 0){
				ret[ii] = true;
				continue;
			}
			PDBAtom aa = r.getC();
			PDBAtom bb = r.getN();
			if(aa == null || bb == null){
				ret[ii] = true;
				continue;
			}
			double dist = aa.distance(bb);
			//System.out.println(dist);
			if(dist > 3.0){//適当
				ret[ii] = true;
			}
		}
		return ret;
	}
	
	public static void calcPhiPsi_CNDist(String infile,String outfile){
		PDBData pdb = PDBData.loadPDBFile(infile);
		try{
			PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile,false),"UTF-8"))); 
			for(String s:pdb.chains.keySet()){

				PDBChain c = pdb.chains.get(s);
				ArrayList<PDBResidue> al = makeFilteredAA(c.residues);
				boolean[] chk = checkChainBreak(al);
				for(int ii = 1;ii < al.size()-1;ii++){
					if(!chk[ii] && !chk[ii-1] ){
						
						double x = phi(al.get(ii),al.get(ii-1));
						double y = psi(al.get(ii),al.get(ii+1));
						double ddist = al.get(ii-1).getC().distance(al.get(ii+1).getN());
						
						if(chk[ii-1]){
							pw.write(">>>");
						}
						pw.write(outfile+"\t"+s+"\t"+al.get(ii).getRepresentativeCode()+"\t"+x+"\t"+y+"\t"+ddist+"\n");
						//System.out.println(al.get(ii).getResidueNumber()+"\t"+x+"\t"+y);
					}else{
						//System.out.println("chainbreak");
					}

				}
			}
			pw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
	}
	
	
	/**
	 * a を底辺とする三角形の高さを返す。
	 * @param a
	 * @param b
	 * @param c 
	 */
	public static double calcTriangleHeight(double a,double b,double c){
		double p = (c*c+a*a-b*b)/2/a;
		double r = c*c-p*p;
		if(r <= 0){
			throw new RuntimeException("??? cannot make triangle.");
		}
		return Math.sqrt(r);
		
	}
	
	
	/**
	 * a, b, c を 3 辺とする三角形の面積を返す。
	 * @param a
	 * @param b
	 * @param c 
	 */	
	public static double calcTriangleSize(double a,double b,double c){
		double r = calcTriangleHeight(a,b,c);
		return r*a/2;	
	}
	
	/**
	 * bcda の四辺の長さと bc の交点から ad の交点に伸びる対角線の長さが分かっている四角形
	 * bc の交点から ad の交点へ延びる対角線と
	 * cd の交点から ab の交点へ延びる対角線の
	 * 交点について、bc の交点から ad の交点へ延びる対角線を bc からみて x:(1-x) に分割するとき
	 * x を返す。
	 * @param b
	 * @param c
	 * @param d
	 * @param a
	 * @param e
	 * @return 
	 */
	
	public static double calcCrossCenterRatio(double b,double c,double d,double a,double e){
		double r = calcTriangleHeight(b,e,a);
		double s = (b*b+e*e-a*a)/(2*b);
		double u = (c*c+e*e-d*d)/(2*e);
		
		double cos1 = s/e;
		double sin1 = r/e;
		
		double q = calcTriangleHeight(e,d,c);
		double sin2 = q/c;
		double cos2 = u/c;
		
		double cos3 = -1*(cos1*cos2-sin1*sin2);
		double sin3 = (cos1*sin2+sin1*cos2);
		double sizall = calcTriangleSize(b,e,a)+calcTriangleSize(e,c,d);
		double f  = (b+c*cos3)*(b+c*cos3)+(c*sin3)*(c*sin3);
		if(f <= 0){
			throw new RuntimeException("??? cannot calc quadlangle size.");
		}
		f = Math.sqrt(f);
		double sizf = calcTriangleSize(b,c,f);
		
		return sizf/sizall;
	}
	
	
	
	
	
	/**
	 * Omega が 180 になるような CA の場所を得る関数だが、
	 * 前の CA-C-N 面と合わせて 180 度回転した方が計算量少ないと思われるので中止
	 * @return 
	 */
	public static double[] getNextOmega180CA(FloatingResidue f1){
		double[] ret = new double[3];
		//2018/02/10
		double pratio = 0.4237513761543342;
		double dist_p_ca = 1.8184933638016645;
		//	1.9902838487766326;
	
		Point3D c = f1.getC().loc;
		Point3D ca = f1.getCA().loc;
		Point3D nn = f1.nextn.loc;
		
		double xx = (c.x-nn.x)*pratio+nn.x;
		double yy = (c.y-nn.y)*pratio+nn.y;
		double zz = (c.z-nn.z)*pratio+nn.z;
		
		double x2 = xx-ca.x;
		double y2 = yy-ca.y;
		double z2 = zz-ca.z;
		
		double llen = x2*x2+y2*y2+z2*z2;
		if(llen == 0){
			ret[0] = xx;
			ret[1] = yy;
			ret[2] = zz;
			return ret;
		}
		llen = Math.sqrt(llen);
		x2 /= llen;
		y2 /= llen;
		z2 /= llen;
		
		x2 *= dist_p_ca;
		y2 *= dist_p_ca;
		z2 *= dist_p_ca;
		
		x2 += xx;
		y2 += yy;
		z2 += zz;
		ret[0] = x2;
		ret[1] = y2;
		ret[2] = z2;
		
		System.out.println(ca.distance(new Point3D(x2,y2,z2)));
		return ret;
		
	}
	
	public static void main(String[] args){
		//main_ggg(args);
		//main_ee(args);
		//if(args[0].equals("-calcphipsi")){
		//	calcPhiPsi_CNDist(args[1],args[2]);
		//}
		
		System.out.println(calcCrossCenterRatio(1,2/(Math.sqrt(3)),1/(Math.sqrt(3))
				,1,1));
	}
	
	
	
	public static void main_ggg(String[] args){
		//内積とかを使う場合とどちらが速いか
		//あまり最適化してないのでフェアでないかもしれないが、内積とかを使ってない方が速かった。
		//For 文の最初に書くと遅くなる感じする。
		//ガベージコレクションとかの都合だろうか。
		double diff = 0;
		double diff2 = 0;
		for(int i = 0;i < 50;i++){
			long astart = System.currentTimeMillis();
			for(int ii = 0;ii < 10000;ii++){
				Point3D p1 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p3 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p5 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				ArrayList<Point3D> al = new ArrayList<>();
				al.add(new Point3D(p2));
				al.add(new Point3D(p5));
				adjustVector3D_Dot(p1,p2,p2b,p3,p4,p4b,al);
			}
			long aend = System.currentTimeMillis();
			System.out.println(aend-astart);

			long bstart = System.currentTimeMillis();
			for(int ii = 0;ii < 5000;ii++){
				
				Point3D p1 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p3 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p5 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				ArrayList<Point3D> al = new ArrayList<>();
				al.add(new Point3D(p2));
				al.add(new Point3D(p5));
				PepProcess.adjustVector3D(p1,p2,p2b,p3,p4,p4b,al);
			}
			long bend = System.currentTimeMillis();
			System.out.println(bend-bstart);
			diff += aend-astart-bend+bstart;
		
			long astart2 = System.currentTimeMillis();
			for(int ii = 0;ii < 10000;ii++){
				Point3D p1 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p2b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p3 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p4b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				Point3D p5 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
				ArrayList<Point3D> al = new ArrayList<>();
				al.add(new Point3D(p2));
				al.add(new Point3D(p5));
				adjustVector3D_Dot(p1,p2,p2b,p3,p4,p4b,al);
			}
			long aend2 = System.currentTimeMillis();
			
			diff2 += aend2-astart2-bend+bstart;
		}
		System.out.println("diff:"+diff);
		System.out.println("diff:"+diff2);
	}
	
	
	public static void main_ee(String[] args){
		//二つの方法で答えが異ならないかのチェック
		for(int ii = 0;ii < 10;ii++){
			
			//Point3D p1 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			//Point3D p2 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			//Point3D p2b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			Point3D p3 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			Point3D p4 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			Point3D p4b = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			Point3D p1 = new Point3D(0,0,0);
			Point3D p2 = new Point3D(0,1,0);
			Point3D p2b = new Point3D(1,0,0);
			//Point3D p3 = new Point3D(0,0,0);
			//Point3D p4 = new Point3D(1,0,0);
			//Point3D p4b = new Point3D(0,1,0);
			
			Point3D p5 = new Point3D(Math.random()*50-25,Math.random()*50-25,Math.random()*50-25);
			
			/*
			Point3D p1 = new Point3D(0,1,1);
			Point3D p2 = new Point3D(0,2,1);
			Point3D p2b = new Point3D(1,1,0);
			
			Point3D p3 = new Point3D(0,1,0);
			Point3D p4 = new Point3D(-1,1,0);
			Point3D p4b = new Point3D(0,2,0);
			
			Point3D p5 = new Point3D(1,1,0);
			*/
			
			ArrayList<Point3D> al = new ArrayList<>();
			al.add(new Point3D(p2));
			al.add(new Point3D(p2));
			al.add(new Point3D(p2b));
			
			ArrayList<Point3D> al2 = new ArrayList<>();
			al2.add(p2);
			al2.add(new Point3D(p2));
			al2.add(new Point3D(p2b));
			//adjustVector3D(p1,p2,p3,p4,al);
			PepProcess.adjustVector3D_Dot(p1,p2,p2b,p3,p4,p4b,al);
			adjustVector3D(p1,p2,p2b,p3,p4,p4b,al2);
			//System.out.println(al2.get(1).z);
			/*
			double ddist1 = p1.distance(p2);
			double ddist2 = p3.distance(p4);
			double ddist3 = al2.get(1).distance(p4);
			if(Math.abs(Math.abs(ddist1-ddist2)-ddist3) >0.001){
				System.out.println(Math.abs(ddist1-ddist2)-ddist3);
			}else{
				System.out.println("OK");
			}*/
			if(Math.abs(al.get(0).distance(al2.get(0))) > 0.00001){
				System.out.println(Math.abs(al.get(1).distance(al2.get(1))));
			}
			
			
		}
		
	}
	
	
	public static void main_(String[] args){
		PDBData pdb = PDBData.loadPDBFile("C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\candidates\\1T3Y.pdb");
		for(String s:pdb.chains.keySet()){
			PDBChain c = pdb.chains.get(s);
			//大体あっているがちょっと違う。
			//個別に出してくれるソフトと比べるべし
			ArrayList<PDBResidue> al = makeFilteredAA(c.residues);
			boolean[] chk = checkChainBreak(al);
			for(int ii = 1;ii < al.size()-1;ii++){
				if(!chk[ii]){
					double x = phi(al.get(ii),al.get(ii-1));
					double y = psi(al.get(ii),al.get(ii+1));
					System.out.println(al.get(ii).getResidueNumber()+"\t"+x+"\t"+y);
				}else{
					//System.out.println("chainbreak");
				}
				
			}
		}
		
	}
	
}
class AtomListComparator implements Comparator<PDBAtom>{
	@SuppressWarnings("unchecked")
	public int compare(PDBAtom arg1, PDBAtom arg2){
		return arg1.getAltCode().compareTo(arg2.getAltCode());
	}
	
}



class AtomDistance{
	PDBAtom target = null;
	double distance = 1000;
	double distance2 = 1000;
	double distance3 = 1000;
	AtomDistance(PDBAtom t,double d,double d2){
		distance = d;
		distance2 = d2;
		target = t; 
	}
}

 class AtomDistanceComparator implements Comparator<AtomDistance>{
	@SuppressWarnings("unchecked")
	public int compare(AtomDistance arg1, AtomDistance arg2){
		
		if(arg1.distance < arg2.distance ){
			return -1;
		}
		if(arg1.distance == arg2.distance ){
			return 0;
		}
			return 1;
	}
	
}
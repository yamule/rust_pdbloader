/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;

/**
 *
 * @author yamule
 */
public class LatticeWorld {
	AtomInLattice latticeBuffer[][][] = null;
	
	//あらかじめ格子に含まれる点の距離を計算しておけばどうかと思ったが、普通に距離を計算した方が早かった。
	//abs 等での判定をせねばならないからだろう。
	//ちなみに次に遅かったのが Math.pow
	//double[][][] distanceBuff = null;
	
	
	int bufferSize = 45;//バッファから外れる時に全原子をチェックしてバッファが更新される
	double resolution = 0.5;//Vertex 間の距離
	int offsetx = 0;
	int offsety = 0;
	int offsetz = 0;
	Location3DInt bufferCenter = new Location3DInt();
	ArrayList<AtomInLattice> atoms = new ArrayList<>();
	
	double repulsivePenalty = 1.0;
	double covalentBondScore = 1.0;
	double covalentBondPenalty = 1.0;
	double maxDist = 20;
	double[][][] distBuff = new double[40][40][40];
	
	
	LatticeWorld(int bs){
		bufferSize = bs;
		latticeBuffer = new AtomInLattice[bs][bs][bs];
	}
	
	public double calcDistance(int x,int y,int z,int x2,int y2,int z2){
		double d =(x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2);
		if(d <= 0){
			return 0;
		}
		return Math.sqrt(d)*resolution;
	}
	
	public void addAtom(AtomInLattice a){
		atoms.add(a);
	}
	
	public int[] getNearestVertexPosition(double xx,double yy,double zz){
		int[] ret = new int[3];
		ret[0] = (int)((xx+0.5)/resolution);
		ret[1] = (int)((yy+0.5)/resolution);
		ret[2] = (int)((zz+0.5)/resolution);
		return ret;
		
	}
	
	/*
	ロードする。
	全てのATOM に一意のString ID を与えて HASHMap で管理する
	ラティスのバーテックス一つに付き一個だけに収まるように変更する関数
	次にどの原子を移動させるか決定する関数。
	ランダムに変更する関数。
	ローカルスコア関数。ちゃんとラティスの距離を実際の距離に直すように。
	グローバルスコア関数。
	結果出力関数。
	*/
	
	
	
	
	
	
	
	
	/**
	 * BufferLattice 内でなく、グローバル Lattice 内の位置にある ATOM を得る。
	 * @param xx
	 * @param yy
	 * @param zz
	 * @return 
	 */
	public AtomInLattice getAtomIn(int xx,int yy, int zz){
		if(isInsideOfBuffer(xx,yy,zz)){
			return latticeBuffer[xx-offsetx][yy-offsety][zz-offsetz];
		}else{
			moveToCenter(xx,yy,zz);
			return latticeBuffer[xx-offsetx][yy-offsety][zz-offsetz];
		}
	}
	
	//シグモイド関数とか後で作ろう
	public double repulsivePenalty(AtomInLattice a,AtomInLattice b){
		double ddist = a.distance(b);
		if(a.repulsiveDist/2+b.repulsiveDist/2 > ddist){
			return -1.0;
		}
		return 0;
	}
	
	//シグモイド関数とか後で作ろう
	public double covalentBondScore(AtomInLattice a,AtomInLattice b){
		double ddist = a.distance(b);
		if((a.covalentBondDist/2+b.covalentBondDist/2)/2 > ddist){
			return -1.0;
		}else if((a.covalentBondDist/2+b.covalentBondDist/2)*1.2 < ddist){
			return (a.covalentBondDist/2+b.covalentBondDist/2)*1.2 - ddist;
		}
		return 1.0;
	}
	
	
	
	public ArrayList<AtomInLattice> getAtomsInDistance(int xx,int yy, int zz,int dist){
		ArrayList<AtomInLattice> ret = new ArrayList<>();
		for(int sx = xx-dist;sx <= xx+dist;sx++){
			for(int sy = yy-dist;sy <= yy+dist;sy++){
				for(int sz = zz-dist;sz <= zz+dist;sz++){
					if(isInsideOfBuffer(sx,sy,sz)){
						ret.add(latticeBuffer[sx-offsetx][sy-offsety][sz-offsetz]);
					}else{
						moveToCenter(xx,yy,zz);
						ret.add(latticeBuffer[sx-offsetx][sy-offsety][sz-offsetz]);
					}
				}
			}
		}
		
		return ret;
	}
	
	
	
	/**
	 * center から margin 分を調べて、バッファ内に収まるかチェックする 。
	 * 収まらない場合 False
	 * @param centerx
	 * @param centery
	 * @param centerz
	 * @param margin 
	 */
	public boolean checkBuffer(int centerx,int centery,int centerz,int margin){
		if(centerx-margin < offsetx || centery -margin < offsety || centerz -offsetz < offsetz ){
			//moveToCenter(centerx,centery,centerz);
			return false;
		}
		if(centerx+margin > offsetx+bufferSize-1 || centery +margin > offsety+bufferSize-1 || centerz +offsetz < offsetz ){
			//moveToCenter(centerx,centery,centerz);
			return false;
		}
		return true;
		
	}
	/**
	 * バッファバウンダリからどれくらいの距離があるか最も近い距離を返す。
	 * バウンダリ内にあるかないかは考えていない。
	 * @param xx
	 * @param yy
	 * @param zz
	 * @return 
	 */
	public int distFromBoundary(int xx,int yy,int zz){
		int ax = Math.abs(xx-bufferCenter.x);
		int ay = Math.abs(yy-bufferCenter.y);
		int az = Math.abs(zz-bufferCenter.z);
		return bufferSize - Math.max(ax,Math.max(ay,az));
	}
	public void moveToCenter(int cx,int cy,int cz){
		int bs = (bufferSize-1)/2;
		offsetx = cx-bs;
		offsety = cy-bs;
		offsetz = cz-bs;
		bufferCenter.set(cx,cy,cz);
		for(int ii = 0;ii < bufferSize;ii++){
			for(int jj = 0;jj < bufferSize;jj++){
				for(int kk = 0;kk < bufferSize;kk++){
					latticeBuffer[ii][jj][kk] = null;
				}
			}
		}
		for(AtomInLattice a:atoms){
			if(isInsideOfBuffer(a)){
				latticeBuffer[a.x-offsetx][a.y-offsety][a.z-offsetz] = a;
			}
		}
	}
	/**
	 * バッファ内に収まるか。
	 * 収まる場合 True
	 * @param a
	 * @return 
	 */
	public boolean isInsideOfBuffer(int xx,int yy,int zz){
		return (xx-offsetx)*(offsetx+bufferSize-1-xx) > 0 && 
		(yy-offsety)*(offsety+bufferSize-1-yy) > 0 && 
		(zz-offsetz)*(offsetz+bufferSize-1-zz) > 0 ;
	}
	public boolean isInsideOfBuffer(AtomInLattice a){
		return isInsideOfBuffer(a.x,a.y,a.z);
	}
	public static void main(String[] args){
		
	}
}

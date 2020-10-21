/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;
import java.util.HashSet;

/**
 *
 * @author yamule
 */
public class AtomInLattice {
	int x = 0;
	int y = 0;
	int z = 0;
	int index = 0;//atoms 内でのインデクス
	
		
	PDBAtom pdbAtom;
	ArrayList<DistanceConstraint> covalentBond = new ArrayList<>();
	
	LatticeWorld world;
	int repulsiveDist_Int = 2;
	double repulsiveDist = 3.0;
	//https://en.wikipedia.org/wiki/Van_der_Waals_radius 参照
	//面倒なので 1.5 A＊2 として考える 
	
	double covalentBondDist = 1.5;
	//https://ja.wikipedia.org/wiki/%E5%85%B1%E6%9C%89%E7%B5%90%E5%90%88%E5%8D%8A%E5%BE%84 参照
	//面倒なので全部 1.5 で考える
	AtomInLattice(LatticeWorld lw,PDBAtom d){
		x = (int)((d.loc.x+0.5)/lw.resolution);
		y = (int)((d.loc.y+0.5)/lw.resolution);
		z = (int)((d.loc.z+0.5)/lw.resolution);
		pdbAtom = d;
		world = lw;
		repulsiveDist_Int = (int)(repulsiveDist/lw.resolution+0.5);
		
	}
	public void setIndex(int i){
		index = i;
	}
	public double basicScore(){
		int dist = repulsiveDist_Int+1;
		int xx = x;
		int yy = y;
		int zz = z;
		double ret = 0;
		for(int sx = xx-dist;sx <= xx+dist;sx++){
			for(int sy = yy-dist;sy <= yy+dist;sy++){
				for(int sz = zz-dist;sz <= zz+dist;sz++){
					AtomInLattice a = null;
					if(world.isInsideOfBuffer(sx,sy,sz)){
						a = world.latticeBuffer[sx-world.offsetx][sy-world.offsety][sz-world.offsetz];
					}else{
						world.moveToCenter(sx,sy,sz);
						a = world.latticeBuffer[sx-world.offsetx][sy-world.offsety][sz-world.offsetz];
					}
					if(a != null){
						if(!covalentBond.contains(a)){
							ret += world.repulsivePenalty(a,this);
						}
					}
				}
			}
		}
		for(DistanceConstraint a: covalentBond){
			ret += a.calcScore();
		}
		
		return ret;
	}
	
	
	public void addCovalentBond(AtomInLattice a){
		//DistanceConstraint(AtomInLattice a1,AtomInLattice a2,double mi, double ma,double sc,double pe){
		covalentBond.add(new DistanceConstraint(this,a
		,covalentBondDist-world.resolution*1.01
		,covalentBondDist+world.resolution*1.01
		,world.covalentBondScore
		,world.covalentBondPenalty
		));
	}
	
	
	
	public double distance(AtomInLattice a){
		double ret = (x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)+(z-a.z)*(z-a.z);
		if(ret == 0){
			return ret;
		}
		return Math.sqrt(ret)*world.resolution;
	}
}

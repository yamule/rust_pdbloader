/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.ArrayList;

/**
 *
 * @author kimidori
 */
public class DockingRefine {
	ChainBuilder cb = new ChainBuilder();
	ArrayList<PDBResidue> allres = new ArrayList<>();
	ArrayList<ChainState> chains = new ArrayList<>();
	ArrayList<FloatingResidue> allres_f = new ArrayList<>();
	ArrayList<ArrayList<Point3D>> domains = new ArrayList<>();
	Threading th = new Threading();
	DockingRefine(ArrayList<ChainState> p1){
		
		//ArrayList<FloatingResidue> pp2 = ChainBuilder.changeToFloating(p2);
		for(ChainState c:p1){
			allres_f.addAll(c.residues);
		}
		chains.addAll(p1);
		allres.addAll(allres_f);
		cb.prepare(allres_f);
	}
	public static double sum(double[] d){
		double ret = 0;
		for(double dd:d){
			ret+=dd;
		}
		return ret;
	}
	public static Point3D center(ArrayList<Point3D> al){
		Point3D p = new Point3D(0,0,0);
		for(Point3D pp:al){
			p.x += pp.x;
			p.y += pp.y;
			p.z += pp.z;
		}
		
		p.x /= al.size();
		p.y /= al.size();
		p.z /= al.size();
		return p;
	}
	public void randomDock(int ite){
		double prevscore = sum(FuzzyDecisionTreeScoring_generator.calcResidueScores(allres,cb.scoring));
		
		for(int ii = 0;ii < ite;ii++){
			if(ii == 0){
				for(FloatingResidue f:allres_f){
					f.saveLoc();
				}
			}
			double deg = 3.0*Math.random() - 1.5;
			double rotx = 0.3*Math.random() - 0.15;
			double roty = 0.3*Math.random() - 0.15;
			double rotz = 0.3*Math.random() - 0.15;
			if(ii%2 == 0){
				deg *= 3;
				rotx *= 3;
				roty *= 3;
				rotz *= 3;
			}
			for(ArrayList<Point3D> p:domains){
				Point3D c = center(p);
				RandomDocker.moves(p,new Point3D(Math.random()*deg,Math.random()*deg,Math.random()*deg));
				RandomDocker.rotates(p,c,rotx,roty,rotz);
			}
			for(ChainState fal:chains){
				th.fitToPrev(fal.groups);
			}
			double s = sum(FuzzyDecisionTreeScoring_generator.calcResidueScores(allres,cb.scoring));
			
			
			if(s > prevscore){
				for(FloatingResidue f:allres_f){
					f.saveLoc();
				}
				prevscore = s;
			}else{
				for(FloatingResidue f:allres_f){
					f.restoreLoc();
				}
			}
		}
	}
	public static void main(String[] args){
		//ここから
		//Floating の Array を Threading からフラグとともに Import するように
	}
}

class ChainState{
	ArrayList<FloatingResidue> residues = new ArrayList<>();
	ArrayList<ArrayList<FloatingResidue>> groups = new ArrayList<>();
	ArrayList<ArrayList<Point3D>> groups_point = new ArrayList<>();//これを domains の代わりに渡すように	
}

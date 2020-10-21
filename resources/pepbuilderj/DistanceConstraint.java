/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

/**
 *
 * @author yamule
 */
public class DistanceConstraint {
	double min = 1.0;
	double max = 10.0;
	double penalty = 0;//must be positive
	double score = 0;
	AtomInLattice a = null;
	AtomInLattice b = null;
	DistanceConstraint(AtomInLattice a1,AtomInLattice a2,double mi, double ma,double sc,double pe){
		a = a1;
		b = a2;
		min = mi;
		max = ma;
		score = sc;
		penalty = pe;
	}
	public double calcScore(double precalc){
		double ret = 0;
		if((min-precalc)*(max-precalc) < 0){
			ret += score;
		}else{
			ret -= penalty;
		}
		return ret;
	}
	public double calcScore(){
		return calcScore(a.distance(b));
	}
	
}

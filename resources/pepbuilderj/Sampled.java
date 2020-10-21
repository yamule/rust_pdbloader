/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.Comparator;

/**
 *
 * @author kimidori
 */
public interface Sampled {
	public double getProb();
}

class SampleComparator implements Comparator<Sampled>{
	@SuppressWarnings("unchecked")
	public int compare(Sampled arg1, Sampled arg2){
		
		if(arg1.getProb() < arg2.getProb()){
			return -1;
		}
		if(arg1.getProb() == arg2.getProb()){
			return 0;
		}
			return 1;
	}
	
}

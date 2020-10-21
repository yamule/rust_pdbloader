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
 *//**
 * 元のインデックスと値を持っておいて、ソート後のランキングに使う。
 * @author kimidori
 */
public class VSorter{
	double val = 0;
	int index = 0;
	VSorter(double v,int i){
		index = i;
		val = v;
	}
}

class VComparator implements Comparator<VSorter>{
	@SuppressWarnings("unchecked")
	public int compare(VSorter arg1, VSorter arg2){
		if(arg1.val < arg2.val ){
			return -1;
		}
		if(arg1.val == arg2.val ){
			return 0;
		}
			return 1;
	}
}
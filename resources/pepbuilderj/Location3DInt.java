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
public class Location3DInt {
	int x = 0;
	int y = 0;
	int z = 0;
	Location3DInt(){
	}
	Location3DInt(int xx,int yy, int zz){
		x = xx;
		y = yy;
		z = zz;
	}
	public void set(int xx,int yy, int zz){
		x = xx;
		y = yy;
		z = zz;
	}
	
}

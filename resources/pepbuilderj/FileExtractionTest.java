/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.util.*;
import java.util.regex.*;
import java.io.*;

public class FileExtractionTest{
	
	public static ArrayList<File> getAllFiles(String dirname){
		return getAllFiles(new File(dirname),null);
	}
	public static ArrayList<File> getAllFiles(String dirname,String pat){
		return getAllFiles(new File(dirname),Pattern.compile(pat));
	}
	public static ArrayList<File> getAllFiles(File dir,Pattern pat){
		ArrayList<File> ret = new ArrayList<File>();
		if(dir.isDirectory()){
			File[] list = dir.listFiles();
			for(File f:list){
				if(f.getName().indexOf(".") != 0){//no system file
					if(f.isDirectory()){
						ret.addAll(getAllFiles(f,pat));
					}else{
						if(pat == null){
							ret.add(f);
						}else{
							if(pat.matcher(f.getName()).find()){
								ret.add(f);
							}
						}
					}
				}
			}
		}
		return ret;
	}
	
	public static void main(String[] args){
		ArrayList<File> files =getAllFiles("test","2");
		for(File f:files){
			System.out.println(f.getPath());
		}
	}
}
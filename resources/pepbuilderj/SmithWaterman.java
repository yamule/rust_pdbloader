/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pepbuilderj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * snippet for SmithWaterman alignment
 * Smith, Temple F., and Michael S. Waterman. 
 * "Identification of common molecular subsequences."
 * Journal of molecular biology 147.1 (1981): 195-197.
 * 
 * author: yamule (https://github.com/yamule/smithwaterman)
 * usage ===
 *	SmithWaterman sw = new SmithWaterman();
 *	SWResult res = sw.align("EEEEMDQNNSLPPYAGGTWRYII","IIIIMDQNNSPPYAQGGTWRYEE");
 *	System.out.println("score: "+res.score);
 *	System.out.println(SmithWaterman.listToString(res.qseq));
 *	System.out.println(SmithWaterman.listToString(res.sseq));
 *		
 * output ===
 * score: 73
 * ----EEEEMDQNNSLPPYA-GGTWRYII--
 * IIII----MDQNNS-PPYAQGGTWRY--EE
 * 
 * @author kimidori
 */
public class SmithWaterman {
	public static int TYPE_MATCH = 0;
	public static int TYPE_GAPINCOL = 1;
	public static int TYPE_GAPINROW = 2;
	SWCell[][] dpMat;
	ScoringMatrix smat = new BLOSUM62();
	double penalO = 10;
	double penalE = 0.5;
	int[] freezed_q = null;
	int[] freezed_s = null;
	boolean freezed = false;
	double freezingScore = 10000;
	
	
	public SWResult align(String qq,String ss){
		String qq2 = qq.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
		String ss2 = ss.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
		
		String name1 = "seq1";
		String name2= "seq2";
		Pattern npat = Pattern.compile(">([^\\s]+)");
		Matcher mat = npat.matcher(qq);
		if(mat.find()){
			name1 = mat.group(1);
		}
		 mat = npat.matcher(ss);
		if(mat.find()){
			name2 = mat.group(1);
		}
		
		Sequence q = new Sequence();
		Sequence s = new Sequence();
		q.name = name1;
		s.name = name2;
		q.add(qq2);
		s.add(ss2);
		smat.prepare(q.seq,s.seq);
		return align(q,s);
	}
	
	public SWResult align_with_PSSM(String qq,String pssm){
		
		PSSMScoring psss = PSSMScoring.loadPSSM(pssm);
		return align_with_PSSM(qq,psss);
	}
	
	public static double calcPSSMScore_nogap(ArrayList<Character> q,ArrayList<Character> s,PSSMScoring pssm){
		double ret = 0;
		int qpos = 0;
		int spos = 0;
		for(int ii = 0;ii < q.size();ii++){
			char qq = q.get(ii);
			char ss = s.get(ii);
			if(qq != '-' && ss != '-'){
				ret += pssm.getScoreOfAA(qq,spos);
			}
			if(qq != '-'){
				qpos++;
			}
			if(ss != '-'){
				spos++;
			}
		}
		return ret;
	}
	
	
	
	public SWResult align_with_PSSM(String qq,PSSMScoring psss){
		String qq2 = qq.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","").replaceAll("[^A-Za-z]","");
		String ss = psss.pssm.getFasta();
		smat = psss;
		String ss2 = ss.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
		
		String name1 = "seq1";
		String name2= "seq2";
		Pattern npat = Pattern.compile(">([^\\s]+)");
		Matcher mat = npat.matcher(qq);
		if(mat.find()){
			name1 = mat.group(1);
		}
		mat = npat.matcher(ss);
		if(mat.find()){
			name2 = mat.group(1);
		}
		
		Sequence q = new Sequence();
		Sequence s = new Sequence();
		q.name = name1;
		s.name = name2;
		q.add(qq2);
		s.add(ss2);
		smat.prepare(q.seq,s.seq);
		return align(q,s);
	}
	
	
	/*
	
	public void freezeAlignmentAt(int q, int s){
	}
	*/
	
	
	public void prepareDPMat(ArrayList<Character> q,ArrayList<Character> s){
		int w = q.size()+1;
		int h = s.size()+1;
		
		dpMat = new SWCell[w][h];
		
		for(int xx = 0;xx < w;xx++){
			for(int yy = 0;yy < h;yy++){
				dpMat[xx][yy] = new SWCell();
			}
		}
		for(int xx = 0;xx < w;xx++){
			
			for(int ii = 0;ii < 3;ii++){
				dpMat[xx][0].setScoreAt(ii,0);
			}
		}
		for(int yy = 0;yy < h;yy++){
			for(int ii = 0;ii < 3;ii++){
				dpMat[0][yy].setScoreAt(ii,0);
			}
		}
		for(int xx = 1;xx <  w;xx++){
			for(int yy = 1;yy <  h;yy++){
				double fm = dpMat[xx-1][yy-1].getScoreAt(TYPE_MATCH) + smat.getScore((xx-1), (yy-1));
				double fc = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINCOL) + smat.getScore((xx-1), (yy-1));
				double fr = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINROW) + smat.getScore((xx-1), (yy-1));
				SWCell currentcell = dpMat[xx][yy];
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				if(fm >= fr){
					if(fm >= fc){
						currentcell.setScoreAt(TYPE_MATCH,fm);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_MATCH);
					}else{
						currentcell.setScoreAt(TYPE_MATCH,fc);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
					}
				}else{
					if(fr > fc){
						currentcell.setScoreAt(TYPE_MATCH,fr);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINROW);
			
					}else{
						currentcell.setScoreAt(TYPE_MATCH,fc);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
					}
				}
				
				
				//Gap in col
				fm = dpMat[xx][yy-1].getScoreAt(TYPE_MATCH)-this.penalO;
				fc = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINCOL)-this.penalE;
				fr = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINROW)-this.penalO;// not used
				
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				
				if(fm >= fc){
					currentcell.setScoreAt(TYPE_GAPINCOL, fm);
					currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_MATCH);
				}else{
					currentcell.setScoreAt(TYPE_GAPINCOL, fc);
					currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_GAPINCOL);
				}
				
				//Gap in row
				fm = dpMat[xx-1][yy].getScoreAt(TYPE_MATCH)-this.penalO;
				fc = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINCOL)-this.penalO;// not used
				fr = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINROW)-this.penalE;
				
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				if(fm > fr){
					currentcell.setScoreAt(TYPE_GAPINROW, fm);
					currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_MATCH);
				}else{
					currentcell.setScoreAt(TYPE_GAPINROW, fr);
					currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_GAPINROW);
				}
				if(freezed){
					
					if(freezed_q.length != dpMat.length || freezed_s.length != dpMat[0].length){
						throw new RuntimeException("Please melt all freezed pair when you reuse the object.");
					}
					if(freezed_q[xx] > -1){
						if(freezed_q[xx] == yy){
							currentcell.setScoreAt(TYPE_MATCH, freezingScore+xx+yy);
							currentcell.setScoreAt(TYPE_GAPINCOL, -1);
							currentcell.setScoreAt(TYPE_GAPINROW, -1);
						}else{
							currentcell.setScoreAt(TYPE_MATCH, -1);
							currentcell.setScoreAt(TYPE_GAPINCOL,-1);
							currentcell.setScoreAt(TYPE_GAPINROW,-1);
						}
					}
					if(freezed_s[yy] > -1){
						if(freezed_s[yy] == xx){//必要ないはず
							currentcell.setScoreAt(TYPE_MATCH, freezingScore+xx+yy);
							currentcell.setScoreAt(TYPE_GAPINCOL, -1);
							currentcell.setScoreAt(TYPE_GAPINROW, -1);
						}else{
							currentcell.setScoreAt(TYPE_MATCH, -1);
							currentcell.setScoreAt(TYPE_GAPINCOL,-1);
							currentcell.setScoreAt(TYPE_GAPINROW,-1);
						}
					}
				}
				
				
			}
		}
		
	}
	public SWResult align(Sequence qq,Sequence ss){
		return align(qq,ss,false);
	}
	
	public void meltAll(){
		freezed = false;
		freezed_q = null;
		freezed_s = null;
	}
	public void prepareFreezing(int qlen,int slen,double fscore){
		freezed_q = new int[qlen+1];
		freezed_s = new int[slen+1];
		freezed = true;
		for(int ii = 0;ii < freezed_q.length;ii++){
			freezed_q[ii] = -1;
		}
		for(int ii = 0;ii < freezed_s.length;ii++){
			freezed_s[ii] = -1;
		}
		freezingScore = fscore;
	}
	public void freeze(int qpos,int spos){
		freezed_q[qpos+1] = spos+1;
		freezed_s[spos+1] = qpos+1;	
	}
	
	public void melt(int qpos){
		if(freezed_q[qpos+1] != -1){
			//System.out.println(qpos+";"+freezed_q[qpos+1]);
			freezed_s[freezed_q[qpos+1]] = -1;
			freezed_q[qpos+1] = -1;
		}
	}
	
	
	public SWResult align(Sequence qq,Sequence ss,boolean prepared){
		//Prepare letters for alignment.
		//Empty objects are removed and unknown letters are changed to 'X'.
		//Please see the filter method for exact process.
		ArrayList<Character> q = smat.filter(qq.seq);
		ArrayList<Character> s = smat.filter(ss.seq);
		if(!prepared){
			prepareDPMat(q,s);
		}
		
		int maxpos[] = getMaxScorePos();
		int cx = maxpos[0];
		int cy = maxpos[1];
		int sx = cx;
		int sy = cy;
		int lx = cx;
		int ly = cy;
		if(cx == -1){//very short sequences which do not have positive value match between two.
			return new SWResult();
		}
		
		
		
		//backtrackpart
		double cs = dpMat[cx][cy].getScoreAt(TYPE_MATCH);
		double maxscore = cs;
		int ct = TYPE_MATCH;
		int cpt = dpMat[cx][cy].getPrevTypeAt(TYPE_MATCH);
		ArrayList<Character> qchar = new ArrayList<>();
		ArrayList<Character> schar = new ArrayList<>();
		qchar.add(q.get(cx-1));
		schar.add(s.get(cy-1));
		while(cs > 0){
			if(ct == TYPE_MATCH){
				cx--;
				cy--;
			}else if(ct == TYPE_GAPINCOL){
				cy--;
			}else if(ct == TYPE_GAPINROW){
				cx--;
			}
			ct = cpt;
			cs = dpMat[cx][cy].getScoreAt(ct);
			cpt = dpMat[cx][cy].getPrevTypeAt(ct);
			if(cs <= 0){
				break;
			}
			if(ct == TYPE_MATCH){
				qchar.add(q.get(cx-1));
				schar.add(s.get(cy-1));
				lx = cx;
				ly = cy;
			}else if(ct == TYPE_GAPINCOL){
				qchar.add('-');
				schar.add(s.get(cy-1));
				ly = cy;
			}else if(ct == TYPE_GAPINROW){
				qchar.add(q.get(cx-1));
				schar.add('-');
				lx = cx;
			}
		}
		
		
		//add unaligned nterminal part
		for(int xx = lx-1;xx > 0;xx--){
			qchar.add(q.get(xx-1));
			schar.add('-');
		}
		for(int yy = ly-1;yy > 0;yy--){
			schar.add(s.get(yy-1));
			qchar.add('-');
		}
		
		Collections.reverse(qchar);
		Collections.reverse(schar);
		
		
		//add unaligned cterminal part
		for(int xx = sx+1;xx < dpMat.length;xx++){
			qchar.add(q.get(xx-1));
			schar.add('-');
		}
		for(int yy = sy+1;yy < dpMat[0].length;yy++){
			schar.add(s.get(yy-1));
			qchar.add('-');
		}
		
		
		return new SWResult(qchar,schar,maxscore);
		
	}
	
	
	
	
	/**
	 * Returns x y position of the cell which has maximum scores.
	 * If there are multiple cells with the same score, the cell which has lowest x and lowest y is selected.
	 * @return 
	 */
	public int[] getMaxScorePos(){
		int w = dpMat.length;
		int h = dpMat[0].length;
		int ret[] = new int[2];
		ret[0] = -1;
		ret[1] = -1;
		double maxscore = 0;
		for(int xx = 0;xx < w;xx++){
			for(int yy = 0;yy < h;yy++){
				for(int i = 0;i < 3;i++){
					double sc = dpMat[xx][yy].getScoreAt(i);
					if(maxscore < sc){
						maxscore = sc;
						ret[0] = xx;
						ret[1] = yy;
					}
				}
			}	
		}
		return ret;
	}
	
	
	
	/**
	 * Changes ArrayList<Character> to String.
	 * @param al
	 * @return 
	 */
	public static String listToString(ArrayList<Character> al){
		StringBuffer ret = new StringBuffer();
		for(Character c:al){
			ret.append(c);
		}
		return ret.toString();
	}
	
	public static ArrayList<Character> toCharList(String s){
		ArrayList<Character> ret = new ArrayList<>();
		for(char c:s.toCharArray()){
			ret.add(c);
		}
		return ret;
	}
	
	public static void main(String[] args){
		if(args.length > 1){
			for(int ii = 0;ii < args.length;ii++){
				System.out.println(ii+" "+args[ii]);
			}
			if(args[1].equals("-pssm")){
				String testfas = args[0];
				String testpssm = args[2];
				SmithWaterman sw = new SmithWaterman();
				ArrayList<Sequence> seq1 = Sequence.loadFasta(testfas);
				SWResult res = sw.align_with_PSSM(seq1.get(0).getSeq(),testpssm);
				System.out.println("score: "+res.score);
				System.out.println(SmithWaterman.listToString(res.qseq));
				System.out.println(SmithWaterman.listToString(res.sseq));
			}else{
				ArrayList<Sequence> seq1 = Sequence.loadFasta(args[0]);
				ArrayList<Sequence> seq2 = Sequence.loadFasta(args[1]);

				SmithWaterman sw = new SmithWaterman();
				SWResult res = sw.align(seq1.get(0),seq2.get(0));
				System.out.println("# score: "+res.score);
				System.out.println(">"+seq1.get(0).name);
				System.out.println(listToString(res.qseq));
				System.out.println(">"+seq2.get(0).name);
				System.out.println(listToString(res.sseq));
			}
		}else{
			System.err.println("usage: java SmithWaterman <fastafile path> <fastafile path>");
			System.err.println("usage: java SmithWaterman <fastafile path> -pssm <asciipssmfile path>");
			
			/*
			// for debug
			SmithWaterman sw = new SmithWaterman();
			SWResult res = sw.align("MTPPPPGRAAPSAPRARVPGPPARLGLPLRLRLLLLLWAAAASAQGHLRSGPRIFAVWKG","PGPGNPSPMSLSPAWPGHPDQPLPREQMTSPAPPRIITSATADPEGTETALAGDTSDGLA");
			System.out.println("score: "+res.score);
			System.out.println(SmithWaterman.listToString(res.qseq));
			System.out.println(SmithWaterman.listToString(res.sseq));
			*/
			/*
			String testpssm = "C:\\dummy\\private\\bioo\\seminar_mine\\20170830\\hemo_test.pssm";
			String testfas = "C:\\dummy\\private\\bioo\\seminar_mine\\20170830\\hemoglobin_b.fas";
			SmithWaterman sw = new SmithWaterman();
			ArrayList<Sequence> seq1 = Sequence.loadFasta(testfas);
			SWResult res = sw.align_with_PSSM(seq1.get(0).getSeq(),testpssm);
			System.out.println("score: "+res.score);
			System.out.println(SmithWaterman.listToString(res.qseq));
			System.out.println(SmithWaterman.listToString(res.sseq));
			*/
		}
		/*
		BLOSUM62 bl = new BLOSUM62();
		SmithWaterman sw = new SmithWaterman();
		ArrayList<Sequence> seqs = Sequence.loadFasta("C:\\Users\\kimidori\\Desktop\\TBM\\sw\\testfas.txt");
		for(Sequence s:seqs){
			System.out.println(">"+s.name+" "+s.desc);
			for(Character c:s.seq){
				System.out.print(String.valueOf(c));
			}
			System.out.println();
		}
		SWResult res = sw.align(seqs.get(0),seqs.get(1));
		System.out.println(res.score);
		System.out.println(seqs.get(0).name);
		System.out.println(listToString(res.qseq));
		System.out.println(seqs.get(1).name);
		System.out.println(listToString(res.sseq));
		*/
	}
}



/**
 * Stores result of SmithWaterman Alignment.
 * @author kimidori
 */
class SWResult{
	ArrayList<Character> qseq = new ArrayList<>();
	ArrayList<Character> sseq = new ArrayList<>();
	double score = -1;
	SWResult(){
	}
	SWResult(ArrayList<Character> qq,ArrayList<Character> ss,double s){
		qseq.addAll(qq);
		sseq.addAll(ss);
		score = s;
	}
}



class SWCell{
	double score[] = new double[3];
	int prevType[] = new int[3];
	public void setScoreAt(int type,double sc){
		score[type] = sc;
	}
	public double getScoreAt(int type){
		return score[type];
	}
	public int getPrevTypeAt(int type){
		return prevType[type];
	}
	public void setPrevTypeAt(int ctype,int ptype){
		prevType[ctype] = ptype;
	}
	
	
}

class BLOSUM62 implements ScoringMatrix{
	
	
	String lines[] = {
		//https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ",
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ",
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ",
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ",
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ",
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ",
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ",
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ",
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ",
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ",
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ",
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ",
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ",
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ",
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ",
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ",
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ",
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ",
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ",
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ",
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ",
"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 "
	};
	
	
	HashSet<Character> acceptable = new HashSet<>();
	int[][] scoreMat = new int[128][128];
	ArrayList<Character> seqA = null;
	ArrayList<Character> seqB = null;
	
	BLOSUM62(){
		for(int ii = 0;ii < 128;ii++){
			for(int jj = 0;jj < 128;jj++){
				scoreMat[ii][jj] = 10000;
			}	
		}
		ArrayList<String> head = splitWithSpace(lines[0]);
		for(int ii = 1;ii < lines.length;ii++){
			ArrayList<String> pt = splitWithSpace(lines[ii]);
			String s = pt.get(0);
			char sc = s.charAt(0);
			acceptable.add(sc);
			for(int jj = 1;jj < pt.size();jj++){
				String q = head.get(jj-1);
				char qc = q.charAt(0);
				acceptable.add(qc);
				if(scoreMat[sc][qc] < 100){
					System.out.println(sc+"-"+qc+" duplicate pair?");
				}
				scoreMat[sc][qc] = Integer.parseInt(pt.get(jj));
			}
		}
		/*
		ArrayList<Character> aa = new ArrayList<>(acceptable);
		for(Character k:aa){
			for(Character q:aa){
				int sc3 = scoreMat[q][k];
				int sc4 = scoreMat[k][q];
				if(sc3 == sc4){
					System.out.println(k+"-"+q+" OK "+sc3);
				}else{
					System.out.println(k+q+" is different, "+sc3+","+sc4+",");
				}
			}
		}
		*/
		
	}
	public void prepare(ArrayList<Character> sa,ArrayList<Character> sb){
		seqA = new ArrayList<>(sa);
		seqB = new ArrayList<>(sb);
		
	}
	
	/**
	 * Change letters which are not compatible with this matrix.
	 * @param al 
	 */
	public ArrayList<Character> filter(ArrayList<Character> al){
		ArrayList<Character> ret = new ArrayList<>();
		for(int ii = 0;ii < al.size();ii++){
			Character a = al.get(ii);
			if(a == null || a == Character.MIN_VALUE){
			}else{
				if(a.equals('-') || a.equals('.')){
				}else{
					if(acceptable.contains(a)){
						ret.add(a);
					}else{
						ret.add('X');
					}
				}
			}
		}
		return ret;
	}
	
	
	
	//For discrepancy between open java & oracle java
	public ArrayList<String> splitWithSpace(String str){
		ArrayList<String> ret = new ArrayList<>();
		String head[] = str.replaceAll("^[\\s]+","").replaceAll("[\\s]+$","").split("[\\s]+");
		for(String s:head){
			if(s.length() > 0){
				ret.add(s);
			}
		}
		return ret;
	}
	
	public double getScore(int qpos,int spos){
		
		return scoreMat[seqA.get(qpos)][seqB.get(spos)];
	}
}



interface ScoringMatrix{
	public double getScore(int q,int s);
	public void prepare(ArrayList<Character> s1,ArrayList<Character> s2);
	public ArrayList<Character> filter(ArrayList<Character> al);
}

class Sequence{
	String name;
	String desc;
	ArrayList<Character> seq = new ArrayList<Character>();
	
	
	public static ArrayList<Sequence> PDBToSeq(PDBData d){
		ArrayList<Sequence> ret = new ArrayList<>();
		for(String ss:d.chains.keySet()){
			PDBChain cc = d.chains.get(ss);
			Sequence s = new Sequence();
			s.desc = "";
			s.name = ss;
			String e = cc.getSequence();
			for(char c:e.toCharArray()){
				s.seq.add(c);
			}
			ret.add(s);
		}
		return ret;
	}
	public static ArrayList<Sequence> loadFasta(String filename){
		ArrayList<Sequence> ret = new ArrayList<>();
		try{
			BufferedReader br = new  BufferedReader(new FileReader(new File(filename)));
			String ln = null;
			Sequence currentseq = new Sequence();
			ret.add(currentseq);
			Pattern pat1 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]*)");
			Pattern pat2 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]+)[\\s]+([^\\r\\n]+)");
			while((ln = br.readLine()) != null){
				Matcher mat = pat1.matcher(ln);
				if(mat.find()){
					String n = mat.group(1);
					String d = "";
					Matcher mat2 = pat2.matcher(ln);
					if(mat2.find()){
						d = mat2.group(2);
					}
					if(ret.size() == 1 && currentseq.seq.size() == 0){
					}else{
						currentseq = new Sequence();
						ret.add(currentseq);
					}
					currentseq.name = n;
					currentseq.desc = d;
				}else{
					currentseq.add(ln);
				}
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	public void add(String s){
		String[] pt = s.replaceAll("[\\r\\n]","").split("");
		for(String pp:pt){
			if(pp.length() == 1){
				seq.add(pp.charAt(0));
			}else{
				if(pp.length() == 0){
					
				}else{
					throw new RuntimeException("java implementation error.");
				}
			}
		}
	}
	public String getSeq(){
		StringBuffer sb = new StringBuffer();
		for(Character c:seq){
			sb.append(c);
		}
		return sb.toString();
	}
}

class PSSMScoring implements ScoringMatrix{
	HashSet<Character> acceptable = new HashSet<>();
	PSSMData pssm = null;
	ArrayList<Character> seqA = null;
	char colhead[] = "ARNDCQEGHILKMFPSTWYV".toCharArray();
	PSSMScoring(){
		for(char c:colhead){
			acceptable.add(c);
		}
		acceptable.add('X');
	}
	
	public static  PSSMScoring merge(ArrayList<PSSMScoring> al){
		PSSMScoring ret = new PSSMScoring();
		ret.acceptable.addAll(al.get(0).acceptable);
		ArrayList<PSSMData> pal = new ArrayList<>();
		for(PSSMScoring p:al){
			pal.add(p.pssm);
			if(p.seqA != null){
				ret.seqA.addAll(p.seqA);
			}
		}
		ret.pssm = PSSMData.merge(pal);
		return ret;
	}
	public ArrayList<PSSMScoring> splitAt(int i){
		ArrayList<PSSMScoring> ret = new ArrayList<>();
		ArrayList<PSSMData> sp = pssm.splitAt(i);
		PSSMScoring r1 = new PSSMScoring();
		PSSMScoring r2 = new PSSMScoring();
		r1.acceptable.addAll(this.acceptable);
		r2.acceptable.addAll(this.acceptable);
		r1.pssm = sp.get(0);
		r2.pssm = sp.get(1);
		ret.add(r1);
		ret.add(r2);
		return ret;
	}
	
	
	public void prepare(ArrayList<Character> q,ArrayList<Character> s){
		//setPSSM(p);
		seqA = q;
	}
	
	public double getScore(int qpos,int spos){
		if(pssm.char_index_map[seqA.get(qpos)] == -1){
			System.out.println(seqA.get(qpos));
		}
		return pssm.scores.get(spos).get(pssm.char_index_map[seqA.get(qpos)]);
	}
	
	public double getScoreOfAA(char qq,int spos){
		return pssm.scores.get(spos).get(pssm.char_index_map[qq]);
	}
	
	public double getScoreAt(int posrow,char cc){
		return pssm.scores.get(posrow).get(pssm.char_index_map[cc]);
	}
	
	public static PSSMScoring loadPSSM(String filename){
		PSSMScoring ret = new PSSMScoring();
		ret.setPSSM(PSSMData.load(filename));
		return ret;
	}
	
	public static PSSMScoring prepare(PSSMData dat){
		PSSMScoring ret = new PSSMScoring();
		ret.setPSSM(dat);
		return ret;
	}
	
	/**
	 * 
	 * @param pssmm 
	 */
	public void setPSSM(PSSMData pssmm){
		if(this.pssm == pssmm){
			return;
		}
		this.pssm = pssmm;
		if(pssmm.char_index_map['X'] == -1){
			//異常なアミノ酸が来た場合の値を入れる
			int lastindex = 0;
			for(int ii = 0; ii < pssmm.char_index_map.length;ii++){
				lastindex = Math.max(lastindex, pssmm.char_index_map[ii]);
			}
			int xindex  = lastindex+1;
			//System.err.println("score of X char was set on index "+xindex);
			
			for(int ii = 0;ii < pssmm.scores.size();ii++){
				ArrayList<Double> dd = pssmm.scores.get(ii);
				double sum = 0;
				double count = 0;
				for(Character cc:acceptable){
					if(cc == 'X'){
					}else{
						sum += dd.get(pssmm.char_index_map[cc]);
						count += 1;
					}
				}
				double av = sum/count;
				if(dd.size() > xindex){
					dd.set(xindex,av);
				}else{
					if(dd.size() != xindex){
						throw new RuntimeException("?? pssm ascii format error?");
					}
					dd.add(av);
				}
			}
		}
		
		
		
	}
	
	/**
	 * Change letters which are not compatible with this matrix.
	 * @param al 
	 */
	public ArrayList<Character> filter(ArrayList<Character> al){
		ArrayList<Character> ret = new ArrayList<>();
		for(int ii = 0;ii < al.size();ii++){
			Character a = al.get(ii);
			if(a == null || a == Character.MIN_VALUE){
			}else{
				if(a.equals('-') || a.equals('.')){
				}else{
					if(acceptable.contains(a)){
						ret.add(a);
					}else{
						ret.add('X');
					}
				}
			}
		}
		return ret;
	}
}
class PSSMData{
	String fileName = "";
	ArrayList<Character> seq = new ArrayList<>();
	//String fastaSeq = "";
	int colNum = 0;
	ArrayList<ArrayList<Double>> scores = new ArrayList<>();
	int[] char_index_map = new int[128];//Char に対応する文字が scores に詰められている配列上で何番目にあるか。
	PSSMData(){
		//char colhead[] = "ARNDCQEGHILKMFPSTWYVX".toCharArray();
		char colhead[] = "ARNDCQEGHILKMFPSTWYV".toCharArray();
		colNum = colhead.length;
		for(int ii = 0;ii < 128;ii++){
			char_index_map[ii] = -1;
		}
		for(int ii = 0;ii < colhead.length;ii++){
			char_index_map[colhead[ii]] = ii;
		}
	}
	public String getFasta(){
		StringBuffer ret = new StringBuffer();
		for(Character c:seq){
			ret.append(c);
			
		}
		return ret.toString();
	}
	
	public static PSSMData merge(ArrayList<PSSMData> al){
		PSSMData ret = new PSSMData();
		PSSMData f = al.get(0);
		ret.fileName = f.fileName;
		ret.colNum = f.colNum;
		ret.char_index_map = new int[f.char_index_map.length];
        System.arraycopy(f.char_index_map,0,ret.char_index_map,0,ret.char_index_map.length);
        //元は同じ PSSM からきていることを前提にしており整合性をチェックしていない。
		for(PSSMData p:al){
			ret.scores.addAll(p.scores);
			ret.seq.addAll(p.seq);
		}
		
		return ret;
	}
	
	public static double calcLambda(ArrayList<ArrayList<Double>> al ){
		
		return 1.0;
	}
	
	
	public ArrayList<PSSMData> splitAt(int i){
		ArrayList<PSSMData> ret = new ArrayList<>();
		PSSMData r1 = new PSSMData();
		PSSMData r2 = new PSSMData();
		r1.fileName = this.fileName;
		r2.fileName = this.fileName;
		r1.colNum = this.colNum;
		r2.colNum = this.colNum;
		
        System.arraycopy(this.char_index_map,0,r1.char_index_map,0,r1.char_index_map.length);
        System.arraycopy(this.char_index_map,0,r2.char_index_map,0,r2.char_index_map.length);
		for(int ii = 0;ii < scores.size();ii++){
			if(ii < i){
				r1.scores.add(scores.get(ii));
				r1.seq.add(this.seq.get(ii));
			}else{
				r2.scores.add(scores.get(ii));
				r2.seq.add(this.seq.get(ii));
			}
		}
		ret.add(r1);
		ret.add(r2);
		return ret;
	}
	public void addResidue(char r,HashMap<Character,Double> score){
		ArrayList<Double> dal = new ArrayList<>();
		for(int ii = 0;ii < colNum;ii++){
			dal.add(0.0);
		}
		for(Character c:score.keySet()){
			if(char_index_map[c] > -1){
				dal.set(char_index_map[c],score.get(c));
			}
		}
		seq.add(r);
		scores.add(dal);
	}
	
	public static PSSMData load(String filename){
		PSSMData ret = new PSSMData();
		ret.fileName = filename;
		try{
			Pattern headpat = Pattern.compile("A[\\s]+R[\\s]+N[\\s]+D[\\s]+C[\\s]+Q[\\s]+E[\\s]+G[\\s]+H[\\s]+I[\\s]+L[\\s]+K[\\s]+M[\\s]+F[\\s]+P[\\s]+S[\\s]+T[\\s]+W[\\s]+Y[\\s]+V");
			BufferedReader br = new  BufferedReader(new FileReader(new File(filename)));
			String ln = null;
			
			Pattern pat1 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]*)");
			Pattern pat2 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]+)[\\s]+([^\\r\\n]+)");
			boolean headflag = false;
			while((ln = br.readLine()) != null){
				if(headpat.matcher(ln).find()){
					headflag = true;
					break;
				}
			}
			if(!headflag){
				throw new RuntimeException("Cannot find header. "+headpat.pattern());
			}
			while((ln = br.readLine()) != null){
				
				String[] ptt = ln.toUpperCase().replaceAll("^[\\s]+","").replaceAll("[\r\n]","").split("[\\s]+");
				if(ptt.length < 2){
					continue;
				}
				if(ptt[0].equals("K") && ptt[1].equals("LAMBDA")){
					break;
				}
				ArrayList<Double> dal = new ArrayList<>();
				ret.scores.add(dal);
				for(int ii = 2;ii < ptt.length;ii++){
					try{
						dal.add(Double.parseDouble(ptt[ii]));
					}catch(Exception exx){
						exx.printStackTrace();
						break;
					}
				}
				
				ret.seq.add(ptt[1].toCharArray()[0]);
			}
			StringBuffer sb = new StringBuffer();
			for(Character c:ret.seq){
				sb.append(c);
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}

}
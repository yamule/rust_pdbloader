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
//blast+2.7.1
//c++/src/algo/blast/core/blast_stat.c
//c++/src/algo/blast/core/blast_psi_priv.c
//c++/src/algo/blast/core/ncbi_math.c
public class PSSMCalc {
	public static double kPSINearIdentical = 0.94;
	public static double kPSIIdentical = 1.0;
	public static  int kQueryIndex = 0;
	public static double kEpsilon = 0.0001;
	public static int kPSIScaleFactor = 200;
	public static double kPositScalingPercent = 0.05;
	public static int kPositScalingNumIterations = 10;
	public static double BLAST_KARLIN_LAMBDA0_DEFAULT = 0.5;
	
	public static final double  BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT = 1.0e-5;
	public static final int BLAST_KARLIN_LAMBDA_ITER_DEFAULT=17;

	public static int BLAST_SCORE_MIN = -32768;
	public static int BLAST_SCORE_MAX = 32768;
	public static double[] background_frequency  = new double[128];
	public static int[] aa_to_index = new int[128];//Char に対応する文字が scores に詰められている配列上で何番目にあるか。
	public static char[] index_to_aa = new char[128];
	public static int nb_aatype = 0;
	public static char[] aa_type;
	static{
		//char colhead[] = "ARNDCQEGHILKMFPSTWYVX".toCharArray();
		aa_type = "ARNDCQEGHILKMFPSTWYV".toCharArray();
		nb_aatype = aa_type.length;
		for(int ii = 0;ii < 128;ii++){
			aa_to_index[ii] = -1;
			background_frequency[ii] = -1;
		}
		for(int ii = 0;ii < aa_type.length;ii++){
			aa_to_index[aa_type[ii]] = ii;
			index_to_aa[ii] = aa_type[ii];
		}
		
		//blast_stat.c
		/* amino acid background frequencies from Robinson and Robinson */
		double[] _background_freq  = new double[128];
		_background_freq['A']=78.05;
		_background_freq['C']=19.25;
		_background_freq['D']=53.64;
		_background_freq['E']=62.95;
		_background_freq['F']=38.56;
		_background_freq['G']=73.77;
		_background_freq['H']=21.99;
		_background_freq['I']=51.42;
		_background_freq['K']=57.44;
		_background_freq['L']=90.19;
		_background_freq['M']=22.43;
		_background_freq['N']=44.87;
		_background_freq['P']=52.03;
		_background_freq['Q']=42.64;
		_background_freq['R']=51.29;
		_background_freq['S']=71.20;
		_background_freq['T']=58.41;
		_background_freq['V']=64.41;
		_background_freq['W']=13.30;
		_background_freq['Y']=32.16;
		
		for(int ii = 0;ii < aa_type.length;ii++){
			background_frequency[aa_to_index[aa_type[ii]]] = _background_freq[aa_type[ii]]/1000.0; 
		}
		
	}
	
	public static int[][] convertFrequencyToPSSM(double[][] f,int[] aatoindex ){
		PSSMCalc cal = new PSSMCalc();
		BlastScoreBlk scoreblk = new BlastScoreBlk();
		scoreblk.std_prob =PSSMCalc.background_frequency;
		
		double[][] freq = new double[f.length][20];
		for(int ii = 0;ii < f.length;ii++){
			double d[] = f[ii];
			double sum = 0;
			for(char cc:aa_type){
				
				freq[ii][PSSMCalc.aa_to_index[cc]] = d[aatoindex[cc]];
				sum += d[aatoindex[cc]]; 
			}
			if(Math.abs(1.0-sum) > 0.01){
				for(char cc:aa_type){
					freq[ii][PSSMCalc.aa_to_index[cc]] /= sum;
				}	
			}
		}
		freq = cal.addPseudoFrequency(freq,PSSMCalc.background_frequency);
		PSIInternalPssmData ps  = new PSIInternalPssmData(freq,scoreblk.std_prob,scoreblk.ideal_lambda);
		cal.PSIScaleMatrix(ps,scoreblk);
		return ps.pssm;
	}
	
	
	
	
	
	public double NlmKarlinLambdaNR(double[] probs
			,int d, int low, int high, double lambda0,
                  double tolx, int itmax, int maxNewton){
		double x0, x, a = 0, b = 1;
		double f = 4;  
		int isNewton = 0;

		x0 = Math.exp( -lambda0 );
		x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;

		for(int k = 0; k < itmax; k++ ) { /* all iteration indices k */
			double g, fold = f;
			int wasNewton = isNewton; 

			isNewton  = 0;            
			g = 0;
			f = probs[low-low];
			//ホーナー法で [0-1]に存在する解を探索する
			for(int i = low + d; i < 0; i += d ) {
				g = x * g + f;
				f = f * x + probs[i-low];
			}
			g = x * g + f;
			f = f * x + probs[0-low] - 1;//この -1 は何
			for(int i = d; i <= high; i += d ) {
				g = x * g + f;
				f = f * x + probs[i-low];
			}

			if( f > 0 ) {
				a = x; /* move the left endpoint */
			} else if( f < 0 ) {
				b = x; /* move the right endpoint */
			} else { /* f == 0 */
				break; /* x is an exact solution */
			}
			if( b - a < 2 * a * ( 1 - b ) * tolx ) {
				x = (a + b) / 2;
				break;
			}

			if( k >= maxNewton ||( wasNewton == 1 && Math.abs( f ) > .9 * Math.abs(fold) ) ||g >= 0){
				x = (a + b)/2;
			} else {
				double p = - f/g;
				double y = x + p;
				if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
					x = (a + b)/2;
				} else { /* The proposed iterate is in (a,b). Accept it. */
					isNewton = 1;
					x = y;
					if( Math.abs( p ) < tolx * x * (1-x) ){
						break;
					}
				}
			} /* else try a Newton step. */
		}
		return -Math.log(x)/d;
	}
	
	//最大公約数を出すものだが大体 1 になりそう
	public static int BLAST_Gcd(int a, int b){
	   int   c;
	   b = Math.abs(b);
	   if (b > a){
			c=a;
			a=b;
			b=c;
	   }

	   while (b != 0) {
		  c = a%b;
		  a = b;
		  b = c;
	   }
	   return a;
	}
	
	public void updateLambdaK(int[][] pssm, BlastScoreBlk sbp){
		Blast_ScoreFreq score_freqs = Blast_ScoreFreq.PSIComputeScoreProbabilities(pssm,sbp.std_prob);
		sbp.kbp_lambda = Blast_KarlinLambdaNR(score_freqs, BLAST_KARLIN_LAMBDA0_DEFAULT);
		if (sbp.kbp_lambda < 0.){
			throw new RuntimeException("??? error in lambda calculation.");
		}
	}

	
	public double Blast_KarlinLambdaNR(Blast_ScoreFreq sfp, double initialLambdaGuess){
	   int  low;        /* Lowest score (must be negative)  */
	   int  high;       /* Highest score (must be positive) */
	   int   d;
	   double[]  sprob;
	   double   returnValue;

	   low = sfp.obs_min;
	   high = sfp.obs_max;
	   if (sfp.score_average >= 0.) {   /* Expected score must be negative */
		  return -1.0;
	   }
	   sprob = sfp.sprob;
	   /* Find greatest common divisor of all scores */
	   //for (i = 1, d = -low; i <= high-low && d > 1; ++i) {
	   //   if (sprob[i+low] != 0.0) {
	   //      d = BLAST_Gcd(d, i);
	   //   }
	   //}
	   d = 1;
	   returnValue =NlmKarlinLambdaNR( sprob, d, low, high,
							   initialLambdaGuess,
							   BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
						 20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT);
	   return returnValue;
	}



	public boolean PSIScaleMatrix(PSIInternalPssmData internal_pssm,	BlastScoreBlk sbp){
		boolean first_time = true;
		int index = 0;	 /* loop index */
		int[][] scaled_pssm = null;//kPSIScaleFactor が掛けられた値
		int[][] pssm = null;
		double factor;
		double factor_low = 1.0;
		double factor_high = 1.0;
		double ideal_lambda = 0.0;	  /* ideal value of ungapped lambda for
										   underlying scoring matrix */
		double new_lambda = 0.0;		/* Karlin-Altschul parameter calculated 
										   from scaled matrix*/
		boolean too_high = true;
		scaled_pssm = internal_pssm.scaled_pssm;
		pssm = internal_pssm.pssm;
		ideal_lambda = sbp.ideal_lambda;

		factor = 1.0;
		for ( ; ; ) {
			int i = 0;
			int j = 0;

			for (i = 0; i < internal_pssm.ncol(); i++) {
				for (j = 0; j < internal_pssm.nrow(); j++) {
					if (scaled_pssm[i][j] != BLAST_SCORE_MIN) {
						pssm[i][j] = BLAST_Nint(factor*scaled_pssm[i][j]/kPSIScaleFactor);
					} else {
						pssm[i][j] = BLAST_SCORE_MIN;
					}
				}
			}
			updateLambdaK(pssm,sbp);

			new_lambda = sbp.kbp_lambda;

			if (new_lambda > ideal_lambda) {
				if (first_time) {
					factor_high = 1.0 + kPositScalingPercent;
					factor = factor_high;
					factor_low = 1.0;
					too_high = true;
					first_time = false;
				} else {
					if (too_high == false) {
						break;
					}
					factor_high += (factor_high - 1.0);
					factor = factor_high;
				}
			} else if (new_lambda > 0) {
				if (first_time) {
					factor_high = 1.0;
					factor_low = 1.0 - kPositScalingPercent;
					factor = factor_low;
					too_high = false;
					first_time = false;
				} else {
					if (too_high == true) {
						break;
					}
					factor_low += (factor_low - 1.0);
					factor = factor_low;
				}
			} else {
				//return PSIERR_POSITIVEAVGSCORE;
				return false;
			}
		}

		/* Binary search for kPositScalingNumIterations times */
		for (index = 0; index < kPositScalingNumIterations; index++) {
			int i = 0;
			int j = 0;
			factor = (factor_high + factor_low)/2;

			for (i = 0; i < internal_pssm.ncol(); i++) {
				for (j = 0; j < internal_pssm.nrow(); j++) {
					if (scaled_pssm[i][j] != BLAST_SCORE_MIN) {
						pssm[i][j] = 
							BLAST_Nint(factor*scaled_pssm[i][j]/kPSIScaleFactor);
					} else {
						pssm[i][j] = BLAST_SCORE_MIN;
					}
				}
			}

			updateLambdaK(pssm,sbp);
			new_lambda = sbp.kbp_lambda;
			if (new_lambda > ideal_lambda) {
				factor_low = factor;
			} else {
				factor_high = factor;
			}
		}
		return true;
	//	return PSI_SUCCESS;
	}

	public static int BLAST_Nint(double x){
		x += (x >= 0. ? 0.5 : -0.5);
		return (int)x;
	}

	public double[][] addPseudoFrequency(double[][] f,double[] std_prob){
		double[][] ret = new double[f.length][f[0].length];
		/* As specified in 2001 paper, underlying matrix frequency 
                   ratios are used here */
		double kBeta = 1;//参考になる値が無い。デフォルトは 30 だが MSA のエントロピーを基に計算される
		double observations = 30;//Effective observation であり MSA から計算される 
		for(int ii = 0;ii < f.length;ii++){
			double rsum = 0;
			for(int r = 0;r < f[ii].length;r++){
				double pseudo = 0;
				//f[ii][i] = seq_weights->match_weights[p][i]  とおいているが。。。
                for (int i = 0; i < f[ii].length; i++) {
                        pseudo += (f[ii][i] *
                                   Blosum62Freq.freq_ratios[r][i]);
               }
                pseudo *= kBeta;

                double numerator =
                    (observations * f[ii][r] / 
                     std_prob[r]) 
                    + pseudo;

                double denominator = observations + kBeta;

                double qOverPEstimate = numerator/denominator;

                /* Note artificial multiplication by standard probability
                 * to normalize */
                ret[ii][r] = qOverPEstimate *std_prob[r];
				rsum += ret[ii][r];
			}
			//これは PSIBLAST には見当たらないが・・・。
			for(int r = 0;r < f[ii].length;r++){
				 ret[ii][r] /= rsum;
			}
		}
		return ret;
	}
	
	
	public static void main(String[] args){
		PSSMCalc cal = new PSSMCalc();
		BlastScoreBlk scoreblk = new BlastScoreBlk();
		scoreblk.std_prob =PSSMCalc.background_frequency;
		
		PSSMData testloaded = PSSMData.load("C:\\dummy\\private\\bioo\\seminar_mine\\20170830\\hemo_test.pssm");
		double[][] freq = new double[testloaded.scores.size()][20];
		char colhead[] = "ARNDCQEGHILKMFPSTWYV".toCharArray();
		
		
		
		for(int ii = 0;ii < testloaded.scores.size();ii++){
			ArrayList<Double> d=testloaded.scores.get(ii);
			for(char cc:colhead){
				freq[ii][PSSMCalc.aa_to_index[cc]] = d.get(testloaded.char_index_map[cc]+20)/100.0;
			}
		}
		freq = cal.addPseudoFrequency(freq,PSSMCalc.background_frequency);
		PSIInternalPssmData ps  = new PSIInternalPssmData(freq,scoreblk.std_prob,scoreblk.ideal_lambda);
		cal.PSIScaleMatrix(ps,scoreblk);
		for(int ii = 0;ii < ps.pssm.length;ii++){
			for(int jj = 0;jj < ps.pssm[0].length;jj++){
				//
				//こんな値が出る
//				-130	-122	-115	-108	-101	-93	-86	-79	-71	-64	-57	-50	-42	-35	-28	-21	-13	-6	1	9	
//-130	-122	-115	-108	-101	-93	-86	-79	-71	-64	-57	-50	-42	-35	-28	-21	-13	-6	1	9	

				System.out.print((int)(ps.pssm[ii][jj])+"\t");
			}
			System.out.println("");
			
		}
	}
}

class PSIInternalPssmData{
	double lambda;
	int[][] pssm;
	int[][] scaled_pssm;
	double[][] freq_ratio;
	PSIInternalPssmData(double[][] f,double[] std_probs,double ideal_lambda){
		
		
		
		
		
		
		
		setFreqRatio(f,std_probs,ideal_lambda);
		
	}
	public int nrow(){
		return pssm[0].length;
	}
	public int ncol(){
		return pssm.length;
	}
	
	/**
	 * internal_pssm.freq を基に internal_pssm.scaled_pssm を計算する
	 * @param internal_pssm
	 * @param sbp
	 * @param std_probs
	 * @return 
	 */
	
	public void setFreqRatio(double[][] f,double[] std_probs,double ideal_lambda){
		freq_ratio = f;
		int ncol = freq_ratio.length;
		int asiz = freq_ratio[0].length;
		scaled_pssm = new int[ncol][asiz];
		pssm = new int[ncol][asiz];
		
		for (int i = 0; i < ncol; i++) {

			/* True if all frequencies in column i are zero */
			for (int j = 0; j < 20; j++) {
				scaled_pssm[i][j] = 0;
				double qOverPEstimate = 0.0;
				if(std_probs[j] < PSSMCalc.kEpsilon){
						System.err.println("std_probs of "+j+" is smaller than kEpsilon. ");
						System.err.println("set as "+PSSMCalc.kEpsilon);
						std_probs[j] = PSSMCalc.kEpsilon;
				}
				qOverPEstimate = 
						freq_ratio[i][j] / std_probs[j];


				if (qOverPEstimate == 0.0) {
					scaled_pssm[i][j] = PSSMCalc.BLAST_SCORE_MIN;
				} else {
					double tmp = Math.log(qOverPEstimate)/ideal_lambda;
					scaled_pssm[i][j] = Math.min(Math.max(PSSMCalc.BLAST_SCORE_MIN,(int)PSSMCalc.BLAST_Nint(PSSMCalc.kPSIScaleFactor * tmp)
					),PSSMCalc.BLAST_SCORE_MAX);
				}
			}
		}

	}

}
class BlastScoreBlk{
	double ideal_lambda=0.321;
	double kbp_lambda;
	double[] std_prob;
	/*
	blast_stat.c 
	デフォルトの Lambda
	
Lambda      K        H
   0.321    0.139    0.431 
Gapped
Lambda      K        H
   0.267   0.0410    0.140 
	この三番目の値が Lambda。
	static array_of_8 blosum62_values[BLOSUM62_VALUES_MAX] = {
    {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660},
    {11, 2, (double) INT2_MAX, 0.297, 0.082, 0.27, 1.1, -10, 0.641766, 12.673800, 12.757600},
    {10, 2, (double) INT2_MAX, 0.291, 0.075, 0.23, 1.3, -15, 0.649362, 16.474000, 16.602600},
    {9, 2, (double) INT2_MAX, 0.279, 0.058, 0.19, 1.5, -19, 0.659245, 22.751900, 22.950000},
    {8, 2, (double) INT2_MAX, 0.264, 0.045, 0.15, 1.8, -26, 0.672692, 35.483800, 35.821300},
    {7, 2, (double) INT2_MAX, 0.239, 0.027, 0.10, 2.5, -46, 0.702056, 61.238300, 61.886000},
    {6, 2, (double) INT2_MAX, 0.201, 0.012, 0.061, 3.3, -58, 0.740802, 140.417000, 141.882000},
    {13, 1, (double) INT2_MAX, 0.292, 0.071, 0.23, 1.2, -11, 0.647715, 19.506300, 19.893100},
    {12, 1, (double) INT2_MAX, 0.283, 0.059, 0.19, 1.5, -19, 0.656391, 27.856200, 28.469900},
    {11, 1, (double) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200},
    {10, 1, (double) INT2_MAX, 0.243, 0.024, 0.10, 2.5, -44, 0.693267, 83.178700, 85.065600},
    {9, 1, (double) INT2_MAX, 0.206, 0.010, 0.052, 4.0, -87, 0.731887, 210.333000, 214.842000},
	}; 
	*/
	
}

class Blast_ScoreFreq{
	public double sprob[];//スコアの Probability
	public double score_average;
	public int obs_min = 0;
	public int obs_max = 0;
	public void setMinMax(int mi,int ma){
		obs_min = mi;
		obs_max = ma;
	}
	public void setScoreProbAt(int i,double d){
		sprob[i-obs_min] = d;
	}
	public double getScoreProbAt(int i){
		return sprob[i-obs_min];
	}
	
	
	public static Blast_ScoreFreq PSIComputeScoreProbabilities(int[][] pssm,double[] std_probs){
		int alphabet_size = 0;                    /* number of elements populated*/
		int p = 0;                                /* index on positions */
		int r = 0;                                /* index on residues */
		int s = 0;                                  /* index on scores */
		int min_score = PSSMCalc.BLAST_SCORE_MAX;            /* minimum score in pssm */
		int max_score = PSSMCalc.BLAST_SCORE_MIN;            /* maximum score in pssm */
		Blast_ScoreFreq score_freqs = new Blast_ScoreFreq();        /* score frequencies */


		//effective_length = _PSISequenceLengthWithoutX(query, query_length);
		int query_length = pssm.length;
		int effective_length = query_length;

		/* Get the minimum and maximum scores */
		for (p = 0; p < query_length; p++) {
			for (r = 0; r < PSSMCalc.nb_aatype; r++) {
				int kScore = pssm[p][r];
				if (kScore <= PSSMCalc.BLAST_SCORE_MIN || kScore >= PSSMCalc.BLAST_SCORE_MAX) {
					continue;
				}
				max_score = Math.max(kScore, max_score);
				min_score = Math.min(kScore, min_score);
			}
		}


		score_freqs.obs_min = min_score;
		score_freqs.obs_max = max_score;
		score_freqs.sprob = new double[score_freqs.obs_max-score_freqs.obs_min+1];
		for (p = 0; p < query_length; p++) {
			for (r = 0; r <  PSSMCalc.nb_aatype; r++) {
				int kScore = pssm[p][r];
				if (kScore <= PSSMCalc.BLAST_SCORE_MIN || kScore >= PSSMCalc.BLAST_SCORE_MAX) {
					continue;
				}
				/* Increment the weight for the score in position p, residue r */
				score_freqs.setScoreProbAt(kScore,score_freqs.getScoreProbAt(kScore)
						+std_probs[r]);
			}
		}

		for (s = min_score; s <= max_score; s++) {
			score_freqs.setScoreProbAt(s,score_freqs.getScoreProbAt(s)/ effective_length);
			score_freqs.score_average += (s * score_freqs.getScoreProbAt(s));
		}

		return score_freqs;
	}

}

class Blosum62Freq{
	//underlying matrix frequency ratios
	//source code says it is described in 2001 paper
	static  double BLOSUM62_FREQRATIOS[][] =
	{{0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
	  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
	 {0.00000000e+00, 3.90294070e+00, 5.64459671e-01, 8.67987664e-01,
	  5.44605275e-01, 7.41264113e-01, 4.64893827e-01, 1.05686961e+00,
	  5.69364849e-01, 6.32481035e-01, 7.75390239e-01, 6.01945975e-01,
	  7.23150342e-01, 5.88307640e-01, 7.54121369e-01, 7.56803943e-01,
	  6.12698600e-01, 1.47210399e+00, 9.84401956e-01, 9.36458396e-01,
	  4.16548781e-01, 7.50000000e-01, 5.42611869e-01, 7.47274948e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.14377313e-01},
	 {0.00000000e+00, 5.64459671e-01, 4.43758048e+00, 3.45226274e-01,
	  4.74290926e+00, 1.33503378e+00, 3.24101420e-01, 7.38524318e-01,
	  9.25449581e-01, 3.33981361e-01, 8.54849426e-01, 2.97257620e-01,
	  4.04640322e-01, 4.07083696e+00, 5.53838329e-01, 9.44103648e-01,
	  7.02873767e-01, 1.05798620e+00, 8.26250098e-01, 3.51280513e-01,
	  2.52855433e-01, 7.50000000e-01, 4.09444638e-01, 1.18382127e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.12208474e-01},
	 {0.00000000e+00, 8.67987664e-01, 3.45226274e-01, 1.95765857e+01,
	  3.01454345e-01, 2.85934574e-01, 4.38990118e-01, 4.20387870e-01,
	  3.55049505e-01, 6.53458801e-01, 3.49128465e-01, 6.42275633e-01,
	  6.11354340e-01, 3.97802620e-01, 3.79562691e-01, 3.65781531e-01,
	  3.08939296e-01, 7.38415701e-01, 7.40551692e-01, 7.55844055e-01,
	  4.49983903e-01, 7.50000000e-01, 4.34203398e-01, 3.16819526e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.46828489e-01},
	 {0.00000000e+00, 5.44605275e-01, 4.74290926e+00, 3.01454345e-01,
	  7.39792738e+00, 1.68781075e+00, 2.98969081e-01, 6.34301019e-01,
	  6.78558839e-01, 3.39015407e-01, 7.84090406e-01, 2.86613046e-01,
	  3.46454634e-01, 1.55385281e+00, 5.98716826e-01, 8.97081129e-01,
	  5.73200024e-01, 9.13504624e-01, 6.94789868e-01, 3.36500142e-01,
	  2.32102315e-01, 7.50000000e-01, 3.45683565e-01, 1.38195506e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.07946931e-01},
	 {0.00000000e+00, 7.41264113e-01, 1.33503378e+00, 2.85934574e-01,
	  1.68781075e+00, 5.46952608e+00, 3.30743991e-01, 4.81267655e-01,
	  9.60040718e-01, 3.30522558e-01, 1.30827885e+00, 3.72873704e-01,
	  5.00342289e-01, 9.11298183e-01, 6.79202587e-01, 1.90173784e+00,
	  9.60797602e-01, 9.50357185e-01, 7.41425610e-01, 4.28943130e-01,
	  3.74300212e-01, 7.50000000e-01, 4.96467354e-01, 4.08949895e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.55631838e-01},
	 {0.00000000e+00, 4.64893827e-01, 3.24101420e-01, 4.38990118e-01,
	  2.98969081e-01, 3.30743991e-01, 8.12879702e+00, 3.40640908e-01,
	  6.51990521e-01, 9.45769883e-01, 3.44043119e-01, 1.15459749e+00,
	  1.00437163e+00, 3.54288952e-01, 2.87444758e-01, 3.33972402e-01,
	  3.80726330e-01, 4.39973597e-01, 4.81693683e-01, 7.45089738e-01,
	  1.37437942e+00, 7.50000000e-01, 2.76938063e+00, 3.31992746e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.06958025e+00},
	 {0.00000000e+00, 1.05686961e+00, 7.38524318e-01, 4.20387870e-01,
	  6.34301019e-01, 4.81267655e-01, 3.40640908e-01, 6.87630691e+00,
	  4.92966576e-01, 2.75009722e-01, 5.88871736e-01, 2.84504012e-01,
	  3.95486600e-01, 8.63711406e-01, 4.77385507e-01, 5.38649627e-01,
	  4.49983999e-01, 9.03596525e-01, 5.79271582e-01, 3.36954912e-01,
	  4.21690355e-01, 7.50000000e-01, 3.48714366e-01, 5.03463109e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.80638726e-01},
	 {0.00000000e+00, 5.69364849e-01, 9.25449581e-01, 3.55049505e-01,
	  6.78558839e-01, 9.60040718e-01, 6.51990521e-01, 4.92966576e-01,
	  1.35059997e+01, 3.26288125e-01, 7.78887490e-01, 3.80675486e-01,
	  5.84132623e-01, 1.22200067e+00, 4.72879831e-01, 1.16798104e+00,
	  9.17048021e-01, 7.36731740e-01, 5.57503254e-01, 3.39447442e-01,
	  4.44088955e-01, 7.50000000e-01, 1.79790413e+00, 1.04047242e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.58533474e-01},
	 {0.00000000e+00, 6.32481035e-01, 3.33981361e-01, 6.53458801e-01,
	  3.39015407e-01, 3.30522558e-01, 9.45769883e-01, 2.75009722e-01,
	  3.26288125e-01, 3.99792994e+00, 3.96372934e-01, 1.69443475e+00,
	  1.47774450e+00, 3.27934752e-01, 3.84662860e-01, 3.82937802e-01,
	  3.54751311e-01, 4.43163582e-01, 7.79816110e-01, 2.41751209e+00,
	  4.08874390e-01, 7.50000000e-01, 6.30388931e-01, 3.50796872e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.63222650e+00},
	 {0.00000000e+00, 7.75390239e-01, 8.54849426e-01, 3.49128465e-01,
	  7.84090406e-01, 1.30827885e+00, 3.44043119e-01, 5.88871736e-01,
	  7.78887490e-01, 3.96372934e-01, 4.76433717e+00, 4.28270363e-01,
	  6.25302816e-01, 9.39841129e-01, 7.03774479e-01, 1.55432308e+00,
	  2.07680867e+00, 9.31919141e-01, 7.92905803e-01, 4.56542720e-01,
	  3.58930071e-01, 7.50000000e-01, 5.32179333e-01, 1.40344922e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.15284382e-01},
	 {0.00000000e+00, 6.01945975e-01, 2.97257620e-01, 6.42275633e-01,
	  2.86613046e-01, 3.72873704e-01, 1.15459749e+00, 2.84504012e-01,
	  3.80675486e-01, 1.69443475e+00, 4.28270363e-01, 3.79662137e+00,
	  1.99429557e+00, 3.10043276e-01, 3.71121724e-01, 4.77325586e-01,
	  4.73919278e-01, 4.28893743e-01, 6.60328975e-01, 1.31423573e+00,
	  5.68037074e-01, 7.50000000e-01, 6.92059423e-01, 4.13275887e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.94078574e+00},
	 {0.00000000e+00, 7.23150342e-01, 4.04640322e-01, 6.11354340e-01,
	  3.46454634e-01, 5.00342289e-01, 1.00437163e+00, 3.95486600e-01,
	  5.84132623e-01, 1.47774450e+00, 6.25302816e-01, 1.99429557e+00,
	  6.48145121e+00, 4.74529655e-01, 4.23898024e-01, 8.64250293e-01,
	  6.22623369e-01, 5.98558924e-01, 7.93801616e-01, 1.26893679e+00,
	  6.10296214e-01, 7.50000000e-01, 7.08364628e-01, 6.41102583e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.78399892e+00},
	 {0.00000000e+00, 5.88307640e-01, 4.07083696e+00, 3.97802620e-01,
	  1.55385281e+00, 9.11298183e-01, 3.54288952e-01, 8.63711406e-01,
	  1.22200067e+00, 3.27934752e-01, 9.39841129e-01, 3.10043276e-01,
	  4.74529655e-01, 7.09409488e+00, 4.99932836e-01, 1.00058442e+00,
	  8.58630478e-01, 1.23152924e+00, 9.84152635e-01, 3.69033853e-01,
	  2.77782896e-01, 7.50000000e-01, 4.86030806e-01, 9.45834265e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.17327197e-01},
	 {0.00000000e+00, 7.54121369e-01, 5.53838329e-01, 3.79562691e-01,
	  5.98716826e-01, 6.79202587e-01, 2.87444758e-01, 4.77385507e-01,
	  4.72879831e-01, 3.84662860e-01, 7.03774479e-01, 3.71121724e-01,
	  4.23898024e-01, 4.99932836e-01, 1.28375437e+01, 6.41280589e-01,
	  4.81534905e-01, 7.55503259e-01, 6.88897122e-01, 4.43082984e-01,
	  2.81833164e-01, 7.50000000e-01, 3.63521119e-01, 6.64534287e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.76634549e-01},
	 {0.00000000e+00, 7.56803943e-01, 9.44103648e-01, 3.65781531e-01,
	  8.97081129e-01, 1.90173784e+00, 3.33972402e-01, 5.38649627e-01,
	  1.16798104e+00, 3.82937802e-01, 1.55432308e+00, 4.77325586e-01,
	  8.64250293e-01, 1.00058442e+00, 6.41280589e-01, 6.24442175e+00,
	  1.40579606e+00, 9.65555228e-01, 7.91320741e-01, 4.66777931e-01,
	  5.09360272e-01, 7.50000000e-01, 6.11094097e-01, 3.58149606e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.38898727e-01},
	 {0.00000000e+00, 6.12698600e-01, 7.02873767e-01, 3.08939296e-01,
	  5.73200024e-01, 9.60797602e-01, 3.80726330e-01, 4.49983999e-01,
	  9.17048021e-01, 3.54751311e-01, 2.07680867e+00, 4.73919278e-01,
	  6.22623369e-01, 8.58630478e-01, 4.81534905e-01, 1.40579606e+00,
	  6.66557707e+00, 7.67165633e-01, 6.77754679e-01, 4.20072316e-01,
	  3.95102106e-01, 7.50000000e-01, 5.55965425e-01, 1.13292384e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.25403989e-01},
	 {0.00000000e+00, 1.47210399e+00, 1.05798620e+00, 7.38415701e-01,
	  9.13504624e-01, 9.50357185e-01, 4.39973597e-01, 9.03596525e-01,
	  7.36731740e-01, 4.43163582e-01, 9.31919141e-01, 4.28893743e-01,
	  5.98558924e-01, 1.23152924e+00, 7.55503259e-01, 9.65555228e-01,
	  7.67165633e-01, 3.84284741e+00, 1.61392097e+00, 5.65223766e-01,
	  3.85303035e-01, 7.50000000e-01, 5.57520051e-01, 9.56235816e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.34703235e-01},
	 {0.00000000e+00, 9.84401956e-01, 8.26250098e-01, 7.40551692e-01,
	  6.94789868e-01, 7.41425610e-01, 4.81693683e-01, 5.79271582e-01,
	  5.57503254e-01, 7.79816110e-01, 7.92905803e-01, 6.60328975e-01,
	  7.93801616e-01, 9.84152635e-01, 6.88897122e-01, 7.91320741e-01,
	  6.77754679e-01, 1.61392097e+00, 4.83210516e+00, 9.80943005e-01,
	  4.30934144e-01, 7.50000000e-01, 5.73156574e-01, 7.60725140e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.08974203e-01},
	 {0.00000000e+00, 9.36458396e-01, 3.51280513e-01, 7.55844055e-01,
	  3.36500142e-01, 4.28943130e-01, 7.45089738e-01, 3.36954912e-01,
	  3.39447442e-01, 2.41751209e+00, 4.56542720e-01, 1.31423573e+00,
	  1.26893679e+00, 3.69033853e-01, 4.43082984e-01, 4.66777931e-01,
	  4.20072316e-01, 5.65223766e-01, 9.80943005e-01, 3.69215640e+00,
	  3.74456332e-01, 7.50000000e-01, 6.58038693e-01, 4.43577702e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.76339815e+00},
	 {0.00000000e+00, 4.16548781e-01, 2.52855433e-01, 4.49983903e-01,
	  2.32102315e-01, 3.74300212e-01, 1.37437942e+00, 4.21690355e-01,
	  4.44088955e-01, 4.08874390e-01, 3.58930071e-01, 5.68037074e-01,
	  6.10296214e-01, 2.77782896e-01, 2.81833164e-01, 5.09360272e-01,
	  3.95102106e-01, 3.85303035e-01, 4.30934144e-01, 3.74456332e-01,
	  3.81077833e+01, 7.50000000e-01, 2.10980812e+00, 4.26541694e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 5.03239261e-01},
	 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
	 {0.00000000e+00, 5.42611869e-01, 4.09444638e-01, 4.34203398e-01,
	  3.45683565e-01, 4.96467354e-01, 2.76938063e+00, 3.48714366e-01,
	  1.79790413e+00, 6.30388931e-01, 5.32179333e-01, 6.92059423e-01,
	  7.08364628e-01, 4.86030806e-01, 3.63521119e-01, 6.11094097e-01,
	  5.55965425e-01, 5.57520051e-01, 5.73156574e-01, 6.58038693e-01,
	  2.10980812e+00, 7.50000000e-01, 9.83220341e+00, 5.40805192e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.66952325e-01},
	 {0.00000000e+00, 7.47274948e-01, 1.18382127e+00, 3.16819526e-01,
	  1.38195506e+00, 4.08949895e+00, 3.31992746e-01, 5.03463109e-01,
	  1.04047242e+00, 3.50796872e-01, 1.40344922e+00, 4.13275887e-01,
	  6.41102583e-01, 9.45834265e-01, 6.64534287e-01, 3.58149606e+00,
	  1.13292384e+00, 9.56235816e-01, 7.60725140e-01, 4.43577702e-01,
	  4.26541694e-01, 7.50000000e-01, 5.40805192e-01, 3.89300249e+00,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.87839626e-01},
	 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
	 {0.00000000e+00, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
	  2.50000000e-01, 1.33300000e+00, 2.50000000e-01, 2.50000000e-01},
	 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
	 {0.00000000e+00, 6.14377313e-01, 3.12208474e-01, 6.46828489e-01,
	  3.07946931e-01, 3.55631838e-01, 1.06958025e+00, 2.80638726e-01,
	  3.58533474e-01, 2.63222650e+00, 4.15284382e-01, 2.94078574e+00,
	  1.78399892e+00, 3.17327197e-01, 3.76634549e-01, 4.38898727e-01,
	  4.25403989e-01, 4.34703235e-01, 7.08974203e-01, 1.76339815e+00,
	  5.03239261e-01, 7.50000000e-01, 6.66952325e-01, 3.87839626e-01,
	  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.81516607e+00}};

	//from blast_encoding.c
	static char NCBISTDAA_TO_AMINOACID[] = {
	'-','A','B','C','D','E','F','G','H','I','K','L','M',
	'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
	'O', 'J'};

	static double[][] freq_ratios = new double[20][20];
	static{
		for(int ii = 0;ii < NCBISTDAA_TO_AMINOACID.length;ii++){
			if(PSSMCalc.aa_to_index[NCBISTDAA_TO_AMINOACID[ii]] > -1){
				for(int jj = 0;jj < NCBISTDAA_TO_AMINOACID.length;jj++){
					if(PSSMCalc.aa_to_index[NCBISTDAA_TO_AMINOACID[jj]] > -1){
						freq_ratios[PSSMCalc.aa_to_index[NCBISTDAA_TO_AMINOACID[ii]]]
								[PSSMCalc.aa_to_index[NCBISTDAA_TO_AMINOACID[jj]]] 
								= BLOSUM62_FREQRATIOS[ii][jj];
							//	System.out.println(NCBISTDAA_TO_AMINOACID[ii]+";"+NCBISTDAA_TO_AMINOACID[jj]
							//			+"\t"+BLOSUM62_FREQRATIOS[ii][jj]);
					}

				}
			}
		}
	}


}
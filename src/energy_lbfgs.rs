

#[allow(unused_imports)]
use std::f64;
#[allow(unused_imports)]
use std::collections::{HashSet,HashMap};
#[allow(unused_imports)]
use rand::prelude::*;
#[allow(unused_imports)]
use std::fs::File;
use super::pp_energy;
#[allow(unused_imports)]
use super::charmm_based_energy;
use super::geometry::Vector3D;
#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};

pub const FPRIME_EPSILON:f64 = 1.0e-8;//数値微分するときに使用する移動量
pub const INNERCHECK_EPSILON:f64 = 1.0e-8;//RHO の分母がこれ以下になった場合失敗とみなして GoldenRUle だけで何とかしようとする
pub const LOGISTICLOSS_STOP:f64 = 1.0e-8;//完全に分離できる場合にオーバーフローを起こすので、予測の失敗がこれ未満になった場合収束とみなす
pub const GOLDEN_LOOPCHECK:usize = 100;//GOLDEN_NUMBER^(LOOP 回⁺１) されるはず。756 で Inf になる。下限は /0.5 されるが、どこまで行けるのかチェックしてない。不明。
pub const LOSSCALC_SORT:bool = false;//小さい数から足していくことで桁あふれを最小限にする

pub const DEFAULT_MAX_ITER:usize = 1000;


pub const GOLDEN_NUMBER:f64 =1.6180339887;//excel で計算
pub const GOLDEN_BREAK_THRESHOLD:f64 = 0.0001;
pub const GOLDEN_BREAK_THRESHOLD_X:f64 = 0.0001;


pub const X:usize = 0;
pub const Y:usize = 1;
pub const Z:usize = 2;


pub trait LinearRegressionPenalty{
	fn calc_penalty(&self,weights:&Vec<f64>)->f64;
}

pub struct PenaltyRedge{
	pub alpha:f64
}
impl LinearRegressionPenalty for PenaltyRedge{
	fn calc_penalty(&self,betas: &Vec<f64>)->f64{
		let mut ret:f64 = 0.0;
		for i in 1..betas.len(){
			ret += betas[i]*betas[i];
		}
		return ret*self.alpha;
	}
}


pub struct PenaltyLasso{
	pub alpha:f64
}
impl LinearRegressionPenalty for PenaltyLasso{
	fn calc_penalty(&self,betas: &Vec<f64>)->f64{
		let mut ret:f64 = 0.0;
		for i in 1..betas.len(){
			ret += betas[i].abs();
		}
		return ret*self.alpha;
	}
}


pub struct PenaltyElasticNet{
	pub alpha:f64,
	pub ratio_l1:f64
}
impl LinearRegressionPenalty for PenaltyElasticNet{
	fn calc_penalty(&self,betas: &Vec<f64>)->f64{
		let mut ret_l1:f64 = 0.0;
		let mut ret_l2:f64 = 0.0;
		for i in 1..betas.len(){
			ret_l1 += betas[i].abs();
			ret_l2 += betas[i]*betas[i];
		}
		return self.alpha*((1.0-self.ratio_l1)/2.0*ret_l2+self.ratio_l1*ret_l1);
	}
}






pub enum LinearRegressionFunctionType{
	BinomialLogistic,
	LinearRMSE,
}


//[row][col]
pub fn matrix_multi(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;y[0].len()];x.len()];
	let rnum = x.len();
	let cnum = y[0].len();
	if x[0].len() != y.len(){
		panic!("The num col in x & numrow in y must be equal x:{},{};y:{},{}",x.len(),x[0].len(),y.len(),y[0].len());
	}
	for ii in 0..rnum{
		for jj in 0..cnum{
			for rr in 0..x[0].len(){
				let mut xn:f64 = 0.0;
				let mut yn:f64 = 0.0;
				if x[ii].len() > rr{
					xn = x[ii][rr];
				} 
				if y.len() > rr{
					yn = y[rr][jj];
				} 
				ret[ii][jj] += xn*yn;
			}
		}
	}
	return ret;
}

pub fn matrix_t(x:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;x.len()];x[0].len()];
	for rr in 0..x.len(){
		for cc in 0..x[0].len(){
			ret[cc][rr] = x[rr][cc];
		}
	}
	return ret;
}
pub fn matrix_scalarmulti(x:&Vec<Vec<f64>>,y:f64)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;x[0].len()];x.len()];
	for rr in 0..x.len(){
		for cc in 0..x[0].len(){
			ret[rr][cc] = x[rr][cc]*y;
		}
	}
	return ret;
}
pub fn matrix_add(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;x[0].len()];x.len()];
	
	if x.len() != y.len() || x[0].len() != y[0].len(){
		panic!("shapes must be the same x:{},{};y:{},{}",x.len(),x[0].len(),y.len(),y[0].len());
	}
	for rr in 0..x.len(){
		for cc in 0..x[0].len(){
			ret[rr][cc] = x[rr][cc]+y[rr][cc];
		}
	}
	return ret;
}

pub fn ave(v1:&Vec<f64>)->f64{
    if v1.len() == 0{
        return 0.0;
    }
    return v1.iter().fold(0.0,|s,m|{s+m})/(v1.len() as f64);
}

pub fn var(v1:&Vec<f64>)->f64{
    if v1.len() == 0{
        return f64::NAN;
    }
    let a:f64 = ave(v1);
    let dff = v1.iter().fold(0.0,|s,m|{s+(m-a)*(m-a)})/(v1.len() as f64);
    return dff;
}

pub fn ave_w(v1:&Vec<f64>,weight:&Vec<f64>)->f64{
    if v1.len() == 0{
        return 0.0;
    }
	let wsum:f64 = weight.iter().fold(0.0,|s,m|{s+m});
	if wsum == 0.0{
		return 0.0;
	}
    return v1.iter().zip(weight.iter()).fold(0.0,|s,m|{s+m.0*m.1})/wsum;
}

pub fn var_w(v1:&Vec<f64>,weight:&Vec<f64>)->f64{
    if v1.len() == 0{
        return f64::NAN;
    }
	let wsum:f64 = weight.iter().fold(0.0,|s,m|{s+m});
	if wsum == 0.0{
		return f64::NAN;
	}
    let a:f64 = ave_w(v1,weight);
    let dff = v1.iter().zip(weight.iter()).fold(0.0,|s,m|{s+(m.0-a)*(m.0-a)*m.1})/wsum;
    return dff;
}



pub fn pearson_correl(v1:&Vec<f64>,v2:&Vec<f64>)->f64{
    let a1:f64 = ave(v1);
    let a2:f64 = ave(v2);
    if v1.len()*v2.len() == 0{
        return f64::NAN;
    }
    if v1.len() != v2.len(){
        panic!("The lengthes of vectors must be the same.");
    }
    let mut cov:f64 = 0.0;
    for ii in 0..v1.len(){
        cov += (a1-v1[ii])*(a2-v2[ii]);
    }
    
    let va1:f64 = var(v1);
    let va2:f64 = var(v2);
    if va1.is_nan() || va2.is_nan() {
        return f64::NAN;
    }
    return cov/(va1.sqrt()*va2.sqrt())/(v1.len() as f64);
}


pub fn matrix_diagx(r:usize,x:f64)-> Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;r];r];
	for rr in 0..r{
		for cc in 0..r{
			if rr == cc{
				ret[rr][cc] = x;
			}
		}
	}
	return ret;
}

pub fn array_multi(x:&Vec<f64>,y:&Vec<f64>,res:&mut Vec<f64>){
	assert_eq!(x.len(),y.len());
	assert_eq!(x.len(),res.len());
	for ii in 0..x.len(){
		res[ii] = x[ii]*y[ii];
	}
}

pub fn array_inner(x:&Vec<f64>,y:&Vec<f64>)->f64{
	let mut res = 0.0;
	for ii in 0..x.len(){
		res += x[ii]*y[ii];
	}
	return res;
}

//xstep 検討する各段階における beta の変化量 To do 増加か減少かの確認
//graddiff 検討する各段階における gradient の変化量
//h_zero 開始時の行列。大抵 1 の直行行列
//prevgrad 直前のステップの gradient
pub fn calc_next_gradient_lbfgs_straight(xstep:&Vec<Vec<f64>>,graddiff:&Vec<Vec<f64>>
,h_zero:&Vec<f64>,prevgrad:&Vec<f64>)
->Vec<f64>{
	let vlen:usize = xstep[0].len();//係数の数
	let slen:usize = xstep.len();//計算に使用するステップ数
	let mut vt:Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0;vlen];vlen];slen];
	let mut rho:Vec<f64> = vec![0.0;slen];
	
	for ii in 0..slen{
		let lb = array_inner(&xstep[ii],&graddiff[ii]);
		rho[ii] = 1.0/lb;
	}
	
	for ii in 0..slen{
		//一次元配列と行列ではカラムと行の扱いが逆になる
		let vm = matrix_multi(&matrix_t(&vec![graddiff[ii].clone()])
		,&vec![xstep[ii].clone()]);
		
		vt[ii] = matrix_t(&matrix_add(&matrix_diagx(vlen,1.0)
		,&matrix_scalarmulti(&vm,-1.0*rho[ii])));
		
	}
	

	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;vlen];vlen];
	for ii in 0..(slen+1){
		let mut kres:Vec<Vec<f64>>=matrix_diagx(vlen,1.0);
		
		if ii < slen{
			for xx in 0..(slen-ii){
				kres = matrix_multi(&kres,&vt[slen-1-xx]);
			}
		}
		let mut ctt:Vec<Vec<f64>>=vec![vec![0.0;vlen];vlen];
		if ii == 0{
			for jj in 0..vlen{
				ctt[jj][jj] = h_zero[jj];
			}
		}else{
			ctt = matrix_scalarmulti(&matrix_multi(
			&matrix_t(&vec![xstep[ii-1].clone()])
			,&vec![xstep[ii-1].clone()])
			,rho[ii-1]);
		}
		kres = matrix_multi(&kres,&ctt);

		if ii < slen{
	
			for xx in ii..slen{
				kres = matrix_multi(&kres,&matrix_t(&vt[xx]));
			}
			
		}
		ret = matrix_add(&ret,&kres);
	}
	let mut tmp:Vec<Vec<f64>> =  matrix_t(&matrix_multi(&ret,&matrix_t(&vec![prevgrad.clone()])));

	return tmp.remove(0);
}




//省メモリ出来るトリックをインプリした奴
pub fn calc_next_gradient_lbfgs(xstep:&Vec<Vec<f64>>
,graddiff:&Vec<Vec<f64>>,prevgrad:&Vec<f64>)
->Vec<f64>{
	let vlen:usize = xstep[0].len();//係数の数
	let slen:usize = xstep.len();//計算に使用するステップ数
	let mut rho:Vec<f64> = vec![0.0;slen];
	let mut q:Vec<f64> = prevgrad.clone();
	let mut alpha:Vec<f64> = vec![0.0;slen];
	for ii in 0..slen{
		let lb = array_inner(&xstep[ii],&graddiff[ii]);
		rho[ii] = 1.0/lb;
	}
	for i in 0..slen{
		let pos = slen-1-i;
		alpha[pos] = array_inner(&xstep[pos],&q)*rho[pos];
		for jj in 0..vlen{
			q[jj] = q[jj]-graddiff[pos][jj]*alpha[pos];
		}
	}
	/*
	let h_zero:Vec<f64> = vec![array_inner(&graddiff[slen-1],&xstep[slen-1])
	/array_inner(&graddiff[slen-1],&graddiff[slen-1]);vlen];
	*/
	let h_zero:Vec<f64> =  vec![1.0;vlen];
	let mut ret:Vec<f64> = vec![0.0;vlen];
	for jj in 0..vlen{
		//ret[jj] = -1.0*q[jj]*h_zero[jj];
		ret[jj] = q[jj]*h_zero[jj];
	}

	for ii in 0..(slen){
		let beta = rho[ii]*array_inner(&graddiff[ii],&ret);
		for jj in 0..vlen{
			ret[jj] = ret[jj]+xstep[ii][jj]*(alpha[ii]-beta);
		}
	}

	return ret;
}




fn calc_enegy_with_beta(
	betas:&Vec<(f64,usize,usize)>
	,buff_array:&mut Vec<f64>//前の値を置いておく配列
	,energyset:&mut pp_energy::PPEnergySet
	,atom_level_energy:&mut Vec<f64>//これはダミーであってよいがいつか消すかも
)->f64{
	for ii in 0..betas.len(){
		let prevval:f64;
		if betas[ii].2 == X{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_x();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_x(prevval+betas[ii].0);
		}else if betas[ii].2 == Y{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_y();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_y(prevval+betas[ii].0);
		}else{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_z();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_z(prevval+betas[ii].0);
		}
		buff_array[ii] = prevval;
	}
	energyset.update_distance();
	let ret = energyset.calc_energy(atom_level_energy);

	for ii in 0..betas.len(){
		if betas[ii].2 == X{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_x(buff_array[ii]);
		}else if betas[ii].2 == Y{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_y(buff_array[ii]);
		}else{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_z(buff_array[ii]);
		}
	}
	energyset.update_distance();
	return ret;
}

fn calc_fprime(betas:&Vec<(f64,usize,usize)>
	,energyset:&mut pp_energy::PPEnergySet
	,atom_level_energy:&mut Vec<f64>//これはダミーであってよいがいつか消すかも
	,_penalty_func:&Option<Box<dyn LinearRegressionPenalty>>
	,betas_buff:&mut Vec<f64>
){
	
	
	energyset.update_distance();
	let f0:f64 = energyset.calc_energy(atom_level_energy);
	for ii in 0..betas.len(){
		let prevval:f64;
		if betas[ii].2 == X{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_x();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_x(prevval+betas[ii].0+FPRIME_EPSILON);
		}else if betas[ii].2 == Y{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_y();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_y(prevval+betas[ii].0+FPRIME_EPSILON);
		}else{
			prevval = energyset.evoef2_env.md_envset.atoms[betas[ii].1].get_z();
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_z(prevval+betas[ii].0+FPRIME_EPSILON);
		}
		energyset.update_distance_one(betas[ii].1);
		let ff:f64 = energyset.calc_energy(atom_level_energy);

		if betas[ii].2 == X{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_x(prevval+betas[ii].0-FPRIME_EPSILON);
		}else if betas[ii].2 == Y{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_y(prevval+betas[ii].0-FPRIME_EPSILON);
		}else{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_z(prevval+betas[ii].0-FPRIME_EPSILON);
		}
		energyset.update_distance_one(betas[ii].1);
		let ff2:f64 = energyset.calc_energy(atom_level_energy);

		if ff > f0 && ff2 > f0{//どっちに動いても損失が大きくなる場合は微分はゼロにする
			betas_buff[ii] = 0.0;
		}else{
			if ff < f0 && ff2 < f0 && ff2 < ff {
				betas_buff[ii] = (f0-ff2)/FPRIME_EPSILON;
			}else{
				betas_buff[ii] = (ff-f0)/FPRIME_EPSILON;
			}
		}
		
		if betas[ii].2 == X{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_x(prevval);
		}else if betas[ii].2 == Y{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_y(prevval);
		}else{
			energyset.evoef2_env.md_envset.atoms[betas[ii].1].set_z(prevval);
		}
		energyset.update_distance_one(betas[ii].1);
	}
}


//beta を掛けて足しただけの値を得る
pub fn calc_beta_res(betas:&Vec<f64>,x:&Vec<f64>)->f64{
	assert_eq!(betas.len(),x.len()+1);
	let mut res:f64 = betas[0];
	for ii in 0..x.len(){
		res += x[ii]*betas[ii+1];
	}
	return res;
}


pub struct LBFGSResult{
	pub betas:Vec<f64>,
	pub min_loss:f64,
	pub converged:bool
}



//係数を返す。第 0 は定数
pub fn run_lbfgs(
 energyset:&mut pp_energy::PPEnergySet
,atom_indices:Vec<usize>
,iter_max:usize
,stop_threshold:f64
,stop_threshold_tolerance:usize
,penalty_func:&Option<Box<dyn LinearRegressionPenalty>>
)->LBFGSResult{
	let num_betas:usize = atom_indices.len()*3;
	let num_atoms:usize = energyset.evoef2_env.md_envset.atoms.len();
	let mut atom_level_energy:Vec<f64> = vec![0.0;num_atoms];

	let mut betas:Vec<(f64,usize,usize)> = vec![(0.0,0,0);num_betas];
	let mut betas_buff:Vec<f64> = vec![0.0;num_betas];
	for (ii,bb) in atom_indices.iter().enumerate(){
		betas[ii*3+X] = (0.0,*bb,X);
		betas[ii*3+Y] = (0.0,*bb,Y);
		betas[ii*3+Z] = (0.0,*bb,Z);
	}

	let mut minimumloss_betas:Vec<f64> = vec![0.0;num_betas];

	let mut minimumloss:f64 = f64::INFINITY;
	
	let mut betas_tmp:Vec<(f64,usize,usize)> = betas.clone();
	let mut buff_array:Vec<f64> = vec![0.0;num_betas];

	let mut direc_current:Vec<f64> = vec![0.0;num_betas];
	let mut grad_prev:Vec<f64> = vec![0.0;num_betas];

	let mut lthresholdcount = 0_usize;
	let lbfgs_m:usize = 10;
	let mut xstep:Vec<Vec<f64>> = vec![];

	let mut graddiff:Vec<Vec<f64>> = vec![];
	let mut prevfailed:bool = false;
	let mut convergedflag = false;

	let mut preddiff_sum_prev:f64 = 100000.0;

	//println!("original:{}",calc_enegy_with_beta(&betas,&mut buff_array, md_envset, energyset,&mut atom_level_energy));

	for ii in 0..direc_current.len(){
		direc_current[ii] = 0.0;
	}
	'tloop: for icou in 0..iter_max{
		let preddiff_sum:f64;
		
		if icou == 0{
			calc_fprime(&betas,energyset,&mut atom_level_energy,penalty_func,&mut betas_buff);
			for ii in 0..betas_buff.len(){
				grad_prev[ii] = betas_buff[ii];//初期値どうすればいいのか
				direc_current[ii] = betas_buff[ii]*-1.0;
			}
		}else{
			if prevfailed{
				calc_fprime(&betas,energyset,&mut atom_level_energy,penalty_func,&mut betas_buff);
				if xstep.len() < lbfgs_m{
					for ii in 0..betas_buff.len(){
						direc_current[ii] = -1.0*betas_buff[ii];
					}
				}else{
					for ii in 0..betas_buff.len(){
						betas_buff[ii] = -1.0*betas_buff[ii];
					}
					direc_current = calc_next_gradient_lbfgs(&xstep,&graddiff,&betas_buff);
				}
			}else{
				for ii in 0..grad_prev.len(){
					grad_prev[ii] = -1.0*grad_prev[ii];
				}
				direc_current = calc_next_gradient_lbfgs(&xstep,&graddiff,&grad_prev);
				for ii in 0..grad_prev.len(){
					grad_prev[ii] = -1.0*grad_prev[ii];
				}
			}
		}

		//println!("{} {}",preddiff_sum,calc_binomial_logistic_loss(&betas,&samples));

		//黄金探索。もっときれいに書けそうだが。。。
		let mut x1:f64 = 0.0;
		let mut x2:f64 = FPRIME_EPSILON;
		
		let mut x4:f64 = x2*GOLDEN_NUMBER+x2;
		let mut x3:f64 = (x4*GOLDEN_NUMBER+x1)/(GOLDEN_NUMBER+1.0);

		//println!("{}",x3-x4+1.0);

		let mut y1 = calc_enegy_with_beta(&betas,&mut buff_array, energyset,&mut atom_level_energy);
		
		for ii in 0..direc_current.len(){
			betas_tmp[ii].0 = x2*direc_current[ii]+betas[ii].0;
		}
		let mut y2 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
		
		for ii in 0..direc_current.len(){
			betas_tmp[ii].0 = x3*direc_current[ii]+betas[ii].0;
		}
		let mut y3:f64;
		
		for ii in 0..direc_current.len(){
			betas_tmp[ii].0 = x4*direc_current[ii]+betas[ii].0;
		}
		let mut y4 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
		
		let mut ccount:usize = 0;
		let mut prestepfailed:bool = false;

		while y1 < y2{//x2 がサンプリング点の内で最小であることを確認する
			ccount+=1;
			if ccount > GOLDEN_LOOPCHECK{//失敗。何か間違えているのだろうか・・・。
				//eprintln!("{} {}---???",x4,x2);
				prestepfailed  = true;
				break;
			}
			x2 = x2*0.5;
			x4 = x2*GOLDEN_NUMBER+x2;
			
			for ii in 0..direc_current.len(){
				betas_tmp[ii].0 = x2*direc_current[ii]+betas[ii].0;
			}
			y2 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
			
			for ii in 0..direc_current.len(){
				betas_tmp[ii].0 = x4*direc_current[ii]+betas[ii].0;
			}
			y4 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
			
		}

		let mut ccount:usize = 0;
		while y4 <= y2{//最小値をまたぐように x4 を設定する
			ccount+=1;
			if ccount > GOLDEN_LOOPCHECK{
				//多分収束している
				//eprintln!("{} {} already converged?",x4,x2);
				
				prestepfailed  = true;
				break;
			}
			x2 = x4;y2 = y4;x4 = x2+x2*GOLDEN_NUMBER;
			for ii in 0..direc_current.len(){
				betas_tmp[ii].0 = x4*direc_current[ii]+betas[ii].0;
			}
			y4 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
			
		}
		
		x3 = (x4*GOLDEN_NUMBER+x1)/(GOLDEN_NUMBER+1.0);

		for ii in 0..direc_current.len(){betas_tmp[ii].0 = x3*direc_current[ii]+betas[ii].0;}
		y3 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);

		if !prestepfailed{
			for _ in 0..10000{
				if y2 > y3{
					
					x1 = x2;
					y1 = y2;
					x2 = x3;
					y2 = y3;

					x3 = (x4*GOLDEN_NUMBER+x1)/(GOLDEN_NUMBER+1.0);
					for ii in 0..direc_current.len(){betas_tmp[ii].0 = x3*direc_current[ii]+betas[ii].0;}
					y3 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
					
				}else{
					if (y1-y2).abs().max((y1-y3).abs()) < GOLDEN_BREAK_THRESHOLD
					&& (y4-y3).abs().max((y4-y2).abs()) < GOLDEN_BREAK_THRESHOLD
					&& (y2-y3).abs() < GOLDEN_BREAK_THRESHOLD{
						break;
					}
					
					x4 = x3;
					y4 = y3;
					x3 = x2;
					y3 = y2;

					x2 = (x1*GOLDEN_NUMBER+x4)/(GOLDEN_NUMBER+1.0);
					for ii in 0..direc_current.len(){betas_tmp[ii].0 = x2*direc_current[ii]+betas[ii].0;}
					y2 = calc_enegy_with_beta(&betas_tmp,&mut buff_array, energyset,&mut atom_level_energy);
					
				}
				//微妙な差で エネルギーが変化する場合無限ループに入るので Break する
				if x4-x1 < GOLDEN_BREAK_THRESHOLD_X{
					break;
				}
			}
		}
		
		if y2 > y3{
			for ii in 0..betas.len(){
				betas_tmp[ii].0 = direc_current[ii]*x3+betas[ii].0;
			}
		}else{
			for ii in 0..betas.len(){
				betas_tmp[ii].0 = direc_current[ii]*x2+betas[ii].0;
			}
		}
		
		//atom num とかで割っていたが、stop_threshold は呼び出し側で制御できるのでやめた
		if  minimumloss-y2.min(y3) < stop_threshold{
			lthresholdcount += 1;
		}else{
			lthresholdcount = 0;
		}
		if y1.is_nan() || y2.is_nan() || y3.is_nan() || y4.is_nan() || minimumloss.is_nan(){
			panic!("LBFGS FAILED!");
		}
		//if icou%100 == 99{
			eprintln!("{} loss: {}->{}",icou,minimumloss,y2.min(y3));
		//}
		if y2.min(y3) < minimumloss {
			minimumloss = y2.min(y3);
			for ii in 0..betas.len(){
				minimumloss_betas[ii] = betas_tmp[ii].0;
			}
		}
		
		if lthresholdcount >= stop_threshold_tolerance{
			//eprintln!("converged!");
			convergedflag = true;
			break 'tloop;
		}

		calc_fprime(&betas_tmp,energyset,&mut atom_level_energy,penalty_func,&mut betas_buff);
		let mut innercheck:f64 = 0.0;
		for ii in 0..betas.len(){
			innercheck += (betas_tmp[ii].0-betas[ii].0)*(betas_buff[ii]-grad_prev[ii]);
		}

		if innercheck.abs() > INNERCHECK_EPSILON {
			let mut xstep_m:Vec<f64>;
			let mut graddiff_m:Vec<f64>;
			if xstep.len() >= lbfgs_m{
				xstep_m = xstep.remove(0);
				graddiff_m = graddiff.remove(0);
			}else{
				xstep_m = vec![0.0;betas.len()];
				graddiff_m = vec![0.0;betas.len()];
			}
			for ii in 0..betas.len(){
				xstep_m[ii] = betas_tmp[ii].0 -betas[ii].0;//上で足しているので冗長
				graddiff_m[ii] =  betas_buff[ii]-grad_prev[ii];
				grad_prev[ii] = betas_buff[ii];
				betas[ii] = betas_tmp[ii];
			}
			graddiff.push(graddiff_m);
			xstep.push(xstep_m);
			prevfailed = false;
		}else{
			//eprintln!("failed!");//二回連続 Fail したらロジック上改善しないかも
			prevfailed = true;
		}

		preddiff_sum =  calc_enegy_with_beta(&betas,&mut buff_array, energyset,&mut atom_level_energy);
		if false{
			eprintln!("{}->{}",preddiff_sum_prev,preddiff_sum);
		}
		//eprintln!("chk:{}->{}",preddiff_sum_prev,preddiff_sum);
		preddiff_sum_prev = preddiff_sum;
		
	}
	if !convergedflag{
		//eprintln!("not converged!");
	}
	return LBFGSResult{betas:minimumloss_betas,min_loss:minimumloss,converged:convergedflag};
}


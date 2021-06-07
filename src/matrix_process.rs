



#[allow(unused_imports)]
use rand::SeedableRng;
#[allow(unused_imports)]
use rand::rngs::StdRng;
#[allow(unused_imports)]
use rand::Rng;
#[allow(dead_code)]
const ERROR_PANIC_THRESHOLD:f64 = 0.0000000001;//Test の時に使う誤差の許容度


//[row][col]
pub fn matrix_multi(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
    if x[0].len() != y.len(){
		panic!("The num col in x & numrow in y must be equal x:{},{};y:{},{}",x.len(),x[0].len(),y.len(),y[0].len());
    }
    
    return matrix_multi_slice(x,y,x.len(),y[0].len(),y.len());
}

pub fn matrix_t(x:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
    return matrix_t_slice(x,x.len(), x[0].len());
}

//x の row の数、y の col の数、合成される要素の数を渡す。
pub fn matrix_multi_slice(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>,xlen:usize,ylen:usize,plen:usize)->Vec<Vec<f64>>{
	
	let rnum = xlen;
	let cnum = ylen;
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;cnum];rnum];
	for ii in 0..rnum{
		for jj in 0..cnum{
			for rr in 0..plen{
				ret[ii][jj] += x[ii][rr]*y[rr][jj];
			}
		}
	}
	return ret;
}

pub fn matrix_t_slice(x:&Vec<Vec<f64>>,rlen:usize,clen:usize)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;rlen];clen];
	for rr in 0..rlen{
		for cc in 0..clen{
			ret[cc][rr] = x[rr][cc];
		}
	}
	return ret;
}

pub fn matrix_multi_old(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
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

pub fn matrix_t_old(x:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
    let mut ret:Vec<Vec<f64>> = vec![vec![0.0;x.len()];x[0].len()];
	for rr in 0..x.len(){
		for cc in 0..x[0].len(){
			ret[cc][rr] = x[rr][cc];
		}
	}
    return ret;
}

pub fn matrix_dot(a:&Vec<Vec<f64>>,b:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
    assert_eq!(a.len(),b.len());
    assert_eq!(a[0].len(),b[0].len());
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;a.len()];b[0].len()];
	for rr in 0..a.len(){
		for cc in 0..b[0].len(){
			ret[rr][cc] = a[rr][cc]*b[rr][cc];
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

pub fn matrix_subtract(x:&Vec<Vec<f64>>,y:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
	let mut ret:Vec<Vec<f64>> = vec![vec![0.0;x[0].len()];x.len()];
	if x.len() != y.len() || x[0].len() != y[0].len(){
		panic!("shapes must be the same x:{},{};y:{},{}",x.len(),x[0].len(),y.len(),y[0].len());
	}
	for rr in 0..x.len(){
		for cc in 0..x[0].len(){
			ret[rr][cc] = x[rr][cc]-y[rr][cc];
		}
	}
	return ret;
}




//逆行列の計算に当該行が使えるかチェックする使えない場合使えそうな行とスワップする。
//使える行が無い場合 -1 が返る。
pub fn rowcheck_inv(px:&mut Vec<Vec<f64>>,pr:&mut Vec<Vec<f64>>,rcou:usize)->i64{
    if px[rcou][rcou] != 0.0{
        return rcou as i64;
    }
    let mut rz:i64 = -1;
    let n:usize = px.len();
    for rr in (rcou+1)..n{
        if px[rr][rcou] != 0.0{
            rz = rr as i64;
            break;
        }
    }
    if rz < 0_i64{
        return rz;
    }
    for rr in 0..n{
        let zpx = px[rcou][rr];
        let zpr = pr[rcou][rr];
        px[rcou][rr] = px[rz as usize][rr];
        pr[rcou][rr] = pr[rz as usize][rr];
        px[rz as usize][rr] = -1.0*zpx;
        pr[rz as usize][rr] = -1.0*zpr;
    }
    return rz;
}
pub fn matrix_inv(x:&Vec<Vec<f64>>)->Option<Vec<Vec<f64>>>{
    assert_eq!(x.len(),x[0].len());
    let n:usize = x.len();
    let mut ret:Vec<Vec<f64>> = vec![vec![0.0;n];n];
    let mut px:Vec<Vec<f64>> = x.clone();
    for ii in 0..n{
        ret[ii][ii] = 1.0;
    }
    for rr in 0..n{
        let chk:i64 = rowcheck_inv(&mut px,&mut ret,rr);
        if chk < 0{
            return None;
        }
        for kk in 0..n{
            let scale:f64 = px[kk][rr]/px[rr][rr];
            if kk == rr{
                continue;
            }
            for r2 in rr..n{
                px[kk][r2] -= px[rr][r2]*scale;
            }
            for r2 in 0..n{
                ret[kk][r2] -= ret[rr][r2]*scale;
            }
        }
        
        let ss:f64 = px[rr][rr];
        for r2 in rr..n{
            px[rr][r2] /= ss;
        }
        for r2 in 0..n{
            ret[rr][r2] /= ss;
        }
    }

    /*
    for rr in 0..n{
        for cc in 0..n{
            if rr == cc{
                assert_eq!(px[rr][cc],1.0);
            }else{
                assert_eq!(px[rr][cc],0.0);
            }
        }
    }
    */
    return Some(ret);
}


pub fn matrix_print(vv:&Vec<Vec<f64>>){
    let rnum:usize = vv.len();
    let cnum:usize = vv[0].len();

    for i in 0..rnum{
        for j in 0..cnum{
            if j > 0{
                print!("\t");
            }
            print!("{}",vv[i][j]);
        }
        print!("\n");
    }
}

pub fn l1sqrt(v:&Vec<Vec<f64>>)->f64{
    let ncol:usize = v.len();
    let nrow:usize = v[0].len();
    let mut ret:f64 = 0.0;
    for cc in 0..ncol{
        for rr in 0..nrow{
            ret += v[cc][rr]*v[cc][rr];
        }
    }
    if ret == 0.0{
        return 0.0;
    }
    return ret.sqrt();
}

pub fn eye(u:usize)->Vec<Vec<f64>>{
    let mut ret:Vec<Vec<f64>> = vec![vec![0.0;u];u];
    for uu in 0..u{
        ret[uu][uu] = 1.0;
    }
    return ret;
}

//
//http://takashiijiri.com/study/miscs/QRfactorization.htm
pub fn house_qr_decomp(mat:&Vec<Vec<f64>>)->(Vec<Vec<f64>>,Vec<Vec<f64>>){
    let ncol:usize = mat[0].len();
    let nrow:usize = mat.len();
    
    let mut qmat:Vec<Vec<f64>> = eye(nrow);
    let mut rmat:Vec<Vec<f64>> = mat.clone();
    for cc in 0..(ncol-1){
        let xx:Vec<Vec<f64>> = rmat[cc..nrow].iter().map(|m|vec![m[cc]]).collect();
        let mut yy:Vec<Vec<f64>> = vec![vec![0.0];xx.len()];
        yy[0][0] = l1sqrt(&xx);
        if yy[0][0] == 0.0{
            //この処理は正しいのか？ToDo 調査。
            continue;
        }
        let mut uu:Vec<Vec<f64>> = matrix_subtract(&xx,&yy);
        let su:f64 = l1sqrt(&uu);
        if su == 0.0{
            continue;
        }

        uu = matrix_scalarmulti(&uu,1.0/su);
        let hmat:Vec<Vec<f64>> = matrix_add(
            &eye(uu.len()),
            &matrix_scalarmulti(
                &matrix_multi(
                    &uu
                    ,&matrix_t(&uu)
                )
            ,-2.0)
        );
        let mut hmat_c:Vec<Vec<f64>> = eye(nrow);
        for rrr in cc..nrow{
            for ccc in cc..nrow{
                hmat_c[rrr][ccc] = hmat[rrr-cc][ccc-cc];
            }
        }
        rmat = matrix_multi(&hmat_c,&rmat);
        qmat = matrix_multi(&hmat_c,&qmat);
    }
    return (matrix_t(&qmat),rmat);
}


#[test]
fn matrix_invtest(){
    //https://mathtrain.jp/inversematrix
    let a:Vec<Vec<f64>> = vec![
        vec![1.0 ,1.0,-1.0],
        vec![-2.0,0.0, 1.0],
        vec![0.0 ,2.0, 1.0]
    ];
    let inva :Vec<Vec<f64>> = matrix_inv(&a).unwrap();
    let b:Vec<Vec<f64>> = vec![
        vec![-1.0/2.0 ,-3.0/4.0,1.0/4.0],
        vec![1.0/2.0,1.0/4.0, 1.0/4.0 ],
        vec![-1.0 ,-1.0/2.0, 1.0/2.0 ]
    ];
    assert_eq!(inva,b);


    let a:Vec<Vec<f64>> = vec![
        vec![1.0 ,1.0,-1.0],
        vec![-2.0,0.0, 1.0],
        vec![-2.0,0.0, 1.0],
    ];
    
    assert_eq!(matrix_inv(&a),None);


    let a:Vec<Vec<f64>> = vec![
        vec![1.0 ,1.0,-1.0],
        vec![-1.0,-1.0,2.0],
        vec![-2.0,0.0, 1.0],
    ];
    let b:Vec<Vec<f64>> = vec![
        vec![0.5 ,0.5,-0.5],
        vec![1.5,0.5,0.5],
        vec![1.0,1.0, 0.0],
    ];
    assert_eq!(matrix_inv(&a).unwrap(),b);

    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    for _ in 0..100{  
        let m:usize = rgen.gen_range(2..15);
        let mut vv:Vec<Vec<f64>> = vec![vec![0.0;m];m];
        for i in 0..m{
            for j in 0..m{
                vv[i][j] = rgen.gen_range(-100.0..100.0);
            }
        }
        let res_:Option<Vec<Vec<f64>>> = matrix_inv(&vv);
        if let None = res_{
            matrix_print(&vv);
            std::process::exit(0);   
        }else{
            let res:Vec<Vec<f64>> = res_.unwrap();
            let imat:Vec<Vec<f64>> = matrix_multi(&res,&vv);
            for i in 0..m{
                for j in 0..m{
                    if i == j{
                        if (imat[i][j] - 1.0).abs()  > ERROR_PANIC_THRESHOLD{
                            eprintln!("Large absolute error was found.");
                            eprintln!("source===");
                            matrix_print(&vv);
                            eprintln!("result===");
                            matrix_print(&res);
                            eprintln!("imat===");
                            matrix_print(&imat);
                            panic!();
                        }
                    }else{
                        if imat[i][j].abs()  > ERROR_PANIC_THRESHOLD{
                            eprintln!("Large absolute error was found.");
                            eprintln!("source===");
                            matrix_print(&vv);
                            eprintln!("result===");
                            matrix_print(&res);
                            eprintln!("imat===");
                            matrix_print(&imat);
                            
                            panic!();
                        }
                    }
                }
            }
            //println!("{:?}",matrix_multi(&vv, &(res.unwrap())));
        }
    }
}

#[test]
fn qrdecomptest(){
    let chkk = vec![
        vec![12.0,-51.0,4.0],
        vec![6.0,167.0,-68.0],
        vec![-4.0,24.0,-41.0]
    ];
    let res = house_qr_decomp(&chkk);
    println!("{:?}\n{:?}",res.0,res.1);
}


#[test]
fn qrdecomptest2(){
    let chkk = vec![
        vec![3.0, -1.5, -8.0],
        vec![0.0, -0.25, 2.0],
        vec![0.0, 0.1875, 1.0]
    ];
    let res = house_qr_decomp(&chkk);
    println!("{:?}\n{:?}",res.0,res.1);
}


#[test]
fn matrix_slicetest(){

    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    for _ in 0..100{  
        let l:usize = rgen.gen_range(2..15);
        let m:usize = rgen.gen_range(2..15);
        let n:usize = rgen.gen_range(2..15);
        let mut vv:Vec<Vec<f64>> = vec![vec![0.0;n];l];
        let mut vv2:Vec<Vec<f64>> = vec![vec![0.0;m];n];
        let mut vv_b:Vec<Vec<f64>> = vec![vec![0.0;n*2];l*2];
        let mut vv2_b:Vec<Vec<f64>> = vec![vec![0.0;m*3];n];
        for v in vv_b.iter_mut(){
            for i in v.iter_mut(){
                *i  = rgen.gen_range(-100.0..100.0);
            }    
        }

        for v in vv2_b.iter_mut(){
            for i in v.iter_mut(){
                *i  = rgen.gen_range(-100.0..100.0);
            }    
        }

        for i in 0..m.max(l){
            for j in 0..n{
                if i < l{
                    vv[i][j] = rgen.gen_range(-100.0..100.0);
                    vv_b[i][j] = vv[i][j];
                }
                if i < m{
                    vv2[j][i] = rgen.gen_range(-100.0..100.0);
                    vv2_b[j][i] = vv2[j][i];
                }
            }
        }
        assert_eq!(matrix_multi(&vv,&vv2),matrix_multi_old(&vv,&vv2));
        assert_eq!(matrix_multi(&vv,&vv2),matrix_multi_slice(&vv_b,&vv2_b,l,m,n));
        assert_eq!(matrix_t(&vv),matrix_t_old(&vv));
        assert_eq!(matrix_t_slice(&vv_b,l,n),matrix_t(&vv));
        
    }


}


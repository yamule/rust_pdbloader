
use std::thread;
use std::sync::mpsc;
use super::pdbdata;
use super::geometry::Vector3D;
#[allow(dead_code,unused_imports)]
use super::sequence_alignment;


#[allow(dead_code,unused_imports)]
use rand::SeedableRng;
#[allow(dead_code,unused_imports)]
use rand::rngs::StdRng;
#[allow(dead_code,unused_imports)]
use rand::Rng;

#[allow(dead_code,unused_imports)]
use super::debug_env;

use std::collections::{HashMap,HashSet};

#[allow(dead_code,unused_imports)]
const ERROR_THRESHOLD:f64 = 0.000000001;
//Kabsh のアルゴリズムとか SVD とかが良く分かんなかったので QR DECOMP で ROTATION MATRIX を出している。


#[allow(non_camel_case_types)]
pub enum DistanceAlignmentType{
    SEED,
    DP
}



pub fn get_blosum62_positive_pair()->HashMap<String,HashSet<String>>{
    let mut lin:Vec<String> = Vec::new();
    //https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
    //#  Matrix made by matblas from blosum62.iij
    //#  * column uses minimum score
    //#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    //#  Blocks Database = /data/blocks_5.0/blocks.dat
    //#  Cluster Percentage: >= 62
    //#  Entropy =   0.6979, Expected =  -0.5209
    lin.push("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *".to_string());
    lin.push("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4".to_string());
    lin.push("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4".to_string());
    lin.push("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4".to_string());
    lin.push("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4".to_string());
    lin.push("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4".to_string());
    lin.push("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4".to_string());
    lin.push("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4".to_string());
    lin.push("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4".to_string());
    lin.push("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4".to_string());
    lin.push("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4".to_string());
    lin.push("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4".to_string());
    lin.push("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4".to_string());
    lin.push("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4".to_string());
    lin.push("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4".to_string());
    lin.push("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4".to_string());
    lin.push("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4".to_string());
    lin.push("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4".to_string());
    lin.push("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4".to_string());
    lin.push("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4".to_string());
    lin.push("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4".to_string());
    lin.push("B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4".to_string());
    lin.push("Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4".to_string());
    lin.push("X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4".to_string());
    lin.push("* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1".to_string());
    let mut lincount:i64 = -1;
    let mut ret:HashMap<String,HashSet<String>> = HashMap::new();
    let mut string_to_index:HashMap<String,usize>  = HashMap::new();
    let mut index_to_string:Vec<String>  = vec![];
    for (_ii,line) in lin.into_iter().enumerate() {
        let bs = line.trim();
        let ptt:Vec<String> = bs.split_whitespace().map(|m| m.to_string()).collect();
        let lcc:char = ptt[0].chars().next().unwrap();
        if lcc == '#'{
            continue;
        }
        lincount+=1;
        if lincount == 0{
            for (_jj,pp) in ptt.into_iter().enumerate(){
                if string_to_index.contains_key(&pp){
                   panic!("{} was already found.",pp);
                }
                let csiz = string_to_index.keys().len();
                string_to_index.insert(pp.clone(),csiz);
                index_to_string.push(pp.clone());
                ret.insert(pp,HashSet::new());
            }
        }else{

            for ll in 1..ptt.len(){
                let ss = match ptt[ll].parse::<f32>(){
                    Ok(x)=>{
                        x
                    },
                    _=>{
                        eprintln!("{} can not be parsed! zero was assigned",ptt[ll]);
                        0.0
                    }
                };
                if ss > 0.0{
                    ret.get_mut(&ptt[0]).unwrap().insert(index_to_string[ll-1].clone());
                }
            }
        }
    }
    println!("{:?}",ret);
    return ret;
}


//sqrt する前の distance を返す。
pub fn distance2(x1:f64,y1:f64,z1:f64,x2:f64,y2:f64,z2:f64)->f64{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

pub fn estimate_aligned_atom_distance(target:&Vec<Vec<f64>>,template:&Vec<Vec<f64>>,x_to_y_indexmap:&Vec<i64>,xi:usize,yi:usize)->f64{
    let mut lcou:usize = 0;
    let mut dist_e:f64 = 0.0;
    for (xxi,yyi_) in x_to_y_indexmap.iter().enumerate(){
        if *yyi_ < 0{
            continue;
        }
        let yyi = *yyi_ as usize;
        if (xi as i64 - xxi as i64)*(yi as i64 - yyi as i64) <= 0{
            //前後関係が逆である場合並んでいるとみなさない
            continue;
        }
        dist_e += (target[xi][xxi]-template[yi][yyi]).abs()/2.0;
        lcou += 1;
    }
    if lcou == 0{
        return 1000.0;
    }
    return dist_e/(lcou as f64);
}


pub fn calc_dist_tmscore_gdtts(target:&Vec<Vec<f64>>,template:&Vec<Vec<f64>>,x_to_y_indexmap:&Vec<i64>,tm_cutoff:f64,tm_d0:f64,tm_lnorm:usize)->(f64,f64){
    let mut gdtts:Vec<usize> = vec![0;4];
    let mut tm_sum:f64 = 0.0;
    let d02:f64 = tm_d0*tm_d0;
    let tm_cutoff2:f64 = tm_cutoff*tm_cutoff;
    for (xi,yi_) in x_to_y_indexmap.iter().enumerate(){
        if *yi_ < 0{
            continue;
        }
        let yi = *yi_ as usize;
        let dist_e = estimate_aligned_atom_distance(target, template, x_to_y_indexmap, xi, yi);
        let dist_e2 = dist_e*dist_e;
        /*
        GDT_TS
        GDT_TS - GlobalDistanceTest_TotalScore
        GDT_TS = (GDT_P1 + GDT_P2 + GDT_P4 + GDT_P8)/4,
        where GDT_Pn denotes percent of residues under distance cutoff <= nÅ
        GDT_HA
        GDT_HA - GDT High Accuracy
        GDT_HA = (GDT_P0.5 + GDT_P1 + GDT_P2 + GDT_P4)/4,
        where GDT_Pn denotes percent of residues under distance cutoff <= nÅ 
        */

        if dist_e2 < 64.0{
            gdtts[3] += 1;
            if dist_e2 < 16.0{
                gdtts[2] += 1;
                if dist_e2 < 4.0{
                    gdtts[1] += 1;
                    if dist_e2 < 1.0{
                        gdtts[0] += 1;
                    }
                }
            }
        }
        
        if dist_e2 <= tm_cutoff2{
            tm_sum += 1.0/(1.0+dist_e2/d02);
        }

        /*
         Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
        if (Lnorm<=19)                    //update 15-->19
            d0=0.168;                   //update 0.5-->0.168
        else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
        D0_MIN=d0+0.8;              //this should be moved to above
        d0=D0_MIN;                  //update: best for search    
        */
        
    }
    tm_sum /= tm_lnorm as f64;
    let mut gsum :f64 = 0.0;
    for ii in 0..gdtts.len(){
        gsum += gdtts[ii] as f64/4.0;
    }
    return (tm_sum,gsum);
}


pub fn align_distance(xvec:&Vec<Vec<f64>>,yvec:&Vec<Vec<f64>>
    ,query_aa:&Vec<String>
    ,template_aa:&Vec<String>
    ,seed_length:usize
    ,tmscore_break:f64,dist_cutoff:f64
    ,aligntype:DistanceAlignmentType)->Option<DistanceAlignmentResult>{

    let lnorm:usize = xvec.len().min(yvec.len());
    let x_length:usize = xvec.len();
    let t_length:usize = yvec.len();
    let d0:f64 = if lnorm <=19{
        0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    
    let (index_map,_scores):(Vec<i64>,(f64,f64)) = 
    match aligntype{
        DistanceAlignmentType::SEED => {
            align_seed(&xvec,&yvec,query_aa,template_aa,100,d0_search,d0,lnorm,seed_length,tmscore_break)
        },
        DistanceAlignmentType::DP => {
            align_dp(&xvec,&yvec,query_aa,template_aa,d0_search,d0,lnorm)
        }
    };


    let mut ret = DistanceAlignmentResult::new();
    
    //short の方のスコア計算
    let lnorm:usize = x_length;
    let d0:f64 = if lnorm <=19{
        0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    let scores_chain1 = calc_dist_tmscore_gdtts(xvec,yvec, &index_map,d0_search,d0,lnorm);
    

    //long の方のスコア計算
    let lnorm:usize = t_length;
    let d0:f64 = if lnorm <=19{0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    let scores_chain2 = calc_dist_tmscore_gdtts(xvec,yvec, &index_map,d0_search,d0,lnorm);

    ret.tmscore_chain1 = scores_chain1.0;
    ret.tmscore_chain2 = scores_chain2.0;

    assert_eq!(scores_chain1.1,scores_chain2.1);
    ret.gdt_ts = scores_chain1.1;

    ret.dist_cutoff = dist_cutoff;
    
    
    let mut rmsd = 0.0;
    let mut aligned_num:usize = 0;
    for (xi,yi_) in index_map.iter().enumerate(){
        if *yi_ < 0{
            continue;
        }
        let yi = *yi_ as usize;
        let dist_e = estimate_aligned_atom_distance(xvec,yvec, &index_map, xi, yi);
        rmsd += dist_e*dist_e;
        aligned_num += 1;
    }
    if aligned_num > 0{

        ret.rmsd_aligned = (rmsd/aligned_num as f64).sqrt();
    }
    let mut xstr:Vec<String> = vec![];
    let mut ystr:Vec<String> = vec![];
    let mut lastcount_x:i64 = -1;
    let mut lastcount_y:i64 = -1;
    let mut num_aligned:usize = 0;
    for ii in 0..index_map.len(){
        if index_map[ii] < 0{
            continue;
        }
        while lastcount_x < ii as i64 -1{
            lastcount_x += 1;
            xstr.push(query_aa[lastcount_x as usize].clone());
            ystr.push("-".to_string());
        }
        while lastcount_y < index_map[ii] as i64 -1{
            lastcount_y += 1;
            ystr.push(template_aa[lastcount_y as usize].clone());
            xstr.push("-".to_string());
        }
        lastcount_y += 1;
        lastcount_x += 1;
        num_aligned += 1;
        xstr.push(query_aa[lastcount_x as usize].clone());
        ystr.push(template_aa[lastcount_y as usize].clone());
    }
    ret.index_map = index_map;
    while lastcount_x < x_length as i64 -1{
        lastcount_x += 1;
        xstr.push(query_aa[lastcount_x as usize].clone());
        ystr.push("-".to_string());
    }
    while lastcount_y < t_length as i64 -1{
        lastcount_y += 1;
        ystr.push(template_aa[lastcount_y as usize].clone());
        xstr.push("-".to_string());
    }
    ret.aligned_chain1 = xstr;
    ret.aligned_chain2 = ystr;
    ret.num_aligned = num_aligned as i64;
    return Some(ret);
}





#[derive(Debug)]
pub struct DistanceAlignmentResult{
    pub tmscore_chain1:f64,
    pub tmscore_chain2:f64,
    pub gdt_ts:f64,
    pub index_map:Vec<i64>,
    
    pub dist_cutoff:f64,
    pub rmsd_aligned:f64,
    pub num_aligned:i64,
    pub aligned_chain1:Vec<String>,
    pub aligned_chain2:Vec<String>,

}


impl DistanceAlignmentResult{
    pub fn new()->DistanceAlignmentResult{
        return DistanceAlignmentResult{
            tmscore_chain1:0.0,
            tmscore_chain2:0.0,
            gdt_ts:0.0,
            index_map:vec![],
            dist_cutoff:0.0,
            rmsd_aligned:0.0,
            num_aligned:0,
            aligned_chain1:vec![],
            aligned_chain2:vec![]
        };
    }
}

#[allow(dead_code)]
fn filt_short_segments(indexmap:&mut Vec<(i64,f64)>){
    let ilen:usize = indexmap.len();
    let mut ii:usize = 0;
    while ii < ilen{
        let mut prevy:i64 = indexmap[ii].0;
        let mut prevx:usize = ii;
        let mut llength:usize = 1;
        let mut jj:usize = ii+1;
        while jj < ilen{
            if indexmap[jj].0 < 0{
                jj += 1;
                continue;
            }
            if indexmap[jj].0-prevy != (jj-prevx) as i64{
                break;
            }else{
                prevx = jj;
                prevy = indexmap[jj].0;
                llength+=1;
            }
            jj += 1;
        }
        if llength < 3{
            for pp in ii..jj{
                indexmap[pp].0 = -1;
                indexmap[pp].1 = 1000000.0;//tekitou
            }
        }
        ii = jj;
    }
}

//マッピングに前後の矛盾がある場合に矛盾のないようにする
//f64 には大きい方が悪い値が入っている RMSD とか。
fn filt_discrepancy(indexmap:&mut Vec<(i64,f64)>){
    let plen:usize = indexmap.len();
    let mut penalty:Vec<(usize,f64)> = vec![(0,0.0);plen];
    loop{

        for ii in 0..plen{
            penalty[ii].0 = ii;
            penalty[ii].1 = 0.0;
        }
        for ii in 0..plen{
            if indexmap[ii].0 < 0{
                continue;
            }
            for jj in ii+1..plen{
                if indexmap[jj].0 < 0{
                    continue;
                }
                if indexmap[ii].0 >= indexmap[jj].0{
                    penalty[ii].1 += 1.0;
                    penalty[jj].1 += 1.0;
                }
            }
        }
        penalty.sort_by(|a,b|a.1.partial_cmp(&b.1).unwrap());
        penalty.reverse();
        if penalty[0].1 == 0.0{
            break;
        }

        let mut candidates:Vec<usize> = vec![penalty[0].0];
        for ii in 1..plen{
            if penalty[ii].1 == penalty[0].1{
                candidates.push(penalty[ii].0);
            }
        }
        if candidates.len() == 1{
            indexmap[candidates[0]] = (-1,indexmap[candidates[0]].1);
        }else{
            let mut mmax = indexmap[candidates[0]].1;
            let mut maxi = candidates[0];
            for ii in 1..candidates.len(){
                if mmax < indexmap[candidates[ii]].1{
                    mmax = indexmap[candidates[ii]].1;
                    maxi = candidates[ii];
                }
            }
            indexmap[maxi]  = (-1,indexmap[maxi].1);
        }
    }
}

//x と y とを並べて使われている残基にフラグを立てる。
//アラインメントのスワップを考えないので一意になるはず。
pub fn get_bool(index_map:&Vec<i64>,ysiz:usize)->Vec<bool>{
    let mut ret:Vec<bool> = vec![false;index_map.len()+ysiz];
    let xsiz = index_map.len();
    for (xi,yi_) in index_map.iter().enumerate(){
        if *yi_ < 0{
            continue;
        }
        ret[xi] = true;
        ret[xsiz+(*yi_ as usize)] = true;
    }
    return ret;
}


pub fn align_seed(
    query:&Vec<Vec<f64>>
    ,template:&Vec<Vec<f64>>
    ,query_aa:&Vec<String>
    ,template_aa:&Vec<String>
    ,num_iter_max:usize
    ,dist_cutoff:f64
    ,tm_d0:f64
    ,lnorm:usize
    ,seed_length:usize
    ,tmscore_break:f64)->(Vec<i64>,(f64,f64)){
    assert_eq!(query.len(),query_aa.len());
    assert_eq!(template.len(),template_aa.len());
    let fragment_length:Vec<usize> = vec![seed_length];
    let q_length:usize =query[0].len(); 
    let t_length:usize =template[0].len(); 
    let check_step:usize = 1;
    let mut max_score:(f64,f64) = (-1.0,-1.0);
    let mut max_map:Vec<i64> = vec![-1;q_length];

    let mut x_to_y_indexmap:Vec<i64> = vec![-1;q_length];
    let aa_filter:HashMap<String,HashSet<String>> = get_blosum62_positive_pair();
    for ff in fragment_length.iter(){
        let mut start:usize=0;
        let mut x_to_y_indexmap_first:Vec<(usize,i64)> = vec![(0,-1);*ff];
        loop{
            println!("{}/{}",start,q_length-ff);
            for jj in 0..(t_length-ff+1){
                for ii in x_to_y_indexmap.iter_mut(){
                    *ii = -1;
                }
                let mut last_q:usize = 0;
                let mut filter_check:bool = true;
                for ii in start..(start+ff){
                    if ii >= q_length{
                        break;
                    }
                    x_to_y_indexmap[ii as usize] = jj as i64 + ii as i64 - start as i64;
                    if !aa_filter.get(&query_aa[ii as usize]).unwrap().contains(&template_aa[(jj as i64 + ii as i64 - start as i64) as usize]){
                        filter_check = false;
                    }
                    x_to_y_indexmap_first[(ii - start) as usize] = (ii as usize,jj as i64 + ii as i64 - start as i64);
                    last_q = (ii - start) as usize;
                }
                if !filter_check{
                    continue;
                }
                //if true{
                // fragment 内だけで Estimate して閾値以上だとスキップするかと思ったけど
                // 多分ノイズが大きすぎて使えないと思う
                //}

                'outiter:for itt in 0..num_iter_max{
                    for ii in x_to_y_indexmap_first.iter(){
                        if ii.1 < 0{
                            continue;
                        }
                        let dist_e = estimate_aligned_atom_distance(query, template, &x_to_y_indexmap, ii.0, ii.1 as usize);
                        
                        if dist_e >= dist_cutoff{
                            //println!("{}",dist_e);
                            //seed の点の距離がカットオフ距離を超えるとこれ以上にシードとして良い場所があるとみなす。
                            break 'outiter;
                        }
                    }

                    let scores = calc_dist_tmscore_gdtts(query,template, &x_to_y_indexmap,dist_cutoff,tm_d0,lnorm);
                    
                    if scores.0 > max_score.0{
                        for ii in 0..q_length{
                            max_map[ii] = x_to_y_indexmap[ii];
                        }
                        max_score = scores;
                    }

                    if itt == num_iter_max-1{
                        break;
                    }

                    if itt > 0 && scores.0 < tmscore_break{
                        break;
                    } 
                    let mut x_to_y_indexmap_t:Vec<(i64,f64)> = vec![(-1,dist_cutoff);q_length];

                    
                    for ii in 0..x_to_y_indexmap_first[0].0{
                        for jj in 0..(x_to_y_indexmap_first[0].1 as usize){
                            let dist_e = estimate_aligned_atom_distance(query, template, &x_to_y_indexmap, ii, jj);
                            if dist_e < x_to_y_indexmap_t[ii].1{
                                x_to_y_indexmap_t[ii] = (jj as i64,dist_e);
                            }
                        }
                    }
                    
                    for ii in (x_to_y_indexmap_first[last_q].0)..q_length{
                        for jj in (x_to_y_indexmap_first[last_q].1 as usize)..t_length{
                            let dist_e = estimate_aligned_atom_distance(query, template, &x_to_y_indexmap, ii, jj);
                            if dist_e < x_to_y_indexmap_t[ii].1{
                                x_to_y_indexmap_t[ii] = (jj as i64,dist_e);
                            }
                        }
                    }


                    filt_discrepancy(&mut x_to_y_indexmap_t);
                    //filt_short_segments(&mut x_to_y_indexmap_t);
                    let mut lflag = false;//マッピングが変化してないと converge したとみなす。
                    for ii in 0..x_to_y_indexmap_t.len(){
                        if x_to_y_indexmap[ii] != x_to_y_indexmap_t[ii].0{
                            lflag = true;
                            x_to_y_indexmap[ii] = x_to_y_indexmap_t[ii].0;
                        }
                    }
                    
                    if !lflag{
                        break;
                    }
                    /*
                    let booll:Vec<bool> = get_bool(&x_to_y_indexmap,t_length);
                    for cc in checked.iter(){
                        if *cc == booll{

                            break 'outiter;
                        }
                    }
                    
                    checked.push(booll);
                    */
                }
            }
            
            if q_length <= start+ff{
                break;
            }else{
                start += check_step;
            }
        }

    }
    return (max_map,max_score);
}

//3.4～12.0 は Contact として、その距離に応じたスコアを与える
//Query が 12.0 未満で Template が 12.0 以上の場合ペナルティを与える
//Template が 12.0 未満で Query が 12.0 以上の場合はスコア 0 で False Negative に対応する
pub fn align_dp(
    query_:&Vec<Vec<f64>>
    ,template_:&Vec<Vec<f64>>
    ,query_aa:&Vec<String>
    ,template_aa:&Vec<String>
    ,dist_cutoff:f64
    ,tm_d0:f64
    ,lnorm:usize)->(Vec<i64>,(f64,f64)){
    assert_eq!(query_.len(),query_aa.len());
    assert_eq!(template_.len(),template_aa.len());



    let q_length:usize =query_[0].len(); 
    let t_length:usize =template_[0].len(); 
    let num_threads:usize = 3;
    let mut lastend:usize = 0;


    
    let (tx, rx):(std::sync::mpsc::Sender<Vec<(usize,usize,f32)>>
                ,std::sync::mpsc::Receiver<Vec<(usize,usize,f32)>>)
     = mpsc::channel();

    let mut start_end:Vec<(usize,usize)> = vec![];
    for tt in 0..num_threads{
        let start:usize = lastend;
        let mut end:usize = (q_length as f64 /(num_threads as f64)*(tt as f64+1.0)) as usize;
        if tt == num_threads-1{
            end = q_length;
        }
        if end > q_length{
            end = q_length;
        }
        lastend = end;
        start_end.push((start,end));
    }
    for (start,end) in start_end{
        
        let tx = tx.clone();
        let query=query_.clone();
        let template=template_.clone();
        let mut estimated_scores:Vec<(usize,usize,f32)> = vec![];
        let _handle = thread::spawn(move || {
            let pm = sequence_alignment::PositionSpecificMatrix::new();
            let mut salign:sequence_alignment::SequenceAlignment = sequence_alignment::SequenceAlignment::new(
            Box::new(pm),-1.0,0.0
            ,sequence_alignment::ALIGN_LOCAL);
            let contact_threshold:f64 = 12.0;
            let mut seq1_dummy:sequence_alignment::SeqData = sequence_alignment::SeqData::new();
            let mut seq2_dummy:sequence_alignment::SeqData = sequence_alignment::SeqData::new();
            seq1_dummy.seq = vec!["X".to_string();q_length];
            seq2_dummy.seq = vec!["X".to_string();t_length];
            salign.prepare(&seq1_dummy,&seq2_dummy);
            for xx_ in start..end{
                
                let xx:usize = end-xx_+start-1;//こうしないとなぜか半分以降が遅い
                //let xx:usize = xx_;
                for yy in 0..t_length{
                    let mut psc:f32 = 0.0;
                    
                    let overlapcount:usize = xx.min(yy)+(t_length-yy-1).min(q_length-xx-1);
                    if overlapcount < 15{
                        //前後の位置関係が合う残基が少ない場合並べない
                        //10 残基 2 秒ほど縮まる

                        psc = -10000000.0;
                    }else{
                        
                        if xx < q_length-1 && yy < t_length-1{
                            for ii in (xx+1)..q_length{
                                for jj in (yy+1)..t_length{
                                    let mut pscore = 0.0;
                                    if query[xx][ii] < contact_threshold && template[yy][jj] > contact_threshold{
                                        pscore = -1.0;
                                    }
                                    if query[xx][ii] < contact_threshold && template[yy][jj] < contact_threshold{
                                        pscore = 1.0;
                                    }
                                    salign.scoring_matrix.set_score(ii-xx-1,jj-yy-1,pscore as f32);
                                }
                            }
                            //ToDo: align_pertial みたいなの作る
                            //→OK
                            //PARTIAL_GLOBAL という一方だけギャップペナルティを入れるタイプを作る？
                            //->使わない。むしろ近くの残基ほど無視した方が良い可能性がある 
                            
                            let res = salign.align_partial(&seq1_dummy,&seq2_dummy,false,Some((q_length-xx-1,t_length-yy-1)),true);
                            //println!("::{} {} {} {} {:?}",start,end,q_length-xx-1,t_length-yy-1,res);
                            psc += res.2;
                        }

                        if xx > 0 && yy > 0{
                                for ii in 0..xx{
                                    for jj in 0..yy{
                                        let mut pscore = 0.0;
                                        if query[xx][ii] < contact_threshold && template[yy][jj] > contact_threshold{
                                            pscore = -1.0;
                                        }
                                        if query[xx][ii] < contact_threshold && template[yy][jj] < contact_threshold{
                                            pscore = 1.0;
                                        }
                                        salign.scoring_matrix.set_score(ii,jj,pscore as f32);
                                    }
                                }
                                           
                                let res = salign.align_partial(&seq1_dummy,&seq2_dummy,false,Some((xx,yy)),true);
                                //println!("{} {} {} {} {:?}",start,end,xx,yy,res);
                                psc += res.2;
                            }
                        }
                        estimated_scores.push((xx,yy,psc));
                    }
                }
                tx.send(estimated_scores)
                
        });
    }

    let mut estimated_scores:Vec<Vec<f32>> = vec![vec![0.0_f32;t_length];q_length];
    
    for _ in 0..num_threads {
        let  escores:Vec<(usize,usize,f32)>= rx.recv().unwrap();
        for ee in escores.into_iter(){
            estimated_scores[ee.0][ee.1] = ee.2;
        }
    }

    let pm = sequence_alignment::PositionSpecificMatrix::new();
    let mut salign:sequence_alignment::SequenceAlignment = sequence_alignment::SequenceAlignment::new(
        Box::new(pm),-1.0,0.0
        ,sequence_alignment::ALIGN_LOCAL);
    let mut seq1_dummy:sequence_alignment::SeqData = sequence_alignment::SeqData::new();
    let mut seq2_dummy:sequence_alignment::SeqData = sequence_alignment::SeqData::new();
    seq1_dummy.seq = query_aa.iter().map(|m|m.to_string()).collect();
    seq2_dummy.seq = template_aa.iter().map(|m|m.to_string()).collect();
    salign.prepare(&seq1_dummy,&seq2_dummy);
    for ii in 0..q_length{
        for jj in 0..t_length{
            salign.scoring_matrix.set_score(ii,jj,estimated_scores[ii][jj]);
        }
    }

    let rres = salign.align(&seq1_dummy,&seq2_dummy,true);
    let mut x_to_y_indexmap = vec![-1;q_length];
    let mut xcou:i64 = -1;
    let mut ycou:i64 = -1;
    for (r1,r2) in rres.0.iter().zip(rres.1.iter()){
        if r1 != "-"{
            xcou += 1;
        }
        if r2 != "-"{
            ycou += 1;
        }
        if r1 != "-" && r2 != "-"{
            x_to_y_indexmap[xcou as usize] = ycou;
        }
    }
    
    let scores = calc_dist_tmscore_gdtts(query_,template_, &x_to_y_indexmap,dist_cutoff,tm_d0,lnorm);
    return (x_to_y_indexmap,scores);
}




pub fn dist_aligntest(){
    
    let pdb_orig = pdbdata::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"6lu7_A.pdb").as_str());
    let pdb = pdbdata::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"2gx4_A.pdb").as_str());
    let orig_aa = pdb_orig.get_aa_sequences();
    let pdb_aa = pdb_orig.get_aa_sequences();
    
    let mut residues_a:Vec<&pdbdata::PDBResidue> = vec![];
    
    for rr in pdb_orig.chains[0].residues.iter(){
        residues_a.push(rr);
    }
    
    
    let mut residues_b:Vec<&pdbdata::PDBResidue> = vec![];
    
    for rr in pdb.chains[0].residues.iter(){
        residues_b.push(rr);
    }

    let mut dist_a:Vec<Vec<f64>> = vec![vec![0.0;residues_a.len()];residues_a.len()];
    let mut dist_b:Vec<Vec<f64>> = vec![vec![0.0;residues_b.len()];residues_b.len()];

    for ii in 0..residues_a.len(){
        for jj in 0..residues_a.len(){
            dist_a[ii][jj] = residues_a[ii].get_CA().unwrap().distance(residues_a[jj].get_CA().unwrap());
        }    
    }
    
    for ii in 0..residues_b.len(){
        for jj in 0..residues_b.len(){
            dist_b[ii][jj] = residues_b[ii].get_CA().unwrap().distance(residues_b[jj].get_CA().unwrap());
        }    
    }
    //let res_:Option<DistanceAlignmentResult> = align_distance(&dist_a,&dist_b,&orig_aa[0].1,&pdb_aa[0].1,3,0.3,8.0,DistanceAlignmentType::SEED);
    //let res = res_.unwrap();
    //println!("{:?}",res);
    let res_:Option<DistanceAlignmentResult> = align_distance(&dist_a,&dist_b,&orig_aa[0].1,&pdb_aa[0].1,3,0.3,8.0,DistanceAlignmentType::DP);
    let res = res_.unwrap();
    println!("{:?}",res);

}

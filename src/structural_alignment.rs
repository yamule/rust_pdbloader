
use super::matrix_process;
use super::pdbdata;
#[allow(unused_imports)]
use super::mmcif_process;
#[allow(dead_code,unused_imports)]
use super::process_3d;
use super::geometry::Vector3D;
use super::sequence_alignment;


#[allow(dead_code,unused_imports)]
use rand::SeedableRng;
#[allow(dead_code,unused_imports)]
use rand::rngs::StdRng;
#[allow(dead_code,unused_imports)]
use rand::Rng;

#[allow(dead_code,unused_imports)]
use super::debug_env;

use std::collections::HashMap;
use std::sync::Mutex;

#[allow(dead_code,unused_imports)]
const ERROR_THRESHOLD:f64 = 0.000000001;
//Kabsh のアルゴリズムとか SVD とかが良く分かんなかったので QR DECOMP で ROTATION MATRIX を出している。


#[allow(non_camel_case_types)]
#[derive(Clone)]
pub enum AlignmentType{
    RAW_FULLRMSD,
    RAW,
    SW,
    MAXIMUM
}

lazy_static! {
    static ref AA_3_TO_1:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
}

pub fn prepare_static(){
    if AA_3_TO_1.lock().unwrap().len() != 0{
        return;
    }
    AA_3_TO_1.lock().unwrap().insert("ALA".to_string(),"A".to_string());
    AA_3_TO_1.lock().unwrap().insert("ARG".to_string(),"R".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASN".to_string(),"N".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASP".to_string(),"D".to_string());
    AA_3_TO_1.lock().unwrap().insert("CYS".to_string(),"C".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLN".to_string(),"Q".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLU".to_string(),"E".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLY".to_string(),"G".to_string());
    
    AA_3_TO_1.lock().unwrap().insert("HIS".to_string(),"H".to_string());
    AA_3_TO_1.lock().unwrap().insert("HSD".to_string(),"H".to_string());
    AA_3_TO_1.lock().unwrap().insert("HSC".to_string(),"H".to_string());

    AA_3_TO_1.lock().unwrap().insert("ILE".to_string(),"I".to_string());
    AA_3_TO_1.lock().unwrap().insert("LEU".to_string(),"L".to_string());
    AA_3_TO_1.lock().unwrap().insert("LYS".to_string(),"K".to_string());
    AA_3_TO_1.lock().unwrap().insert("MET".to_string(),"M".to_string());
    AA_3_TO_1.lock().unwrap().insert("PHE".to_string(),"F".to_string());
    AA_3_TO_1.lock().unwrap().insert("PRO".to_string(),"P".to_string());
    AA_3_TO_1.lock().unwrap().insert("SER".to_string(),"S".to_string());
    AA_3_TO_1.lock().unwrap().insert("THR".to_string(),"T".to_string());
    AA_3_TO_1.lock().unwrap().insert("TRP".to_string(),"W".to_string());
    AA_3_TO_1.lock().unwrap().insert("TYR".to_string(),"Y".to_string());
    AA_3_TO_1.lock().unwrap().insert("VAL".to_string(),"V".to_string());
}

//二つのアラインされた構造の座標を与えて、q->t のインデクスのマップとその距離を返す。
//一次元が座標である場合と二次元が座標である場合があるので注意
pub fn get_query_tmp_indexmap(qpos:&Vec<Vec<f64>>,tpos:&Vec<Vec<f64>>,distance_cutoff:f64)
-> Vec<(i64,f64)>{
    let distance_cutoff2:f64 = distance_cutoff*distance_cutoff;//sqrt の計算コストを節約するため
    let q_length:usize = qpos.len();
    let t_length:usize = tpos.len();
    let mut x_to_y_indexmap_t:Vec<(i64,f64)> = vec![(-1,distance_cutoff2);q_length];

    for ii in 0..q_length{
        for jj in 0..t_length{
            let ddis2 = distance2(
                qpos[ii][0],qpos[ii][1],qpos[ii][2],
                tpos[jj][0],tpos[jj][1],tpos[jj][2]
            );
            if ddis2 < x_to_y_indexmap_t[ii].1{
                x_to_y_indexmap_t[ii] = (jj as i64,ddis2);
            }
        }
    }
    for tt in x_to_y_indexmap_t.iter_mut(){
        tt.1 = tt.1.sqrt();
    }
    return x_to_y_indexmap_t;
}

//sqrt する前の distance を返す。
pub fn distance2(x1:f64,y1:f64,z1:f64,x2:f64,y2:f64,z2:f64)->f64{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}


pub fn get_rotation_matrix(a:&Vec<Vec<f64>>,b:&Vec<Vec<f64>>,pnum:usize)->Option<Vec<Vec<f64>>>{
    let mut pmat:Vec<Vec<f64>> = matrix_process::matrix_multi_slice(&b,&matrix_process::matrix_t_slice(&a,3,pnum),3,3,pnum);
    let aat_:Option<Vec<Vec<f64>>> = matrix_process::matrix_inv(&matrix_process::matrix_multi_slice(&a,&matrix_process::matrix_t_slice(&a,3,pnum),3,3,pnum));
    if let None = aat_{
        return None;
    }
    let aat = aat_.unwrap();
    pmat = matrix_process::matrix_multi(&pmat,&aat);

    let qr:(Vec<Vec<f64>>,Vec<Vec<f64>>) = matrix_process::house_qr_decomp(&pmat);
    return Some(qr.0);
}


pub fn calc_tmscore_gdtts(target:&Vec<Vec<f64>>,template:&Vec<Vec<f64>>,x_to_y_indexmap:&Vec<i64>,tm_cutoff:f64,tm_d0:f64,tm_lnorm:usize)->(f64,f64){
    let mut gdtts:Vec<usize> = vec![0;4];
    let mut tm_sum:f64 = 0.0;
    let d02:f64 = tm_d0*tm_d0;
    let tm_cutoff2:f64 = tm_cutoff*tm_cutoff;
    for (xi,yi_) in x_to_y_indexmap.iter().enumerate(){
        if *yi_ < 0{
            continue;
        }
        let yi = *yi_ as usize;
        let dist2 = distance2(
            target[0][xi],target[1][xi],target[2][xi],
            template[0][yi],template[1][yi],template[2][yi]
        );
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

        if dist2 < 64.0{
            gdtts[3] += 1;
            if dist2 < 16.0{
                gdtts[2] += 1;
                if dist2 < 4.0{
                    gdtts[1] += 1;
                    if dist2 < 1.0{
                        gdtts[0] += 1;
                    }
                }
            }
        }
        
        if dist2 <= tm_cutoff2{
            tm_sum += 1.0/(1.0+dist2/d02);
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


pub fn align_pdb(chain1:&Vec<&pdbdata::PDBComp>,chain2:&Vec<&pdbdata::PDBComp>,alitype:AlignmentType,tmscore_break:f64)->Option<StructuralAlignmentResult>{
    prepare_static();
    let mut xvec:Vec<Vec<f64>> = vec![];
    let mut yvec:Vec<Vec<f64>> = vec![];
    let mut unrelated:bool = true;
    match alitype{
        AlignmentType::SW| AlignmentType::RAW| AlignmentType::RAW_FULLRMSD=>{
            unrelated = false;
            let mut xseq:Vec<(String,(f64,f64,f64))> = vec![];//回転した後に近い残基を取るだけなので配列上の位置情報は捨ててしまう
            let mut yseq:Vec<(String,(f64,f64,f64))> = vec![];  
            for rr in chain1.iter(){
                let mut caatom:Option<&pdbdata::PDBAtom> = None;
                for aa in rr.iter_atoms(){
                    if aa.atom_code == "CA"{
                        caatom = Some(aa);
                    }
                }
                if let Some(aa) = caatom{
                    xseq.push((rr.comp_id.clone(),(aa.get_x(),aa.get_y(),aa.get_z())));
                }
            }
            
            for rr in chain2.iter(){
                let mut caatom:Option<&pdbdata::PDBAtom> = None;
                for aa in rr.iter_atoms(){
                    if aa.atom_code == "CA"{
                        caatom = Some(aa);
                    }
                }
                if let Some(aa) = caatom{
                    yseq.push((rr.comp_id.clone(),(aa.get_x(),aa.get_y(),aa.get_z())));
                }
            }
            match alitype{
                AlignmentType::RAW| AlignmentType::RAW_FULLRMSD=>{
                    if yseq.len() != xseq.len(){
                        panic!("Sequence lengths must be equal with this mode!");
                    }
                    for ii in 0..xseq.len(){
                        xvec.push(vec![(xseq[ii].1).0,(xseq[ii].1).1,(xseq[ii].1).2]);
                        yvec.push(vec![(yseq[ii].1).0,(yseq[ii].1).1,(yseq[ii].1).2]);
                    }
                },
                AlignmentType::SW =>{
                   

                    let mut ss1 = sequence_alignment::SeqData::new();
                    let mut ss2 = sequence_alignment::SeqData::new();
                    let mut sw = sequence_alignment::SequenceAlignment::new(Box::new(sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,sequence_alignment::ALIGN_LOCAL);

                    for ss in xseq.iter(){
                        ss1.seq.push(AA_3_TO_1.lock().unwrap().get(&ss.0).unwrap_or(&("X".to_string())).clone());
                    }
                    for ss in yseq.iter(){
                        ss2.seq.push(AA_3_TO_1.lock().unwrap().get(&ss.0).unwrap_or(&("X".to_string())).clone());
                    }
                    
                    let res = sw.align(&ss1,&ss2,true);
                    let alen:usize = res.0.len();
                    let mut xpos:usize = 0;
                    let mut ypos:usize = 0;
                    let mut mcou:usize = 0;
                    for ii in 0..alen{
                        if res.0[ii] != "-" && res.1[ii] != "-"{
                            xvec.push(vec![(xseq[xpos].1).0,(xseq[xpos].1).1,(xseq[xpos].1).2]);
                            yvec.push(vec![(yseq[ypos].1).0,(yseq[ypos].1).1,(yseq[ypos].1).2]);
                            mcou += 1;
                        }
                        if res.0[ii] != "-"{
                            xpos += 1;        
                        }
                        if res.1[ii] != "-"{
                            ypos += 1;
                        }
                    } 
                    if mcou < 10{
                        panic!("Regions aligned by sequence alignment was too short! {}\n Try -type max!",mcou);
                    }
                },
                _=>{}
            }
            
            

            //sw 結果
            //let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
            //let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
            //parameter_set4final による

        },
        AlignmentType::MAXIMUM =>{
            for rr in chain1.iter(){
                let mut caatom:Option<&pdbdata::PDBAtom> = None;
                for aa in rr.iter_atoms(){
                    if aa.atom_code == "CA"{
                        caatom = Some(aa);
                    }
                }
                if let Some(aa) = caatom{
                    xvec.push(vec![aa.get_x(),aa.get_y(),aa.get_z()]);
                }
            }
            
            for rr in chain2.iter(){
                let mut caatom:Option<&pdbdata::PDBAtom> = None;
                for aa in rr.iter_atoms(){
                    if aa.atom_code == "CA"{
                        caatom = Some(aa);
                    }
                }
                
                if let Some(aa) = caatom{
                    yvec.push(vec![aa.get_x(),aa.get_y(),aa.get_z()]);
                }
            }
        }
    }

    let lnorm:usize = xvec.len().min(yvec.len());
    let d0:f64 = if lnorm <=19{
        0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    
    let xvec = matrix_process::matrix_t(&xvec);
    let yvec = matrix_process::matrix_t(&yvec);
    let mat:Option<(Vec<Vec<f64>>,(f64,f64))> = if !unrelated{
        align(&xvec,&yvec,100,d0_search,d0,lnorm)
    }else{//ToDo パラメータわたすようにする
        align_unrelated(&xvec,&yvec,100,d0_search,d0,lnorm,4,3.0,tmscore_break,0.2)
    };

    
        /*
         Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
        if (Lnorm<=19)                    //update 15-->19
            d0=0.168;                   //update 0.5-->0.168
        else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
        D0_MIN=d0+0.8;              //this should be moved to above
        d0=D0_MIN;                  //update: best for search    

        */


    
    if let None = mat{
        panic!("Structures couldn't be aligned!");
    };
    let rotmat:Vec<Vec<f64>> = mat.unwrap().0;
    let mut xpos_all:Vec<Vec<f64>> = vec![];
    let mut xpos_index:Vec<usize> = vec![];
    let mut ypos_all:Vec<Vec<f64>> = vec![];  
    let mut ypos_index:Vec<usize> = vec![];
    for (rii,rr) in chain1.iter().enumerate(){
        let mut caatom:Option<&pdbdata::PDBAtom> = None;
        for aa in rr.iter_atoms(){
            if aa.atom_code == "CA"{
                caatom = Some(aa);
            }
        }
        if let Some(aa) = caatom{
            xpos_all.push(vec![aa.get_x(),aa.get_y(),aa.get_z(),1.0]);
            xpos_index.push(rii);
        }else{
            eprintln!("Chain1: {} {} {} does not have CA atom!",rr.get_name(),rr.get_seq_id(),rr.get_ins_code());
        }
    }

    for (rii,rr) in chain2.iter().enumerate(){
        let mut caatom:Option<&pdbdata::PDBAtom> = None;
        for aa in rr.iter_atoms(){
            if aa.atom_code == "CA"{
                caatom = Some(aa);
            }
        }
        if let Some(aa) = caatom{
            ypos_all.push(vec![aa.get_x(),aa.get_y(),aa.get_z()]);
            ypos_index.push(rii);
        }else{
            eprintln!("Chain2: {} {} {} does not have CA atom!",rr.get_name(),rr.get_seq_id(),rr.get_ins_code());
        }
    }
    //最適な回転行列が計算できたので
    //最初の SW アラインメントとかで並ばなかった残基についてもアラインできる可能性を考慮し
    //全体を回転させて近くに CA が来るかどうか見る。
    //最初の SW で並んだ残基が並ばなくなる可能性もある。
    let dist_cutoff:f64 = 8.0;
    let x_length = xpos_all.len();
    let y_length = ypos_all.len();
    let mut x_to_y_indexmap_t:Vec<(i64,f64)> = vec![(-1,dist_cutoff*dist_cutoff);x_length];
    let xvec_aligned = matrix_process::matrix_multi(&rotmat,&matrix_process::matrix_t(&xpos_all));//chain1 を chain2 に合うように回転と移動
    let ypos_all = matrix_process::matrix_t(&ypos_all);
    let mut rmsd:f64 = 0.0;
    let mut rcount:usize = 0;
    
    match alitype{
        AlignmentType::RAW_FULLRMSD=>{
            for ii in 0..x_length{
                let ddis2 = distance2(
                xvec_aligned[0][ii],xvec_aligned[1][ii],xvec_aligned[2][ii],
                ypos_all[0][ii],ypos_all[1][ii],ypos_all[2][ii]
                );
                x_to_y_indexmap_t[ii] = (ii as i64,ddis2);
                rcount += 1;
                rmsd += ddis2;
                
            }
        },
        _=>{
            for ii in 0..x_length{
                for jj in 0..y_length{
                    let ddis2 = distance2(
                    xvec_aligned[0][ii],xvec_aligned[1][ii],xvec_aligned[2][ii],
                    ypos_all[0][jj],ypos_all[1][jj],ypos_all[2][jj]
                    );
                    if ddis2 < x_to_y_indexmap_t[ii].1{
                        x_to_y_indexmap_t[ii] = (jj as i64,ddis2);
                    }
                }
            }
            
            filt_discrepancy(&mut x_to_y_indexmap_t);
            for ii in 0..x_length{
                if x_to_y_indexmap_t[ii].0 < 0{
                    continue;
                }
                rcount += 1;
                if x_to_y_indexmap_t[ii].1 > 0.0{
                    rmsd += x_to_y_indexmap_t[ii].1;
                }
            }
            if rcount == 0{
                panic!("Structures could not be aligned.");
            }
        }
    }
    rmsd = (rmsd/(rcount as f64)).sqrt();

    let mut ret = StructuralAlignmentResult::new();
    let indexmap:Vec<i64> = x_to_y_indexmap_t.iter().map(|m|m.0).collect();

    
    //short の方のスコア計算
    let lnorm:usize = x_length;
    let d0:f64 = if lnorm <=19{
        0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    let scores_chain1 = calc_tmscore_gdtts(&xvec_aligned,&ypos_all, &indexmap,d0_search,d0,lnorm);
    

    //long の方のスコア計算
    let lnorm:usize = y_length;
    let d0:f64 = if lnorm <=19{0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    let scores_chain2 = calc_tmscore_gdtts(&xvec_aligned,&ypos_all, &indexmap,d0_search,d0,lnorm);

    ret.tmscore_chain1 = scores_chain1.0;
    ret.tmscore_chain2 = scores_chain2.0;

    assert_eq!(scores_chain1.1,scores_chain2.1);
    ret.gdt_ts = scores_chain1.1;
    ret.transform_matrix = rotmat;
    ret.dist_cutoff = dist_cutoff;
    ret.rmsd_aligned = rmsd;
    let mut xstr:Vec<String> = vec![];
    let mut ystr:Vec<String> = vec![];
    let mut lastcount_x:i64 = -1;
    let mut lastcount_y:i64 = -1;
    let mut num_aligned:usize = 0;
    for ii in 0..x_to_y_indexmap_t.len(){
        if x_to_y_indexmap_t[ii].0 < 0{
            continue;
        }
        while lastcount_x < ii as i64 -1{
            lastcount_x += 1;
            xstr.push(AA_3_TO_1.lock().unwrap().get(chain1[xpos_index[lastcount_x as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
            ystr.push("-".to_string());
        }
        while lastcount_y < x_to_y_indexmap_t[ii].0 as i64 -1{
            lastcount_y += 1;
            ystr.push(AA_3_TO_1.lock().unwrap().get(chain2[ypos_index[lastcount_y as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
            xstr.push("-".to_string());
        }
        
        lastcount_y += 1;
        lastcount_x += 1;
        num_aligned += 1;
        xstr.push(AA_3_TO_1.lock().unwrap().get(chain1[xpos_index[lastcount_x as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
        ystr.push(AA_3_TO_1.lock().unwrap().get(chain2[ypos_index[lastcount_y as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
    }
    
    while lastcount_x < xpos_index.len() as i64 -1{
        lastcount_x += 1;
        xstr.push(AA_3_TO_1.lock().unwrap().get(chain1[xpos_index[lastcount_x as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
        ystr.push("-".to_string());
    }
    while lastcount_y < ypos_index.len() as i64 -1{
        lastcount_y += 1;
        ystr.push(AA_3_TO_1.lock().unwrap().get(chain2[ypos_index[lastcount_y as usize]].get_comp_id()).unwrap_or(&("X".to_owned())).clone());
        xstr.push("-".to_string());
    }
    ret.aligned_chain1 = xstr;
    ret.aligned_chain2 = ystr;
    if let AlignmentType::RAW_FULLRMSD =  alitype{
        ret.num_aligned = -1;
    }else{
        ret.num_aligned = num_aligned as i64;
    }
    return Some(ret);
}

#[derive(Debug)]
pub struct StructuralAlignmentResult{
    pub tmscore_chain1:f64,
    pub tmscore_chain2:f64,
    pub gdt_ts:f64,
    pub transform_matrix:Vec<Vec<f64>>,
    
    pub dist_cutoff:f64,
    pub rmsd_aligned:f64,
    pub num_aligned:i64,
    pub aligned_chain1:Vec<String>,
    pub aligned_chain2:Vec<String>,

}


impl StructuralAlignmentResult{
    pub fn new()->StructuralAlignmentResult{
        return StructuralAlignmentResult{
            tmscore_chain1:0.0,
            tmscore_chain2:0.0,
            gdt_ts:0.0,
            transform_matrix:matrix_process::eye(4),

            dist_cutoff:0.0,
            rmsd_aligned:0.0,
            num_aligned:0,
            aligned_chain1:vec![],
            aligned_chain2:vec![]
        };
    }
}


//既にある程度どの原子とどの原子が並ぶか分かっている場合に使用するモード
//対応する原子が query と template として与えられる Vec の同じインデクス上にあるべき
//MISSING な場合相手の原子も削除しておく
//(transformation matrix ,(tmscore,gdtts))
//を返す。
//transformation matrix は rotation + 移動
//x
//y
//z
//1
//に左からかけると y になるように。
pub fn align(query:&Vec<Vec<f64>>,template:&Vec<Vec<f64>>,num_iter_max:usize,dist_cutoff:f64,tm_d0:f64,lnorm:usize)->Option<(Vec<Vec<f64>>,(f64,f64))>{
    assert_eq!(query[0].len(),template[0].len());
    let dist_cutoff2:f64 = dist_cutoff*dist_cutoff;
    //let fragment_length:Vec<usize> = vec![3,5,10,15,20];
    let fragment_length:Vec<usize> = vec![3,5];
    let q_length:usize =query[0].len(); 
    let t_length:usize =template[0].len(); 
    let indexmap_dummy:Vec<i64> = (0..q_length).into_iter().map(|m|m as i64).collect();

    let check_step:usize = 1;
    let mut max_score:(f64,f64) = (-1.0,-1.0);

    //rotationmatrix, xcenter, ycenter
    let mut max_matrix:(Vec<Vec<f64>>,(f64,f64,f64),(f64,f64,f64)) = (vec![],(0.0,0.0,0.0),(0.0,0.0,0.0));
    let mut x_to_y_indexmap:Vec<i64> = vec![-1;q_length];
    let mut qmatbuff:Vec<Vec<f64>> = query.clone();
    let mut tmatbuff:Vec<Vec<f64>> = template.clone();
    let calc_center = |xx:&Vec<Vec<f64>>,anum:usize|->(f64,f64,f64){
        let mut ret:(f64,f64,f64) = (0.0,0.0,0.0);
        for i in 0..anum{
            ret.0 += xx[0][i];
            ret.1 += xx[1][i];
            ret.2 += xx[2][i];
        }
        ret.0 /= anum as f64;
        ret.1 /= anum as f64;
        ret.2 /= anum as f64;
        return ret
    };

    for ff in fragment_length.iter(){
        let mut start:usize=0;
        loop{
            for ii in x_to_y_indexmap.iter_mut(){
                *ii = -1;
            }
            for ii in start..(start+ff){
                if ii >= q_length{
                    break;
                }
                x_to_y_indexmap[ii as usize] = ii as i64;
            }

            for itt in 0..num_iter_max{
                let mut num_align:usize = 0;
                for i in x_to_y_indexmap.iter(){
                    if *i > -1{
                        let iu = *i as usize;
                        qmatbuff[0][num_align] = query[0][iu];
                        qmatbuff[1][num_align] = query[1][iu];
                        qmatbuff[2][num_align] = query[2][iu];

                        tmatbuff[0][num_align] = template[0][iu];
                        tmatbuff[1][num_align] = template[1][iu];
                        tmatbuff[2][num_align] = template[2][iu];

                        num_align += 1;
                    }
                }
                let xcenter:(f64,f64,f64) = calc_center(&qmatbuff,num_align);
                let ycenter:(f64,f64,f64) = calc_center(&tmatbuff,num_align);
                for ii in 0..num_align{
                    qmatbuff[0][ii] -= xcenter.0;
                    qmatbuff[1][ii] -= xcenter.1;
                    qmatbuff[2][ii] -= xcenter.2;
                    
                    tmatbuff[0][ii] -= ycenter.0;
                    tmatbuff[1][ii] -= ycenter.1;
                    tmatbuff[2][ii] -= ycenter.2;
                }

                let rotmat_ = get_rotation_matrix(&qmatbuff,&tmatbuff,num_align);
                if let None = rotmat_{
                }else{
                    let rotmat = rotmat_.unwrap();
                    for ii in 0..q_length{
                        let mres = matrix_process::matrix_multi(&rotmat
                            ,&vec![
                             vec![query[0][ii]-xcenter.0]
                            ,vec![query[1][ii]-xcenter.1]
                            ,vec![query[2][ii]-xcenter.2]
                            ]);
                        qmatbuff[0][ii] = mres[0][0];
                        qmatbuff[1][ii] = mres[1][0];
                        qmatbuff[2][ii] = mres[2][0];
                    }

                    for ii in 0..t_length{
                        tmatbuff[0][ii] = template[0][ii]-ycenter.0;
                        tmatbuff[1][ii] = template[1][ii]-ycenter.1;
                        tmatbuff[2][ii] = template[2][ii]-ycenter.2;
                    }

                    let scores = calc_tmscore_gdtts(&qmatbuff,&tmatbuff, &indexmap_dummy,dist_cutoff,tm_d0,lnorm);
                    if scores.0 > max_score.0{
                        max_matrix = (rotmat,xcenter,ycenter);
                        max_score = scores;
                    }
                    
                    if itt == num_iter_max-1{
                        break;
                    }

                    let mut x_to_y_indexmap_t:Vec<i64> = vec![-1;q_length];

                    for ii in 0..q_length{
                        let ddis2 = distance2(
                        qmatbuff[0][ii],qmatbuff[1][ii],qmatbuff[2][ii],
                        tmatbuff[0][ii],tmatbuff[1][ii],tmatbuff[2][ii]
                        );
                        
                        if ddis2 < dist_cutoff2{
                            x_to_y_indexmap_t[ii] = ii as i64;
                        }
                    }
                    let mut lflag = false;//マッピングが変化してないと converge したとみなす。
                    for ii in 0..x_to_y_indexmap_t.len(){
                        if x_to_y_indexmap[ii] != x_to_y_indexmap_t[ii]{
                            lflag = true;
                            x_to_y_indexmap[ii] = x_to_y_indexmap_t[ii];
                        }
                    }

                    if !lflag{
                        break;
                    }
                }

            }
            if q_length <= start+ff{
                break;
            }else{
                start += check_step;
            }
        }

    }


    let mut movx:Vec<Vec<f64>> = matrix_process::eye(4);
    movx[0][3] = -(max_matrix.1).0;
    movx[1][3] = -(max_matrix.1).1;
    movx[2][3] = -(max_matrix.1).2;


    let mut movy:Vec<Vec<f64>> = matrix_process::eye(4);
    movy[0][3] = (max_matrix.2).0;
    movy[1][3] = (max_matrix.2).1;
    movy[2][3] = (max_matrix.2).2;

    let mut m_ret_:Vec<Vec<f64>> = matrix_process::eye(4);
    for ii in 0..3{
        for jj in 0..3{
            m_ret_[ii][jj] = max_matrix.0[ii][jj];
        }
    }

    let m_ret:Vec<Vec<f64>> = matrix_process::matrix_multi(&movy,&matrix_process::matrix_multi(&m_ret_,&movx));

    return Some((m_ret,max_score));

}


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
pub fn filt_discrepancy(indexmap:&mut Vec<(i64,f64)>){
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


//dist cutoff 内にある原子をグルーピングしていく
//クラスターの大きい順に返す
pub fn group_atoms(pos:&Vec<Vec<f64>>,dist_cutoff:f64)->Vec<Vec<usize>>{
    let plen:usize = pos.len();
    let mut indexmaplist:Vec<usize> = (0..plen).into_iter().collect();

    let get_min_index = |target:usize,indexmaplist_:&Vec<usize>|->usize{
        let mut ret:usize = target;
        loop{
            if indexmaplist_[ret] == ret{
                break;
            }
            assert!(indexmaplist_[ret] < ret);
            ret = indexmaplist_[ret];
        }
        return ret;
    };

    for ii in 0..plen{
        for jj in ii+1..plen{
            let ddis =  process_3d::distance(&(pos[ii][0],pos[ii][1],pos[ii][2]), &(pos[jj][0],pos[jj][1],pos[jj][2]));
            if ddis < dist_cutoff{
                let p = get_min_index(ii, &indexmaplist);
                let q = get_min_index(jj, &indexmaplist);
                let pq = p.min(q).min(ii);
                indexmaplist[jj] = pq;
                indexmaplist[ii] = pq;
                indexmaplist[p] = pq;
                indexmaplist[q] = pq;
            }
        }
    }
    let mut groups:HashMap<usize,Vec<usize>> = HashMap::new();
    for ii in 0..plen{
        let mmin = get_min_index(ii, &indexmaplist);
        if !groups.contains_key(&mmin){
            groups.insert(mmin,vec![]);
        }
        groups.get_mut(&mmin).unwrap().push(ii);
    }
    let mut vvec:Vec<Vec<usize>> = groups.into_iter().map(|m|m.1).collect();
    vvec.sort_by(|a,b|a.len().cmp(&b.len()));
    vvec.reverse();
    return vvec;
}


//二つ構造を並べて構造の違いからドメインを見積もり、インデクスの Vec として返す。
pub fn comparative_domain_split(
    query:&Vec<Vec<f64>>
    ,template:&Vec<Vec<f64>>
    ,min_length:usize
    ,max_dist_diff:f64
    ,num_iter_max:usize
    ,seed_length:usize
    ,seed_dist_limit:f64
    ,tmscore_break:f64//数回 ITER かましてこれより TMSCORE が小さい場合 Break
    ,stop_tmscore_ratio:f64//数回 ITER かましてこれより現在TMSCORE/最大 TMSCORE が小さい場合 Break
)->Vec<Vec<usize>>{
    
    
    let mut query_remained:Vec<(usize,Vec<f64>)> = query.iter().enumerate().map(|m|{(m.0,m.1.clone())}).collect();
    let mut template_remained:Vec<(usize,Vec<f64>)> = template.iter().enumerate().map(|m|{(m.0,m.1.clone())}).collect();
    
    let mut query_remained_flag:Vec<bool> = vec![true;query.len()];
    let mut template_remained_flag:Vec<bool> = vec![true;template.len()];
    let mut ret:Vec<Vec<usize>> = vec![];//ドメインに含まれる座標の INDEX
    loop{
        if query_remained.len() < min_length || template_remained.len() < min_length{
            break;
        }
        let query_tmp:Vec<Vec<f64>> = query_remained.iter().map(|m|m.1.clone()).collect();
        let template_tmp:Vec<Vec<f64>> = template_remained.iter().map(|m|m.1.clone()).collect();
        
        let lnorm:usize = template_tmp.len().min(query_tmp.len());
        let tm_d0:f64 = if lnorm <=19{
            0.5
        }else{
            1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
        };
        let dist_cutoff:f64 = tm_d0.max(4.5).min(8.0);
        

        let pres = align_unrelated(
        &matrix_process::matrix_t(&query_tmp)
        ,&matrix_process::matrix_t(&template_tmp)
        ,num_iter_max
        ,dist_cutoff
        ,tm_d0
        ,lnorm
        ,seed_length
        ,seed_dist_limit
        ,tmscore_break
        ,stop_tmscore_ratio);
        println!("{:?}",pres);
        if let None = pres{
            break;
        }
        let mat:Vec<Vec<f64>> = pres.unwrap().0;
        let mut query_rotated:Vec<Vec<f64>> = vec![];
        for qq in query_tmp.iter(){
            let pos:Vec<Vec<f64>> = matrix_process::matrix_multi(&mat,&vec![vec![qq[0]],vec![qq[1]],vec![qq[2]],vec![1.0]]);
            query_rotated.push(vec![pos[0][0],pos[1][0],pos[2][0]]);
        }

        let qt_index_map_z:Vec<(i64,f64)> = get_query_tmp_indexmap(&query_rotated,&template_tmp,max_dist_diff);

        //一番大きなグループ以外はマスクする 1 にはグループサイズ
        let group_filtered = |qt_index_map_:&Vec<(i64,f64)>,gpos:&Vec<Vec<f64>>| -> (Vec<(i64,f64)>,usize){
            let mut qs:Vec<(usize,Vec<f64>)> = vec![];
            for (qii,qq) in qt_index_map_.iter().enumerate(){
                if qq.0 > -1{
                    qs.push((qii,gpos[qii].clone()));
                }
            }
            let ppos:Vec<Vec<f64>> = qs.iter().map(|m|m.1.clone()).collect();//メモリ再確保は冗長だが気にしない
            let groups = group_atoms(&ppos,11.0);//cutoff は適当。
            let mut ret:Vec<(i64,f64)> = vec![(-1,max_dist_diff);qt_index_map_.len()];
            for gg_ in groups[0].iter(){
                let gg = qs[*gg_].0;
                ret[gg] = qt_index_map_[gg];
            }
            
            filt_discrepancy(&mut ret);
            let mut ccou:usize = 0;
            for rr in ret.iter(){
                if rr.0 > -1{
                    ccou += 1;
                }
            }
            return (ret,ccou);
        };
        let (mut qt_index_map,mut max_group_size):(Vec<(i64,f64)>,usize) =group_filtered(&qt_index_map_z,&query_rotated);
        if max_group_size < min_length{
            break;
        }
        loop{
            let mut query_tmp2:Vec<Vec<f64>> = vec![];
            let mut template_tmp2:Vec<Vec<f64>> = vec![];
            for (tii,tt) in qt_index_map.iter().enumerate(){
                if tt.0 > -1{
                    query_tmp2.push(query_tmp[tii].clone());
                    template_tmp2.push(template_tmp[tt.0 as usize].clone());
                }
            }
            if query_tmp2.len() < min_length || template_tmp2.len() < min_length{
                break;
            }
            //matrix を得るだけなので
            let pres = align(&matrix_process::matrix_t(&query_tmp2),&matrix_process::matrix_t(&template_tmp2),num_iter_max,dist_cutoff,tm_d0,lnorm);
            if let None = pres{
                break;
            }

            let mut query_rotated:Vec<Vec<f64>> = vec![];
            //align は全原子でやる
            for qq in query_tmp.iter(){
                let pos:Vec<Vec<f64>> = matrix_process::matrix_multi(&mat,&vec![vec![qq[0]],vec![qq[1]],vec![qq[2]],vec![1.0]]);
                query_rotated.push(vec![pos[0][0],pos[1][0],pos[2][0]]);
            }
            let qt_index_map_t_:Vec<(i64,f64)> = get_query_tmp_indexmap(&query_rotated,&template_tmp,max_dist_diff);
            let (qt_index_map_t,group_size_t):(Vec<(i64,f64)>,usize) = group_filtered(&qt_index_map_t_,&query_rotated);
            if group_size_t < max_group_size{
                break;
            }
            let mut updated:bool = false;
            let qlen:usize  = qt_index_map.len();
            for ii in 0..qlen{
                if qt_index_map[ii].0 < 0{
                    continue;
                }else{
                    if qt_index_map[ii].0 != qt_index_map_t[ii].0{
                        updated = true;
                        qt_index_map[ii] = qt_index_map_t[ii]
                    }
                }
            }
            
            //無限ループを避けるために前回と同じサイズである場合停止する
            if group_size_t == max_group_size{
                break;
            }
            max_group_size = group_size_t;
            if !updated{
                break;
            }
        }
        let mut dom:Vec<usize> = vec![];
        for (qii,qq) in qt_index_map.iter().enumerate(){
            if qq.0 > -1{
                query_remained_flag[query_remained[qii].0] = false;
                dom.push(query_remained[qii].0);
                template_remained_flag[template_remained[qq.0 as usize].0] = false;
            }
        }
        if dom.len() > min_length{
            ret.push(dom);
        }
        query_remained.clear();
        for (qii,qq) in query.iter().enumerate(){
            if query_remained_flag[qii]{
                query_remained.push((qii,qq.clone()));
            }
        }
        template_remained.clear();
        for (tii,tt) in template.iter().enumerate(){
            if template_remained_flag[tii]{
                template_remained.push((tii,tt.clone()));
            }
        }
    }
    return ret;
}




// まったくアラインの仕方に手掛かりがない構造同士を並べる
//transformation matrix は rotation + 移動
//x
//y
//z
//1
//に左からかけると y になるように。
//seed_length: この長さの配列が完全にアラインできると考える
//seed_dist_limit: seed の先頭と最後の点の距離の差の絶対値がこれを超えるとアラインできないと考える
//tmscore_break: これより低いスコアが出るとアイテレーションしても改善しないと考えてアイテレーションループを終了する。ただし数回はアイテレーションする
//stop_tmscore_ratio: 現在の maximum の tmscore*この ratio 未満の値が出た時アイテレーションループを終了する
pub fn align_unrelated(query:&Vec<Vec<f64>>,template:&Vec<Vec<f64>>,num_iter_max:usize
    ,dist_cutoff:f64
    ,tm_d0:f64
    ,lnorm:usize
    ,seed_length:usize
    ,seed_dist_limit:f64
,tmscore_break:f64
,stop_tmscore_ratio:f64)->Option<(Vec<Vec<f64>>,(f64,f64))>{
    
    let dist_cutoff2:f64 = dist_cutoff*dist_cutoff;
    //let fragment_length:Vec<usize> = vec![3,5,10,15,20];
    //let fragment_length:Vec<usize> = vec![3,5];
    let fragment_length:Vec<usize> = vec![seed_length];
    let q_length:usize =query[0].len(); 
    let y_length:usize =template[0].len(); 

    let check_step:usize = 1;
    let mut max_score:(f64,f64) = (-1.0,-1.0);

    //rotationmatrix, xcenter, ycenter
    let mut max_matrix:(Vec<Vec<f64>>,(f64,f64,f64),(f64,f64,f64)) = (vec![],(0.0,0.0,0.0),(0.0,0.0,0.0));
    let mut x_to_y_indexmap:Vec<i64> = vec![-1;q_length];
    let mut qmatbuff:Vec<Vec<f64>> = query.clone();
    let mut tmatbuff:Vec<Vec<f64>> = template.clone();
    let calc_center = |xx:&Vec<Vec<f64>>,anum:usize|->(f64,f64,f64){
        let mut ret:(f64,f64,f64) = (0.0,0.0,0.0);
        for i in 0..anum{
            ret.0 += xx[0][i];
            ret.1 += xx[1][i];
            ret.2 += xx[2][i];
        }
        ret.0 /= anum as f64;
        ret.1 /= anum as f64;
        ret.2 /= anum as f64;
        return ret
    };

    for ff in fragment_length.iter(){
        let mut start:usize=0;
        let mut x_to_y_indexmap_first:Vec<(usize,i64)> = vec![(0,-1);*ff];
        loop{
            for jj in 0..(y_length-ff+1){
                for ii in x_to_y_indexmap.iter_mut(){
                    *ii = -1;
                }
                for ii in start..(start+ff){
                    if ii >= q_length{
                        break;
                    }
                    x_to_y_indexmap[ii as usize] = jj as i64 + ii as i64 - start as i64;
                    x_to_y_indexmap_first[(ii - start) as usize] = (ii as usize,jj as i64 + ii as i64 - start as i64);
                }
                if true{

                    let xi = x_to_y_indexmap_first[0].0;
                    let yi = x_to_y_indexmap_first[0].1 as usize;
                    let xi2 = x_to_y_indexmap_first[0+(*ff-1)].0;
                    let yi2 = x_to_y_indexmap_first[0+(*ff-1)].1 as usize;
                    let pdist1 = distance2(query[0][xi],query[1][xi],query[2][xi],query[0][xi2],query[1][xi2],query[2][xi2])+0.00001;
                    let pdist2 = distance2(template[0][yi],template[1][yi],template[2][yi],template[0][yi2],template[1][yi2],template[2][yi2])+0.00001;
                    
                    //println!("{}",(pdist1.sqrt()-pdist2.sqrt()).abs());
                    if (pdist1.sqrt()-pdist2.sqrt()).abs() > seed_dist_limit{
                        continue;
                    }
                    
                }

                'outiter:for itt in 0..num_iter_max{
                    let mut num_align:usize = 0;
                    for (qii,i) in x_to_y_indexmap.iter().enumerate(){
                        if *i > -1{
                            let iu = *i as usize;
                            qmatbuff[0][num_align] = query[0][qii];
                            qmatbuff[1][num_align] = query[1][qii];
                            qmatbuff[2][num_align] = query[2][qii];

                            tmatbuff[0][num_align] = template[0][iu];
                            tmatbuff[1][num_align] = template[1][iu];
                            tmatbuff[2][num_align] = template[2][iu];

                            num_align += 1;
                        }
                    }

                    let xcenter:(f64,f64,f64) = calc_center(&qmatbuff,num_align);
                    let ycenter:(f64,f64,f64) = calc_center(&tmatbuff,num_align);
                    for ii in 0..num_align{
                        qmatbuff[0][ii] -= xcenter.0;
                        qmatbuff[1][ii] -= xcenter.1;
                        qmatbuff[2][ii] -= xcenter.2;
                        
                        tmatbuff[0][ii] -= ycenter.0;
                        tmatbuff[1][ii] -= ycenter.1;
                        tmatbuff[2][ii] -= ycenter.2;
                    }

                    let rotmat_ = get_rotation_matrix(&qmatbuff,&tmatbuff,num_align);
                    if let None = rotmat_{
                    }else{
                        let rotmat = rotmat_.unwrap();
                        for ii in x_to_y_indexmap_first.iter(){
                            if ii.1 < 0{
                                continue;
                            }
                            let mres = matrix_process::matrix_multi(&rotmat
                                ,&vec![
                                vec![query[0][ii.0]-xcenter.0]
                                ,vec![query[1][ii.0]-xcenter.1]
                                ,vec![query[2][ii.0]-xcenter.2]
                                ]);
                            let dist2 = distance2(mres[0][0],mres[1][0],mres[2][0],template[0][ii.1 as usize]-ycenter.0,template[1][ii.1 as usize]-ycenter.1,template[2][ii.1 as usize]-ycenter.2);
                            if dist2 >= dist_cutoff2{
                                //seed の点の距離がカットオフ距離を超えるとこれ以上にシードとして良い場所があるとみなす。
                                break 'outiter;
                            }
                        }

                        for ii in 0..q_length{
                            let mres = matrix_process::matrix_multi(&rotmat
                                ,&vec![
                                vec![query[0][ii]-xcenter.0]
                                ,vec![query[1][ii]-xcenter.1]
                                ,vec![query[2][ii]-xcenter.2]
                                ]);
                            qmatbuff[0][ii] = mres[0][0];
                            qmatbuff[1][ii] = mres[1][0];
                            qmatbuff[2][ii] = mres[2][0];
                        }

                        for ii in 0..y_length{
                            tmatbuff[0][ii] = template[0][ii]-ycenter.0;
                            tmatbuff[1][ii] = template[1][ii]-ycenter.1;
                            tmatbuff[2][ii] = template[2][ii]-ycenter.2;
                        }

                        let scores = calc_tmscore_gdtts(&qmatbuff,&tmatbuff, &x_to_y_indexmap,dist_cutoff,tm_d0,lnorm);
                        
                        if scores.0 > max_score.0{
                            max_matrix = (rotmat,xcenter,ycenter);
                            max_score = scores;
                        }
                        
                        if itt == num_iter_max-1{
                            break;
                        }
                        
                        if itt > 10{//適当 これくらい iteration したら並ぶ場合は並ぶだろうという判断
                            if scores.0 < tmscore_break{
                                break;
                            } 

                            if scores.0 < max_score.0*stop_tmscore_ratio{
                                break;
                            }
                        }

                        let mut x_to_y_indexmap_t:Vec<(i64,f64)> = vec![(-1,dist_cutoff2);q_length];

                        for ii in 0..q_length{
                            for jj in 0..y_length{
                                let ddis2 = distance2(
                                qmatbuff[0][ii],qmatbuff[1][ii],qmatbuff[2][ii],
                                tmatbuff[0][jj],tmatbuff[1][jj],tmatbuff[2][jj]
                                );
                                if ddis2 < x_to_y_indexmap_t[ii].1{
                                    x_to_y_indexmap_t[ii] = (jj as i64,ddis2);
                                }
                            }
                        }
                        filt_discrepancy(&mut x_to_y_indexmap_t);
                        filt_short_segments(&mut x_to_y_indexmap_t);
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
                    }

                }
            }
            
            if q_length <= start+ff{
                break;
            }else{
                start += check_step;
            }
        }

    }


    let mut movx:Vec<Vec<f64>> = matrix_process::eye(4);
    movx[0][3] = -(max_matrix.1).0;
    movx[1][3] = -(max_matrix.1).1;
    movx[2][3] = -(max_matrix.1).2;


    let mut movy:Vec<Vec<f64>> = matrix_process::eye(4);
    movy[0][3] = (max_matrix.2).0;
    movy[1][3] = (max_matrix.2).1;
    movy[2][3] = (max_matrix.2).2;

    let mut m_ret_:Vec<Vec<f64>> = matrix_process::eye(4);
    for ii in 0..3{
        for jj in 0..3{
            m_ret_[ii][jj] = max_matrix.0[ii][jj];
        }
    }

    let m_ret:Vec<Vec<f64>> = matrix_process::matrix_multi(&movy,&matrix_process::matrix_multi(&m_ret_,&movx));

    return Some((m_ret,max_score));

}


#[test]

fn aligntest(){
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    let mut errorsum:f64 = 0.0;
    for _ in 0..100{
        let mut avec_:Vec<Vec<f64>> = vec![];
        let mut bvec_:Vec<Vec<f64>> = vec![];
        for _ in 0..100{
            avec_.push(vec![
            rgen.gen_range(-100.0,100.0)
            ,rgen.gen_range(-100.0,100.0)
            ,rgen.gen_range(-100.0,100.0)
            ]);
        }
        let avec = matrix_process::matrix_t(&avec_);

        let radi:f64 = rgen.gen_range(-3.14,3.14);
        let ssin:f64 = radi.sin();
        let ccos:f64 =  radi.cos();
        
        let radi:f64 = rgen.gen_range(-3.14,3.14);
        let ssin2:f64 = radi.sin();
        let ccos2:f64 =  radi.cos();
        
        let v1 = vec![
            vec![ccos,-ssin,0.0]
            ,vec![ssin,ccos,0.0]
            ,vec![0.0,0.0,1.0]
        ];
        let v2 = vec![
            vec![1.0,0.0,0.0]
            ,vec![0.0,ccos2,-ssin2]
            ,vec![0.0,ssin2,ccos2]
        ];
        let rotmat = matrix_process::matrix_multi(&v1,&v2);

        for ii in 0..avec[0].len(){
            let mres = matrix_process::matrix_multi(&rotmat,&vec![vec![avec[0][ii]],vec![avec[1][ii]],vec![avec[2][ii]]]);
            bvec_.push(vec![mres[0][0],mres[1][0],mres[2][0]]);
        }
        let bvec = matrix_process::matrix_t(&bvec_);


        let res = get_rotation_matrix(&avec,&bvec,avec[0].len()).unwrap();
        let rlen:usize = res.len();
        let clen:usize = res[0].len();
        for rr in 0..rlen{
            for cc in 0..clen{
                errorsum += (res[rr][cc] - rotmat[rr][cc]).abs();
                if (res[rr][cc] - rotmat[rr][cc]).abs() > ERROR_THRESHOLD{
                    println!("{:?}",rotmat);
                    panic!("{:?}",res);
                }
            }
        }
    }
    println!("errorsum: {}",errorsum);

}

#[test]
fn domain_splittest(){
    let mut pdb = mmcif_process::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"2gx4_A.pdb").as_str());
    let mut cas:Vec<Vec<f64>> = vec![];
    for cc in pdb.get_mut_model_at(0).get_mut_entity_at(0).iter_mut_asyms(){
        for rr in cc.iter_mut_comps(){
            for aa in rr.iter_mut_atoms(){
                if aa.atom_code == "CA"{
                    cas.push(vec![aa.get_x(),aa.get_y(),aa.get_z()]);
                }
            }
        }
    }
    let rnum:usize = cas.len();
    let s0:usize = rnum/3;
    let s1:usize = rnum/3*2;
    let mut moved = cas.clone();
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(10);

    let m1:(f64,f64,f64) = (rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0));
    let r1:(f64,f64,f64) = (rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0));
    for ii in 0..s0{
        let orig = moved[ii].clone();
        moved[ii][0] = orig[0]*r1.0.cos()-orig[1]*r1.0.sin();
        moved[ii][1] = orig[1]*r1.0.cos()+orig[0]*r1.0.sin();
        moved[ii][0] += m1.0;
        moved[ii][1] += m1.1;
        moved[ii][2] += m1.2;
    }
    
    let m1:(f64,f64,f64) = (rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0));
    let r1:(f64,f64,f64) = (rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0),rgen.gen_range(-20.0,20.0));
    for ii in s0..s1{
        let orig = moved[ii].clone();
        moved[ii][1] = orig[1]*r1.0.cos()-orig[2]*r1.0.sin();
        moved[ii][2] = orig[2]*r1.0.cos()+orig[1]*r1.0.sin();
        moved[ii][0] += m1.0;
        moved[ii][1] += m1.1;
        moved[ii][2] += m1.2;
    }
    
    let mut res:Vec<Vec<usize>> = comparative_domain_split(
    &cas
    ,&moved
    ,10
    ,1.0
    ,100
    ,4
    ,3.0
    ,0.2
    ,0.2);
    for rr in res.iter_mut(){
        rr.sort();
    }
    res.sort_by(|a,b|a[0].cmp(&b[0]));
    assert!(res[0][0] == 0);
    assert!(res[1][0] == s0);
    assert!(res[2][0] == s1);

}
#[test]
fn pdbaligntest(){
    
    let pdb_orig = mmcif_process::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"6lu7_A.pdb").as_str());
    let mut pdb = mmcif_process::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"2gx4_A.pdb").as_str());

    let mut residues_a:Vec<&pdbdata::PDBComp> = vec![];
    for cc in pdb_orig.get_all_asyms().iter(){
        for rr in cc.iter_comps(){
            residues_a.push(rr);
        }
    }
    
    let mut residues_b:Vec<&pdbdata::PDBComp> = vec![];
    for cc in pdb.get_all_asyms().iter(){
        for rr in cc.iter_comps(){
            residues_b.push(rr);
        }
    }
    let res_:Option<StructuralAlignmentResult> = align_pdb(&residues_b,&residues_a,AlignmentType::SW,0.0);
    let res = res_.unwrap();
    //pdb.save("test/testrot.pdb");
    for cc in pdb.get_mut_all_asyms(){
        for rr in cc.iter_mut_comps(){
            for aa in rr.iter_mut_atoms(){
                let mres = matrix_process::matrix_multi(&res.transform_matrix,&vec![vec![aa.get_x()],vec![aa.get_y()],vec![aa.get_z()],vec![1.0]]);
                aa.set_xyz(mres[0][0],mres[1][0],mres[2][0]);
            }
        }
    }
    println!("{:?}",res);
    //かきさし
    pdb.save("test/testrot2.pdb");
    pdb_orig.save("test/testrot_orig.pdb");
}

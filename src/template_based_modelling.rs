#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use std::fs;
use regex::Regex;

use std::f64::consts::PI;

use rand::seq::SliceRandom;
use std::collections::{HashMap,HashSet};
use super::misc_util::*;
use super::charmm_based_energy;
use super::backbone_sample;
use super::side_chain_sample;
use super::chain_builder;
use super::charmm_param;
use super::pdbdata;
#[allow(unused_imports)]
use super::mmcif_process;
use super::sequence_alignment;
use super::structural_alignment;
use super::matrix_process;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
#[allow(unused_imports)]
use rand_distr::{Normal, Distribution};

use super::geometry::Vector3D;
use super::geometry::Point3D;


use super::peptide_backbone_dihedral_energy;
use super::distance_energy;
#[allow(unused_imports)]
use super::energy_lbfgs;

use super::evoef2_energy;
use super::pp_energy::PPEnergyWeight;
use super::pp_energy::PPEnergySet;
use super::pp_energy;
use super::pp_energy_mc;
use super::process_3d;


pub struct RefinementParam{
    pub num_iter:Vec<usize>,
    pub inv_edge:Vec<u64>,
    pub inv_dist:Vec<f64>,
    pub wdec_factor:Vec<f64>,
    pub num_target_atoms:Vec<usize>,
    pub mov_dist:Vec<f64>,
    pub refine_retry:Vec<usize>,
    pub lbfgs:Vec<usize>,
    pub refine_loop_threshold:Vec<f64>,
    pub group_rotate:Vec<f64>,

    pub sorter:bool,
    pub shuffle:bool,
    pub backbone_cb_only:bool,
    pub accept_bound:(f64,f64),
    pub steps_bound_check:usize,
    pub use_template:bool
}
impl RefinementParam{
    pub fn generate_lines(&self)->Vec<RefinementParamLine>{
        let mut ret:Vec<RefinementParamLine> = vec![];
        for ni in self.num_iter.iter(){
            for ie in self.inv_edge.iter(){
                for id in self.inv_dist.iter(){
                    for df in self.wdec_factor.iter(){
                        for nt in self.num_target_atoms.iter(){
                            for md in self.mov_dist.iter(){
                                for ra in self.refine_retry.iter(){
                                    for lb in self.lbfgs.iter(){
                                        for rl in self.refine_loop_threshold.iter(){
                                            for rg in self.group_rotate.iter(){
                                                ret.push(
                                                    RefinementParamLine{
                                                        num_target_atoms:*nt,
                                                        inv_edge:*ie,
                                                        inv_dist:*id,
                                                        dec_factor:*df,
                                                        mov_dist:*md,
                                                        num_iter:*ni,
                                                        refine_retry:*ra,
                                                        lbfgs:*lb,
                                                        refine_loop_threshold:*rl,
                                                        group_rotate:*rg
                                                    }
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return ret;
    }
}
#[derive(Debug)]
pub struct RefinementParamLine{
    num_target_atoms:usize,
    inv_edge:u64,
    inv_dist:f64,
    dec_factor:f64,
    mov_dist:f64,
    num_iter:usize,
    lbfgs:usize,
    refine_retry:usize,
    refine_loop_threshold:f64,
    group_rotate:f64,
}
pub fn try_mapping(fastafiles:&Vec<String>,backbonedir:&str,rotamerdir:&str,random_seed:Option<u64>,outfilename_:&str){
    let mut rgen:StdRng;//乱数生成器
    //Dihedral が低いので値確かめ
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }

    let mut chain_hm:HashMap<String,Vec<pdbdata::PDBComp>> =HashMap::new();
    let mut flag_string:Vec<String> = vec![];//chain residuename residue_index mapped_or_not(bool)
    
    let outfilename = outfilename_.to_string();
    let bset = backbone_sample::BackboneSet::new(backbonedir);
    let sset = side_chain_sample::SideChainSet::new(rotamerdir);

    let mut residue_count:usize = 0;
    for fastafile in fastafiles.iter(){
        let seqs:Vec<sequence_alignment::SeqData> = sequence_alignment::SeqData::load_fasta(fastafile,false);

        let mut template_file:String = "".to_owned();
        let mut template_chain:String = "".to_owned();
        let filex:Regex = Regex::new("file=([^\\s]+)").unwrap();
        let chainx:Regex = Regex::new("chain=([^\\s]+)").unwrap();
        let startx:Regex = Regex::new("start=([^\\s]+)").unwrap();
        let mut template_seq:Vec<String> = vec![];
        let mut query_seq:Vec<String> = vec![];
        let mut query_start_:Option<i64> = None;
        let mut query_chain_name:String = "A".to_owned();
        for ss in seqs.into_iter(){
            if ss.name == "template"{
                if let Some(x) = filex.captures(&ss.desc){
                    template_file = x.get(1).unwrap().as_str().to_owned();
                }else{
                    panic!("Template must have 'file' section!");
                }
                if let Some(x) = chainx.captures(&ss.desc){
                    template_chain = x.get(1).unwrap().as_str().to_owned();
                }else{
                    panic!("Template must have 'chain' section!");
                }
                template_seq = ss.seq;
            }else if ss.name == "query"{
                if let Some(x) = startx.captures(&ss.desc){
                    query_start_ = Some(x.get(1).unwrap().as_str().parse::<i64>().unwrap());
                }
                if let Some(x) = chainx.captures(&ss.desc){
                    query_chain_name = x.get(1).unwrap().as_str().to_owned();
                }
                query_seq = ss.seq;
            }else{
                panic!("The alignment file must contauins only 'query' and 'template'.");
            }
        }
        
        let pdbb:pdbdata::PDBEntry = mmcif_process::load_pdb(&template_file);
        let mut targetchain:Option<pdbdata::PDBAsym> =None;
        for mut cc in pdbb.get_model_at(0).get_entity_at(0).iter_mut_asyms(){
            cc.remove_alt(None);
            if cc.chain_name == template_chain{
                targetchain = Some(*cc);
                break;
            }
        }
        if let None = targetchain{
            panic!("The chain {} was not found.",template_chain);
        }
        if !chain_hm.contains_key(&query_chain_name){
            chain_hm.insert(query_chain_name.clone(),vec![]);
        }
        let query_start:i64;
        if let None = query_start_{
            if chain_hm.get(&query_chain_name).unwrap().len() == 0{
                query_start = 1;
            }else{
                let p = chain_hm.get(&query_chain_name).unwrap();
                query_start = p[p.len()-1].get_seq_id()+1;
            }
        }else{
            query_start = query_start_.unwrap();
        }

        let ress_and_flag:Vec<(pdbdata::PDBComp,bool)> = chain_builder::build_from_alignment(&query_seq
            ,&template_seq,&(targetchain.unwrap().iter_comps().map(|m|{m}).collect()),&bset,&sset);
            
        let mut floating_mostclose:HashMap<usize,(f64,f64,f64)> = HashMap::new();
        let llen = ress_and_flag.len();
        for ii in 0..llen{
            if !ress_and_flag[ii].1{
                for kk in 0..llen{
                    if ii >= kk && ress_and_flag[ii-kk].1{
                        floating_mostclose.insert(ii,ress_and_flag[ii-kk].0.get_CA().unwrap().get_xyz());
                        break;
                    }
                    if ii+kk < llen && ress_and_flag[ii+kk].1{
                        floating_mostclose.insert(ii,ress_and_flag[ii+kk].0.get_CA().unwrap().get_xyz());
                        break;
                    }
                }
            }
        }

        for (ii,(mut rr,ff)) in ress_and_flag.into_iter().enumerate(){
            rr.set_seq_id(query_start+ii as i64);
            if !ff{
                let ca:&pdbdata::PDBAtom = rr.get_CA().unwrap();
                let mut pbase = floating_mostclose.get(&ii).unwrap_or(&(rgen.gen_range(-6.0,6.0),rgen.gen_range(-6.0,6.0),rgen.gen_range(-6.0,6.0))).clone();
                pbase.0 -= ca.get_x();
                pbase.1 -= ca.get_y();
                pbase.2 -= ca.get_z();
                let xx:f64 = rgen.gen_range(-3.0,3.0)+pbase.0;
                let yy:f64 = rgen.gen_range(-3.0,3.0)+pbase.1;
                let zz:f64 = rgen.gen_range(-3.0,3.0)+pbase.2;
                for aa in rr.iter_mut_atoms(){
                    let p:(f64,f64,f64) = aa.get_xyz();
                    aa.set_xyz(xx+p.0,yy+p.1,zz+p.2);
                }
            }
            flag_string.push(format!("chain:{}\tresidue_name:{}\tresidue_number:{}\tplaced:{}",query_chain_name,rr.get_comp_id(),rr.get_seq_id(),ff));
            chain_hm.get_mut(&query_chain_name).unwrap().push(rr);
            residue_count += 1;
        }
    }

    let mut pdb = pdbdata::PDBEntry::prepare_base();
    pdb.entry_id = "XXXX".to_owned();
    let mut chainnames:Vec<String> = chain_hm.iter().map(|m|m.0.to_string()).collect();
    chainnames.sort();
    for cc in chainnames.into_iter(){
        let mut chain:pdbdata::PDBAsym = pdbdata::PDBAsym::new(&cc);
        let mut cvec = chain_hm.get_mut(&cc).unwrap();
        while cvec.len() > 0{
            chain.add_comp(cvec.remove(0));
        }
        pdb.get_mut_model_at(0).get_entity_at(0).add_asym(chain);
    }
    pdb.save(&outfilename);
    write_to_file(&(outfilename+".flag"),flag_string);
}


pub fn load_missing_and_freezed_atoms(flag_file:&str,atoms:&Vec<charmm_based_energy::MDAtom>)->(HashSet<usize>,HashSet<usize>){
    let file = File::open(flag_file).unwrap();
    let reader = BufReader::new(file);
    let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
    let mut chain_resnum_to_resname:HashMap<String,String> =HashMap::new();//ずれのチェック用
    let mut missing_chain_resnum:HashSet<String> = HashSet::new();
    let mut freezed_chain_resnum:HashSet<String> = HashSet::new();
    for (_lcount,line) in lines.iter().enumerate() {
        if start_with(line,"#"){
            continue;
        }
        let hs = line_to_hash(line);
        let chain_resnum = hs.get("chain").unwrap().to_string()+"#"+hs.get("residue_number").unwrap();
        if hs.get("missing").unwrap_or(&("".to_owned())) == "true" || hs.get("placed").unwrap_or(&("".to_owned())) == "false"{
            missing_chain_resnum.insert(chain_resnum.clone());
            chain_resnum_to_resname.insert(chain_resnum.clone(),hs.get("residue_name").unwrap().clone());
        }
        if hs.get("freezed").unwrap_or(&("".to_owned())) == "true"{
            freezed_chain_resnum.insert(chain_resnum.clone());
            chain_resnum_to_resname.insert(chain_resnum,hs.get("residue_name").unwrap().clone());
        }
    }
    let mut ret_fre:HashSet<usize> = HashSet::new();    
    let mut ret_miss:HashSet<usize> = HashSet::new();    
    
    for (aii,aa) in atoms.iter().enumerate(){
        let chain_resnum:String = aa.chain_name.to_string()+"#"+aa.residue_number.to_string().as_str();
        if missing_chain_resnum.contains(&chain_resnum){
            ret_miss.insert(aii);
            if chain_resnum_to_resname.get(&chain_resnum).unwrap() != &aa.residue_name{
                if chain_resnum_to_resname.get(&chain_resnum).unwrap() == "HIS" && &aa.residue_name == "HSD"{
                }else{
                    eprintln!("Discrepancy: {} {} {}. is it ok?",chain_resnum,chain_resnum_to_resname.get(&chain_resnum).unwrap(),aa.residue_name);
                }
            }
        } 
        if freezed_chain_resnum.contains(&chain_resnum){
            ret_fre.insert(aii);
            if chain_resnum_to_resname.get(&chain_resnum).unwrap() != &aa.residue_name{
                if chain_resnum_to_resname.get(&chain_resnum).unwrap() == "HIS" && &aa.residue_name == "HSD"{
                }else{
                    eprintln!("Discrepancy: {} {} {}. is it ok?",chain_resnum,chain_resnum_to_resname.get(&chain_resnum).unwrap(),aa.residue_name);
                }
            }
        } 
    }
    return (ret_miss,ret_fre);
}


pub fn load_group_atoms(flag_file:&str,atoms:&Vec<charmm_based_energy::MDAtom>)->Vec<Vec<usize>>{
    let file = File::open(flag_file).unwrap();
    let reader = BufReader::new(file);
    let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
    let mut chain_resnum_to_resname:HashMap<String,String> =HashMap::new();//ずれのチェック用
    let mut chain_resnum_to_group:HashMap<String,Vec<usize>> = HashMap::new();//chain_resnum_string,groupid
    let mut group_atomid:HashMap<usize,HashSet<usize>> = HashMap::new();//groupid,atom

    for (_lcount,line) in lines.iter().enumerate() {
        if start_with(line,"#"){
            continue;
        }
        let hs = line_to_hash(line);
        let chain_resnum = hs.get("chain").unwrap().to_string()+"#"+hs.get("residue_number").unwrap();
        //if !hs.contains_key("group"){
        //    panic!("group file must contains group section.");
        //}
        
        chain_resnum_to_resname.insert(chain_resnum.clone(),hs.get("residue_name").unwrap().clone());
        
        if !chain_resnum_to_group.contains_key(&chain_resnum){
            chain_resnum_to_group.insert(chain_resnum.clone(),vec![]);
        }
        if hs.contains_key("group"){
            let groups:Vec<&str> = hs.get("group").unwrap().split(",").collect();
            for gg in groups.into_iter(){
                if gg.len() > 0{
                    let u:usize = gg.parse::<usize>().unwrap_or_else(|e|panic!("{}",e));
                    chain_resnum_to_group.get_mut(&chain_resnum).unwrap().push(u);
                }
            }
        }
    }
    for (aii,aa) in atoms.iter().enumerate(){
        let chain_resnum:String = aa.chain_name.to_string()+"#"+aa.residue_number.to_string().as_str();
        if chain_resnum_to_group.contains_key(&chain_resnum){
            if chain_resnum_to_resname.get(&chain_resnum).unwrap() != &aa.residue_name{
                if chain_resnum_to_resname.get(&chain_resnum).unwrap() == "HIS" && &aa.residue_name == "HSD"{

                }else{
                    eprintln!("Discrepancy: {} {} {}. is it ok?",chain_resnum,chain_resnum_to_resname.get(&chain_resnum).unwrap(),aa.residue_name);
                }
            }
            for uu in chain_resnum_to_group.get(&chain_resnum).unwrap().iter(){
                if !group_atomid.contains_key(uu){
                    group_atomid.insert(*uu,HashSet::new());
                }
                group_atomid.get_mut(uu).unwrap().insert(aii);
            }
        } 
    }
    let mut ret:Vec<Vec<usize>> = vec![];
    for gg in group_atomid.into_iter(){
        let mut gv:Vec<usize> = gg.1.into_iter().collect();
        gv.sort();
        ret.push(gv);
    }
    return ret;
}
pub fn make_residue_group(atoms:&Vec<charmm_based_energy::MDAtom>)->Vec<Vec<usize>>{
    let mut hss:HashMap<String,Vec<usize>> = HashMap::new();
    for (aii,aa) in atoms.iter().enumerate(){
        let chain_resnum:String = aa.chain_name.to_string()+"#"+aa.residue_number.to_string().as_str();
        if ! hss.contains_key(chain_resnum.as_str()){
            hss.insert(chain_resnum.clone(),vec![]);
        }
        hss.get_mut(&chain_resnum).unwrap().push(aii);
    }
    let mut ret:Vec<Vec<usize>> =  hss.into_iter().map(|m|m.1).collect();
    for vv in ret.iter_mut(){
        let vlen:usize = vv.len();
        let mut n:i64 = -1;
        let mut c:i64 = -1;
        let mut np:i64 = -1;
        let mut cp:i64 = -1;
        for vvv in vv.iter().enumerate(){
            if atoms[*(vvv.1)].atom_name == "C"{
                cp = vvv.0 as i64;
                c = *vvv.1 as i64;
            }
            if atoms[*(vvv.1)].atom_name == "N"{
                np = vvv.0 as i64;
                n = *vvv.1 as i64;

            }
        }
        if c > -1{
            let en = vv[vlen-1];
            vv[vlen-1] = c as usize;
            vv[cp as usize] = en;
        }
        if n > -1{
            let st = vv[0];
            vv[0] = n as usize;
            if np as usize == vlen-1{
                vv[cp as usize] = st;
            }else{
                vv[np as usize] = st;
            }
        }
        let mut chk:HashSet<usize> = HashSet::new();
        for v in vv.iter(){
            if chk.contains(v){
                panic!();
            }
            chk.insert(*v);
        }
    }
    return ret;
}

pub fn get_sidechain_atoms(atoms:&Vec<charmm_based_energy::MDAtom>)->HashSet<usize>{
    let mut ret:HashSet<usize> = HashSet::new();    
    for (aii,aa) in atoms.iter().enumerate(){
        if aa.atom_name == "CA"
        || aa.atom_name == "C"
        || aa.atom_name == "N"{
            
        }else{
            ret.insert(aii);
        }
    }
    return ret;
}
pub fn get_backbone_atoms(atoms:&Vec<charmm_based_energy::MDAtom>)->HashSet<usize>{
    let mut ret:HashSet<usize> = HashSet::new();    
    for (aii,aa) in atoms.iter().enumerate(){
        if aa.atom_name == "CA"
        || aa.atom_name == "C"
        || aa.atom_name == "N"{
            ret.insert(aii);
        }
    }
    return ret;
}

pub fn fit_to_chain_template(sub_energyset:&mut PPEnergySet,sub_chain:&Vec<Vec<usize>>,template_n_ca_c:&Vec<(Vec<f64>,Vec<f64>,Vec<f64>)>){
    let mut q_cas:Vec<Vec<f64>> = vec![];
    let mut q_idx:Vec<usize> = vec![];
    for (avii,avv) in sub_chain.iter().enumerate(){
        for aii in avv.iter(){
            if &sub_energyset.get_atom(*aii).atom_name == "CA"{
                let xyz = sub_energyset.get_atom(*aii).get_xyz();
                q_cas.push(vec![xyz.0,xyz.1,xyz.2]);
                q_idx.push(avii);
                break;
            }
        } 
    }
    let tvec:Vec<Vec<f64>> = template_n_ca_c.iter().map(|m|m.1.clone()).collect();
    let lnorm:usize = q_cas.len().min(tvec.len());
    let d0:f64 = if lnorm <=19{
        0.5
    }else{
        1.24*((lnorm as f64)-15.0).powf(1.0/3.0)-1.8
    };
    let d0_search:f64 = d0.max(4.5).min(8.0);
    
    //temp の方を動かす
    
    let res = structural_alignment::align_unrelated(
        &matrix_process::matrix_t(&tvec)
        ,&matrix_process::matrix_t(&q_cas)
        ,100
        ,d0_search
        ,d0
        , lnorm,4,4.0,0.2,0.3
        );
    if let None = res{
        return;
    }
    if (res.as_ref().unwrap().1).0 < 0.5{
        return;
    }
    let rotmat:Vec<Vec<f64>> = res.unwrap().0;
    let mut ttvec:Vec<((f64,f64,f64),(f64,f64,f64),(f64,f64,f64))> = vec![];
    for tvv in template_n_ca_c.iter(){
        let n = matrix_process::matrix_multi(&rotmat,&vec![vec![tvv.0[0]],vec![tvv.0[1]],vec![tvv.0[2]],vec![1.0]]);
        let ca = matrix_process::matrix_multi(&rotmat,&vec![vec![tvv.1[0]],vec![tvv.1[1]],vec![tvv.1[2]],vec![1.0]]);
        let c = matrix_process::matrix_multi(&rotmat,&vec![vec![tvv.2[0]],vec![tvv.2[1]],vec![tvv.2[2]],vec![1.0]]);
        ttvec.push(((n[0][0],n[1][0],n[2][0]),(ca[0][0],ca[1][0],ca[2][0]),(c[0][0],c[1][0],c[2][0])));
    }
    let tlen = ttvec.len();
    let qlen = q_idx.len();
    let mut qdiff:Vec<(i64,f64)> = vec![(-1,8.0);qlen];
    for qq in 0..qlen{
        let mut mindist:f64 = 8.0;
        let mut minid:i64 = -1;
        let qt:(f64,f64,f64) = (q_cas[qq][0],q_cas[qq][1],q_cas[qq][2]);
        for tt in 0..tlen{
            let dist = process_3d::distance(&qt,&(ttvec[tt].1));
            if dist < mindist{
                mindist = dist;
                minid = tt as i64;
            }
        }
        qdiff[qq].0 = minid;
        qdiff[qq].1 = mindist;
    }
    structural_alignment::filt_discrepancy(&mut qdiff);
    let avv:&Vec<Vec<usize>> = sub_chain;
    for qq in 0..qlen{
        if qdiff[qq].0 > -1{
            let mut n:i64 = -1;
            let mut ca:i64 = -1;
            let mut c:i64 = -1;
            for z in avv[q_idx[qq]].iter(){
                if sub_energyset.get_atom(*z).atom_name == "N"{
                    n = *z as i64;
                }
                if sub_energyset.get_atom(*z).atom_name == "CA"{
                    ca = *z as i64;
                }
                if sub_energyset.get_atom(*z).atom_name == "C"{
                    c = *z as i64;
                }
            }
            
            if n < 0 || ca < 0 || c < 0{
                continue;
            }
            let pn:Point3D = Point3D::from_vector3d(sub_energyset.get_atom(n as usize));
            let pca:Point3D = Point3D::from_vector3d(sub_energyset.get_atom(ca as usize));
            let pc:Point3D = Point3D::from_vector3d(sub_energyset.get_atom(c as usize));
            
            let ptn:Point3D = Point3D::from_tuple(&ttvec[qdiff[qq].0 as usize].0);
            let ptca:Point3D = Point3D::from_tuple(&ttvec[qdiff[qq].0 as usize].1);
            let ptc:Point3D = Point3D::from_tuple(&ttvec[qdiff[qq].0 as usize].2);
            
            let mut rvec:Vec<&mut dyn Vector3D> = vec![];
            let hss:HashSet<usize> = avv[q_idx[qq]].iter().map(|m|*m).collect();
            for (zii,z) in sub_energyset.evoef2_env.md_envset.atoms.iter_mut().enumerate(){
                if hss.contains(&zii){
                    rvec.push(z);
                }
            }
            process_3d::fit_to_vector(&pn,&pca,&pc,&ptn,&ptca,&ptc,&mut rvec);
            //たしか三角形の起点を合わせるので
            let mut cdiff:(f64,f64,f64) = ptca.get_xyz();
            let caa = sub_energyset.get_atom(ca as usize).get_xyz();
            cdiff.0 -= caa.0;
            cdiff.1 -= caa.1;
            cdiff.2 -= caa.2;
            sub_energyset.get_atom_mut(c as usize).set_xyz(ptc.get_x(),ptc.get_y(),ptc.get_z());        
            sub_energyset.get_atom_mut(n as usize).set_xyz(ptn.get_x(),ptn.get_y(),ptn.get_z());        
            
            if cdiff.0.abs() > 0.01 ||  cdiff.1.abs() > 0.01 ||  cdiff.2.abs() > 0.01{
                for aii in avv[q_idx[qq]].iter(){
                    if n as usize == *aii || c as usize == *aii{
                        continue;
                    }
                    let att_p = sub_energyset.get_atom(*aii).get_xyz();
                    sub_energyset.get_atom_mut(*aii).set_xyz(att_p.0+cdiff.0,att_p.1+cdiff.1,att_p.2+cdiff.2);
                }
            }

            /*
            //隣の残基がインサーション扱いになっている場合 C と N は元の場所に置こうとしたが、Group で動かしているとズレるのでやめた
            if qq > 0 && qdiff[qq-1].0 < 0{
                sub_energyset.get_atom_mut(n as usize).set_xyz(pn.get_x(),pn.get_y(),pn.get_z());
            }
            if qq < qlen-1 && qdiff[qq+1].0 < 0{
                sub_energyset.get_atom_mut(c as usize).set_xyz(pc.get_x(),pc.get_y(),pc.get_z());
            }
            */
        }
    }
}

//オレオレフォーマットの
//c1:xxx\tr1:xxx\ti1:xxx\tc2:xxx\tr2:xxx\ti2:xxx\tvalue:xxx
//みたいなファイルを読み込んで
//chainname#residuennumber#insersioncode
//chainname#residuennumber#insersioncode
//value(lower is better)
//の Vec で返す。
pub fn load_contact_file(filename:&str)->Vec<(String,String,f64)>{
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut ret:Vec<(String,String,f64)> = vec![];
    let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
    let mut positive_flag:bool = false;
    for (_lcount,line) in lines.iter().enumerate() {
        if start_with(line,"#"){
            continue;
        }
        let hs = line_to_hash(line);
        let c1:String = hs.get("c1").unwrap_or(&("".to_owned())).to_string();
        let c2:String = hs.get("c2").unwrap_or(&("".to_owned())).to_string();
        let r1 = hs.get("r1").unwrap_or_else(||panic!("Line must have r1."));
        let r2 = hs.get("r2").unwrap_or_else(||panic!("Line must have r2."));
        let i1:String = hs.get("i1").unwrap_or(&("".to_owned())).to_string();
        let i2:String = hs.get("i2").unwrap_or(&("".to_owned())).to_string();
        let val:f64 = hs.get("value").unwrap_or_else(||panic!("Line must have value.")).parse::<f64>().unwrap();
        if val > 0.0{
            positive_flag = true;
        }
        ret.push((
            c1+"#"+r1+"#"+&i1
            ,c2+"#"+r2+"#"+&i2
            ,val
        ));

    }
    if positive_flag{
        eprintln!("Warning: Positive value was found in contact file. (The value was considered as energy thus lower is better.)");
    }
    return ret;
}

pub fn get_md_residue_key_from_atom(atom:&charmm_based_energy::MDAtom)->(String,i64){
    return (atom.chain_name.clone(),atom.residue_index_in_chain);
}

pub fn get_md_residue_key(chain_id:&str,residue_inde_in_chain:i64)->(String,i64){
    return (chain_id.to_string(),residue_inde_in_chain);
}

pub struct MissingBuilderParam{
    pub param1:Vec<RefinementParam>,
    pub param2:Vec<RefinementParam>,
    pub phipsi_top_x:usize,
    pub num_structures_step1:usize,
    pub num_structures_step2:usize,
    pub fix_omega:bool
}

pub fn merge_structure(
    basestructure:&Vec<pdbdata::PDBComp>
    ,samplestructure:&Vec<pdbdata::PDBComp>
    ,chain_name:&str
    ,topfile:&str
    ,paramfile:&str
    ,evoefdir:&str
    ,backbone_dihedral_angle_file:&str
    ,cb_distance_file:&str
    ,residue_contact_file:&str
    ,params:&Vec<RefinementParam>
    ,random_seed:Option<u64>
    ,outfile:&str
    ,top_x:usize){
    
    let mut qresidues:Vec<&pdbdata::PDBComp> = vec![];
    let mut tresidues:Vec<&pdbdata::PDBComp> = vec![];

    for rr in samplestructure.iter(){
        qresidues.push(rr); 
    }
    
    for rr in basestructure.iter(){
        tresidues.push(rr); 
    }

    let num_residues:usize = qresidues.len();

    if tresidues.len() != qresidues.len() {
        panic!("The two structures must be same!");
    }
    let res_ = structural_alignment::align_pdb(&qresidues,&tresidues,structural_alignment::AlignmentType::SW.clone(),0.2);
    let res:structural_alignment::StructuralAlignmentResult = res_.unwrap();
    let aligned_length:usize = res.aligned_chain1.len();
    let mut sampcount:usize = 0;
    let mut basecount:usize = 0;
    let mut residue_aligned:HashSet<usize> = HashSet::new();
    for xx in 0..aligned_length{
        if res.aligned_chain1[xx] != "-"
        &&  res.aligned_chain2[xx] != "-"{
            if sampcount == basecount{ 
                if res.aligned_chain1[xx] != res.aligned_chain2[xx]{
                    panic!("!?????");
                }
                residue_aligned.insert(basecount);
            }
        }
        if res.aligned_chain1[xx] != "-"{
            sampcount +=1;
        }
        
        if res.aligned_chain2[xx] != "-"{
            basecount +=1;
        }
    }
    let mut unaligned_region:Vec<Vec<usize>> = vec![];
    if !residue_aligned.contains(&0){
        unaligned_region.push(vec![0]);
    }
    for aa in 1..num_residues{
        if !residue_aligned.contains(&aa){
            if !residue_aligned.contains(&(aa-1)){
                unaligned_region.last_mut().unwrap().push(aa);
            }else{
                unaligned_region.push(vec![aa]);
            }
        }
    }
    if unaligned_region.len() == 0{
        panic!("The two structures were completely aligned.");
    }
    let mut top_x_structures:Vec<(f64,Vec<(String,(String,i64,String),pdbdata::PDBAtom)>)> = vec![];
    for (str_count,vv) in unaligned_region.iter().enumerate(){
        let mut pdbb:pdbdata::PDBEntry = pdbdata::PDBEntry::prepare_base();
        let mut chain:pdbdata::PDBAsym = pdbdata::PDBAsym::new(chain_name);
        let swapper:HashSet<usize> = vv.iter().map(|m|*m).collect();
        for rr in 0..num_residues{
            if swapper.contains(&rr){
                chain.add_comp(samplestructure[rr].get_copy_wo_parents());
            }else{
                chain.add_comp(basestructure[rr].get_copy_wo_parents());
            }
        }
        pdbb.get_mut_model_at(0).get_mut_entity_at(0).add_asym(chain);
        let res:pp_energy::PPEnergySet = refinement(pdbb
            ,None
            ,None
            ,None
            ,0
            ,topfile
            ,paramfile
            ,evoefdir
            ,backbone_dihedral_angle_file
            ,cb_distance_file
            , residue_contact_file
            ,params
            ,None
            ,None
            ,random_seed.clone()
            ,0
            ,outfile);
        let mdatomnum = res.evoef2_env.md_envset.atoms.len();
        let mut atom_level_energy:Vec<f64> = vec![0.0;mdatomnum];
        let ene:f64 = res.calc_energy(&mut atom_level_energy);
        if top_x > 0{
            if top_x_structures.len() >= top_x{
                top_x_structures.sort_by(|a,b| (a.0).partial_cmp(&b.0).unwrap());
                if ene > top_x_structures.last().unwrap().0{
                    top_x_structures.pop();
                }
            }
            if top_x_structures.len() < top_x{
                let mut vvec:Vec<(String,(String,i64,String),pdbdata::PDBAtom)> = vec![]; 
                for mm in 0..mdatomnum{
                    let  (chainid,(resname,resnum,altcode),mut att) = res.get_atom(mm).to_pdbatom();
                    att.temp_factor = atom_level_energy[mm].max(0.0).min(999.0);
                    vvec.push((chainid,(resname,resnum,altcode),att));
                }
            }
        }else{
            let mut lines:Vec<String> = vec![];
            for mm in 0..mdatomnum{
                let  (chainid,(resname,resnum,altcode),mut att) = res.get_atom(mm).to_pdbatom();
                att.temp_factor = atom_level_energy[mm].max(0.0).min(999.0);
                lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
            }
            write_to_file(format!("{}.{}.pdb",outfile,str_count).as_str(),lines);
        }
        
    }
    if top_x > 0{
        for (str_count,strs) in top_x_structures.into_iter(){
            let mut lines:Vec<String> = vec![];
            for (chainid,(resname,resnum,altcode),att) in strs.into_iter(){
                lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
            }
            write_to_file(format!("{}.{}.pdb",outfile,str_count).as_str(),lines);
        }
    }
}

pub fn docking(pdbb:pdbdata::PDBEntry
    ,groupfile:Option<String>
    ,topfile:&str
    ,paramfile:&str
    ,evoefdir:&str
    ,backbone_dihedral_angle_file:&str
    ,cb_distance_file:&str
    ,residue_contact_file:&str
    ,params:&Vec<RefinementParam>
    ,num_structures_step2:usize
    ,resolution:Option<f64>
    ,random_seed:Option<u64>
    )->PPEnergySet{
    let mut rgen:StdRng;//乱数生成器
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }
    let parr = charmm_param::CHARMMParam::load_chamm19(topfile,paramfile);
    let (mut md_envset,md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&pdbb.get_all_asyms().iter().collect(),&parr,true);
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(
        &backbone_dihedral_angle_file);
    let (torsion,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,false,true),false);

    let mut distbe:Vec<distance_energy::AtomBinnedDistanceEnergy> = vec![];
    if cb_distance_file.len() > 0{
        let dist_bvvec:Vec<(String,usize,String,usize,distance_energy::AtomBinnedDistanceEnergy)> = distance_energy::AtomBinnedDistanceEnergy::load_name_mapped(cb_distance_file);
        distbe = distance_energy::AtomBinnedDistanceEnergy::assign_atom_ids(&md_envset,dist_bvvec);
    }
    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
    println!("energy:{:?}",scoresum);
   
    //dihed が重複して計算されているのでウェイトをつける。今は適当に 0.5
    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&md_varset.dihedvec,&torsion,&omegas_general);
    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.5;
    }

    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,evoefdir,false),
        backbone_energy_omega:omegas_general,
        backbone_energy_phi_psi:torsion,
        atom_distance_energy:vec![],
        atom_binned_distance_energy:distbe,
        atom_contact_energy:vec![],
        weights:PPEnergyWeight::new()
    };
    ppenergyset.weights.dihed_weight_charmm = weight_dihed;
    ppenergyset.weights.charmm_electro = 0.0;
    ppenergyset.weights.charmm_lj = 1.0;
    
    ppenergyset.weights.backbone_energy_phi_psi = 1.0;
    ppenergyset.weights.backbone_energy_omega = 1.0;


    let sorted_atoms:HashMap<String,Vec<Vec<usize>>> = pp_energy_mc::make_sorted_residue_group(&ppenergyset.evoef2_env.md_envset.atoms);
    let mut groupatoms:Vec<Vec<usize>> = if let Some(x) = groupfile.as_ref(){
        load_group_atoms(&x,&ppenergyset.evoef2_env.md_envset.atoms)
    }else{
        panic!("This mode needs group!");
    };

    groupatoms.sort_by(|a,b|a.len().cmp(&b.len()));
    groupatoms.reverse();
    
    let mut residue_label_to_indexinchain:HashMap<String,(String,i64)> = HashMap::new();
    for aa in ppenergyset.evoef2_env.md_envset.atoms.iter(){
        let label:String = aa.chain_name.clone()+"#"+&(aa.residue_number.to_string())+"#"+&aa.residue_ins_code;
        if !residue_label_to_indexinchain.contains_key(&label){
            //println!("{:?}",(aa.chain_name.clone(),aa.residue_index_in_chain));
            residue_label_to_indexinchain.insert(label,(aa.chain_name.clone(),aa.residue_index_in_chain));
        }else{
            if residue_label_to_indexinchain.get(&label).unwrap() != &(aa.chain_name.clone(),aa.residue_index_in_chain){
                eprintln!("{} is assigned to {:?} but {:?} was going to assin. Some function might not run successfully.",
                label,
                residue_label_to_indexinchain.get(&label).unwrap(),
                (aa.chain_name.clone(),aa.residue_index_in_chain),
                );
            }
        }
    }

    //ある原子の distance ではなく Residue に存在する全原子のうちのいずれかが 6.0 angstrom 内である場合に Energy を加える。よってマイナスである必要がある
    //6.0 は evoef2 の attr の cutoff distance
    //chainname+"#"+residue_number#insertion_code
    let mut residue_contact_:Vec<(String,String,f64)> = if residue_contact_file.len() > 0{
        load_contact_file(residue_contact_file)
    }else{
        vec![]
    };
    residue_contact_.sort_by(|a,b| a.2.partial_cmp(&b.2).unwrap_or_else(||panic!("Error in sort contact.")));
    
    let mut residue_contact:Vec<(String,i64,String,i64,f64)> = vec![];
    for (r1_,r2_,val) in residue_contact_.iter(){
        if !residue_label_to_indexinchain.contains_key(r1_.as_str()){
            println!("{} is not found in pdb file.",r1_);
            continue;
        }
        if !residue_label_to_indexinchain.contains_key(r2_.as_str()){
            println!("{} is not found in pdb file.",r2_);
            continue;
        }
        let r1:&(String,i64) = residue_label_to_indexinchain.get(r1_.as_str()).unwrap();
        let r2:&(String,i64) = residue_label_to_indexinchain.get(r2_.as_str()).unwrap();
        residue_contact.push((r1.0.clone(),r1.1,r2.0.clone(),r2.1,*val));
    }
    for rr in residue_contact.iter(){
        let (c1,r1,c2,r2,v) = rr;
        let mut atoms_a:Vec<usize> = vec![];
        for ss in sorted_atoms.get(c1).unwrap_or_else(||panic!("error processing contact.")).iter(){
            if ppenergyset.get_atom(ss[0]).residue_index_in_chain == *r1{
                atoms_a.append(&mut ss.clone());
            }
        }
        let mut atoms_b:Vec<usize> = vec![];
        for ss in sorted_atoms.get(c2).unwrap_or_else(||panic!("error processing contact.")).iter(){
            if ppenergyset.get_atom(ss[0]).residue_index_in_chain == *r2{
                atoms_b.append(&mut ss.clone());
            }
        }

        if atoms_a.len() > 0 && atoms_b.len() > 0{
            ppenergyset.atom_contact_energy.push(pp_energy::AtomContactEnergy{
                atoms_a:atoms_a,
                atoms_b:atoms_b,
                threshold:6.0,//ここ CB の時は 8.0 にする
                linear_penalty:0.1,
                energy:*v
            });
        }else{
            eprintln!("???contact {:?} was not found. a {} b {}", (c1,r1,c2,r2,v) ,atoms_a.len(),atoms_b.len());
        }
    }
    let glen = groupatoms.len();
    if glen == 1{
        panic!("This mode needs more than 2 groups.");
    }
    let mut base:HashSet<usize> = groupatoms[0].iter().map(|m|*m).collect();
    
    for gg in 1..glen{
        let mover:HashSet<usize> = groupatoms[gg].iter().map(|m|*m).collect();
        zatsuna_docking(
            &mut ppenergyset
            ,& base
            ,& mover
            ,resolution
            ,&mut rgen
            ,params
            ,num_structures_step2
        );
        for mm in mover.into_iter(){
            base.insert(mm);
        }
        let mut lines:Vec<String> = vec![];
        for (_,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
            let  (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
            lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
        }
        write_to_file(format!("test.{}.pdb",gg).as_str(),lines);
    }

    return ppenergyset;
}
pub fn zerofill(arr:&mut Vec<f64>){
    for aa in arr.iter_mut(){
        *aa = 0.0;
    }
}

pub fn zatsuna_docking(ppenergyset:&mut PPEnergySet
    ,base:&HashSet<usize>
    ,mover:&HashSet<usize>
    ,full_assign_resolution:Option<f64>
    ,rgen:&mut StdRng
    ,params:&Vec<RefinementParam>//全体最適化するためのパラメータ
    ,step2_trial:usize//作成した仮構造のうち何個を最適化するか
){
    let anum:usize = ppenergyset.get_num_atoms();

    let mut sparse_atoms:Vec<usize> = vec![];
    for aa in base.iter(){
        if mover.contains(aa){
            panic!("{} is found in both base and mover.",aa);
        }
        sparse_atoms.push(*aa);
    }
    for aa in mover.iter(){
        sparse_atoms.push(*aa);
    }
    sparse_atoms.sort();

    let subatomnum:usize = sparse_atoms.len();
    let mut sub_to_orig:Vec<usize> = vec![0;subatomnum];
    let (mut sub_energyset,submapper):(PPEnergySet,Vec<i64>) = ppenergyset.make_set_for_subenv(&sparse_atoms);
    
    let mut fixflag_sub:HashSet<usize> = HashSet::new();
    let mut sparse_base:Vec<(usize,(f64,f64,f64))> = vec![];
    let mut sparse_mover:Vec<(usize,(f64,f64,f64))> = vec![];
    for ii in 0..anum{
        if submapper[ii] > -1{
            sub_to_orig[submapper[ii] as usize] = ii;
            if base.contains(&ii){
                sparse_base.push((submapper[ii] as usize,ppenergyset.get_atom(ii).get_xyz()));//後で中心座標基準に変更する
                fixflag_sub.insert(submapper[ii] as usize);
            }
            if mover.contains(&ii){
                sparse_mover.push((submapper[ii] as usize,ppenergyset.get_atom(ii).get_xyz()));//後で中心座標基準に変更する
            }
        }
    }
    //中心点を取ってその点をゼロにする。また中心点の座標を返す。
    let calc_center = |tupp:&mut Vec<(usize,(f64,f64,f64))>|->(f64,f64,f64){
        
        let mut ret:(f64,f64,f64) = (0.0,0.0,0.0);
        for ss in tupp.iter(){
            ret.0 += (ss.1).0;
            ret.1 += (ss.1).1;
            ret.2 += (ss.1).2;
        }
        ret.0 /= tupp.len() as f64;
        ret.1 /= tupp.len() as f64;
        ret.2 /= tupp.len() as f64;
        for ss in tupp.iter_mut(){
            (ss.1).0 -= ret.0;
            (ss.1).1 -= ret.1;
            (ss.1).2 -= ret.2;
        }
        return ret;
    };

    let rotate_xyz = |vv:&(f64,f64,f64),rotx:f64,roty:f64,rotz:f64| -> (f64,f64,f64){
        let mut xxx = vv.0;
        let mut yyy = vv.1;
        let mut zzz = vv.2;

        let ppy = rotx.sin()*zzz+rotx.cos()*yyy;
        let ppz = rotx.cos()*zzz-rotx.sin()*yyy;

        yyy = ppy;
        zzz = ppz;

        let ppz = roty.sin()*xxx+roty.cos()*zzz;
        let ppx = roty.cos()*xxx-roty.sin()*zzz;

        zzz = ppz;
        xxx = ppx;

        let ppy = rotz.sin()*xxx+rotz.cos()*yyy;
        let ppx = rotz.cos()*xxx-rotz.sin()*yyy;

        yyy = ppy;
        xxx = ppx;
        return (xxx,yyy,zzz); 
    };

    let base_center = calc_center(&mut sparse_base);
    
    for ss in 0..subatomnum{
        let mut xyz = sub_energyset.get_atom(ss).get_xyz();
        xyz.0 -= base_center.0;
        xyz.1 -= base_center.1;
        xyz.2 -= base_center.2;
        sub_energyset.get_atom_mut(ss).set_xyz(xyz.0,xyz.1,xyz.2);
    }
    let mut step1_topx:Vec<(f64,Vec<(f64,f64,f64)>)> = vec![];
    let mut sub_atom_energy:Vec<f64> = vec![0.0;sub_energyset.get_num_atoms()];
    let num_baseatom:usize = sparse_base.len();
    let num_movatom:usize = sparse_mover.len();
    
    
    if let Some(rr) = full_assign_resolution{
        let _mover_center = calc_center(&mut sparse_mover);//sparse_mover は中心が 0.0 に移されている
        
        let mover_centered_orig:Vec<(f64,f64,f64)> = sparse_mover.iter().map(|m|m.1.clone()).collect();
        let rot_num:usize = (PI*2.0/rr).round() as usize +1;
        
        for xx in 0..rot_num{
            let rotx:f64 = xx as f64 *rr;
            for yy in 0..rot_num{
                let roty:f64 = yy as f64 *rr;
                for zz in 0..rot_num{
                    let rotz:f64 = zz as f64 *rr;
                    let mut mover_pcent:Vec<(f64,f64,f64)> = vec![];
                    for vv in mover_centered_orig.iter(){
                        mover_pcent.push(rotate_xyz(vv,rotx,roty,rotz));
                    }
                    
                    //x 軸中心と y 軸中心だけに回転させる。y 軸中心は x が 0 とか 360 度付近だと細かすぎることになるがまあ気にしない。
                    for px in 0..rot_num{
                        let pxx:f64 = px as f64 *rr;
                        for py in 0..rot_num{
                            let pyy:f64 = py as f64 *rr;
                            let direc_:(f64,f64,f64) = rotate_xyz(&(0.0,1.0,0.0),pxx,pyy,0.0);
                            let direc:(f64,f64,f64) = process_3d::standardize(direc_.0, direc_.1, direc_.2);
                            let mut mmv:Vec<(f64,f64,f64)> = mover_pcent.clone();
                            let mut ddist:Vec<Vec<f64>> = vec![vec![100000000.0;num_movatom];num_baseatom];
                            let trynum:usize = 1000;
                            for _ in 0..trynum{
                                for mm in mmv.iter_mut(){
                                    mm.0 += direc.0;    
                                    mm.1 += direc.1;    
                                    mm.2 += direc.2;    
                                }
                                let mut updated = false;
                                for mm in 0..num_movatom{
                                    for bb in 0..num_baseatom{
                                        let dis = process_3d::distance(&sparse_base[bb].1,&mmv[mm]);
                                        if dis < ddist[mm][bb]{
                                            ddist[mm][bb] = dis;
                                            updated = true;
                                        }
                                    }
                                }
                                if !updated{
                                    //VDW が効かないところまで移動させる
                                    for mm in mmv.iter_mut(){
                                        mm.0 += direc.0*5.0;    
                                        mm.1 += direc.1*5.0;    
                                        mm.2 += direc.2*5.0;    
                                    }
                                    break;
                                }
                            }
                            for ss in 0..num_movatom{
                                let idx:usize = sparse_mover[ss].0;
                                sub_energyset.get_atom_mut(idx).set_xyz(mmv[ss].0,mmv[ss].1,mmv[ss].2);
                            }
                            
                            sub_energyset.update_distance();
                            zerofill(&mut sub_atom_energy);
                            let ene = sub_energyset.calc_energy(&mut sub_atom_energy);
                            if step1_topx.len() < step2_trial{
                                step1_topx.push((ene,mmv));
                                step1_topx.sort_by(|a,b|(a.0).partial_cmp(&b.0).unwrap());
                            }else{
                                if step1_topx.last().unwrap().0 > ene{ 
                                    step1_topx[step2_trial-1] = (ene,mmv);
                                }
                            }
                        }
                    }
                }
            }
        }

    }else{
        let mut vvec:Vec<(f64,f64,f64)> = vec![];
        for ss in sparse_mover.iter(){
            let xyz = sub_energyset.get_atom(ss.0).get_xyz();
            vvec.push(xyz);
        }
        step1_topx.push((0.0,vvec));
    }

    let num_subatom:usize = sub_energyset.get_num_atoms();
    let mut tmpatoms_prevpos:Vec<(f64,f64,f64)> = vec![(0.0,0.0,0.0);num_subatom];
    let mut atoms_mover:Vec<(usize,(f64,f64,f64))> = vec![(0,(0.0,0.0,0.0));num_movatom];
    let mut tmpatoms_checkpoint:Vec<(f64,f64,f64)> = vec![(0.0,0.0,0.0);num_subatom];

    sub_energyset.update_distance();
    zerofill(&mut sub_atom_energy);
    let mut min_score:f64 = sub_energyset.calc_energy(&mut sub_atom_energy);
    let mut min_pos:Vec<(f64,f64,f64)> = vec![(0.0,0.0,0.0);num_subatom];
    
    for (_ff,st) in step1_topx.into_iter(){
        for ss in 0..num_movatom{
            let idx:usize = sparse_mover[ss].0;
            sub_energyset.get_atom_mut(idx).set_xyz(st[ss].0,st[ss].1,st[ss].2);
            atoms_mover[ss].0 = idx;
            atoms_mover[ss].1 = st[ss];
        }
        
        for ii in 0..num_subatom{
            tmpatoms_checkpoint[ii].0 = sub_energyset.get_atom(ii).get_x();
            tmpatoms_checkpoint[ii].1 = sub_energyset.get_atom(ii).get_y();
            tmpatoms_checkpoint[ii].2 = sub_energyset.get_atom(ii).get_z();
            min_pos[ii].0 = sub_energyset.get_atom(ii).get_x();
            min_pos[ii].1 = sub_energyset.get_atom(ii).get_y();
            min_pos[ii].2 = sub_energyset.get_atom(ii).get_z();
        }
        for (_pii,para) in params.iter().enumerate(){
            let mut pline:Vec<RefinementParamLine> = para.generate_lines();
            if para.shuffle{
                pline.shuffle(rgen);
            }
            for pp in pline.iter(){
                sub_energyset.update_distance();
                zerofill(&mut sub_atom_energy);
                let mut prev_scoresum:f64 = sub_energyset.calc_energy(&mut sub_atom_energy);
                let mut first_scoresum = prev_scoresum;
                println!("{:?}",pp);
                if prev_scoresum < min_score{
                    min_score = prev_scoresum;
                    for ii in 0..num_subatom{
                        min_pos[ii].0 = sub_energyset.get_atom(ii).get_x();
                        min_pos[ii].1 = sub_energyset.get_atom(ii).get_y();
                        min_pos[ii].2 = sub_energyset.get_atom(ii).get_z();
                    }
                }

                let mut num_accepted:f64 = 1.0;
                let mut num_negative:f64 = 1.0;
                let mut current_iter:usize = 0;
                let mut prev_checked:usize = 0;
                let mut beta:f64 = 1.0;
                loop{
                    current_iter += 1;
                    for ii in 0..num_subatom{
                        tmpatoms_prevpos[ii].0 = sub_energyset.get_atom(ii).get_x();
                        tmpatoms_prevpos[ii].1 = sub_energyset.get_atom(ii).get_y();
                        tmpatoms_prevpos[ii].2 = sub_energyset.get_atom(ii).get_z();
                    }

                    for ii in 0..num_movatom{
                        (atoms_mover[ii].1).0 = sub_energyset.get_atom(atoms_mover[ii].0).get_x();
                        (atoms_mover[ii].1).1 = sub_energyset.get_atom(atoms_mover[ii].0).get_y();
                        (atoms_mover[ii].1).2 = sub_energyset.get_atom(atoms_mover[ii].0).get_z();
                    }
                
                    let center = calc_center(&mut atoms_mover);
                    let rot = pp.group_rotate.abs();
                    let rotx:f64 = rgen.gen_range(rot*-1.0,rot);
                    let roty:f64 = rgen.gen_range(rot*-1.0,rot);
                    let rotz:f64 = rgen.gen_range(rot*-1.0,rot);

                    let mov = pp.mov_dist.abs();
                    let mx:f64 = rgen.gen_range(mov*-1.0,mov);
                    let my:f64 = rgen.gen_range(mov*-1.0,mov);
                    let mz:f64 = rgen.gen_range(mov*-1.0,mov);
                    for ii in 0..num_movatom{
                        let mut xyz = rotate_xyz(&atoms_mover[ii].1,rotx,roty,rotz);
                        xyz.0 += mx+center.0;
                        xyz.1 += my+center.1;
                        xyz.2 += mz+center.2;
                        sub_energyset.get_atom_mut(atoms_mover[ii].0).set_xyz(xyz.0,xyz.1,xyz.2);
                    }

                    sub_energyset.update_distance();
                    zerofill(&mut sub_atom_energy);
                    let scoresum:f64 = sub_energyset.calc_energy(&mut sub_atom_energy);
                    if current_iter%100 == 0{
                        println!("iter:{} {} {} {}",current_iter,scoresum,prev_scoresum,min_score);
                    }
                    if scoresum < prev_scoresum{
                        //accepted
                        prev_scoresum = scoresum;
                        if prev_scoresum < min_score{
                            min_score = prev_scoresum;
                            for ii in 0..num_subatom{
                                min_pos[ii].0 = sub_energyset.get_atom(ii).get_x();
                                min_pos[ii].1 = sub_energyset.get_atom(ii).get_y();
                                min_pos[ii].2 = sub_energyset.get_atom(ii).get_z();
                            }
                        }
                    }else{
                        let dd:f64 = rgen.gen_range(0.0,1.0);
                        let bb:f64 = ((scoresum-prev_scoresum)*-1.0*beta).exp();
                        if prev_scoresum < scoresum{//同値ではない
                            num_negative += 1.0;
                            if para.accept_bound.1 > 0.0 && dd < bb{
                            //if  dd < bb{
                                num_accepted += 1.0;
                                prev_scoresum = scoresum;
                            }else{
                                //rejected
                                for ii in 0..num_subatom{
                                    sub_energyset.get_atom_mut(ii).set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                                }
                                
                                //sub_energyset.update_distance();
                                //println!("*****{} {}",prev_scoresum,sub_energyset.calc_energy(&mut sub_atom_energy));
                            }
                        }
                    }
                    let checkflag =  current_iter%para.steps_bound_check == 0 || current_iter == pp.num_iter;
                    if checkflag{
                        if num_negative > 0.0 && num_accepted/num_negative > para.accept_bound.1{
                            beta *= 5.0;
                            current_iter = prev_checked;
                            for ii in 0..num_subatom{
                                sub_energyset.get_atom_mut(ii).set_xyz(tmpatoms_checkpoint[ii].0,tmpatoms_checkpoint[ii].1,tmpatoms_checkpoint[ii].2);
                            }
                        }else{
                            if num_negative > 0.0 && num_accepted/num_negative < para.accept_bound.0{
                                beta /= 5.0;
                            }
                            prev_checked = current_iter;
                            for ii in 0..num_subatom{
                                tmpatoms_checkpoint[ii].0 = sub_energyset.get_atom(ii).get_x();
                                tmpatoms_checkpoint[ii].1 = sub_energyset.get_atom(ii).get_y();
                                tmpatoms_checkpoint[ii].2 = sub_energyset.get_atom(ii).get_z();
                            }
                        }
                        num_negative = 0.0;
                        num_accepted = 0.0;
                    }
                    if current_iter == pp.num_iter{
                        let ediff = first_scoresum-min_score;
                        if pp.refine_loop_threshold <= 0.0 || pp.refine_loop_threshold > ediff{
                            break;
                        }
                        first_scoresum = min_score;
                        current_iter = 0;
                        prev_checked = 0;
                        for ii in 0..num_subatom{
                            tmpatoms_checkpoint[ii].0 = sub_energyset.get_atom(ii).get_x();
                            tmpatoms_checkpoint[ii].1 = sub_energyset.get_atom(ii).get_y();
                            tmpatoms_checkpoint[ii].2 = sub_energyset.get_atom(ii).get_z();
                        }
                    }
                }
            }
        }
    }
    
    for ii in 0..num_subatom{
        sub_energyset.get_atom_mut(ii).set_xyz(min_pos[ii].0+base_center.0,min_pos[ii].1+base_center.1,min_pos[ii].2+base_center.2);
    }
    charmm_based_energy::CharmmEnv::accept_subenv_atom(&mut ppenergyset.evoef2_env.md_envset
        ,&sub_energyset.evoef2_env.md_envset,&submapper);
}

pub struct RebuildParam{
    pub rebuild_iter:usize,
    pub rebuild_region:usize,
    pub rebuild_num:usize,
}


pub fn refinement(pdbb:pdbdata::PDBEntry
    //querychain, templatefile,templatechain
    ,templatefiles:Option<Vec<(String,String,String)>>
    ,flagfile:Option<String>
    ,groupfile:Option<String>
    ,_num_incremental_build:usize//todo: remove
    ,topfile:&str
    ,paramfile:&str
    ,evoefdir:&str
    ,backbone_dihedral_angle_file:&str
    ,cb_distance_file:&str
    ,residue_contact_file:&str
    ,params:&Vec<RefinementParam>
    ,missing_params:Option<MissingBuilderParam>
    ,rebuild_param_:Option<RebuildParam>
    ,random_seed:Option<u64>
    ,checkpoint:usize
    ,outfile:&str
    )->PPEnergySet{
        
    let mut rgen:StdRng;//乱数生成器
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }

    let rebuild_param = rebuild_param_.unwrap_or(RebuildParam{rebuild_iter:0,rebuild_num:0,rebuild_region:0});

    let parr = charmm_param::CHARMMParam::load_chamm19(topfile,paramfile);
    //QUERY CHAIN NAME, TEMPLATE N CA C ATOM COORDINATES
    //Template の構造を読み込む
    let mut templates:Vec<(String,Vec<(Vec<f64>,Vec<f64>,Vec<f64>)>)> = vec![];
    if let Some(x) = templatefiles{
        for xx in x{
            let pdbb:pdbdata::PDBEntry = mmcif_process::load_pdb(xx.1.as_str());
            for mut cc in pdbb.get_model_at(0).get_entity_at(0).iter_mut_asyms(){
                cc.remove_alt(None);
                if cc.chain_name == xx.2{
                    let mut vvec:Vec<(Vec<f64>,Vec<f64>,Vec<f64>)> = vec![];
                    for rr in cc.iter_mut_comps(){
                        let n = rr.get_N();
                        let ca = rr.get_CA();
                        let c = rr.get_C();
                        if let Some(nn) = n{
                            if let Some(caa) = ca{
                                if let Some(cc) = c{
                                    vvec.push((
                                        vec![ nn.get_x(), nn.get_y(), nn.get_z()],
                                        vec![caa.get_x(),caa.get_y(),caa.get_z()],
                                        vec![ cc.get_x(), cc.get_y(), cc.get_z()]
                                        )
                                    );
                                }
                            }
                        }
                    }
                    templates.push((xx.0.to_string(),vvec));
                }
            }       
        }
    }
    
    let (mut md_envset,md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&pdbb.get_all_asyms().iter().map(|m|m).collect(),&parr,true);
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(
        &backbone_dihedral_angle_file);
    let (torsion,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,false,true),false);

    let mut distbe:Vec<distance_energy::AtomBinnedDistanceEnergy> = vec![];
    if cb_distance_file.len() > 0{
        let dist_bvvec:Vec<(String,usize,String,usize,distance_energy::AtomBinnedDistanceEnergy)> = distance_energy::AtomBinnedDistanceEnergy::load_name_mapped(cb_distance_file);
        distbe = distance_energy::AtomBinnedDistanceEnergy::assign_atom_ids(&md_envset,dist_bvvec);
    }


    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    md_envset.update_distance();
    let scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
    println!("energy:{:?}",scoresum);
   
    //dihed が重複して計算されているのでウェイトをつける。今は適当に 0.5
    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&md_varset.dihedvec,&torsion,&omegas_general);
    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.5;
    }

    //対象の全原子数
    let num_atoms:usize = md_envset.atoms.len();
    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,evoefdir,false),
        backbone_energy_omega:omegas_general,
        backbone_energy_phi_psi:torsion,
        atom_distance_energy:vec![],
        atom_binned_distance_energy:distbe,
        atom_contact_energy:vec![],
        weights:PPEnergyWeight::new()
    };
    ppenergyset.weights.dihed_weight_charmm = weight_dihed;
    ppenergyset.weights.charmm_electro = 0.0;
    ppenergyset.weights.charmm_lj = 1.0;
    
    ppenergyset.weights.backbone_energy_phi_psi = 1.0;
    ppenergyset.weights.backbone_energy_omega = 1.0;

    let mut fixflag:HashSet<usize> = HashSet::new();



    
    //chainname,residueindex in chain
    let mut residue_fixed:HashMap<String,HashSet<i64>> = HashMap::new();
    let mut residue_missing:HashMap<String,HashSet<usize>> = HashMap::new();
    if let Some(x) = flagfile{
        let (missing,freezed) = load_missing_and_freezed_atoms(&x,&ppenergyset.evoef2_env.md_envset.atoms);
        for ii in 0..num_atoms{
            if !residue_fixed.contains_key(&ppenergyset.get_atom(ii).chain_name){
                residue_fixed.insert(ppenergyset.get_atom(ii).chain_name.clone(),HashSet::new());
            }
            if !residue_missing.contains_key(&ppenergyset.get_atom(ii).chain_name){
                residue_missing.insert(ppenergyset.get_atom(ii).chain_name.clone(),HashSet::new());
            }
            if missing.contains(&ii){
                if ppenergyset.get_atom(ii).residue_index_in_chain > -1{
                    residue_missing.get_mut(&ppenergyset.get_atom(ii).chain_name).unwrap().insert(ppenergyset.get_atom(ii).residue_index_in_chain as usize);
                }
            }
            if freezed.contains(&ii){
                fixflag.insert(ii);
                residue_fixed.get_mut(&ppenergyset.get_atom(ii).chain_name).unwrap().insert(ppenergyset.get_atom(ii).residue_index_in_chain);
            }
        }
    };
    

    //fix された Residue が無い場合は作れないので一個をランダムで入れる
    //今使ってないが必要になるかも？？？
    //for (_kk,vv) in residue_fixed.iter_mut(){
    //    if vv.len() == 0{
    //       vv.insert(rgen.gen_range(residue_min.get(_kk).unwrap(),residue_max.get(_kk).unwrap()+1));
    //    }
    //}
    

    let sorted_atoms:HashMap<String,Vec<Vec<usize>>> = pp_energy_mc::make_sorted_residue_group(&ppenergyset.evoef2_env.md_envset.atoms);
    

    let groupatoms:Vec<Vec<usize>> = if let Some(x) = groupfile.as_ref(){
        load_group_atoms(&x,&ppenergyset.evoef2_env.md_envset.atoms)
    }else{
        vec![]
    };

    let mut residue_label_to_indexinchain:HashMap<String,(String,i64)> = HashMap::new();
    for aa in ppenergyset.evoef2_env.md_envset.atoms.iter(){
        let label:String = aa.chain_name.clone()+"#"+&(aa.residue_number.to_string())+"#"+&aa.residue_ins_code;
        if !residue_label_to_indexinchain.contains_key(&label){
            //println!("{:?}",(aa.chain_name.clone(),aa.residue_index_in_chain));
            residue_label_to_indexinchain.insert(label,(aa.chain_name.clone(),aa.residue_index_in_chain));
        }else{
            if residue_label_to_indexinchain.get(&label).unwrap() != &(aa.chain_name.clone(),aa.residue_index_in_chain){
                eprintln!("{} is assigned to {:?} but {:?} was going to assin. Some function might not run successfully.",
                label,
                residue_label_to_indexinchain.get(&label).unwrap(),
                (aa.chain_name.clone(),aa.residue_index_in_chain),
                );
            }
        }
    }

    //ある原子の distance ではなく Residue に存在する全原子のうちのいずれかが 6.0 angstrom 内である場合に Energy を加える。よってマイナスである必要がある
    //6.0 は evoef2 の attr の cutoff distance
    //chainname+"#"+residue_number#insertion_code
    let mut residue_contact_:Vec<(String,String,f64)> = if residue_contact_file.len() > 0{
        load_contact_file(residue_contact_file)
    }else{
        vec![]
    };
    residue_contact_.sort_by(|a,b| a.2.partial_cmp(&b.2).unwrap_or_else(||panic!("Error in sort contact.")));
    
    let mut residue_contact:Vec<(String,i64,String,i64,f64)> = vec![];
    for (r1_,r2_,val) in residue_contact_.iter(){
        if !residue_label_to_indexinchain.contains_key(r1_.as_str()){
            println!("{} is not found in pdb file.",r1_);
            continue;
        }
        if !residue_label_to_indexinchain.contains_key(r2_.as_str()){
            println!("{} is not found in pdb file.",r2_);
            continue;
        }
        let r1:&(String,i64) = residue_label_to_indexinchain.get(r1_.as_str()).unwrap();
        let r2:&(String,i64) = residue_label_to_indexinchain.get(r2_.as_str()).unwrap();
        residue_contact.push((r1.0.clone(),r1.1,r2.0.clone(),r2.1,*val));
    }
    for rr in residue_contact.iter(){
        let (c1,r1,c2,r2,v) = rr;
        let mut atoms_a:Vec<usize> = vec![];
        for ss in sorted_atoms.get(c1).unwrap_or_else(||panic!("error processing contact.")).iter(){
            if ppenergyset.get_atom(ss[0]).residue_index_in_chain == *r1{
                atoms_a.append(&mut ss.clone());
            }
        }
        
        let mut atoms_b:Vec<usize> = vec![];
        for ss in sorted_atoms.get(c2).unwrap_or_else(||panic!("error processing contact.")).iter(){
            if ppenergyset.get_atom(ss[0]).residue_index_in_chain == *r2{
                atoms_b.append(&mut ss.clone());
            }
        }

        if atoms_a.len() > 0 && atoms_b.len() > 0{
            ppenergyset.atom_contact_energy.push(pp_energy::AtomContactEnergy{
                atoms_a:atoms_a,
                atoms_b:atoms_b,
                threshold:6.0,//ここ CB の時は 8.0 にする
                linear_penalty:0.1,
                energy:*v
            });
        }else{
            eprintln!("???contact {:?} was not found. a {} b {}", (c1,r1,c2,r2,v) ,atoms_a.len(),atoms_b.len());
        }
    }


    if residue_missing.len() > 0{
        //Missing 領域の作成
        if let Some(p1) = missing_params.as_ref(){
        build_missing(
            &mut ppenergyset
            , &sorted_atoms
            , &residue_missing
            , &fixflag
            , p1.fix_omega
            , &mut rgen
            , 1.0
            , p1.phipsi_top_x
            , &p1.param1
            , &p1.param2
            ,p1.num_structures_step1
            ,p1.num_structures_step2);
        }
    }


    

    let mut atom_based_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut check_count:usize = 0;
    for rebuilder in 0..(rebuild_param.rebuild_iter+1){
        if rebuilder > 0 && rebuild_param.rebuild_region > 0{
            if let Some(p1) = missing_params.as_ref(){
                /*
                ppenergyset.update_distance();
                for aa in atom_based_energies.iter_mut(){
                    *aa = 0.0;
                }
                ppenergyset.calc_energy(&mut atom_based_energies);
                let mut res_energies:Vec<(String,usize,f64)> = vec![];
                for (kk,vv) in sorted_atoms.iter(){
                    for (vii,vvv) in vv.iter().enumerate(){
                        let mut scc:f64 = 0.0;
                        for aa in vvv.iter(){
                            scc += atom_based_energies[*aa];
                        }
                        res_energies.push((kk.clone(),vii,scc));
                    }
                }
                res_energies.sort_by(|a,b|(a.2).partial_cmp(&b.2).unwrap());
                */
                
                for (kk,vv) in sorted_atoms.iter(){
                    let rnum:usize = vv.len();
                    let mut partial:Vec<(f64,Vec<usize>)> = vec![];
                    
                    ppenergyset.update_distance();
                    zerofill(&mut atom_based_energies);
                    let _scoree = ppenergyset.calc_energy(&mut atom_based_energies);
                    let plen:usize = rebuild_param.rebuild_region;//構築を試みる領域の長さ
                    let num_region:usize = rebuild_param.rebuild_num;//構築を試みる領域の数
                    let loopnum = rnum/plen+1;
                    for _ in 0..loopnum{
                        let s:usize = rgen.gen_range(0,rnum);
                        let mut pp:Vec<usize> = vec![s];
                        for o in 1..(plen/2){
                            let po:i64 = s as i64 - o as i64;
                            if po > -1{
                                pp.push(po as usize);
                            }
                            let po:i64 = s as i64 + o as i64;
                            if po < rnum as i64{
                                pp.push(po as usize);
                            }    
                        }
                        pp.sort();
                        let mut pscore:f64 = 0.0;
                        for ppp in pp.iter(){
                            for aaa in vv[*ppp].iter(){
                                pscore += atom_based_energies[*aaa];
                            }
                        }
                        partial.push((pscore,pp));
                    }
                    
                    partial.sort_by(|a,b|(a.0).partial_cmp(&b.0).unwrap());
                    partial.reverse();

                    let mut prevpos:Vec<(f64,f64,f64)> = vec![(0.0,0.0,0.0);num_atoms];
                    for aa in 0..num_atoms{
                        prevpos[aa] = ppenergyset.get_atom(aa).get_xyz();
                    }
                    let lj_weight = ppenergyset.weights.charmm_lj;
                    ppenergyset.weights.charmm_lj = 0.0;//lj は後のステップで最適化する予定なのでこの段階で考えない。
                    let prev_score = ppenergyset.calc_energy(&mut atom_based_energies);
                    ppenergyset.weights.charmm_lj = lj_weight;//構築するときは LJ も考える
                    for (pii,(_,pt)) in partial.into_iter().enumerate(){
                        if num_region <= pii{
                            break;
                        }
                        if pt[0] != 0 && *pt.last().unwrap() != rnum-1{
                            let ps:usize = pt[0];
                            let pe:usize = *pt.last().unwrap();
                            let mut base:HashSet<usize> = HashSet::new();
                            for ss in 0..ps{
                                for aa in vv[ss].iter(){
                                    base.insert(*aa);
                                }
                            }
                            let mut mover:HashSet<usize> = HashSet::new();
                            for ss in (pe+1)..rnum{
                                for aa in vv[ss].iter(){
                                    mover.insert(*aa);
                                }
                            }
                            if base.len() == 0 || mover.len() == 0{
                                eprintln!("??? must not be zero.");//ないとは思う
                            }else{
                                zatsuna_docking(&mut ppenergyset, &base,&mover,None, &mut rgen, &params,10);//ToDo 値指定できるように
                            }
                            
                        }
                        let partial_missing:HashSet<usize> = pt.into_iter().collect();
                        
                        let mut mi:HashMap<String,HashSet<usize>> = HashMap::new();
                        mi.insert(kk.clone(),partial_missing);
                        ppenergyset.weights.charmm_lj = lj_weight;
                        build_missing(
                        &mut ppenergyset
                        , &sorted_atoms
                        , &mi
                        , &fixflag
                        , p1.fix_omega
                        , &mut rgen
                        , 1.0
                        , p1.phipsi_top_x
                        , &p1.param1
                        , &p1.param2
                        ,p1.num_structures_step1
                        ,p1.num_structures_step2);
                        ppenergyset.update_distance();
                        ppenergyset.weights.charmm_lj = 0.0;
                        let pscore = ppenergyset.calc_energy(&mut atom_based_energies);
                        if pscore > prev_score{
                            for aa in 0..num_atoms{
                                ppenergyset.get_atom_mut(aa).set_xyz(prevpos[aa].0,prevpos[aa].1,prevpos[aa].2);
                            }
                        }
                        ppenergyset.weights.charmm_lj = lj_weight;
                    }
                }
            }
        }
        for (_pii,para) in params.iter().enumerate(){
            let mut pline:Vec<RefinementParamLine> = para.generate_lines();
            if para.shuffle{
                pline.shuffle(&mut rgen);
            }
            
            for prr in pline.into_iter(){
                println!("{:?}",prr);
                let bond_restrictor:Vec<charmm_based_energy::BondVars> = ppenergyset.evoef2_env.charmm_vars.bondvec.clone();
                
                let rss:u64 = rgen.gen_range(0,std::u64::MAX);
                ppenergyset.update_edges(prr.inv_edge.max(4)+1);
                let res_in_chain:HashMap<String,Vec<Vec<usize>>> = pp_energy_mc::make_sorted_residue_group(&ppenergyset.evoef2_env.md_envset.atoms);

                if templates.len() > 0 && para.use_template{
                    for tt in templates.iter(){
                        if !res_in_chain.contains_key(&tt.0){
                            eprintln!("{} was not found in subchains.",tt.0);
                        }else{
                            println!("Trying template assist...");
                            fit_to_chain_template(&mut ppenergyset, res_in_chain.get(&tt.0).unwrap(),&tt.1);
                        }
                    }
                }else if templates.len() > 0 || para.use_template{
                    println!("template assist is off:\ntemplate num:{} flag:{} ",templates.len(),para.use_template);
                }

                ppenergyset.update_distance();

                loop{
                    check_count += 1;

                    for aa in atom_based_energies.iter_mut(){
                        *aa = 0.0;
                    }
                    
                    let prev_energy:f64 = ppenergyset.calc_energy(&mut atom_based_energies);
                    pp_energy_mc::mc_iter_array(&mut ppenergyset
                        ,prr.num_iter
                        ,&fixflag
                        ,Some(rss)
                        ,prr.num_target_atoms
                        ,para.sorter
                        ,prr.mov_dist
                        ,(para.accept_bound.0,para.accept_bound.1,para.steps_bound_check)
                        ,&bond_restrictor
                        ,&vec![]
                        ,(1.1,prr.dec_factor)
                        ,(prr.inv_dist,prr.inv_edge)
                        ,prr.refine_retry
                        ,prr.lbfgs
                        ,if prr.group_rotate > 0.0{Some(( //ToDo 仮置きなので変更
                            &groupatoms
                            ,prr.group_rotate))}else{None}
                    );
                    
                    for aa in atom_based_energies.iter_mut(){
                        *aa = 0.0;
                    }
                    let pres = ppenergyset.calc_energy(&mut atom_based_energies);
                    let ediff:f64 = prev_energy-pres;
                    println!("improve: {}",ediff);
                    if checkpoint > 0 && check_count%checkpoint == 0{
                        let mut lines:Vec<String> = vec![];
                        for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
                            let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
                            att.temp_factor = atom_based_energies[aii].max(0.0).min(999.0);
                            lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                        }
                        write_to_file(format!("{}.{}.pdb",outfile,check_count).as_str(),lines);
                    }
                    if prr.refine_loop_threshold <= 0.0 || prr.refine_loop_threshold > ediff{
                        break;
                    }
                }
            }
        }
    }
    return ppenergyset;
}

#[derive(Clone)]
struct BackboneStat{
    chain_name:String,
    index_in_chain:usize,
    self_n_ca_c_cb:(usize,usize,usize,i64),
    prev_n_ca_c:(i64,i64,i64),
    next_n_ca_c:(i64,i64,i64),
    self_dihed:i64,
    prev_dihed:i64,
    next_dihed:i64,
    self_pomega:i64,
    self_nomega:i64,
    prev_placed:bool,
    next_placed:bool,
    dead:bool
}




//missing_residues は、sorted_atoms.get(ss)[ii][jj] の ii を指定する
pub fn build_missing(ppenergyset:&mut PPEnergySet
    ,sorted_atoms:&HashMap<String,Vec<Vec<usize>>>
    ,missing_residues:&HashMap<String,HashSet<usize>>
    ,fixed_atoms:&HashSet<usize>
    ,fix_omega:bool
    ,rgen:&mut StdRng
    ,sigma:f64//phi-psi omega を選ぶ際の正規分布のパラメータ
    ,restrict_top:usize //phi-psi omega は正規分布を使用するが TopX のものしかとらない
    ,params:&Vec<RefinementParam>//Backbone と CB だけで作成した仮構造を最適化するためのパラメータ
    ,params2:&Vec<RefinementParam>//全体最適化するためのパラメータ
    ,step1_trial:usize//Backbone と CB だけの仮構造を何個作るか
    ,step2_trial:usize//作成した仮構造のうち何個を全体最適化するか
    ){
    
    let normal_distribution:Normal<f64> = Normal::new(0.0,sigma).unwrap();
    let anum:usize = ppenergyset.evoef2_env.md_envset.atoms.len();
    let orig_pos:Vec<(f64,f64,f64)> = ppenergyset.evoef2_env.md_envset.atoms.iter().map(|m|m.get_xyz()).collect();
    let mut n_ca_c_cb:HashMap<String,Vec<(i64,i64,i64,i64)>> = HashMap::new();
    let mut residue_placed:HashMap<String,HashSet<usize>> = HashMap::new();

    //前の残基への n_ca_c への index, n_ca_c への index,次の残基への n_ca_c の index, prevc_phipsi, nextn_phipsi への index
    let mut missing_backbone_stat:HashMap<String,Vec<BackboneStat>> = HashMap::new();

    //atomidx->dihed
    let mut ca_phipsi:HashMap<usize,usize> = HashMap::new();
    let mut startca_omega:HashMap<usize,usize> = HashMap::new();
    let mut endca_omega:HashMap<usize,usize> = HashMap::new();
    
    let mut atoms_notmissing:HashSet<usize> = HashSet::new();


    //sorted residues から n ca c cb のインデクスを取り出す
    let get_n_ca_c_cb = |pp:&PPEnergySet,idxs:&Vec<usize>|->(i64,i64,i64,i64){
        let mut n:i64 = -1;
        let mut ca:i64 = -1;
        let mut c:i64 = -1;
        let mut cb:i64 = -1;
        for aa in idxs.iter(){
            if pp.get_atom(*aa).atom_name == "CA"{
                ca = *aa as i64;
            }else if pp.get_atom(*aa).atom_name == "C"{
                c = *aa as i64;
            }else if pp.get_atom(*aa).atom_name == "N"{
                n = *aa as i64;
            }else if pp.get_atom(*aa).atom_name == "CB"{
                cb = *aa as i64;
            }
        }
        return (n,ca,c,cb);
    };

    //重複は見てない
    for pp in ppenergyset.backbone_energy_phi_psi.iter().enumerate(){
        ca_phipsi.insert(pp.1.atoms.2,pp.0);
    }
    for pp in ppenergyset.backbone_energy_omega.iter().enumerate(){
        startca_omega.insert(pp.1.atoms.0,pp.0);
        endca_omega.insert(pp.1.atoms.3,pp.0);
    }

    for (k,v) in sorted_atoms.iter(){
        let mut vvec:Vec<(i64,i64,i64,i64)> = vec![];//sorted atoms が全部入っている必要がある
        let mut vvec_missing:Vec<(i64,usize,i64)> = vec![];
        let mut placed:Vec<bool> = vec![true;v.len()];
        let mut bvec:Vec<BackboneStat> = vec![];
        for (vii,vv) in v.iter().enumerate(){
            vvec.push(get_n_ca_c_cb(ppenergyset,vv));
            if !missing_residues.contains_key(k){
                for aa in vv.iter(){
                    atoms_notmissing.insert(*aa);
                }
                continue;
            }
            if missing_residues.get(k).unwrap().contains(&vii){
                placed[vii] = false;
                let mut nn:i64 =-1;
                let mut cc:i64 =-1;
                if vii > 0{
                    nn = vii as i64 -1;
                }
                if vii < v.len()-1{
                    cc = vii as i64 +1;
                }
                vvec_missing.push((nn,vii,cc));
            }else{
                for aa in vv.iter(){
                    atoms_notmissing.insert(*aa);
                }
            }
        }

        for vv in vvec_missing.iter(){
            if vvec[vv.1].0 < 0 || vvec[vv.1].1 < 0 || vvec[vv.1].2 < 0{
                continue;
            }

            let mut bs = BackboneStat{
                chain_name:k.clone(),
                index_in_chain:vv.1,
                self_n_ca_c_cb:(vvec[vv.1].0 as usize,vvec[vv.1].1 as usize,vvec[vv.1].2 as usize,vvec[vv.1].3),
                prev_n_ca_c:(-1,-1,-1),
                next_n_ca_c:(-1,-1,-1),
                self_dihed :-1,
                prev_dihed :-1,
                next_dihed :-1,
                self_pomega:-1,
                self_nomega:-1,
                prev_placed:false,
                next_placed:false,
                dead:false,
            };
            let v1u:usize = vvec[vv.1].1 as usize;
            if ca_phipsi.contains_key(&v1u){
                bs.self_dihed = *ca_phipsi.get(&v1u).unwrap() as i64;
            }
            if startca_omega.contains_key(&v1u){
                bs.self_nomega = *startca_omega.get(&v1u).unwrap() as i64;
            }
            if endca_omega.contains_key(&v1u){
                bs.self_pomega = *endca_omega.get(&v1u).unwrap() as i64;
            }

            if vv.0 > -1{
                let vu = vv.0 as usize;  
                bs.prev_n_ca_c = (vvec[vu].0, vvec[vu].1, vvec[vu].2);
                
                if vvec[vu].1 > -1{
                    let vvu:usize = vvec[vu].1 as usize;
                    if ca_phipsi.contains_key(&vvu){
                        bs.prev_dihed = *ca_phipsi.get(&vvu).unwrap() as i64;
                    }
                }
            }
            
            if vv.2 > -1{
                let vu = vv.2 as usize;  
                bs.next_n_ca_c = (vvec[vu].0, vvec[vu].1, vvec[vu].2);
                
                if vvec[vu].1 > -1{
                    let vvu:usize = vvec[vu].1 as usize;
                    if ca_phipsi.contains_key(&vvu){
                        bs.next_dihed = *ca_phipsi.get(&vvu).unwrap() as i64;
                    }
                }
            }
            bvec.push(bs);
        }
        let mut placed_hs:HashSet<usize> = HashSet::new();
        for pp in placed.into_iter().enumerate(){
            if pp.1 {
                placed_hs.insert(pp.0);
            }
        }
        n_ca_c_cb.insert(k.clone(),vvec);
        missing_backbone_stat.insert(k.clone(),bvec);
        residue_placed.insert(k.clone(),placed_hs);
    }
    let mut next_residues_array:Vec<Vec<BackboneStat>> = vec![];
    loop {
        let mut next_residues:Vec<BackboneStat> = vec![];
        for (kk,vv) in missing_backbone_stat.iter_mut(){
            let mut remindex:HashSet<usize> = HashSet::new();
            let clen = n_ca_c_cb.get(kk).unwrap().len();
            for (iii,vvv) in vv.iter_mut().enumerate(){
                if vvv.index_in_chain > 0{
                    if residue_placed.get(kk).unwrap().contains(&(vvv.index_in_chain-1)){
                        vvv.prev_placed = true;
                        remindex.insert(iii);
                    }
                }
                if vvv.index_in_chain < clen-1{
                    if residue_placed.get(kk).unwrap().contains(&(vvv.index_in_chain+1)){
                        vvv.next_placed = true;
                        remindex.insert(iii);
                    }
                }
            }
            for rr in remindex.iter(){
                next_residues.push(vv[*rr].clone());
                vv[*rr].dead = true;
            }
        }
        
        for (_kk,vv) in missing_backbone_stat.iter_mut(){
            loop {
                let vlen = vv.len();
                let mut updated = false;
                for vi in 0..vlen{
                    if vv[vi].dead{
                        vv.swap_remove(vi);
                        updated = true;
                        break;
                    }
                }
                if !updated{
                    break;
                }
            }
        }
        //何も Residue が置かれない場合
        if next_residues.len() == 0{
            //冗長だが最大でも chain 数しか呼び出されないはずなので気にしない
            let mut all_b:Vec<(String,usize)> = vec![];
            for (kk,vv) in missing_backbone_stat.iter_mut(){
                if vv.len() > 0{
                    all_b.push((kk.clone(),vv.len()));
                }
            }
            if all_b.len() == 0{
                eprintln!(";;;???");
                break;
            }
            all_b.shuffle(rgen);
            let rp:usize = rgen.gen_range(0,all_b[0].1);
            next_residues.push(missing_backbone_stat.get(&all_b[0].0).unwrap()[rp].clone());
            if missing_backbone_stat.get_mut(&all_b[0].0).unwrap()[rp].dead{
                panic!("*****????");
            }
            missing_backbone_stat.get_mut(&all_b[0].0).unwrap()[rp].dead = true;
        }
        if next_residues.len() == 0{
            eprintln!("???");
            continue;
        }
        //使用した Backbone の削除
        for (_kk,vv) in missing_backbone_stat.iter_mut(){
            loop {
                let vlen = vv.len();
                let mut updated = false;
                for vi in 0..vlen{
                    if vv[vi].dead{
                        vv.swap_remove(vi);
                        updated = true;
                        break;
                    }
                }
                if !updated{
                    break;
                }
            }
        }
        for bb in next_residues.iter(){
            residue_placed.get_mut(&bb.chain_name).unwrap().insert(bb.index_in_chain);
        }
        next_residues_array.push(next_residues);
        let mut remained:bool = false;
        for (_kk,vv) in missing_backbone_stat.iter_mut(){
            let vlen = vv.len();
            if vlen > 0{
                remained = true;
            }
        }
        if !remained{
            break;
        }
    }
    let mut backbone_tops:Vec<(f64,Vec<(usize,(f64,f64,f64))>)> = vec![];
    let mut backbone_originalpos:Vec<(usize,(f64,f64,f64))> = vec![];
    for ii in 0..anum{
        let aa = ppenergyset.get_atom(ii);
        if aa.atom_name == "CA"
        || aa.atom_name == "C" 
        || aa.atom_name == "CB" 
        || aa.atom_name == "N"{
            backbone_originalpos.push((ii,aa.get_xyz()));
        }
    }
    for _ in 0..step1_trial{
        let stepnum:usize = next_residues_array.len();
        for i in 0..stepnum{
            let znum:usize = next_residues_array[i].len();
            for ii in 0..znum{
                let bb:&BackboneStat = &next_residues_array[i][ii];
                let mut placen:bool = bb.prev_placed;
                let mut fix:bool = false;
                if bb.prev_placed && bb.next_placed{
                    placen = rgen.gen_range(-1.0,1.0) < 0.0;
                }
                if !bb.prev_placed && !bb.next_placed{
                    fix = true;
                }
                if fix{
                    //bond の長さのチェックはするか
                    continue;
                }
                if placen {
                    //http://www.cryst.bbk.ac.uk/PPS95/course/9_quaternary/3_geometry/torsion.html
                    let mut phipsi_self:(f64,f64) = (-50.0/180.0*PI,-60.0/180.0*PI);
                    if bb.self_dihed > -1{
                        phipsi_self = ppenergyset.backbone_energy_phi_psi[bb.self_dihed as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                    }
                    let mut phipsi_prev:(f64,f64) = (PI,PI);
                    if bb.prev_dihed > -1{
                        phipsi_prev = ppenergyset.backbone_energy_phi_psi[bb.prev_dihed as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                    }
                    let mut omega:f64 = PI;
                    if !fix_omega{
                        if bb.self_pomega > -1{
                            omega = ppenergyset.backbone_energy_omega[bb.self_pomega as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                        }
                    }
                    let pos:Vec<Point3D> = build_peptide(ppenergyset.get_atom(bb.prev_n_ca_c.0 as usize)
                    ,ppenergyset.get_atom(bb.prev_n_ca_c.1 as usize)
                    ,ppenergyset.get_atom(bb.prev_n_ca_c.2 as usize)
                    ,phipsi_prev.1,omega,phipsi_self.0);
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.0).set_xyz(pos[0].get_x(),pos[0].get_y(),pos[0].get_z());
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.1).set_xyz(pos[1].get_x(),pos[1].get_y(),pos[1].get_z());
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.2).set_xyz(pos[2].get_x(),pos[2].get_y(),pos[2].get_z());
                    if bb.self_n_ca_c_cb.3 > -1{
                        let cb:Point3D = get_cb_pos(&pos[0],&pos[1],&pos[2]);
                        ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.3 as usize).set_xyz(cb.get_x(),cb.get_y(),cb.get_z());
                    }
                }else{
                    let mut phipsi_self:(f64,f64) = (-50.0/180.0*PI,-60.0/180.0*PI);
                    if bb.self_dihed > -1{
                        phipsi_self = ppenergyset.backbone_energy_phi_psi[bb.self_dihed as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                    }
                    let mut phipsi_next:(f64,f64) = (PI,PI);
                    if bb.next_dihed > -1{
                        phipsi_next = ppenergyset.backbone_energy_phi_psi[bb.next_dihed as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                    }
                    let mut omega:f64 = PI;
                    if !fix_omega{
                        if bb.self_nomega > -1{
                            omega = ppenergyset.backbone_energy_omega[bb.self_nomega as usize].get_random_angle_radian(&normal_distribution, rgen, restrict_top);
                        }
                    }
                    let pos:Vec<Point3D> = build_peptide_backward(ppenergyset.get_atom(bb.next_n_ca_c.0 as usize)
                    ,ppenergyset.get_atom(bb.next_n_ca_c.1 as usize)
                    ,ppenergyset.get_atom(bb.next_n_ca_c.2 as usize)
                    ,phipsi_self.1,omega,phipsi_next.0);
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.0).set_xyz(pos[0].get_x(),pos[0].get_y(),pos[0].get_z());
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.1).set_xyz(pos[1].get_x(),pos[1].get_y(),pos[1].get_z());
                    ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.2).set_xyz(pos[2].get_x(),pos[2].get_y(),pos[2].get_z());
                    if bb.self_n_ca_c_cb.3 > -1{
                        let cb:Point3D = get_cb_pos(&pos[0],&pos[1],&pos[2]);
                        ppenergyset.get_atom_mut(bb.self_n_ca_c_cb.3 as usize).set_xyz(cb.get_x(),cb.get_y(),cb.get_z());
                    }
                }
            }
        }
        let mut sparse_atoms:Vec<usize> = vec![];
        for ii in 0..anum{
            //if atoms_notmissing.contains(&ii){
            //    sparse_atoms.push(ii);
            //}else{
                let aa = ppenergyset.get_atom(ii);
                if aa.atom_name == "CA"
                || aa.atom_name == "C" 
                || aa.atom_name == "CB" 
                || aa.atom_name == "N"{
                    sparse_atoms.push(ii);
                }
            //}
        }
        let subatomnum:usize = sparse_atoms.len();
        let (mut sub_energyset,submapper):(PPEnergySet,Vec<i64>) = ppenergyset.make_set_for_subenv(&sparse_atoms);
        let mut fixflag_sub:HashSet<usize> = HashSet::new();
        for ii in 0..anum{
            if submapper[ii] > -1{
                if fixed_atoms.contains(&ii){
                    fixflag_sub.insert(submapper[ii] as usize);
                }
            }
        }
        for (_pii,para) in params.iter().enumerate(){
            let pline:Vec<RefinementParamLine> = para.generate_lines();
            for (_prrii,prr) in pline.iter().enumerate(){
                println!("{:?}",prr);
                sub_energyset.weights.charmm_angle = 1.0;
                sub_energyset.weights.charmm_electro = 0.0;
                sub_energyset.weights.charmm_lj = 1.0;
                sub_energyset.weights.evoef2_vdw =  vec![0.0;4];
                sub_energyset.weights.evoef2_elec = 0.0;
                sub_energyset.weights.evoef2_hb = vec![0.0;9];
                sub_energyset.weights.evoef2_desolv_polar =  0.0;
                sub_energyset.weights.evoef2_desolv_nonpolar = 0.0;
                sub_energyset.weights.evoef2_ss = 0.0;
                sub_energyset.weights.backbone_energy_phi_psi = 1.0;
                sub_energyset.weights.backbone_energy_omega = 1.0;
                let sub_bond_restrictor:Vec<charmm_based_energy::BondVars> = sub_energyset.evoef2_env.charmm_vars.bondvec.clone();
                let rss:u64 = rgen.gen_range(0,std::u64::MAX);
                sub_energyset.update_edges(prr.inv_edge.max(4)+1);
                sub_energyset.update_distance();
                //Residue レベルのグルーピングがあると面倒なことになるかも
                let mut subatom_based_energies:Vec<f64> = vec![0.0;subatomnum];
                loop{
                    for aa in subatom_based_energies.iter_mut(){
                        *aa = 0.0;
                    }
                    let prev_energy:f64 = sub_energyset.calc_energy(&mut subatom_based_energies);
                    pp_energy_mc::mc_iter_array(&mut sub_energyset
                        ,prr.num_iter
                        ,&fixflag_sub
                        ,Some(rss)
                        ,prr.num_target_atoms
                        ,para.sorter
                        ,prr.mov_dist
                        ,(para.accept_bound.0,para.accept_bound.1,para.steps_bound_check)
                        ,&sub_bond_restrictor
                        ,&vec![]
                        ,(1.1,prr.dec_factor)
                        ,(prr.inv_dist,prr.inv_edge)
                        ,prr.refine_retry
                        ,prr.lbfgs
                        ,None
                    );
                    for aa in subatom_based_energies.iter_mut(){
                        *aa = 0.0;
                    }
                    let ediff = prev_energy - sub_energyset.calc_energy(&mut subatom_based_energies);
                    if prr.refine_loop_threshold <= 0.0 || prr.refine_loop_threshold > ediff{
                        break;
                    }
                }
                for aa in subatom_based_energies.iter_mut(){
                    *aa = 0.0;
                }
                let ene = sub_energyset.calc_energy(&mut subatom_based_energies);
                if backbone_tops.len() < step2_trial{
                    let mut ppos:Vec<(usize,(f64,f64,f64))> = vec![];
                    for (vii,_vv) in backbone_originalpos.iter(){
                        ppos.push((*vii,sub_energyset.get_atom(submapper[*vii] as usize).get_xyz()));
                    }
                    backbone_tops.push((ene,ppos));
                    backbone_tops.sort_by(|a,b|{(a.0).partial_cmp(&b.0).unwrap()});
                }else if ene < backbone_tops[backbone_tops.len()-1].0{
                    let mut ppos:Vec<(usize,(f64,f64,f64))> = vec![];
                    for (vii,_vv) in backbone_originalpos.iter(){
                        ppos.push((*vii,sub_energyset.get_atom(submapper[*vii] as usize).get_xyz()));
                    }
                    let plenn = backbone_tops.len();
                    backbone_tops[plenn-1] = (ene,ppos);
                    backbone_tops.sort_by(|a,b|{(a.0).partial_cmp(&b.0).unwrap()});
                }

                /*
                let mut lines:Vec<String> = vec![];
                for (aii,aa) in sub_energyset.evoef2_env.md_envset.atoms.iter().enumerate(){
                    let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
                    att.temp_factor = subatom_based_energies[aii].max(0.0).min(999.0);
                    lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                }
                write_to_file(format!("test_{}.pdb",backbone_tops.len()).as_str(),lines);
                */

            }
        }
    }
    let mut min_energy:f64 = std::f64::INFINITY;
    let mut min_pos:Vec<(f64,f64,f64)> = vec![(0.0,0.0,0.0);anum];
    let mut atom_based_energies:Vec<f64> = vec![0.0;anum];
    for bb in backbone_tops.into_iter(){
        for aa in 0..anum{
            ppenergyset.get_atom_mut(aa).set_xyz(orig_pos[aa].0,orig_pos[aa].1,orig_pos[aa].2);   
        }
        for (vii,vv) in bb.1.iter(){
            ppenergyset.get_atom_mut(*vii).set_xyz(vv.0,vv.1,vv.2);
        }
        for (kk,rr) in missing_residues.iter(){
            for rrr in rr.iter(){
                //backbone に sidechain をフィットさせる。元の ROTAMER をそのまま使う。
                let vvec:&Vec<usize> = &sorted_atoms.get(kk).unwrap()[*rrr];
                let n_ca_c_cb:(i64,i64,i64,i64) = get_n_ca_c_cb(&ppenergyset,vvec);
                let mut target_atoms:Vec<Point3D> = vec![];
                let mut tempatoms:Vec<&dyn Vector3D> = vec![];
                let mut origatoms:Vec<Point3D> = vec![];
                if n_ca_c_cb.0 > -1{
                    let aa:usize = n_ca_c_cb.0 as usize;
                    tempatoms.push(ppenergyset.get_atom(aa));
                    origatoms.push(Point3D::new(orig_pos[aa].0,orig_pos[aa].1,orig_pos[aa].2));
                }
                if n_ca_c_cb.1 > -1{
                    let aa:usize = n_ca_c_cb.1 as usize;
                    tempatoms.push(ppenergyset.get_atom(aa));
                    origatoms.push(Point3D::new(orig_pos[aa].0,orig_pos[aa].1,orig_pos[aa].2));
                }
                if n_ca_c_cb.2 > -1{
                    let aa:usize = n_ca_c_cb.2 as usize;
                    tempatoms.push(ppenergyset.get_atom(aa));
                    origatoms.push(Point3D::new(orig_pos[aa].0,orig_pos[aa].1,orig_pos[aa].2));
                }
                if tempatoms.len() != 3{
                    eprintln!("???? can not find n ca c.");
                    continue;
                }
                for aa in vvec.iter(){
                    target_atoms.push(Point3D::new(orig_pos[*aa].0,orig_pos[*aa].1,orig_pos[*aa].2));
                }
                //let mover:Vec<&mut dyn Vector3D> = target_atoms.iter_mut().map(|m| &mut (*m)).collect();//build できない。。。
                let mut mover:Vec<&mut dyn Vector3D> = vec![];
                for tt in target_atoms.iter_mut(){
                    mover.push(tt);
                }
                process_3d::fit_to_vector(
                    &origatoms[0],&origatoms[1],&origatoms[2],
                    tempatoms[0],tempatoms[1],tempatoms[2]
                ,&mut mover);
                for (aii,aa) in vvec.iter().enumerate(){
                    ppenergyset.get_atom_mut(*aa).set_xyz(target_atoms[aii].get_x(),target_atoms[aii].get_y(),target_atoms[aii].get_z());
                }
            }
        }
        
        for (vii,vv) in bb.1.iter(){//計算を間違っていなければほとんど同じだろうが一応座標を計算で得たものにもう一度合わせる
            ppenergyset.get_atom_mut(*vii).set_xyz(vv.0,vv.1,vv.2);
        }

        for (_pii,para) in params2.iter().enumerate(){
            let mut pline:Vec<RefinementParamLine> = para.generate_lines();
            if para.shuffle{
                pline.shuffle(rgen);
            }
            for prr in pline.into_iter(){
                println!("{:?}",prr);
                let bond_restrictor:Vec<charmm_based_energy::BondVars> = ppenergyset.evoef2_env.charmm_vars.bondvec.clone();
                
                let rss:u64 = rgen.gen_range(0,std::u64::MAX);
                ppenergyset.update_edges(prr.inv_edge.max(4)+1);
                ppenergyset.update_distance();
                for aa in atom_based_energies.iter_mut(){
                    *aa = 0.0;
                }
                let prev_energy:f64 = ppenergyset.calc_energy(&mut atom_based_energies);
                pp_energy_mc::mc_iter_array(ppenergyset
                    ,prr.num_iter
                    ,&fixed_atoms
                    ,Some(rss)
                    ,prr.num_target_atoms
                    ,para.sorter
                    ,prr.mov_dist
                    ,(para.accept_bound.0,para.accept_bound.1,para.steps_bound_check)
                    ,&bond_restrictor
                    ,&vec![]
                    ,(1.1,prr.dec_factor)
                    ,(prr.inv_dist,prr.inv_edge)
                    ,prr.refine_retry
                    ,prr.lbfgs
                    ,None
                );
                for aa in atom_based_energies.iter_mut(){
                    *aa = 0.0;
                }
                let pres = ppenergyset.calc_energy(&mut atom_based_energies);
                let ediff:f64 = prev_energy-pres;
                println!("improve: {}",ediff);
                if prr.refine_loop_threshold <= 0.0 || prr.refine_loop_threshold > ediff{
                    break;
                }
            }
        }
        
        for aa in atom_based_energies.iter_mut(){
            *aa = 0.0;
        }
        let pres = ppenergyset.calc_energy(&mut atom_based_energies);
        if min_energy > pres{ 
            min_energy = pres;
            for aa in 0..anum{
                min_pos[aa] = ppenergyset.get_atom(aa).get_xyz();
            }
        }
    }
    for aa in 0..anum{
        ppenergyset.get_atom_mut(aa).set_xyz(min_pos[aa].0,min_pos[aa].1,min_pos[aa].2);
    }
}



pub fn calc_energies(pdbfile:&str
    ,topfile:&str
    ,paramfile:&str
    ,evoefdir:&str
    ,backbone_dihedral_angle_file:&str
    ,cb_distance_file:&str
    ,outfile:&str
    ){
    let parr = charmm_param::CHARMMParam::load_chamm19(topfile,paramfile);
    let mut pdbb:pdbdata::PDBEntry = mmcif_process::load_pdb(pdbfile);
    for cc in pdbb.get_model_at(0).get_entity_at(0).iter_mut_asyms(){
        cc.remove_alt(None);
        charmm_based_energy::MDAtom::change_to_charmmnames(&mut cc.iter_mut_comps().map(|m|*m).collect());
    }
    let (md_envset,md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&pdbb.get_all_asyms().iter().map(|m|m).collect(),&parr,true);
    
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(
        &backbone_dihedral_angle_file);
    let (torsion,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,true,false),false);
    
    let mut distbe:Vec<distance_energy::AtomBinnedDistanceEnergy> = vec![];
    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&md_varset.dihedvec,&torsion,&omegas_general);

    if cb_distance_file.len() > 0{
        let dist_bvvec:Vec<(String,usize,String,usize,distance_energy::AtomBinnedDistanceEnergy)> = distance_energy::AtomBinnedDistanceEnergy::load_name_mapped("test/6F3H_B.cbdistbin.dat");
        distbe = distance_energy::AtomBinnedDistanceEnergy::assign_atom_ids(&md_envset,dist_bvvec);
    }

    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.0;
    }
    let num_atoms:usize = md_envset.atoms.len();
    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,evoefdir,false),
        backbone_energy_omega:omegas_general,
        backbone_energy_phi_psi:torsion,
        atom_distance_energy:vec![],
        atom_binned_distance_energy:distbe,
        atom_contact_energy:vec![],
        weights:PPEnergyWeight::new()
    };
    ppenergyset.weights.dihed_weight_charmm = weight_dihed;
    ppenergyset.weights.charmm_electro = 0.0;
    ppenergyset.weights.charmm_lj = 1.0;
    

    let mut lines:Vec<String> = vec![];
    let mut atom_based_energies:Vec<f64> = vec![0.0;num_atoms];
    let pres = ppenergyset.calc_energy(&mut atom_based_energies);
    println!("energy:{}",pres);
    for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
        att.temp_factor = atom_based_energies[aii].max(0.0).min(999.0);
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file(format!("{}.pdb",outfile).as_str(),lines);


    let mut sparse_atoms:Vec<usize> = vec![];
    for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        if aa.atom_name == "CA"
        || aa.atom_name == "CB" 
        || aa.atom_name == "C" 
        || aa.atom_name == "N"{
        sparse_atoms.push(aii);
        }
    }

    let subatomnum:usize = sparse_atoms.len();
    let (mut sub_energyset,_submapper):(PPEnergySet,Vec<i64>) = ppenergyset.make_set_for_subenv(&sparse_atoms);
    
    
    //sub_energyset.weights.dihed_weight_charmm = weight_dihed;
    sub_energyset.weights.charmm_angle = 1.0;
    sub_energyset.weights.charmm_electro = 0.0;
    sub_energyset.weights.charmm_lj = 1.0;
    sub_energyset.weights.evoef2_vdw =  vec![0.0;4];
    
    sub_energyset.weights.evoef2_vdw[evoef2_energy::VDW_INTER_ATT] = 0.0;
    sub_energyset.weights.evoef2_vdw[evoef2_energy::VDW_INTRA_ATT] = 0.0;

    sub_energyset.weights.evoef2_elec = 0.0;
    sub_energyset.weights.evoef2_hb = vec![0.0;9];
    sub_energyset.weights.evoef2_desolv_polar =  0.0;
    sub_energyset.weights.evoef2_desolv_nonpolar = 0.0;
    sub_energyset.weights.evoef2_ss = 0.0;
    sub_energyset.weights.backbone_energy_phi_psi = 1.0;
    sub_energyset.weights.backbone_energy_omega = 1.0;

    sub_energyset.update_distance();
    
    let mut subatom_based_energies:Vec<f64> = vec![0.0;subatomnum];             
    let pres = sub_energyset.calc_energy(&mut subatom_based_energies);
    let mut lines:Vec<String> = vec![];
    for (aii,aa) in sub_energyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
        att.temp_factor = subatom_based_energies[aii].max(0.0).min(999.0);
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file(format!("{}.cb.pdb",outfile).as_str(),lines);
    println!("energy_cb:{}",pres);
}



pub fn build_peptide(n:&dyn Vector3D,ca:&dyn Vector3D,c:&dyn Vector3D,prev_psi_radian:f64,omega_radian:f64,phi_radian:f64)->Vec<Point3D>{
    
    //charmm19 
    let angle_ca_c_n:f64 = 117.5;
    let angle_c_ca_n:f64 = 111.6;
    let angle_c_n_ca:f64 = 120.0; 

    let bond_c_ca:f64 = 1.52;
    let bond_ca_n:f64 = 1.45;
    let bond_c_n :f64 = 1.33;
    

    let mut new_n  = get_zigzag_pos(n,ca,c,angle_ca_c_n/180.0*PI,bond_c_n);
    let mut new_ca = get_zigzag_pos(ca,c,&new_n,angle_c_n_ca/180.0*PI,bond_ca_n);
    let mut new_c  = get_zigzag_pos(c,&new_n,&new_ca,angle_c_ca_n/180.0*PI,bond_c_ca);

    rotate_bond(ca,c,prev_psi_radian-PI,&mut vec![&mut new_n,&mut new_ca,&mut new_c]);
    if omega_radian != PI{
        rotate_bond(c,&new_n.get_copy(),omega_radian-PI,&mut vec![&mut new_ca,&mut new_c]);
    }
    
    rotate_bond(&new_n.get_copy(),&new_ca.get_copy(),phi_radian-PI,&mut vec![&mut new_c]);
    
    return vec![new_n,new_ca,new_c];
}


pub fn get_cb_pos(n:&dyn Vector3D,ca:&dyn Vector3D,c:&dyn Vector3D)->Point3D{
    
    //charmm19 
    let improper_n_ca_c_cb:f64 = 45.0;
    let angle_cb_ca_c:f64 = 109.5; 

    let bond_cb_ca :f64 = 1.52;
    
    let mut ppos:Point3D = get_impr_pos(n,ca,c,angle_cb_ca_c/180.0*PI,bond_cb_ca);    
    rotate_bond(ca,c,improper_n_ca_c_cb/180.0*PI,&mut vec![&mut ppos]);
    
    return ppos;
}


pub fn build_peptide_backward(n:&dyn Vector3D,ca:&dyn Vector3D,c:&dyn Vector3D,psi_radian:f64,omega_radian:f64,next_phi_radian:f64)->Vec<Point3D>{
    
    //charmm19 
    let angle_ca_c_n:f64 = 117.5;
    let angle_c_ca_n:f64 = 111.6;
    let angle_c_n_ca:f64 = 120.0; 

    let bond_c_ca:f64 = 1.52;
    let bond_ca_n:f64 = 1.45;
    let bond_c_n :f64 = 1.33;
    

    let mut new_c  = get_zigzag_pos(c,ca,n,angle_c_n_ca/180.0*PI,bond_c_n);
    let mut new_ca = get_zigzag_pos(ca,n,&new_c,angle_ca_c_n/180.0*PI,bond_c_ca);
    let mut new_n  = get_zigzag_pos(n,&new_c,&new_ca,angle_c_ca_n/180.0*PI,bond_ca_n);

    
    rotate_bond(ca,n,next_phi_radian-PI,&mut vec![&mut new_c,&mut new_ca,&mut new_n]);
    if omega_radian != PI{
        rotate_bond(n,&new_c.get_copy(),omega_radian-PI,&mut vec![&mut new_ca,&mut new_n]);
    }
    rotate_bond(&new_c.get_copy(),&new_ca.get_copy(),psi_radian-PI,&mut vec![&mut new_n]);
    
    return vec![new_n,new_ca,new_c];
}


//v1,v2,v3 上の平面にあって、v3-v2-new のアングルが radian, v3-new の距離が dist の点を返す。
pub fn get_impr_pos(v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D,radian:f64,dist:f64)->Point3D{
    let mut norm0:Point3D = Point3D::new_tup(&process_3d::calc_norm(
        v1.get_x() - v2.get_x(),
        v1.get_y() - v2.get_y(),
        v1.get_z() - v2.get_z(),
        v3.get_x() - v2.get_x(),
        v3.get_y() - v2.get_y(),
        v3.get_z() - v2.get_z()
    ));
    norm0.standardize();
    let mut newpos:Point3D = Point3D::new_tup(&process_3d::subtracted_t(&(v3.get_xyz()),&(v2.get_xyz())));
    process_3d::rotate_3d(&mut vec![&mut newpos],&norm0,radian*-1.0);
    newpos.standardize();
    newpos.multiply(dist);
    newpos.add(&(v2.get_xyz()));
    return newpos;
}


//v1,v2,v3 上の平面にあって、v2-v3-new のアングルが radian, v3-new の距離が dist の点を返す。
pub fn get_zigzag_pos(v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D,radian:f64,dist:f64)->Point3D{
    let mut norm0:Point3D = Point3D::new_tup(&process_3d::calc_norm(
        v1.get_x() - v2.get_x(),
        v1.get_y() - v2.get_y(),
        v1.get_z() - v2.get_z(),
        v3.get_x() - v2.get_x(),
        v3.get_y() - v2.get_y(),
        v3.get_z() - v2.get_z()
    ));
    norm0.standardize();
    let mut newpos:Point3D = Point3D::new_tup(&process_3d::subtracted_t(&(v2.get_xyz()),&(v3.get_xyz())));
    process_3d::rotate_3d(&mut vec![&mut newpos],&norm0,radian);
    newpos.standardize();
    newpos.multiply(dist);
    newpos.add(&(v3.get_xyz()));
    return newpos;
}

//vstart->vend というベクトルを中心に radian 度回転させる
pub fn rotate_bond(vstart:&dyn Vector3D,vend:&dyn Vector3D,radian:f64,points:&mut Vec<&mut dyn Vector3D>){
    let mut v_angle = Point3D::new_tup(&process_3d::subtracted_t(&vend.get_xyz(),&vstart.get_xyz()));
    v_angle.standardize();

    let dff =  vstart.get_xyz();
    for pp in points.iter_mut(){
        let pt:(f64,f64,f64) = pp.get_xyz();
        pp.set_xyz(pt.0-dff.0,pt.1-dff.1,pt.2-dff.2);
    }
    
    process_3d::rotate_3d(points, &v_angle ,radian);
    for pp in points.iter_mut(){
        let pt:(f64,f64,f64) = pp.get_xyz();
        pp.set_xyz(pt.0+dff.0,pt.1+dff.1,pt.2+dff.2);
    }
}




#[test]
fn peptide_build_test(){
    let mut natom = pdbdata::PDBAtom::new();
    natom.set_x(0.0);
    natom.set_y(0.0);
    natom.set_z(0.0);
   
    let mut caatom = pdbdata::PDBAtom::new();
    caatom.set_x(1.0);
    caatom.set_y(0.0);
    caatom.set_z(0.0);
   
    let mut catom = pdbdata::PDBAtom::new();
    catom.set_x(1.0);
    catom.set_y(1.0);
    catom.set_z(0.0);
    let mut atoms:Vec<pdbdata::PDBAtom> = vec![natom,caatom,catom];
    for _ in 0..10{
        let alen:usize = atoms.len();
        let pos:Vec<Point3D> = build_peptide(&atoms[alen-3],&atoms[alen-2],&atoms[alen-1],-50.0/180.0*PI,PI,-60.0/180.0*PI);
        //let pos:Vec<Point3D> = build_peptide(&atoms[alen-3],&atoms[alen-2],&atoms[alen-1],PI,PI,PI);
        
        for pp in pos.iter(){
            let mut atom = pdbdata::PDBAtom::new();
            atom.set_xyz(pp.get_x(),pp.get_y(),pp.get_z());
            atoms.push(atom);
        }
    }
    let mut lines:Vec<String> = vec![];
    for (aii,aa) in atoms.iter_mut().enumerate(){
        if aii%3 == 0{
            aa.atom_code = "N".to_owned();
        }
        if aii%3 == 1{
            aa.atom_code = "CA".to_owned();
        }
        if aii%3 == 2{
            aa.atom_code = "C".to_owned();
        }
        lines.push(aa.get_pdb_atom_line_string("A","UNK",aii as i64/3,""));
    }

    write_to_file("test/pep_built.pdb",lines);
}

#[test]
fn peptide_build_test2(){
    let mut first:Vec<(f64,f64,f64)> = vec![];
    let mut last:Vec<(f64,f64,f64)> = vec![];
    let mut phi_psi_omega:Vec<(f64,f64,f64)> = vec![];
    let mut pdbb:pdbdata::PDBEntry = mmcif_process::load_pdb("example_files/1a4w_part.pdb");
    for cc in pdbb.get_model_at(0).get_entity_at(0).iter_mut_asyms(){
        cc.remove_alt(None);
        let rnum =cc.num_comps();
        let dist_threshold:f64 = 2.0;
        for rr in 0..rnum{
            
            let mut prevca_:Option<&pdbdata::PDBAtom> = if rr == 0 {None}else{cc.get_comp_at(rr-1).get_CA()};
            let mut prevc_:Option<&pdbdata::PDBAtom> = if rr == 0 {None}else{cc.get_comp_at(rr-1).get_C()};
            let currentn_:Option<&pdbdata::PDBAtom> = cc.get_comp_at(rr).get_N();
            let currentca_:Option<&pdbdata::PDBAtom> = cc.get_comp_at(rr).get_CA();
            let currentc_:Option<&pdbdata::PDBAtom> = cc.get_comp_at(rr).get_C();
            let mut nextn_:Option<&pdbdata::PDBAtom> = if rr == rnum-1{None}else{cc.get_comp_at(rr+1).get_N()};
            let mut nextca_:Option<&pdbdata::PDBAtom> = if rr == rnum-1{None}else{cc.get_comp_at(rr+1).get_CA()};
            if rr == 0{
                first.push(currentn_.unwrap().get_xyz());
                first.push(currentca_.unwrap().get_xyz());
                first.push(currentc_.unwrap().get_xyz());
            }
            if rr == rnum-1{
                last.push(currentn_.unwrap().get_xyz());
                last.push(currentca_.unwrap().get_xyz());
                last.push(currentc_.unwrap().get_xyz());
            }
            
            //前の残基とは結合してないとみなす
            if let None = prevc_{
                prevca_ = None;
                prevc_ = None;
            }else if let None = currentn_{
                prevca_ = None;
                prevc_ = None;
            }else{
                if process_3d::distance(&prevc_.clone().unwrap().get_xyz(),&currentn_.clone().unwrap().get_xyz()) > dist_threshold{
                    prevca_ = None;
                    prevc_ = None;
                }
            }

            
            //後ろの残基とは結合してないとみなす
            if let None = nextn_{
                nextn_ = None;
                nextca_ = None;
            }else if let None = currentc_{
                nextn_ = None;
                nextca_ = None;
            }else{
                if process_3d::distance(&currentc_.clone().unwrap().get_xyz(),&nextn_.clone().unwrap().get_xyz()) > dist_threshold{
                    nextn_ = None;
                    nextca_ = None;
                }
            } 


            let mut phi:f64 = 0.0;
            let mut psi:f64 = 0.0;
            let mut _prev_omega:f64 = 0.0;
            let mut next_omega:f64 = 0.0;


            //dihedral angle は charmm とは逆回り（なはず）なので -1.0 を掛ける

            //こんな書き方しかできないだろうか。。。
            let mut phiflag = true;
            let mut psiflag = true;
            let mut pomegaflag = true;
            let mut nomegaflag = true;

            if let None = prevca_{
                pomegaflag = false;
            }
            if let None = prevc_{
                pomegaflag = false;
                phiflag = false;
            }
            if let None = currentn_{
                pomegaflag = false;
                phiflag = false;
                psiflag = false;
            }
            if let None = currentca_{
                pomegaflag = false;
                phiflag = false;
                psiflag = false;
                nomegaflag = false;
            }
            if let None = currentc_{
                phiflag = false;
                psiflag = false;
                nomegaflag = false;
            }
            if let None = nextn_{
                psiflag = false;
                nomegaflag = false;
            }
            if let None = nextca_{
                nomegaflag = false;
            }

            if pomegaflag{
                _prev_omega = charmm_based_energy::calc_dihedral_angle_radian(prevca_.clone().unwrap(),prevc_.clone().unwrap(),currentn_.clone().unwrap(),currentca_.clone().unwrap())*-1.0;
            }
            
            
            if rr == 0{
                phiflag = true;    
            }
            if rr == rnum -1{
                psiflag = true;
                nomegaflag = true;
            }

            if !phiflag || !psiflag || !nomegaflag{
                break;
            }

            if rr != 0 && phiflag{
                phi = charmm_based_energy::calc_dihedral_angle_radian(prevc_.clone().unwrap(),currentn_.clone().unwrap(),currentca_.clone().unwrap(),currentc_.clone().unwrap())*-1.0;
            }
            
            if rr != rnum-1 && phiflag {
                psi = charmm_based_energy::calc_dihedral_angle_radian(currentn_.clone().unwrap(),currentca_.clone().unwrap(),currentc_.clone().unwrap(),nextn_.clone().unwrap())*-1.0;
            }
            if rr != rnum-1 && nomegaflag {
                next_omega = charmm_based_energy::calc_dihedral_angle_radian(currentca_.clone().unwrap(),currentc_.clone().unwrap(),nextn_.clone().unwrap(),nextca_.clone().unwrap())*-1.0;
            }
            phi_psi_omega.push((phi,psi,next_omega));
        }
        break;
    }

    let mut atoms:Vec<pdbdata::PDBAtom> = vec![];

    for ii in 0..phi_psi_omega.len(){
        if ii == 0{
            for aa in 0..3{
                let mut atom = pdbdata::PDBAtom::new();
                atom.set_xyz(first[aa].0,first[aa].1,first[aa].2);
                atoms.push(atom);
            }
            
            let cb:Point3D = get_cb_pos(&atoms[0],&atoms[1],&atoms[2]);
            let mut cbatom = pdbdata::PDBAtom::new();
            cbatom.set_xyz(cb.get_x(),cb.get_y(),cb.get_z());
            atoms.push(cbatom);

            continue;
        }

        let alen:usize = atoms.len();
        let pos:Vec<Point3D> = build_peptide(&atoms[alen-4],&atoms[alen-3],&atoms[alen-2],phi_psi_omega[ii-1].1,phi_psi_omega[ii-1].2,phi_psi_omega[ii].0);
        
        for pp in pos.iter(){
            let mut atom = pdbdata::PDBAtom::new();
            atom.set_xyz(pp.get_x(),pp.get_y(),pp.get_z());
            atoms.push(atom);
        }
        let cb:Point3D = get_cb_pos(&pos[0],&pos[1],&pos[2]);
        let mut cbatom = pdbdata::PDBAtom::new();
        cbatom.set_xyz(cb.get_x(),cb.get_y(),cb.get_z());
        atoms.push(cbatom);
    }
    let mut lines:Vec<String> = vec![];
    for (aii,aa) in atoms.iter_mut().enumerate(){
        if aii%4 == 0{
            aa.atom_code = "N".to_owned();
        }
        if aii%4 == 1{
            aa.atom_code = "CA".to_owned();
        }
        if aii%4 == 2{
            aa.atom_code = "C".to_owned();
        }
        if aii%4 == 3{
            aa.atom_code = "CB".to_owned();
        }
        lines.push(aa.get_pdb_atom_line_string("A","ALA",aii as i64/4,""));
    }
    write_to_file("test/pep_built_phiphiomega.pdb",lines);


    let mut atoms:Vec<pdbdata::PDBAtom> = vec![];
    phi_psi_omega.reverse();
    for ii in 0..phi_psi_omega.len(){
        //println!("{:?}",phi_psi_omega[ii]);
        if ii == 0{
            for aa in 0..3{
                let mut atom = pdbdata::PDBAtom::new();
                atom.set_xyz(last[aa].0,last[aa].1,last[aa].2);
                atoms.push(atom);
            }
            continue;
        }

        let alen:usize = atoms.len();
        let pos:Vec<Point3D> = build_peptide_backward(&atoms[alen-3],&atoms[alen-2],&atoms[alen-1],phi_psi_omega[ii].1,phi_psi_omega[ii].2,phi_psi_omega[ii-1].0);
        for pp in pos.iter(){
            let mut atom = pdbdata::PDBAtom::new();
            atom.set_xyz(pp.get_x(),pp.get_y(),pp.get_z());
            atoms.push(atom);
        }
    }
    atoms.reverse();
    let mut lines:Vec<String> = vec![];
    for (aii,aa) in atoms.iter_mut().enumerate(){
        if aii%3 == 2{
            aa.atom_code = "N".to_owned();
        }
        if aii%3 == 1{
            aa.atom_code = "CA".to_owned();
        }
        if aii%3 == 0{
            aa.atom_code = "C".to_owned();
        }
        lines.push(aa.get_pdb_atom_line_string("A","GLY",aii as i64/3,""));
    }
    write_to_file("test/pep_built_phiphiomega_back.pdb",lines);
   
}
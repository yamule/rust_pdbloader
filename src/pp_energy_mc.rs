extern crate rand;
use super::charmm_based_energy;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
use self::rand::prelude::*;
use super::misc_util::*;
use super::geometry::Vector3D;
use super::geometry::Point3D;
use super::energy_lbfgs;
#[allow(unused_imports)]
use super::pdbdata;
#[allow(unused_imports)]
use super::charmm_param;
#[allow(unused_imports)]
use super::process_3d;
#[allow(unused_imports)]
use super::debug_env;
#[allow(unused_imports)]
use super::backbone_sample;
#[allow(unused_imports)]
use super::side_chain_sample;
#[allow(unused_imports)]
use super::chain_builder;
#[allow(unused_imports)]
use super::sequence_alignment;

use super::peptide_backbone_dihedral_energy;
use super::distance_energy;

use std::collections::HashMap;
use std::collections::HashSet;

use super::evoef2_energy;
use super::pp_energy::PPEnergyWeight;
use super::pp_energy::PPEnergySet;

//axisvec を軸として radian 度 回転する
//axisvec は standardize されている必要がある
//https://ja.wikipedia.org/wiki/%E3%83%AD%E3%83%89%E3%83%AA%E3%82%B2%E3%82%B9%E3%81%AE%E5%9B%9E%E8%BB%A2%E5%85%AC%E5%BC%8F
pub fn rotate_atom_3d(target:&mut Vec<charmm_based_energy::MDAtom>,atom_idx:&Vec<usize>,axisvec:&dyn Vector3D,radian:f64){
    let cc = radian.cos();
    let ss = radian.sin();
    let nx = axisvec.get_x();
    let ny = axisvec.get_y();
    let nz = axisvec.get_z();

    let vv:Vec<Vec<f64>> = vec![
         vec![cc+nx*nx*(1.0-cc)     ,nx*ny*(1.0-cc)-nz*ss   ,nx*nz*(1.0-cc)+ny*ss   ]
        ,vec![ny*nx*(1.0-cc)+nz*ss  ,cc+ny*ny*(1.0-cc)      ,ny*nz*(1.0-cc)-nx*ss   ]
        ,vec![nz*nx*(1.0-cc)-ny*ss  ,nz*ny*(1.0-cc)+nx*ss   ,cc+nz*nz*(1.0-cc)      ]
    ];
    
    for tt in atom_idx.iter(){
        let x = target[*tt].get_x();
        let y = target[*tt].get_y();
        let z = target[*tt].get_z();
        target[*tt].set_xyz(
            x*vv[0][0] + y*vv[0][1] + z*vv[0][2],
            x*vv[1][0] + y*vv[1][1] + z*vv[1][2],
            x*vv[2][0] + y*vv[2][1] + z*vv[2][2]
        );
    }
}


pub fn random_group_rotate(atoms:&mut Vec<charmm_based_energy::MDAtom>
    ,group:&Vec<Vec<usize>>,degree_max:f64,mov_max:f64,num_movgroup:usize,rgen:&mut StdRng){
    for _ in 0..num_movgroup{
        let mut dvec:Point3D = Point3D::new(rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0));
        dvec.standardize();
        let target_group:usize = rgen.gen_range(0,group.len());
        let bb:f64 = rgen.gen_range(-1.0,1.0);
        mov_group(atoms,&group[target_group],degree_max,mov_max,rgen,bb >= 0.0,bb < 0.0);
    }
}

//group として与えられた index の原子を回転および移動する
pub fn mov_group(
    atoms:&mut Vec<charmm_based_energy::MDAtom>
    ,group:&Vec<usize>
    ,degree_max:f64
    ,mov_max:f64
    ,rgen:&mut StdRng
    ,fix_start:bool
    ,fix_end:bool){
    
    if degree_max > 0.0{
        let mut dvec:Point3D = Point3D::new(rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0));
        dvec.standardize();
        let startpos:(f64,f64,f64) = atoms[group[0]].get_xyz();
        let endpos:(f64,f64,f64) = atoms[group[group.len()-1]].get_xyz(); 
        rotate_atom_3d(atoms,&group,&dvec,rgen.gen_range(-degree_max,degree_max));
        if fix_start{
            let tpos:(f64,f64,f64) = atoms[group[0]].get_xyz();
            let pxx = startpos.0-tpos.0;
            let pxy = startpos.1-tpos.1;
            let pxz = startpos.2-tpos.2;
            for aa in group.iter(){
                let xyz = atoms[*aa].get_xyz();
                atoms[*aa].set_xyz(xyz.0+pxx,xyz.1+pxy,xyz.2+pxz);
            }
        }else if fix_end{
            let tpos:(f64,f64,f64) = atoms[group[group.len()-1]].get_xyz();
            let pxx = endpos.0-tpos.0;
            let pxy = endpos.1-tpos.1;
            let pxz = endpos.2-tpos.2;
            for aa in group.iter(){
                let xyz = atoms[*aa].get_xyz();
                atoms[*aa].set_xyz(xyz.0+pxx,xyz.1+pxy,xyz.2+pxz);
            }
        }else{//ToDo 中心点を取る
            if rgen.gen_range(-1.0,1.0) < 0.0{
                let tpos:(f64,f64,f64) = atoms[group[0]].get_xyz();
                let pxx = startpos.0-tpos.0;
                let pxy = startpos.1-tpos.1;
                let pxz = startpos.2-tpos.2;
                for aa in group.iter(){
                    let xyz = atoms[*aa].get_xyz();
                    atoms[*aa].set_xyz(xyz.0+pxx,xyz.1+pxy,xyz.2+pxz);
                }
            }else{
                let tpos:(f64,f64,f64) = atoms[group[group.len()-1]].get_xyz();
                let pxx = endpos.0-tpos.0;
                let pxy = endpos.1-tpos.1;
                let pxz = endpos.2-tpos.2;
                for aa in group.iter(){
                    let xyz = atoms[*aa].get_xyz();
                    atoms[*aa].set_xyz(xyz.0+pxx,xyz.1+pxy,xyz.2+pxz);
                }
            }

        }
    }
    if mov_max > 0.0{
        let mut pvec:Point3D = Point3D::new(rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0),rgen.gen_range(-1.0,1.0));
        pvec.standardize();

        let prat:f64 = rgen.gen_range(-mov_max,mov_max);
        let pxx = prat*pvec.get_x();
        let pxy = prat*pvec.get_y();
        let pxz = prat*pvec.get_z();
        for aa in group.iter(){
            let xyz = atoms[*aa].get_xyz();
            atoms[*aa].set_xyz(xyz.0+pxx,xyz.1+pxy,xyz.2+pxz);
        }
    }

}


pub fn get_connected_group(
    atomnum:usize
    ,atoms_group:&Vec<Vec<usize>>
    ,bonds:&Vec<charmm_based_energy::BondVars>
    )->HashMap<usize,Vec<usize>>{
    let mut atom_to_group:Vec<usize> = (0..atomnum).into_iter().collect();

    let get_smallest = |a:usize,mapper:&Vec<usize>|->usize{ 
        let mut c = a;
        while c != mapper[c]{
            c = mapper[c];
        }
        c
    };

    for gg in atoms_group.iter(){
        let mut sorter:Vec<&usize> = gg.iter().collect();
        sorter.sort();
        let mut gmin:usize = *sorter[0];
        for ss in sorter.iter(){
            let ssm = get_smallest(**ss,&atom_to_group);
            if ssm < gmin{
                gmin = ssm;
            }
        }
        for ss in sorter.iter(){
            let ssm = get_smallest(**ss,&atom_to_group);
            atom_to_group[**ss] = gmin;
            atom_to_group[ssm] = gmin;
        }
    }
    for bb in bonds.iter(){
        let ssm0 = get_smallest(bb.atoms.0,&atom_to_group);
        let ssm1 = get_smallest(bb.atoms.1,&atom_to_group);
        let gmin = ssm0.min(ssm1);
        atom_to_group[ssm1] = gmin;
        atom_to_group[ssm0] = gmin;
        atom_to_group[bb.atoms.0] = gmin;
        atom_to_group[bb.atoms.1] = gmin;

    }
    
    for aa in 0..atomnum{
        if atom_to_group[aa] != aa{
            assert!(atom_to_group[aa] < aa);
            let ssm = get_smallest(aa,&atom_to_group);
            atom_to_group[aa] = ssm;
        }
    }
    let mut ret:HashMap<usize,Vec<usize>> = HashMap::new();
    for (aii,aa) in atom_to_group.iter().enumerate(){
        if !ret.contains_key(aa){
            ret.insert(*aa,vec![]);
        }
        ret.get_mut(aa).unwrap().push(aii);
    }
    return ret;
}


pub fn mc_iter_group(
    energyset:&mut PPEnergySet
    ,iternum:usize
    ,atoms_fixed:&HashSet<usize>
    ,atoms_group:&Vec<Vec<usize>>
    ,random_seed:Option<u64>
    ,degree_max:f64
    ,mov_max:f64
    ,num_movgroup:usize
    ,acceptance_bound_:(f64,f64,usize)
    ,bond_restrictor:&Vec<charmm_based_energy::BondVars>
    ,bond_attractor:&Vec<charmm_based_energy::BondVars>
    ){
    let mut rgen:StdRng;
    let mut ffloating:Vec<usize> = vec![];
    for ii in 0..energyset.evoef2_env.md_envset.atoms.len(){
        if !atoms_fixed.contains(&ii){
            ffloating.push(ii);
        }
    }
    let checkstep:usize = acceptance_bound_.2;//checstep ごとに採択率を計算し、改善が無い場合の採択率が lowerbound と upperbound の間に収まるようにする
    let accept_lowerbound:f64 = acceptance_bound_.0;
    let accept_upperbound:f64 = acceptance_bound_.1;
    let _groupnum:usize = atoms_group.len();
    //Dihedral が低いので値確かめ
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }
    
    let mut tmpatoms_prevpos:Vec<(f64,f64,f64)> =vec![];
    let mut tmpatoms_checkpoint:Vec<(f64,f64,f64)> =vec![];

    for aa in energyset.evoef2_env.md_envset.atoms.iter(){
        tmpatoms_prevpos.push((aa.get_x(),aa.get_y(),aa.get_z()));
        tmpatoms_checkpoint.push((aa.get_x(),aa.get_y(),aa.get_z()));
    }

    let num_atoms:usize = energyset.evoef2_env.md_envset.atoms.len();
 
    let mut checkpoint_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut prev_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut current_energies:Vec<f64> = vec![0.0;num_atoms];
    
    let mut minenergy_pos:Vec<(f64,f64,f64)> =vec![(0.0,0.0,0.0);num_atoms];
    let mut min_energy:f64;
    fix_group_bond(
        &mut energyset.evoef2_env.md_envset.atoms
    ,atoms_group
    ,bond_restrictor
    ,bond_attractor
    ,&ffloating
    ,&vec![]
    ,0.2
    ,atoms_group.len()*10
    );
    energyset.update_distance();
   

    let mut prev_scoresum = energyset.calc_energy(&mut current_energies);
    for ii in 0..num_atoms{
        prev_energies[ii] = current_energies[ii];
    }

    min_energy = prev_scoresum;
    for ii in 0..num_atoms{
        tmpatoms_checkpoint[ii].0 = energyset.get_atom(ii).get_x();
        tmpatoms_checkpoint[ii].1 = energyset.get_atom(ii).get_y();
        tmpatoms_checkpoint[ii].2 = energyset.get_atom(ii).get_z();
        minenergy_pos[ii].0  = energyset.get_atom(ii).get_x();
        minenergy_pos[ii].1  = energyset.get_atom(ii).get_y();
        minenergy_pos[ii].2  = energyset.get_atom(ii).get_z();
        checkpoint_energies[ii] = current_energies[ii];
    }

    let mut current_iter:usize = 0;
    let mut prev_checked:usize = 0;//最終 iter は checkstep の倍数でないことがあるので current_iter - checkstep でなく前回のチェックポイントを覚えておく
    let mut num_accepted:f64 = 0.0;//スコアが小さかったが採択された回数
    let mut num_negative:f64 = 0.0;//スコアが小さくなり棄却するかどうかの試行が行われた回数

    let mut beta:f64 = 1.0;

    loop{
        current_iter += 1;
        let deg:f64 = rgen.gen_range(0.01,degree_max);
        let mov:f64 = rgen.gen_range(0.01,mov_max);

        let mut idxx:Vec<usize> = (0..num_atoms).into_iter().collect();
        idxx.shuffle(&mut rgen);
        
        for ii in 0.. energyset.evoef2_env.md_envset.atoms.len(){
            tmpatoms_prevpos[ii].0 = energyset.get_atom(ii).get_x();
            tmpatoms_prevpos[ii].1 = energyset.get_atom(ii).get_y();
            tmpatoms_prevpos[ii].2 = energyset.get_atom(ii).get_z();
        }

        random_group_rotate(&mut energyset.evoef2_env.md_envset.atoms
            ,atoms_group
            ,deg
            ,mov
            ,num_movgroup
            ,&mut rgen);
        fix_group_bond(
             &mut energyset.evoef2_env.md_envset.atoms
            ,atoms_group
            ,bond_restrictor
            ,bond_attractor
            ,&ffloating
            ,&vec![]
            ,0.2
            ,atoms_group.len()*10
        );
        energyset.evoef2_env.md_envset.update_distance();
        
        for ii in 0..num_atoms{
            current_energies[ii] = 0.0;
        }
        let scoresum = energyset.calc_energy(&mut current_energies);
        if scoresum < min_energy{
            for ii in 0..num_atoms{
                minenergy_pos[ii].0 = energyset.get_atom(ii).get_x();
                minenergy_pos[ii].1 = energyset.get_atom(ii).get_y();
                minenergy_pos[ii].2 = energyset.get_atom(ii).get_z();
            }
            println!("min changed: {}->{} ",min_energy,scoresum);
            min_energy = scoresum;
        }
        if current_iter%100 == 0{
            println!("iter:{} {} {} {}",current_iter,scoresum,prev_scoresum,min_energy);
        }
        if true{
            if scoresum < prev_scoresum{
                //accepted
                prev_scoresum = scoresum;
                for ii in 0..num_atoms{
                    prev_energies[ii] = current_energies[ii];
                }
            }else{
                let dd:f64 = rgen.gen_range(0.0,1.0);
                let bb:f64 = ((scoresum-prev_scoresum)*-1.0*beta).exp();
                if prev_scoresum < scoresum{
                    num_negative += 1.0;
                    if dd < bb{
                        num_accepted += 1.0;
                        prev_scoresum = scoresum;
                        for ii in 0..num_atoms{
                            prev_energies[ii] = current_energies[ii];
                        }
                    }else{
                        for ii in 0..num_atoms{
                            energyset.get_atom_mut(ii).set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                        }
                    }
                }
            }
        }

        let checkflag =  current_iter%checkstep == 0 || current_iter == iternum;
        if checkflag{
            //println!("accepted:{} rejected:{}",num_accepted,num_negative);
            if num_negative > 0.0 && num_accepted/num_negative > accept_upperbound{
                //採択されすぎなので戻す
                beta *= 5.0;
                current_iter = prev_checked;
                for ii in 0..num_atoms{
                    energyset.get_atom_mut(ii).set_xyz(tmpatoms_checkpoint[ii].0,tmpatoms_checkpoint[ii].1,tmpatoms_checkpoint[ii].2);
                    prev_energies[ii] = checkpoint_energies[ii];
                }
            }else{
                if num_negative > 0.0 && num_accepted/num_negative < accept_lowerbound{
                    //採択されなさすぎなので採択率を上げる
                    beta /= 5.0;
                }
                prev_checked = current_iter;
                for ii in 0..num_atoms{
                    tmpatoms_checkpoint[ii].0 = energyset.get_atom(ii).get_x();
                    tmpatoms_checkpoint[ii].1 = energyset.get_atom(ii).get_y();
                    tmpatoms_checkpoint[ii].2 = energyset.get_atom(ii).get_z();
                    checkpoint_energies[ii] = current_energies[ii];
                }
            }
            num_negative = 0.0;
            num_accepted = 0.0;
        }
        if current_iter >= iternum{
            break;
        }
    }
    
    for ii in 0..num_atoms{
        energyset.get_atom_mut(ii).set_xyz(minenergy_pos[ii].0,minenergy_pos[ii].1,minenergy_pos[ii].2);
    }
    energyset.update_distance();
}

/*
chain ごとに残基順に並んだ原子の ID をグループ化して返す
*/
pub fn make_sorted_residue_group(atoms:&Vec<charmm_based_energy::MDAtom>)->HashMap<String,Vec<Vec<usize>>>{
    let mut ret:HashMap<String,Vec<Vec<usize>>> = HashMap::new();
    let mut hss:HashMap<(String,i64),Vec<usize>> = HashMap::new();

    for (aii,aa) in atoms.iter().enumerate(){
        let chain_resnum:(String,i64) = (aa.chain_name.to_string(),aa.residue_index_in_chain);
        if ! hss.contains_key(&chain_resnum){
            hss.insert(chain_resnum.clone(),vec![]);
        }
        if !ret.contains_key(&aa.chain_name){
            ret.insert(aa.chain_name.clone(),vec![]);
        }
        hss.get_mut(&chain_resnum).unwrap().push(aii);
    }
    let mut sorted:Vec<(String,i64)> = hss.iter().map(|m|m.0.clone()).collect();
    
    //chain は後で分解するので residue だけで sort する
    sorted.sort_by(|a,b|{a.1.cmp(&b.1)});
    for kk in sorted.iter_mut(){
        let mut vv = hss.remove(&kk).unwrap();
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
        let cname:&str = &atoms[vv[0]].chain_name;
        ret.get_mut(cname).unwrap().push(vv);
    }
    return ret;
}


//num_moveatom 個の原子をランダムに移動し、MH 基準で採択する
pub fn mc_iter_array(
    energyset:&mut PPEnergySet
    ,iternum:usize
    ,atoms_fixed:&HashSet<usize>
    ,random_seed:Option<u64>
    ,num_moveatom:usize
    ,movatom_energysort:bool
    ,movement_max:f64
    ,acceptance_bound_:(f64,f64,usize)
    ,bond_restrictor:&Vec<charmm_based_energy::BondVars>
    ,bond_upper_bound:&Vec<charmm_based_energy::BondVars>
    ,selection_weight_factor:(f64,f64)//ループごとに 0 が掛け算され、選ばれるごとに 1 が掛け算される
    ,involution_factor:(f64,u64)
    ,refine_retry:usize
    ,try_refine_lbfgs:usize
    ,group_rotate:Option<(&Vec<Vec<usize>>,f64)>//
    )->(f64,i64){//minenergy, lastupdate
        assert!(selection_weight_factor.0 >= 1.0);
        assert!(selection_weight_factor.1 <= 1.0);
    let mut rgen:StdRng;//乱数生成器
    let mut lastupdate:i64 = -1;
    let checkstep:usize = acceptance_bound_.2;//checstep ごとに採択率を計算し、改善が無い場合の採択率が lowerbound と upperbound の間に収まるようにする
    let accept_lowerbound:f64 = acceptance_bound_.0;
    let accept_upperbound:f64 = acceptance_bound_.1;

    //Dihedral が低いので値確かめ
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }
    let mut tmpatoms_prevpos:Vec<(f64,f64,f64)> =vec![];
    let mut tmpatoms_checkpoint:Vec<(f64,f64,f64)> =vec![];

    for aa in energyset.evoef2_env.md_envset.atoms.iter(){
        tmpatoms_prevpos.push((aa.get_x(),aa.get_y(),aa.get_z()));
        tmpatoms_checkpoint.push((aa.get_x(),aa.get_y(),aa.get_z()));
    }

    //peptide atoms というのと ligand atoms というのを作ると良さそうだ
    let (groups_,rot_max,mov_max):(Vec<Vec<usize>>,f64,f64) = 
    if let Some(x) = group_rotate {
        if x.1 > -0.5{
            //(make_sorted_residue_group(&energyset.evoef2_env.md_envset.atoms),x.1,movement_max)
            (x.0.clone(),x.1,movement_max)
        }else{
            (vec![],-1.0,-1.0)
        }
    }else{
        (vec![],-1.0,-1.0)
    };

    
    let groupmode =  rot_max > -0.5;

    let mut groups:Vec<Vec<usize>> = vec![];
    for gg in groups_.into_iter(){
        let mut fflag:bool = false;
        for ggg in gg.iter(){
            if atoms_fixed.contains(ggg){
                fflag = true;
                break;
            }
        }
        if !fflag{
            groups.push(gg);
        }
    }
    if groupmode && groups.len() == 0{
        eprintln!("All groups are fixed!");
        return (std::f64::INFINITY,0);
    }
    let num_atoms:usize = energyset.evoef2_env.md_envset.atoms.len();

    assert!(num_atoms > atoms_fixed.len());
    
    let mut checkpoint_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut prev_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut current_energies:Vec<f64> = vec![0.0;num_atoms];
    let mut mov_weight:Vec<f64> = vec![1.0;num_atoms];

    let mut selection_weight:Vec<f64> = if movatom_energysort {vec![1.0;num_atoms]}else{vec![]};
    
    let ignores:Vec<usize> = vec![];
    let mut floatings:Vec<usize> = vec![];
    //ToDo multi chain 対応//loop 中にも同じコードがあるので注意
    for ff in 0..energyset.evoef2_env.md_envset.atoms.len(){
        if atoms_fixed.contains(&ff){
            continue;
        }
        floatings.push(ff);
    }
    if floatings.len() == energyset.evoef2_env.md_envset.atoms.len(){
        let remover:usize = rgen.gen_range(0,energyset.evoef2_env.md_envset.atoms.len());
        floatings.remove(remover);
    }
    
    if groups.len() > 0{
        println!("fx floating can not be used with group.");
    }
    if bond_restrictor.len()+bond_upper_bound.len() > 0{
        if groups.len() == 0{
        assign_floating_atoms(
            &mut energyset.evoef2_env.md_envset.atoms
            ,bond_restrictor
            ,bond_upper_bound
            ,&floatings
            ,&ignores,0.2);
        }
    }
    energyset.update_distance();

    
    let mut minenergy_pos:Vec<(f64,f64,f64)> =vec![(0.0,0.0,0.0);num_atoms];
    let mut min_energy:f64;
    
    let mut prev_scoresum = energyset.calc_energy(&mut current_energies);
    for ii in 0..num_atoms{
        prev_energies[ii] = current_energies[ii];
    }
    
    min_energy = prev_scoresum;
    for ii in 0..num_atoms{
        tmpatoms_checkpoint[ii].0 = energyset.get_atom(ii).get_x();
        tmpatoms_checkpoint[ii].1 = energyset.get_atom(ii).get_y();
        tmpatoms_checkpoint[ii].2 = energyset.get_atom(ii).get_z();
        minenergy_pos[ii].0  = energyset.get_atom(ii).get_x();
        minenergy_pos[ii].1  = energyset.get_atom(ii).get_y();
        minenergy_pos[ii].2  = energyset.get_atom(ii).get_z();
        checkpoint_energies[ii] = current_energies[ii];
    }

    let mut current_iter:usize = 0;
    let mut prev_checked:usize = 0;//最終 iter は checkstep の倍数でないことがあるので current_iter - checkstep でなく前回のチェックポイントを覚えておく
    let mut num_accepted:f64 = 0.0;//スコアが小さかったが採択された回数
    let mut num_negative:f64 = 0.0;//スコアが小さくなり棄却するかどうかの試行が行われた回数

    let mut beta:f64 = 1.0;
    let wmax:f64 = 1000000.0;
    //groupatoms, fix_n, fix_c
    let mut moved_groups:Vec<(Vec<usize>,bool,bool)> = vec![];
    let mut movedgroup_prev_energies:Vec<f64> = vec![];
    loop{
        current_iter += 1;
        //let mov:f64 = rgen.gen_range(0.01,movement_max);
        let mov:f64 = movement_max;
        let mmin:f64 = prev_energies.iter().fold(std::f64::INFINITY,|s,m|s.min(*m));
        
        for aii in 0..num_atoms{
            mov_weight[aii] = 0.0;
        }
        if movatom_energysort{
            for aii in 0..num_atoms{
                selection_weight[aii] *= selection_weight_factor.0;
                if selection_weight[aii] > wmax  {
                    selection_weight[aii] = wmax;
                }
            }
            let mut sorter:Vec<(usize,f64)> = prev_energies.iter().enumerate().map(|m|(m.0
                , (*m.1-mmin)*selection_weight[m.0])).collect();
            sorter.sort_by(|a,b|a.1.partial_cmp(&b.1).unwrap_or_else(||panic!("error in sort!")));
            sorter.reverse();
            let mut scou:usize = 0;
            for ii in 0.. sorter.len(){
                if atoms_fixed.contains(&sorter[ii].0){
                    mov_weight[sorter[ii].0] = 0.0;
                    selection_weight[sorter[ii].0] = 1.0;
                }else{
                    if scou < num_moveatom{
                        let sii = sorter[ii].0;
                        if try_refine_lbfgs > 0{
                            println!("{}/{}:{}{}/{}",current_iter,ii,energyset.get_atom(sii).residue_name,energyset.get_atom(sii).residue_number,energyset.get_atom(sii).atom_name);
                            println!("weight:{} energy:{}",selection_weight[sorter[ii].0], sorter[ii].1);
                        }
                        //mov_weight[sorter[ii].0] = 1.0*(num_moveatom as f64-scou as f64)/(num_moveatom as f64);
                        mov_weight[sorter[ii].0] = 1.0;
                        for pp in 0..num_atoms{
                            if energyset.evoef2_env.md_envset.dist[sii][pp] < involution_factor.0
                            || energyset.evoef2_env.md_envset.num_edges[sii][pp] < involution_factor.1 {
                                mov_weight[pp] = 1.0;
                            }
                        }
                        scou += 1;
                        selection_weight[sii] *= selection_weight_factor.1;
                    }else{
                        break;
                    }
                }
            }
            for ii in 0..num_atoms{
                if atoms_fixed.contains(&ii){
                    mov_weight[ii] = 0.0;
                }
            }
        }else{
            let mut idxx:Vec<usize> = (0..num_atoms).into_iter().collect();
            idxx.shuffle(&mut rgen);
            let mut scou:usize = 0;
            for ii in 0.. idxx.len(){
                let sii = idxx[ii];
                if atoms_fixed.contains(&sii){
                    mov_weight[sii] = 0.0;
                }else{
                    if scou < num_moveatom{
                        mov_weight[sii] = 1.0;
                        for pp in 0..num_atoms{
                            if energyset.evoef2_env.md_envset.dist[sii][pp] < involution_factor.0
                            || energyset.evoef2_env.md_envset.num_edges[sii][pp] < involution_factor.1 {
                                mov_weight[pp] = 1.0;
                            }
                        }
                        scou += 1;
                    }else{
                        break;
                    }
                }
            }
            for ii in 0..num_atoms{
                if atoms_fixed.contains(&ii){
                    mov_weight[ii] = 0.0;
                }
            }
        }
        
        for ii in 0.. energyset.evoef2_env.md_envset.atoms.len(){
            tmpatoms_prevpos[ii].0 = energyset.get_atom(ii).get_x();
            tmpatoms_prevpos[ii].1 = energyset.get_atom(ii).get_y();
            tmpatoms_prevpos[ii].2 = energyset.get_atom(ii).get_z();
        }

        
        if rot_max > -0.5{//汚いので変更予定
            moved_groups.clear();
            for vv in groups.iter(){
                moved_groups.push((vv.clone(),false,false));
            }
            movedgroup_prev_energies.clear();
            for gg in moved_groups.iter(){
                let mut sc:f64 = 0.0;
                for aa in gg.0.iter(){
                    sc += prev_energies[*aa];
                }
                movedgroup_prev_energies.push(sc);
            }
            for mm in moved_groups.iter(){
                mov_group(&mut energyset.evoef2_env.md_envset.atoms,&mm.0,rot_max,mov_max,&mut rgen,mm.1,mm.2);
            }
        }else{
            random_movement(&mut energyset.evoef2_env.md_envset,mov,num_moveatom,&mut rgen,&mov_weight);
        }
        let ignores:Vec<usize> = vec![];
        let mut floatings:Vec<usize> = vec![];
        //ToDo multi chain 対応
        for ff in 0..energyset.evoef2_env.md_envset.atoms.len(){
            if atoms_fixed.contains(&ff){
                continue;
            }
            floatings.push(ff);
        }
        if floatings.len() == energyset.evoef2_env.md_envset.atoms.len(){
            let remover:usize = rgen.gen_range(0,energyset.evoef2_env.md_envset.atoms.len());
            floatings.remove(remover);
        }
        let fxfloating = |m:&mut Vec<charmm_based_energy::MDAtom>|{
            if groups.len() == 0{
                assign_floating_atoms(
                m
                ,bond_restrictor
                ,bond_upper_bound
                ,&floatings
                ,&ignores,0.2);
            }
        };
        fxfloating(&mut energyset.evoef2_env.md_envset.atoms);
        energyset.update_distance();
        

        if try_refine_lbfgs > 0{
            let mut patom:Vec<usize> = vec![];
            for ff in 0..num_atoms{
                if mov_weight[ff] > 0.0{
                    patom.push(ff);
                }
            }

            let bret = energy_lbfgs::run_lbfgs(
                energyset
                ,patom.clone()
                ,try_refine_lbfgs
                ,0.001
                ,5
                ,&None
            );
            let patomlen:usize = patom.len();
            for ii in 0..patomlen{
                let xyz:(f64,f64,f64) = energyset.get_atom(patom[ii]).get_xyz();
                energyset.get_atom_mut(patom[ii]).set_xyz(
                    bret.betas[ii*3]+xyz.0
                    ,bret.betas[ii*3+1]+xyz.1
                    ,bret.betas[ii*3+2]+xyz.2
                );
            }
            energyset.update_distance();
            println!("{} {}",patomlen,bret.betas.len());
        }

        for ii in 0..num_atoms{
            current_energies[ii] = 0.0;
        }
        let mut scoresum = energyset.calc_energy(&mut current_energies);
        let dloop = refine_retry;
        if dloop > 0{
            let ploop = if rot_max > -0.5{dloop}else{1};
            for pp in 0..ploop{
                if rot_max > -0.5{
                    //グループの再検討
                    let dratio = 1.0/(2.0_f64.powf(pp as f64+1.0));
                    let mut movedgroup_current_energies:Vec<f64> = vec![];
                    for gg in moved_groups.iter(){
                        let mut sc:f64 = 0.0;
                        for aa in gg.0.iter(){
                            sc += current_energies[*aa];
                        }
                        movedgroup_current_energies.push(sc);
                    }
                    for (mii,mm) in moved_groups.iter().enumerate(){
                        if movedgroup_current_energies[mii] > movedgroup_prev_energies[mii]{
                            for aa in mm.0.iter(){
                                energyset.get_atom_mut(*aa).set_xyz(
                                    tmpatoms_prevpos[*aa].0
                                    ,tmpatoms_prevpos[*aa].1
                                    ,tmpatoms_prevpos[*aa].2);
                            }
                            mov_group(&mut energyset.evoef2_env.md_envset.atoms,&mm.0,rot_max*dratio,mov_max*dratio,&mut rgen,mm.1,mm.2);
                        }
                    }
                    
                    fxfloating(&mut energyset.evoef2_env.md_envset.atoms);
                    energyset.update_distance();
                    for ii in 0..num_atoms{
                        current_energies[ii] = 0.0;
                    }
                    scoresum =  energyset.calc_energy(&mut current_energies);
                }else{
                    for dd in 0..dloop{//おそらくこれを入れた方が改善は早い
                        let dratio = 1.0/(2.0_f64.powf(dd as f64+1.0));
                        for ii in 0..num_atoms{
                            if current_energies[ii] > prev_energies[ii]{
                                if atoms_fixed.contains(&ii){
                                    continue;
                                }
                                if rot_max > -0.5{
                                    let tpos = energyset.get_atom(ii).get_xyz();
                                    if dd < dloop-1{
                                        energyset.get_atom_mut(ii).set_xyz(
                                        tpos.0+rgen.gen_range(mov*dratio*-1.0,mov*dratio)
                                        ,tpos.1+rgen.gen_range(mov*dratio*-1.0,mov*dratio)
                                        ,tpos.2+rgen.gen_range(mov*dratio*-1.0,mov*dratio));
                                    }else{
                                        energyset.get_atom_mut(ii).set_xyz(tpos.0,tpos.1,tpos.2);
                                    }
                                }else{
                                    if dd < dloop-1{
                                        //移動後の位置基準で動かしたら駄目だった
                                        energyset.get_atom_mut(ii).set_xyz(
                                        tmpatoms_prevpos[ii].0+rgen.gen_range(mov*dratio*-1.0,mov*dratio)
                                        ,tmpatoms_prevpos[ii].1+rgen.gen_range(mov*dratio*-1.0,mov*dratio)
                                        ,tmpatoms_prevpos[ii].2+rgen.gen_range(mov*dratio*-1.0,mov*dratio));
                                    }else{
                                        energyset.get_atom_mut(ii).set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                                    }
                                }
                            }
                        }
                    
                        fxfloating(&mut energyset.evoef2_env.md_envset.atoms);
                        energyset.update_distance();
                
                        for ii in 0..num_atoms{
                            current_energies[ii] = 0.0;
                        }

                        energyset.update_distance();
                        scoresum =  energyset.calc_energy(&mut current_energies);
                    }
                }
            }
        }
        if scoresum < min_energy{
            for ii in 0..num_atoms{
                minenergy_pos[ii].0 = energyset.get_atom(ii).get_x();
                minenergy_pos[ii].1 = energyset.get_atom(ii).get_y();
                minenergy_pos[ii].2 = energyset.get_atom(ii).get_z();
            }
            println!("min updated: {}->{} ",min_energy,scoresum);
            min_energy = scoresum;
            lastupdate = current_iter as i64;
        }
        if current_iter%100 == 0{
            println!("iter:{} {} {} {}",current_iter,scoresum,prev_scoresum,min_energy);
        }
        if scoresum < prev_scoresum{
            //accepted
            prev_scoresum = scoresum;
            for ii in 0..num_atoms{
                prev_energies[ii] = current_energies[ii];
            }
            
        }else{
            let dd:f64 = rgen.gen_range(0.0,1.0);
            let bb:f64 = ((scoresum-prev_scoresum)*-1.0*beta).exp();
            if prev_scoresum < scoresum{//同値ではない
                num_negative += 1.0;
                if accept_upperbound > 0.0 && dd < bb{
                    num_accepted += 1.0;
                    prev_scoresum = scoresum;
                    for ii in 0..num_atoms{
                        prev_energies[ii] = current_energies[ii];
                    }
                    
                }else{
                    for ii in 0..num_atoms{
                        energyset.get_atom_mut(ii).set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                    }
                }
            }
        }
        if accept_upperbound > 0.0{
            let checkflag =  current_iter%checkstep == 0 || current_iter == iternum;
            if checkflag{
                //println!("accepted:{} rejected:{}",num_accepted,num_negative);
                if num_negative > 0.0 && num_accepted/num_negative > accept_upperbound{
                    //採択されすぎなので戻す
                    beta *= 5.0;
                    current_iter = prev_checked;
                    for ii in 0..num_atoms{
                        energyset.get_atom_mut(ii).set_xyz(tmpatoms_checkpoint[ii].0,tmpatoms_checkpoint[ii].1,tmpatoms_checkpoint[ii].2);
                        prev_energies[ii] = checkpoint_energies[ii];
                    }
                }else{
                    if num_negative > 0.0 && num_accepted/num_negative < accept_lowerbound{
                        //採択されなさすぎなので採択率を上げる
                        beta /= 5.0;
                    }
                    prev_checked = current_iter;
                    for ii in 0..num_atoms{
                        tmpatoms_checkpoint[ii].0 = energyset.get_atom(ii).get_x();
                        tmpatoms_checkpoint[ii].1 = energyset.get_atom(ii).get_y();
                        tmpatoms_checkpoint[ii].2 = energyset.get_atom(ii).get_z();
                        checkpoint_energies[ii] = current_energies[ii];
                    }
                }
                num_negative = 0.0;
                num_accepted = 0.0;
            }
        }
        if current_iter >= iternum{
            break;
        }
    }
    
    for ii in 0..num_atoms{
        energyset.get_atom_mut(ii).set_xyz(minenergy_pos[ii].0,minenergy_pos[ii].1,minenergy_pos[ii].2);
    }
    energyset.update_distance();
    return (min_energy,lastupdate);
}


pub fn get_fixed_point(fl:usize
    ,fixed_idx_attractorflag:&Vec<(usize,bool,f64)>
    ,atoms:&Vec<charmm_based_energy::MDAtom>
    ,ffixed:&HashSet<usize>
    ,tolerance:f64)->((f64,f64,f64),bool){
    let mut next_pos:Vec<(f64,f64,f64)> = vec![];
    let mut skipcount:usize = 0;
    for (fx,attractor,b0) in fixed_idx_attractorflag.iter(){
        if ! ffixed.contains(fx){
            continue;
        }
        let fxpos:(f64,f64,f64) = atoms[*fx].get_xyz();
        let mut flpos:(f64,f64,f64) = atoms[fl].get_xyz();
        
        let ddis:f64 = process_3d::distance(&fxpos,&flpos);
        let mut stayflag:bool = false;
        if ddis == 0.0{
            skipcount += 1;
            stayflag = true;
        }else if *attractor && ddis - tolerance < *b0{
            skipcount += 1;
            stayflag = true;
        }else if (ddis - b0).abs() < tolerance{
            skipcount += 1;
        }
        if !stayflag{
            flpos.0 -= fxpos.0;
            flpos.1 -= fxpos.1;
            flpos.2 -= fxpos.2;

            flpos.0 /= ddis;
            flpos.1 /= ddis;
            flpos.2 /= ddis;
            
            flpos.0 *= b0;
            flpos.1 *= b0;
            flpos.2 *= b0;
            
            flpos.0 += fxpos.0;
            flpos.1 += fxpos.1;
            flpos.2 += fxpos.2;
        }
        next_pos.push(flpos);
    }
    
    if next_pos.len() == skipcount {
        return (atoms[fl].get_xyz(),false);
    }

    let mut x:f64 = 0.0;
    let mut y:f64 = 0.0;
    let mut z:f64 = 0.0;

    for nn in next_pos.iter(){
        x += nn.0;
        y += nn.1;
        z += nn.2;
    }
    
    x /= next_pos.len() as f64;
    y /= next_pos.len() as f64;
    z /= next_pos.len() as f64;

    return ((x,y,z),true);
}


//他 group と結合しているのは、最初の原子か最後の原子である必要がある。
pub fn fix_group_bond(
    atoms:&mut Vec<charmm_based_energy::MDAtom>
    ,group:&Vec<Vec<usize>>
    ,restrictor_:&Vec<charmm_based_energy::BondVars>
    ,attractor_:&Vec<charmm_based_energy::BondVars>
    ,ffloating_:&Vec<usize>
    ,iignore_:&Vec<usize>
    ,tolerance:f64
    ,max_iter:usize
){

    let mut ffloating:HashSet<usize> = HashSet::new();
    let ffloating_all:HashSet<usize> = ffloating_.iter().map(|m|*m).collect();

    let mut terminal_atoms:HashSet<usize> = HashSet::new();
    //groupindex, index_in_vec // floating なグループの最初の ATOM と最後の ATOM だけ入っている
    let mut terminal_to_group:HashMap<usize,Vec<(usize,usize)>> = HashMap::new();
    //atom_index, attractflag, dist // この Vec には TERMINAL 以外の ATOM も入っている
    let mut terminal_to_pair:HashMap<usize,Vec<(usize,bool,f64)>> = HashMap::new();
    
    //let mut atom_to_group:HashMap<usize,Vec<usize>> = HashMap::new();
    'outer:for (gii,gg) in group.iter().enumerate(){
        for aa  in gg.iter(){
            if !ffloating_all.contains(aa){//floating と fixed が混ざっているやつは動かさない
                continue 'outer;
            }
        }
        
        /*for (gii,aa)  in gg.iter().enumerate(){//同一グループ内の BOND をフィルタしようと思ったがやめた。引数の段階から与えない
            if !atom_to_group.contains_key(aa){
                atom_to_group.insert(*aa,vec![]);
            }
            atom_to_group.get_mut(aa).unwrap().push(gii);
        }*/

        terminal_atoms.insert(gg[0]);
        if !terminal_to_group.contains_key(&gg[0]){
            terminal_to_group.insert(gg[0],vec![]);
            terminal_to_pair.insert(gg[0],vec![]);
        }
        terminal_to_group.get_mut(&gg[0]).unwrap().push((gii,0));

        if gg[0] != gg[gg.len()-1]{
            terminal_atoms.insert(gg[gg.len()-1]);
            if !terminal_to_group.contains_key(&gg[gg.len()-1]){
                terminal_to_group.insert(gg[gg.len()-1],vec![]);
                terminal_to_pair.insert(gg[gg.len()-1],vec![]);
            }
            terminal_to_group.get_mut(&gg[gg.len()-1]).unwrap().push((gii,gg.len()-1));
        }
    }
    
    for ff in ffloating_.iter(){
        if terminal_atoms.contains(ff){
            ffloating.insert(*ff);
        }
    }
    let mut iignore:HashSet<usize> = HashSet::new();
    for ii in iignore_.iter(){
        if terminal_atoms.contains(ii){
            iignore.insert(*ii);
        }
    }

    let mut ffixed:HashSet<usize> = HashSet::new();
 
    for tt in 0..atoms.len(){
        if !iignore.contains(&tt) && !ffloating_all.contains(&tt){
            ffixed.insert(tt);
        }
    }

    let mut restrictor:Vec<&charmm_based_energy::BondVars> = vec![];
    for rr in restrictor_.iter(){
        if terminal_atoms.contains(&rr.atoms.0) || terminal_atoms.contains(&rr.atoms.1) {
            restrictor.push(rr);
        }
        if terminal_atoms.contains(&rr.atoms.0){//グループの中間の ATOM と結合することもある
            terminal_to_pair.get_mut(&rr.atoms.0).unwrap().push((rr.atoms.1,false,rr.b0));
        }
        if terminal_atoms.contains(&rr.atoms.1){
            terminal_to_pair.get_mut(&rr.atoms.1).unwrap().push((rr.atoms.0,false,rr.b0));
        }
    }
    
    let mut attractor:Vec<&charmm_based_energy::BondVars> = vec![];
    for aa in attractor_.iter(){
        if terminal_atoms.contains(&aa.atoms.0)
        || terminal_atoms.contains(&aa.atoms.1)
        {
            attractor.push(aa);
        }
        if terminal_atoms.contains(&aa.atoms.0){
            if terminal_atoms.contains(&aa.atoms.1) || ffixed.contains(&aa.atoms.1) {
                terminal_to_pair.get_mut(&aa.atoms.0).unwrap().push((aa.atoms.1,true,aa.b0));
            }
        }
        if terminal_atoms.contains(&aa.atoms.1){
            if terminal_atoms.contains(&aa.atoms.0) || ffixed.contains(&aa.atoms.0) {
                terminal_to_pair.get_mut(&aa.atoms.1).unwrap().push((aa.atoms.0,true,aa.b0));
            }
        }
    }


    let mut updated_groups_:HashSet<usize> = HashSet::new();
    for ff in ffloating.iter(){
        updated_groups_.extend::<Vec<usize>>(terminal_to_group.get(ff).unwrap().iter().map(|m|m.0).collect());
    }
    let mut updated_groups:Vec<usize> = updated_groups_.into_iter().collect();
    updated_groups.sort();
    /* Fixed された Group に間接的にも一つも繋がっていない Group は適当に採択されるようにしようかと思ったがやめた
    必ず間接的にはつながっているべき。
    let mut groups_checked_flag:Vec<bool> = vec![true;group.len()];
    for uu in updated_groups.iter(){
        groups_checked_flag[*uu] = false;
    }
    */
    let mut breakcount:usize = 0;
    'dloop:loop{
        let mut updated_groups_next:Vec<usize> =vec![];
        let mut ffixed_next:HashSet<usize> = HashSet::new();//今のところ一回目のループにしか効かないと思う
        while updated_groups.len() > 0{
            let mut nextindex_:i64 = -1;
            let mut to_fix:(bool,bool) = (false,false);
            'tloop:for ii in 0..4{
                for (uii,uu) in updated_groups.iter().enumerate(){
                    //前方の ATOM の位置が決定できるか
                    let starter:bool = ii < 2;
                    let ok0:bool = if terminal_to_pair.get(&group[*uu][0]).unwrap().len() > 0 {
                        terminal_to_pair.get(&((group[*uu])[0])).unwrap().iter().fold(starter,|s,m| 
                            if ii < 2 {
                                s && ffixed.contains(&(m.0))
                            }else{
                                s || ffixed.contains(&(m.0))
                            }
                        )
                    }else{
                        false
                    };
                    //後方の ATOM の位置が決定できるか
                    let oklast:bool = if terminal_to_pair.get(&group[*uu][group[*uu].len()-1]).unwrap().len() > 0 {
                        terminal_to_pair.get(&((group[*uu])[group[*uu].len()-1])).unwrap().iter().fold(starter,|s,m| 
                            if ii < 2{
                                s && ffixed.contains(&(m.0))
                            }else{
                                s || ffixed.contains(&(m.0))
                            })
                    }else{
                        false
                    };
                    if ok0 && oklast{
                        nextindex_ = uii as i64;
                        to_fix = (ok0,oklast);
                        break 'tloop;
                    }
                    if (ii == 1 || ii ==3) && (ok0 || oklast){
                        nextindex_ = uii as i64;
                        to_fix = (ok0,oklast);
                        break 'tloop;
                    }
                }
            }

            let target_group:usize;
            let mut donmovflag:bool = false;
            if nextindex_ < 0{
                updated_groups.clear();
                break;
            }else{
                target_group = updated_groups.remove(nextindex_ as usize);

                let mut pos0:(f64,f64,f64) = (0.0,0.0,0.0);
                let mut poslast:(f64,f64,f64) = (0.0,0.0,0.0);
                let mut success0:bool = false;
                let mut successlast:bool = false;
                let mut ifixed:usize = 0;

                if to_fix.0{
                    let pres = get_fixed_point(
                        group[target_group][0]
                        , terminal_to_pair.get(&(group[target_group][0])).as_ref().unwrap()
                        , atoms
                        ,&ffixed
                        , tolerance
                    );
                    pos0 = pres.0;
                    success0 = pres.1;
                    ifixed += 1;
                }
                
                if to_fix.1{
                    let pres = get_fixed_point(
                            group[target_group][group[target_group].len()-1]
                            , terminal_to_pair.get(&(group[target_group][group[target_group].len()-1])).as_ref().unwrap()
                            , atoms
                            ,&ffixed
                            , tolerance
                        );
                        
                    poslast = pres.0;
                    successlast = pres.1;
                    ifixed += 1;
                }
                if ifixed == 2 && group[target_group].len() != 1{//両方が FIX できて一点でない場合回転もする
                    if !successlast && !success0{
                        donmovflag = true;
                    }else{
                        let mut nexx_:Vec<Point3D> = vec![];
                        for pp in group[target_group].iter(){
                            nexx_.push(Point3D::new(atoms[*pp].get_x(),atoms[*pp].get_y(),atoms[*pp].get_z()));
                        }
                        let mut nexx:Vec<&mut dyn Vector3D> = vec![];
                        for nn in nexx_.iter_mut(){
                            nexx.push(nn);
                        }
                        process_3d::fit_to_vector_2(
                            &atoms[(group[target_group])[0]]
                            ,&atoms[group[target_group][group[target_group].len()-1]]
                            ,&Point3D::new(pos0.0, pos0.1, pos0.2)
                            ,&Point3D::new(poslast.0, poslast.1, poslast.2)
                            , &mut nexx);

                        //上の関数は開始点が同一座標になるので、中心点を合わせる
                        let center1:(f64,f64,f64) = (
                            pos0.0/2.0+poslast.0/2.0
                            ,pos0.1/2.0+poslast.1/2.0
                            ,pos0.2/2.0+poslast.2/2.0
                        );
                        let pos0b:(f64,f64,f64) = nexx[0].get_xyz();
                        let poslastb:(f64,f64,f64) = nexx[nexx.len()-1].get_xyz();
                        let center2:(f64,f64,f64) = (
                            pos0b.0/2.0+poslastb.0/2.0
                            ,pos0b.1/2.0+poslastb.1/2.0
                            ,pos0b.2/2.0+poslastb.2/2.0
                        );
                        let mov:(f64,f64,f64) = (center1.0-center2.0,center1.1-center2.1,center1.2-center2.2);
                        if process_3d::distance(
                            &(nexx[0].get_x() + mov.0, nexx[0].get_y() + mov.1, nexx[0].get_z() + mov.2)
                            ,&atoms[group[target_group][0]].get_xyz()
                        ) + process_3d::distance(
                            &(nexx[nexx.len()-1].get_x() + mov.0, nexx[nexx.len()-1].get_y() + mov.1, nexx[nexx.len()-1].get_z() + mov.2)
                            ,&atoms[group[target_group][group[target_group].len()-1]].get_xyz()
                        ) < tolerance{
                            //拘束条件が複数あり、いくつかの拘束条件に当てはまらない場合に対応する
                            donmovflag = true;
                        }else {
                            for (pii,pp) in group[target_group].iter().enumerate(){
                                atoms[*pp].set_xyz(
                                    nexx[pii].get_x()+mov.0,
                                    nexx[pii].get_y()+mov.1,
                                    nexx[pii].get_z()+mov.2
                                );
                            }
                        }
                    }
                }else{
                    let mut mov:(f64,f64,f64) = (0.0,0.0,0.0);
                    if to_fix.0{
                        if !success0{
                            donmovflag = true;
                        }else{
                            let bpp:(f64,f64,f64) = atoms[group[target_group][0]].get_xyz();
                            mov = (pos0.0-bpp.0,pos0.1-bpp.1,pos0.2-bpp.2);
                        }
                    }else if to_fix.1{
                        if !successlast{
                            donmovflag = true;
                        }else{
                            let bpp:(f64,f64,f64) = atoms[group[target_group][group[target_group].len()-1]].get_xyz();
                            mov = (poslast.0-bpp.0,poslast.1-bpp.1,poslast.2-bpp.2);
                        }
                    }else{
                        panic!("???");
                    }
    
                    if !donmovflag{
                        for pp in group[target_group].iter(){
                            let p = atoms[*pp].get_xyz();
                            atoms[*pp].set_xyz(
                                p.0+mov.0,
                                p.1+mov.1,
                                p.2+mov.2
                            );
                        }
                    }
                }
                for tt in group[target_group].iter(){
                    if terminal_atoms.contains(tt) && ffloating.contains(tt){
                        ffixed_next.insert(*tt);
                    }
                }
            

                if !donmovflag{
                    let mut dupcheck:HashSet<usize> = updated_groups_next.iter().map(|m|*m).collect();
                    for pp in group[target_group].iter(){
                        if terminal_to_pair.contains_key(pp){
                            let vvec:&Vec<(usize,bool,f64)> = terminal_to_pair.get(pp).unwrap();
                            for vv in vvec.iter(){
                                if terminal_to_group.contains_key(&vv.0){
                                    for tt in terminal_to_group.get(&vv.0).unwrap().iter(){
                                        if !dupcheck.contains(&tt.0) {
                                            if tt.0 == target_group{
                                                continue;
                                            }
                                            updated_groups_next.push(tt.0);
                                            dupcheck.insert(tt.0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            breakcount += 1;
            if breakcount >= max_iter{
                break 'dloop;
            }
        }
        if updated_groups_next.len() == 0{
            break;
        }
        updated_groups = updated_groups_next;
        ffixed.extend(ffixed_next.into_iter());
    }
    //println!("{}/{}******",breakcount,max_iter);
}


pub fn assign_floating_atoms(
    atoms:&mut Vec<charmm_based_energy::MDAtom>
    ,restrictor:&Vec<charmm_based_energy::BondVars>
    ,attractor:&Vec<charmm_based_energy::BondVars>
    ,ffloating_:&Vec<usize>
    ,iignore_:&Vec<usize>
    ,tolerance:f64
){
    let mut ffixed:HashSet<usize> = HashSet::new();
    let mut ffloating:HashSet<usize> = ffloating_.iter().map(|m|*m).collect();
    let iignore:HashSet<usize> = iignore_.iter().map(|m|*m).collect();
    let anum:usize =atoms.len();
    for aa in 0..anum{
        if iignore.contains(&aa){
            continue;
        }
        if ffloating.contains(&aa){
            continue;
        }
        ffixed.insert(aa);
    }

    let mut bond_res_checked:Vec<bool> = vec![false;restrictor.len()];
    let mut bond_att_checked:Vec<bool> = vec![false;attractor.len()];

    loop{
        
        let mut updated:bool = false;
        let mut next_pos:HashMap<usize,Vec<((f64,f64,f64),bool)>> = HashMap::new();
        let bondvec:Vec<(&Vec<charmm_based_energy::BondVars>,bool)> = vec![(attractor,true),(restrictor,false)];
        for (bvec,attrac) in bondvec.into_iter(){
            for (bii,bb) in bvec.iter().enumerate(){
                if (!attrac && bond_res_checked[bii]) || (attrac && bond_att_checked[bii]){
                    continue;
                }
                let mut fx:i64 = -1;
                let mut fl:i64 = -1;
                if ffloating.contains(&bb.atoms.1) && ffixed.contains(&bb.atoms.0) { 
                    fx = bb.atoms.0 as i64;
                    fl = bb.atoms.1 as i64;
                }
                if ffloating.contains(&bb.atoms.0) && ffixed.contains(&bb.atoms.1) { 
                    fx = bb.atoms.1 as i64;
                    fl = bb.atoms.0 as i64;
                }
                if fx < 0{
                    if ffixed.contains(&bb.atoms.0) && ffixed.contains(&bb.atoms.1) {
                        if !attrac {
                            bond_res_checked[bii] = true;
                        }else{
                            bond_att_checked[bii] = true;
                        }
                    }
                    continue;
                }
                let fxpos:(f64,f64,f64) = atoms[fx as usize].get_xyz();
                let mut flpos:(f64,f64,f64) = atoms[fl as usize].get_xyz();
                
                let ddis:f64 = process_3d::distance(&fxpos,&flpos);
                let mut stayflag:bool = false;
                let skipped:bool;
                
                if ddis == 0.0{//ランダムで動かすか？
                    stayflag = true;
                    skipped = true;
                }else if attrac && ddis-tolerance < bb.b0{
                    stayflag = true;
                    skipped = true;
                }else if (ddis - bb.b0).abs() < tolerance{
                    stayflag = true;
                    skipped = true;
                }else{
                    skipped = false;
                }
                if !stayflag{
                    flpos.0 -= fxpos.0;
                    flpos.1 -= fxpos.1;
                    flpos.2 -= fxpos.2;

                    flpos.0 /= ddis;
                    flpos.1 /= ddis;
                    flpos.2 /= ddis;
                    if ddis > bb.b0{
                        flpos.0 *= bb.b0+tolerance;
                        flpos.1 *= bb.b0+tolerance;
                        flpos.2 *= bb.b0+tolerance;
                    }else{
                        flpos.0 *= bb.b0-tolerance;
                        flpos.1 *= bb.b0-tolerance;
                        flpos.2 *= bb.b0-tolerance;

                    }
                    
                    flpos.0 += fxpos.0;
                    flpos.1 += fxpos.1;
                    flpos.2 += fxpos.2;
                }
                if !next_pos.contains_key(&(fl as usize)){
                    next_pos.insert(fl as usize,vec![]);
                }
                next_pos.get_mut(&(fl as usize)).unwrap().push((flpos,skipped));
                
                updated = true;
            }
            for (kk,vv) in next_pos.iter(){
                let skipper:bool = vv.iter().fold(true,|s,m|s && m.1);
                if !skipper{
                    let mut xx:f64 = 0.0;
                    let mut yy:f64 = 0.0;
                    let mut zz:f64 = 0.0;
                    for (vvv,_flag) in vv.iter(){
                        xx += vvv.0;
                        yy += vvv.1;
                        zz += vvv.2;
                    }
                    xx /= vv.len() as f64;
                    yy /= vv.len() as f64;
                    zz /= vv.len() as f64;
                    atoms[*kk].set_xyz(xx,yy,zz);
                }
                ffloating.remove(kk);
                ffixed.insert(*kk);
            }
        }
        
        if !updated{
            break;
        }
    }
}


pub fn atom_nelder_mead(md_envset:&mut charmm_based_energy::CharmmEnv
    ,energyset:&PPEnergySet
    ,target_atom_index:Vec<usize>
    ,iternum:usize
    )->Vec<f64>{
    let num_variables:usize = target_atom_index.len()*3;
    let num_points:usize = num_variables+1;
    let mut ranks:Vec<usize> = (0..num_points).into_iter().collect();
    let mut betas:Vec<(Vec<f64>,f64)> = vec![(vec![0.0;num_variables],0.0);num_points];

    let alpha:f64 = 1.0;
    let gamma:f64 = 2.0;
    let rho:f64 = 0.5;
    let sigma:f64 = 0.5;
    let lambda:f64 = 1.0;//高めの方が良いかも
    
    for (pp,_tt) in target_atom_index.iter().enumerate(){
        betas[pp*3].0[pp*3] = lambda;
        betas[pp*3+1].0[pp*3+1] = lambda;
        betas[pp*3+2].0[pp*3+2] = lambda;
    }

    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let mut calc_energy_tmp = |tmpbeta:&Vec<f64>|{
        let mut prevval:Vec<f64> = vec![0.0;num_variables];
        for (ii,tt) in target_atom_index.iter().enumerate(){
            prevval[ii*3] = md_envset.atoms[*tt].get_x();
            prevval[ii*3+1] = md_envset.atoms[*tt].get_y();
            prevval[ii*3+2] = md_envset.atoms[*tt].get_z();
            md_envset.atoms[*tt].set_xyz(
                tmpbeta[ii*3]+prevval[ii*3],
                tmpbeta[ii*3+1]+prevval[ii*3+1],
                tmpbeta[ii*3+2]+prevval[ii*3+2]
            );
        }
        md_envset.update_distance();
        let ret = energyset.calc_energy(&mut dummy);
        for (ii,tt) in target_atom_index.iter().enumerate(){
            md_envset.atoms[*tt].set_xyz(
                prevval[ii*3],
                prevval[ii*3+1],
                prevval[ii*3+2]
            );
        }
        ret
    };
    for pp in 0..num_points{
        let benergy:f64 = calc_energy_tmp(&betas[pp].0);
        betas[pp].1 = benergy;
    }
    let mut prevenergy:f64 = std::f64::INFINITY;
    let mut no_improve_limit_count:usize = 0;
    let no_improve_limit:usize = 10000;
    for _i in 0..iternum{
        ranks.sort_by(|a,b|betas[*a].1.partial_cmp(&betas[*b].1).unwrap());
        if prevenergy-0.001 <= betas[ranks[0]].1{//改善が数回見られなかった場合 break に変更
            no_improve_limit_count += 1;
        }else{
            no_improve_limit_count = 0;
        }
        if no_improve_limit <= no_improve_limit_count{
            break;
        }
        println!("current: {} {}",ranks[0],betas[ranks[0]].1);
        prevenergy = betas[ranks[0]].1;
        let mut xcenter:Vec<f64> = vec![0.0;num_points];
        for pp in 0..num_points-1{
            for vv in 0..num_variables{
                xcenter[vv] += betas[pp].0[vv];
            }
        }
        for vv in 0..num_variables{
            xcenter[vv] /= num_points as f64 -1.0;
        }

        let mut xnex:Vec<f64> = vec![0.0;num_variables];
        for ii in 0..num_variables{
            xnex[ii] = xcenter[ii]+alpha*(xcenter[ii]-betas[ranks[num_points-1]].0[ii]);
        }

        let snex:f64 = calc_energy_tmp(&xnex);

        if snex >= betas[ranks[0]].1 && snex < betas[ranks[num_points-2]].1 {
            betas[ranks[num_points-1]] = (xnex,snex);
        }else if snex < betas[ranks[0]].1{
            let mut pnex:Vec<f64> = vec![0.0;num_variables];
            for ii in 0..num_variables{
                pnex[ii] = xcenter[ii]+gamma*(betas[ranks[num_points-1]].0[ii]-xcenter[ii]);
            }
            let spnex:f64 = calc_energy_tmp(&pnex);
            if spnex < snex{
                betas[ranks[num_points-1]] = (pnex,spnex);
            }else{
                betas[ranks[num_points-1]] = (xnex,snex);
            }
        }else{
            let mut pnex:Vec<f64> = vec![0.0;num_variables];
            for ii in 0..num_variables{
                pnex[ii] = xcenter[ii]+rho*(betas[ranks[num_points-1]].0[ii]-xcenter[ii]);
            }
            let spnex:f64 = calc_energy_tmp(&pnex);
            if spnex < betas[ranks[num_points-1]].1 {
                betas[ranks[num_points-1]] =  (pnex,spnex);
            }else{
                for pp in 0..num_points{
                    if ranks[0] == pp{
                        continue;
                    }
                    for vv in 0..num_variables{
                        betas[pp].0[vv] = betas[ranks[0]].0[vv]+sigma*(betas[pp].0[vv]-betas[ranks[0]].0[vv]);
                    }
                    betas[pp].1 = calc_energy_tmp(&betas[pp].0);
                }
            }
        }
    }
    return betas[ranks[0]].0.clone();
}

pub fn make_straight_split_group(num_atoms:usize)
->Vec<Vec<usize>>{
    let mut ret:Vec<Vec<usize>> = vec![];
    for ii in 0..num_atoms{
        let mut pret:Vec<usize> = vec![];
        if ii >= num_atoms/2{
            for jj in ii..num_atoms{
                pret.push(jj);
            }
        }else{
            for jj in 0..ii{
                pret.push(jj);
            }
        }
        if pret.len() > 1{
            ret.push(pret);
        }
    }
    return ret;
}
//Peptide 結合からなる Chain について、Bond 回転できるようにグループを作る
pub fn make_peptide_atom_groups(atoms_:&Vec<charmm_based_energy::MDAtom>,bonds:&Vec<&charmm_based_energy::BondVars>)
->Vec<Vec<usize>>{
    let mut chains:HashMap<String,Vec<usize>> = HashMap::new();
    for (aii,aa) in atoms_.iter().enumerate(){
        if !chains.contains_key(&aa.chain_name){
            chains.insert(aa.chain_name.clone(),vec![]);
        }
        chains.get_mut(&aa.chain_name).unwrap().push(aii);
    }
    let num_all_atoms:usize = atoms_.len();
    let mut chainnames:Vec<String> = chains.iter().map(|m|m.0.to_string()).collect();
    chainnames.sort();

    let mut ret:Vec<Vec<usize>> = vec![];

    for chainname in chainnames.iter(){
        let mut atoms:Vec<usize> = chains.get(chainname).unwrap().clone();
        atoms.sort_by(|a,b|{atoms_[*a].residue_index_in_chain.cmp(&atoms_[*b].residue_index_in_chain)});
        let mut n_:i64 = -1;
        let mut ca_:i64 = -1;
        let mut c_:i64 = -1;
        let num_atoms:usize = atoms.len();        
        for pp in vec!["N","CA","C"].into_iter().enumerate(){
            for ii in 0..num_atoms{
                if atoms_[atoms[ii]].atom_name == pp.1{
                    if pp.1 == "N"{
                        n_ = atoms[ii] as i64;
                    }
                    if pp.1 == "CA"{
                        ca_ = atoms[ii] as i64;
                    }
                    if pp.1 == "C"{
                        c_ = atoms[ii] as i64;
                    }
                    break;
                }
            }
        }

        if n_ < 0 && ca_ < 0 && c_ < 0{
            eprintln!("Chain {} does not have backbone atom!",chainname);
            continue;
        }
        let startatom:usize = if n_ > -1{
            n_ as usize
        }else if ca_ > -1{
            ca_ as usize
        }else{
            c_ as usize
        };

       
        let mut num_edges:HashMap<usize,usize> = HashMap::new();
        let mut edges:HashMap<usize,Vec<usize>> = HashMap::new();
        for aa in 0..atoms.len(){
            edges.insert(atoms[aa],vec![]);
            num_edges.insert(atoms[aa],num_atoms);
        }

        let atoms_in_chain:HashSet<usize> = atoms.into_iter().collect();

        for bb in bonds.iter(){
            if atoms_in_chain.contains(&bb.atoms.0)
            &&  atoms_in_chain.contains(&bb.atoms.1){
                edges.get_mut(&bb.atoms.0).unwrap().push(bb.atoms.1);
                edges.get_mut(&bb.atoms.1).unwrap().push(bb.atoms.0);
            }
        }

        num_edges.insert(startatom,0);
        let mut atoms_updated:Vec<usize> = vec![];
        //atomindex edgenumfromstart
        let mut backbones_idx_edgenum:Vec<(usize,usize)> = vec![];
        backbones_idx_edgenum.push((startatom,0));
        //start atom からの距離
        let mut allatoms_dist:Vec<usize> = vec![num_atoms;num_all_atoms];

        let mut backbones_hs:HashSet<usize> = HashSet::new();
        atoms_updated.push(startatom);
        allatoms_dist[startatom] = 0;

        while atoms_updated.len() > 0{
            let current:usize = atoms_updated.remove(0);
            let currentnum:usize = *num_edges.get(&current).unwrap();
            assert_eq!(currentnum,allatoms_dist[current]);//問題ないようなら num_edges は削除
            for nn in edges.get(&current).unwrap().iter(){
                if *num_edges.get(nn).unwrap() == num_atoms{
                    num_edges.insert(*nn,currentnum+1);
                    atoms_updated.push(*nn);
                    if atoms_[*nn].atom_name == "CA"
                    || atoms_[*nn].atom_name == "C"
                    || atoms_[*nn].atom_name == "N"{
                        if !backbones_hs.contains(nn){
                            backbones_idx_edgenum.push((*nn,currentnum+1));
                            backbones_hs.insert(*nn);
                        }
                    }
                    if allatoms_dist[*nn] == num_atoms{
                        allatoms_dist[*nn] = currentnum+1;
                    }
                }
            }
        }
        
        backbones_idx_edgenum.sort_by(|a,b|a.1.cmp(&b.1));
        let mut side_atoms:HashMap<usize,Vec<usize>> =HashMap::new();//sidechain + o とか

        for bb in backbones_idx_edgenum.iter(){
            //backbone に結合している Atom を Vec に入れる
            //Sidechain だけでなく O や H も入っている
            let mut atoms_updated:Vec<usize> = vec![];
            atoms_updated.push(bb.0);
            let mut checked:HashSet<usize> = HashSet::new();
            let mut vvec:Vec<usize> = vec![];//sort されている方を後で使うため
            while atoms_updated.len() > 0{
                let current:usize = atoms_updated.remove(0);
                for nn in edges.get(&current).unwrap().iter(){
                    if atoms_[*nn].atom_name != "CA"
                    && atoms_[*nn].atom_name != "C"
                    && atoms_[*nn].atom_name != "N"
                    && !checked.contains(nn){
                        checked.insert(*nn);
                        vvec.push(*nn);
                        atoms_updated.push(*nn);
                    }
                }
            }
            side_atoms.insert(bb.0,vvec);
        }

        //Backbone の結合を曲げられるようなグループを作る
        let mut groups_all:Vec<Vec<usize>> = vec![];
        let num_backbones:usize = backbones_idx_edgenum.len();
        for ii in 0..num_backbones{
            let mut groupatoms:Vec<usize> = vec![];
            if ii < num_backbones/2{
                for jj in 0..=ii{
                    groupatoms.push(backbones_idx_edgenum[jj].0);
                    if side_atoms.contains_key(&(backbones_idx_edgenum[jj].0)){//Side についている Atom も全部入れる
                        groupatoms.extend(side_atoms.get(&(backbones_idx_edgenum[jj].0)).unwrap().clone());
                    }
                }
            }else{
                for jj in ii..num_backbones{
                    groupatoms.push(backbones_idx_edgenum[jj].0);
                    if side_atoms.contains_key(&(backbones_idx_edgenum[jj].0)){//Side についている Atom も全部入れる
                        groupatoms.extend(side_atoms.get(&(backbones_idx_edgenum[jj].0)).unwrap().clone());
                    }
                }
            }
            if groupatoms.len() > 1{
                groupatoms.sort();
                groups_all.push(groupatoms);
            }
        }
        

        for ii in 0..num_backbones{//sidechain と O や H を入れる
            if side_atoms.contains_key(&(backbones_idx_edgenum[ii].0)){
                let backboneidx:usize = backbones_idx_edgenum[ii].0;
                let mut side_all_hs:HashSet<usize> = side_atoms.get(&backboneidx).unwrap().iter().map(|m|*m).collect();
                side_all_hs.insert(backboneidx);
                
                let cgg_:Vec<usize> = edges.get(&backboneidx).unwrap().clone();
                //atom idx, connected_with_ca
                let mut cgg:Vec<(usize,bool)> = vec![];
                for nn in cgg_.iter(){
                    if atoms_[*nn].atom_name != "CA"
                    && atoms_[*nn].atom_name != "C"
                    && atoms_[*nn].atom_name != "N"{//BACKBONEATOM 方向からの確認
                        cgg.push((*nn,true));
                    }
                }

                let mut sidegroup:Vec<Vec<usize>> = vec![];
                let cgg_:Vec<usize> = side_atoms.get(&backboneidx).unwrap().clone();
                for nn in cgg_.iter(){
                    let kcc:&Vec<usize> = edges.get(nn).unwrap();
                    let base = allatoms_dist[*nn];
                    if kcc.len() == 1{
                        assert!(base > allatoms_dist[kcc[0]]);
                        cgg.push((*nn,false));
                    }
                }

                for (cc,cbb) in cgg.iter(){
                    let mut checked:HashSet<usize> = HashSet::new();
                    let mut atoms_connected_:Vec<usize> = vec![];
                    if *cbb{//CA 方向からの確認
                        atoms_connected_.push(backboneidx);
                        checked.insert(backboneidx);
                    }
                    let mut atoms_updated:Vec<usize> = vec![];
                    atoms_updated.push(*cc);

                    while atoms_updated.len() > 0{
                        let current:usize = atoms_updated.remove(0);
                        if checked.contains(&current){
                            continue;
                        }
                        atoms_connected_.push(current);
                        checked.insert(current);
                        if edges.contains_key(&current){
                            if edges.get(&current).unwrap().len() <= 2{
                                //二股に分かれる場合そこで止める
                                //終端以外は 2 個ついているはず
                                for nn in edges.get(&current).unwrap().iter(){
                                    if atoms_[*nn].atom_name != "CA"
                                    && atoms_[*nn].atom_name != "C"
                                    && atoms_[*nn].atom_name != "N"
                                    && !checked.contains(nn){
                                        atoms_updated.push(*nn);
                                    }
                                }
                            }
                        }
                    }
                    

                    if *cbb{
                            
                        let mut atoms_this_all:Vec<usize> = vec![];
                        atoms_this_all.push(backboneidx);
                        //cc から伸びる結合の先の ATOM を全部取る
                        atoms_updated.push(*cc);
                        checked.clear();
                        while atoms_updated.len() > 0{
                            let current:usize = atoms_updated.remove(0);
                            if checked.contains(&current){
                                continue;
                            }
                            checked.insert(current);
                            atoms_this_all.push(current);
                            if edges.contains_key(&current){
                                for nn in edges.get(&current).unwrap().iter(){
                                    if atoms_[*nn].atom_name != "CA"
                                    && atoms_[*nn].atom_name != "C"
                                    && atoms_[*nn].atom_name != "N"
                                    && !checked.contains(nn){
                                        atoms_updated.push(*nn);
                                    }
                                }
                            }
                        }
                                            
                        //CA 等 base 側からの場合、元のインデクス以上の Atom をグループ化する
                        let cnum:usize = atoms_connected_.len();
                        for ci in 0..cnum{
                            let mut gatoms:Vec<usize> = vec![];
                            let mut skipper:HashSet<usize> = HashSet::new();
                            for cj in 0..ci{//0 の場合 skipper は 0。全部が skipper に入ることはない
                                skipper.insert(atoms_connected_[cj]);
                            }
                            for sj in 0..atoms_this_all.len(){
                                if !skipper.contains(&atoms_this_all[sj]){
                                    gatoms.push(atoms_this_all[sj]);
                                }
                            }
                            gatoms.sort();
                            sidegroup.push(gatoms);
                        }
                    }else{
                        let cnum:usize = atoms_connected_.len();
                        for ci in 0..cnum{
                            let mut gatoms:Vec<usize> = vec![];
                            for sj in 0..=ci{//0 の場合は最後に外れるはずだが、今後変わるかも
                                //CA まで到達する場合は上のブロックでカバーできるはず
                                gatoms.push(atoms_connected_[sj]);
                            }
                            gatoms.sort();
                            sidegroup.push(gatoms);
                        }
                    }
                }
                let mut vec_checked:HashSet<Vec<usize>> = HashSet::new();
                for ss in sidegroup.into_iter(){
                    if ss.len() <= 1{
                        continue;
                    }
                    if !vec_checked.contains(&ss){
                        vec_checked.insert(ss.clone());
                        //backbone から数えるものは上で処理しているので注意
                        groups_all.push(ss);
                    }
                }
            }
        }
        ret.append(&mut groups_all);
    }
    return ret;
}



pub fn random_movement(envv:&mut charmm_based_energy::CharmmEnv,maxval:f64,targetnum:usize
    ,rgen:&mut StdRng
    ,mov_weight:&Vec<f64>){
    let alen = envv.atoms.len();
    let mut indices:Vec<usize> = (0..alen).collect();
    indices.shuffle(rgen);
    let mut mov_count:usize = 0;
    for ii in 0..alen{
        if mov_weight[indices[ii]] <= 0.0{
            continue;
        }
        let mw = mov_weight[indices[ii]];
        let xyz = envv.atoms[indices[ii]].get_xyz();
        envv.atoms[indices[ii]].set_xyz(
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.0,
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.1,
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.2
        );
        mov_count+=1;
        if mov_count >= targetnum{
            break;
        }
    }
}
pub fn copy_xyz_tmp_to_atom(src:&Vec<(f64,f64,f64)>,dest:&mut Vec<charmm_based_energy::MDAtom>){

    for (ii,aa) in src.iter().enumerate(){
        dest[ii].set_x(aa.0);
        dest[ii].set_y(aa.1);
        dest[ii].set_z(aa.2);
    }
}

pub fn copy_xyz_atom_to_tmp(src:& Vec<charmm_based_energy::MDAtom>,dest:&mut Vec<(f64,f64,f64)>){
    for (ii,aa) in src.iter().enumerate(){
        dest[ii].0 = aa.get_x();
        dest[ii].1 = aa.get_y();
        dest[ii].2 = aa.get_z();
    }
}


//残基一個だけのテスト
#[test]
fn test_aa(){
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(&(debug_env::RESOURCE_DIR.to_string()+"/"+"angle_distribution_energy.dat"));
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    //let mut pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("test/test_lys_nohz3.pdb");
    let mut pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("test/lys_mdd.pdb");
    //let mut pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("test/lys_changed.pdb");
    charmm_based_energy::MDAtom::change_to_charmmnames(&mut pdbb.chains[0].residues);
    let (mut md_envset,mut md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&vec![pdbb.chains.remove(0)],&parr,true);
    let (torsion,omegas) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,true,false),false);
    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&mut md_varset.dihedvec,&torsion,&omegas);

    let mut lines:Vec<String> = vec![];
    for aa in md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
        
    println!("energy:{:?}",scoresum);
    

    write_to_file("test/lys_out.pdb",lines);

    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.2;
        println!("{}",ii);
    }

    for dd in md_varset.anglevec.iter(){
        eprintln!("{} {} {}",&md_envset.atoms[dd.atoms.0].atom_name
        ,&md_envset.atoms[dd.atoms.1].atom_name
        ,&md_envset.atoms[dd.atoms.2].atom_name
        );
    }
  
    let anum:usize = md_envset.atoms.len();
    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,debug_env::RESOURCE_DIR,false),
        backbone_energy_omega:omegas,
        backbone_energy_phi_psi:torsion,
        atom_distance_energy:vec![],
        atom_binned_distance_energy:vec![],
        atom_contact_energy:vec![],
        weights:PPEnergyWeight::new()
    };

    let bret = energy_lbfgs::run_lbfgs(
        &mut ppenergyset
        ,(0..anum).collect()
        ,50
        ,0.001
        ,2
        ,&None
    );

    for ii in 0..anum{
        let xyz:(f64,f64,f64) = ppenergyset.get_atom(ii).get_xyz();
        ppenergyset.get_atom_mut(ii).set_xyz(
            bret.betas[ii*3]+xyz.0
            ,bret.betas[ii*3+1]+xyz.1
            ,bret.betas[ii*3+2]+xyz.2
        );
    }

    let mut lines:Vec<String> = vec![];
    for aa in ppenergyset.evoef2_env.md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file((format!("test/lys_mdd.pdb")).as_str(),lines);
}


#[test]
fn subgroup_assign_test(){
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);
    let seq2 = sequence_alignment::SeqData::load_fasta("test/6F3H_B.pdb_d1.fas",false);
    let mut allaa_:Vec<String> = vec![];
    for ss in seq2[0].seq.iter(){
        allaa_.push(ss.clone());
    }

    let mut ress:Vec<pdbdata::PDBResidue> = chain_builder::build_dirty_chain(&chain_builder::convert_aa_1_to_3(&allaa_),&bset,&sset);
    charmm_based_energy::MDAtom::change_to_charmmnames(&mut ress);
    
    let mut chain:pdbdata::PDBChain = pdbdata::PDBChain::new("A");
    for rr in ress.into_iter(){
        chain.add_residue(rr,true);
    }

    let (mut md_envset,mut md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&vec![chain],&parr,true);
    
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(&(debug_env::RESOURCE_DIR.to_string()+"/"+"angle_distribution_energy.dat"));
    let (_,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,true,false),false);
    
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped("test/6F3H_B.dihed.dat");
    let (torsion,_) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,false,true),false);
    
    let dist_vvec:Vec<(String,usize,String,usize,distance_energy::AtomDistanceEnergy)> = distance_energy::AtomDistanceEnergy::load_name_mapped("test/6F3H_B.cbdist.dat");
    let diste:Vec<distance_energy::AtomDistanceEnergy> = distance_energy::AtomDistanceEnergy::assign_atom_ids(&md_envset,dist_vvec);

    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&mut md_varset.dihedvec,&torsion,&omegas_general);

    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let _scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
        
    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.2;
    }
    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,debug_env::RESOURCE_DIR,false),
        backbone_energy_omega:omegas_general,
        backbone_energy_phi_psi:torsion,
        atom_distance_energy:diste,
        atom_binned_distance_energy:vec![],
        atom_contact_energy:vec![],
        weights:PPEnergyWeight::new()
    };
    ppenergyset.weights.dihed_weight_charmm = weight_dihed;


    eprintln!("{}",ppenergyset.calc_energy(&mut dummy));
    
    let mut sparse_atoms:Vec<usize> = vec![];
    for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        //let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        if aa.atom_name == "CA"
         || aa.atom_name == "CB" 
         || aa.atom_name == "H"
         || aa.atom_name == "C" 
         || aa.atom_name == "O" 
         || aa.atom_name == "N"{
            sparse_atoms.push(aii);
         } 
    }

    let (mut sub_energyset,_mapper):(PPEnergySet,Vec<i64>) = ppenergyset.make_set_for_subenv(&sparse_atoms);
    sub_energyset.update_distance();
    sub_energyset.evoef2_env.md_envset.update_edges(&sub_energyset.evoef2_env.charmm_vars.bondvec,5);
    let mut subdummy:Vec<f64> = vec![0.0;sub_energyset.evoef2_env.md_envset.atoms.len()];
    let score1 = sub_energyset.calc_energy(&mut subdummy);

    let mut sparse_atoms:Vec<usize> = vec![];
    for (aii,aa) in sub_energyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        //let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        if aa.atom_name == "CA"
         || aa.atom_name == "CB" 
         || aa.atom_name == "H"
         || aa.atom_name == "C" 
         || aa.atom_name == "O" 
         || aa.atom_name == "N"{
            sparse_atoms.push(aii);
         } 
    }

    let (mut sub_energyset,_mapper) = sub_energyset.make_set_for_subenv(&sparse_atoms);
    sub_energyset.update_distance();
    let mut subdummy2:Vec<f64> = vec![0.0;sub_energyset.evoef2_env.md_envset.atoms.len()];
    let score2 = sub_energyset.calc_energy(&mut subdummy2);
    assert_eq!(score1,score2);
}
pub fn ditribute_atoms_globular(mdenv:&mut charmm_based_energy::CharmmEnv,box_size:f64,random_seed:Option<u64>){
    
    let mut rgen:StdRng;//乱数生成器
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }

    for aa in mdenv.atoms.iter_mut(){
        let xx = rgen.gen_range(box_size*-0.5_f64,box_size*0.5_f64);
        let yy = rgen.gen_range(box_size*-0.5_f64,box_size*0.5_f64);
        let zz = rgen.gen_range(box_size*-0.5_f64,box_size*0.5_f64);
        aa.set_xyz(xx,yy,zz);
    }
}

//#[test]
pub fn build_test(){
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());

    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);
    let seq2 = sequence_alignment::SeqData::load_fasta("test/6F3H_B.pdb_d1.fas",false);
    let mut allaa_:Vec<String> = vec![];
    for ss in seq2[0].seq.iter(){
        allaa_.push(ss.clone());
    }

    let mut ress:Vec<pdbdata::PDBResidue> = chain_builder::build_dirty_chain(&chain_builder::convert_aa_1_to_3(&allaa_),&bset,&sset);
    charmm_based_energy::MDAtom::change_to_charmmnames(&mut ress);
    
    let mut chain:pdbdata::PDBChain = pdbdata::PDBChain::new("A");
    for rr in ress.into_iter(){
        chain.add_residue(rr,true);
    }


    let (mut md_envset,mut md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&vec![chain],&parr,true);
    
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(&(debug_env::RESOURCE_DIR.to_string()+"/"+"angle_distribution_energy.dat"));
    let (_,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,true,false),false);
    
    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped("test/6F3H_B.dihed.dat");
    let (torsion,_) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,false,true),false);
    
    //let dist_vvec:Vec<(String,usize,String,usize,distance_energy::AtomDistanceEnergy)> = distance_energy::AtomDistanceEnergy::load_name_mapped("test/6F3H_B.cbdist.dat");
    //let diste:Vec<distance_energy::AtomDistanceEnergy> = distance_energy::AtomDistanceEnergy::assign_atom_ids(&md_envset,dist_vvec);

    let dist_bvvec:Vec<(String,usize,String,usize,distance_energy::AtomBinnedDistanceEnergy)> = distance_energy::AtomBinnedDistanceEnergy::load_name_mapped("test/6F3H_B.cbdistbin.dat");
    let distbe:Vec<distance_energy::AtomBinnedDistanceEnergy> = distance_energy::AtomBinnedDistanceEnergy::assign_atom_ids(&md_envset,dist_bvvec);

    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&mut md_varset.dihedvec,&torsion,&omegas_general);

    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
        
    println!("energy:{:?}",scoresum);
   
    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.2;
    }
    
    let mut ppenergyset = PPEnergySet{
        evoef2_env:evoef2_energy::EvoEF2Env::new(md_envset,md_varset,debug_env::RESOURCE_DIR,false),
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
    ppenergyset.weights.atom_binned_distance_energy = 1.0;

    eprintln!("{}",ppenergyset.calc_energy(&mut dummy));
    
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

    ditribute_atoms_globular(&mut sub_energyset.evoef2_env.md_envset,50.0,Some(1234_u64));
    sub_energyset.update_distance();
    sub_energyset.update_edges(5);


    let mut cbatoms:Vec<usize> = vec![];
    for (aii,aa) in sub_energyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        if aa.residue_name == "GLY" && aa.atom_name == "CA"{
            cbatoms.push(aii);
            continue;
        }
        if aa.atom_name == "CB"{
            cbatoms.push(aii);
            continue;
        }
    }
    
    let cbatomnum:usize = cbatoms.len();
    let (mut cb_energyset,cbmapper):(PPEnergySet,Vec<i64>) = sub_energyset.make_set_for_subenv(&cbatoms);
    let mut cbdummy:Vec<f64> = vec![0.0;cbatomnum];
    ditribute_atoms_globular(&mut cb_energyset.evoef2_env.md_envset,50.0,Some(1234_u64));
    cb_energyset.update_distance();
    cb_energyset.gen_pseudo_edges(5);
    let mut cb_pseudobond:Vec<charmm_based_energy::BondVars> = vec![];
    
    for ii in 0..(cbatomnum-1){
        if cb_energyset.evoef2_env.md_envset.atoms[ii+1].chain_name == cb_energyset.get_atom(ii).chain_name 
        && cb_energyset.evoef2_env.md_envset.atoms[ii+1].residue_index_in_chain - cb_energyset.get_atom(ii).residue_index_in_chain == 1{
            cb_pseudobond.push(
                charmm_based_energy::BondVars{
                    atoms:(ii,ii+1),
                    kb:0.0,
                    b0:6.3
                }
            );
        }
    }

    let ignores:Vec<usize> = vec![];
    let mut lines:Vec<String> = vec![];
    for aa in sub_energyset.evoef2_env.md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file((format!("test/6F3H_B_d1_built_cb_mc_start1_b.pdb")).as_str(),lines);
    let sub_bond_restrictor:Vec<charmm_based_energy::BondVars> = sub_energyset.evoef2_env.charmm_vars.bondvec.clone();
    assign_floating_atoms(&mut sub_energyset.evoef2_env.md_envset.atoms,&sub_bond_restrictor,&vec![],&((1..subatomnum).into_iter().collect()),&ignores,0.5);
    sub_energyset.update_distance();

    let mut lines:Vec<String> = vec![];
    for aa in sub_energyset.evoef2_env.md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file((format!("test/6F3H_B_d1_built_cb_mc_start_b.pdb")).as_str(),lines);
    let mut subdummy:Vec<f64> = vec![0.0;sub_energyset.evoef2_env.md_envset.atoms.len()];
    let fixedatoms:HashSet<usize> = HashSet::new();

    eprintln!("!!!{}",sub_energyset.calc_energy(&mut subdummy));
    for tt in 0..30{
        for kkk in 0..100{
            if tt == 0 && kkk == 0{
                for ddd in 0..1{
                    if false{
                        
                        mc_iter_group(&mut cb_energyset
                        , 100000
                        , &fixedatoms
                        ,&make_straight_split_group(cbatomnum)
                        , Some(123)
                        , 0.3
                        , 0.5
                        , 1
                        ,(0.3,0.5,200)
                        ,&vec![]
                        ,&cb_pseudobond
                        );

                    }else if false{
                        mc_iter_array(&mut cb_energyset
                            ,100000
                            ,&fixedatoms
                            ,Some(123)
                            ,cbatomnum/(ddd*2+2)
                            ,ddd > 2
                            ,5.0
                            ,(0.3,0.5,200)
                            ,&vec![]
                            ,&cb_pseudobond
                            ,(1.0,1.0)
                            ,(0.0,0)
                            ,0
                            ,0
                            ,None
                        );
                    }else{
                        let mut cb_partial_minset:Option<PPEnergySet> = None;
                        let mut minmapper:Vec<i64> = vec![];
                        let mut minenergy:f64 = 1000000000.0;
                        for seedpos in (0..cbatomnum).step_by(1000){
                            let loopnum = 50;
                            let mut prevatoms:HashSet<usize> = HashSet::new();
                            for cbb in 0..loopnum{
                                let mut cb_partial_atoms:Vec<usize> = vec![];
                                let mut fixedatoms_partial:HashSet<usize> = HashSet::new();
                                for (aii_,_aa) in cb_energyset.evoef2_env.md_envset.atoms.iter().enumerate(){
                                    if cbb < loopnum-1{
                                        if cb_partial_atoms.len() > (cbb+1)*(cbatomnum/loopnum+1){
                                            break;
                                        }
                                    }
                                    let aii1:i64 = seedpos as i64 -aii_ as i64;
                                    let aii2:i64 = aii_ as i64 +seedpos as i64;
                                    if aii1 >= 0 {
                                        let aii:usize = aii1 as usize;
                                        if prevatoms.contains(&aii){
                                            fixedatoms_partial.insert(cb_partial_atoms.len());
                                        }
                                        cb_partial_atoms.push(aii);
                                    }
                                    if aii_ != 0 && aii2 < cbatomnum as i64{
                                        let aii:usize = aii2 as usize;
                                        if prevatoms.contains(&aii){
                                            fixedatoms_partial.insert(cb_partial_atoms.len());
                                        }
                                        cb_partial_atoms.push(aii);
                                    }
                                }
                                cb_partial_atoms.sort();
                                
                                let cb_partial_atomnum:usize = cb_partial_atoms.len();
                                let (mut cb_partial_energyset,cb_partial_mapper):(PPEnergySet,Vec<i64>) = cb_energyset.make_set_for_subenv(&cb_partial_atoms);
                                let mut cb_partial_dummy:Vec<f64> = vec![0.0;cb_partial_atomnum];
                                cb_partial_energyset.update_distance();
                                cb_partial_energyset.gen_pseudo_edges(5);
                                let mut cb_partial_pseudobond:Vec<charmm_based_energy::BondVars> = vec![];
                                
                                for ii in 0..(cb_partial_atomnum-1){
                                    if cb_partial_energyset.get_atom(ii+1).chain_name == cb_partial_energyset.get_atom(ii).chain_name 
                                    && cb_partial_energyset.get_atom(ii+1).residue_index_in_chain - cb_partial_energyset.get_atom(ii).residue_index_in_chain == 1{
                                        cb_partial_pseudobond.push(
                                            charmm_based_energy::BondVars{
                                                atoms:(ii,ii+1),
                                                kb:0.0,
                                                b0:6.3
                                            }
                                        );
                                    }
                                }
                                

                                mc_iter_array(&mut cb_partial_energyset
                                    ,5000
                                    ,&fixedatoms_partial
                                    ,Some(123)
                                    ,100000
                                    ,false
                                    ,5.0
                                    ,(0.3,0.5,200)
                                    ,&vec![]
                                    ,&cb_partial_pseudobond
                                    ,(1.0,1.0)
                                    ,(0.0,0)
                                    ,0
                                    ,0
                                    ,None
                                );
                                
                                mc_iter_group(&mut cb_partial_energyset
                                    , 5000
                                    , &fixedatoms_partial
                                    ,&make_straight_split_group(cb_partial_atomnum)
                                    , Some(123)
                                    , 0.3
                                    , 0.5
                                    , 1
                                    ,(0.3,0.5,200)
                                    ,&vec![]
                                    ,&cb_partial_pseudobond
                                    );

                                mc_iter_array(&mut cb_partial_energyset
                                    ,5000
                                    //,&fixedatoms_partial
                                    ,&HashSet::new()
                                    ,Some(123)
                                    ,100000
                                    ,false
                                    ,5.0
                                    ,(0.3,0.5,200)
                                    ,&vec![]
                                    ,&cb_partial_pseudobond
                                    ,(1.0,1.0)
                                    ,(0.0,0)
                                    ,0
                                    ,0
                                    ,None
                                );
                            
                                mc_iter_group(&mut cb_partial_energyset
                                    , 5000
                                    , &HashSet::new()
                                    ,&make_straight_split_group(cb_partial_atomnum)
                                    , Some(123)
                                    , 0.1
                                    , 0.5
                                    , 1
                                    ,(0.3,0.5,200)
                                    ,&vec![]
                                    ,&cb_partial_pseudobond
                                    );

                                charmm_based_energy::CharmmEnv::accept_subenv_atom(&mut cb_energyset.evoef2_env.md_envset,&cb_partial_energyset.evoef2_env.md_envset,&cb_partial_mapper);
                                for pp in cb_partial_atoms.into_iter(){
                                    prevatoms.insert(pp);
                                }
                                println!("{}/{}================",cbb+1,loopnum+1);
                                let denergy:f64 = cb_partial_energyset.calc_energy(&mut cb_partial_dummy);

                                if prevatoms.len() == cbatomnum || cbb == loopnum-1{
                                    if let None = cb_partial_minset{
                                        cb_partial_minset = Some(cb_partial_energyset);
                                        minenergy = denergy;
                                        minmapper = cb_partial_mapper;
                                    }else if denergy < minenergy{
                                        cb_partial_minset = Some(cb_partial_energyset);
                                        minenergy = denergy;
                                        minmapper = cb_partial_mapper;
                                    }
                                    break;
                                }
                            }
                        }
                        charmm_based_energy::CharmmEnv::accept_subenv_atom(&mut cb_energyset.evoef2_env.md_envset,&cb_partial_minset.unwrap().evoef2_env.md_envset,&minmapper);
                        

                        let mut lines:Vec<String> = vec![];
                        for aa in cb_energyset.evoef2_env.md_envset.atoms.iter(){
                            let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                            lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                        }
                        write_to_file((format!("test/6F3H_B_d1_built_cbs{}_{}_min{}.pdb",tt,kkk,minenergy)).as_str(),lines);
                    }
                    let prevenergy:f64 = cb_energyset.calc_energy(&mut cbdummy);
                    println!("current: {}",prevenergy);
                    let bret = energy_lbfgs::run_lbfgs(
                        &mut cb_energyset
                        ,(0..cbatomnum).collect()
                        ,1000
                        ,0.001
                        ,5
                        ,&None
                    );
                    for ii in 0..cbatomnum{
                        let xyz:(f64,f64,f64) = cb_energyset.get_atom(ii).get_xyz();
                        cb_energyset.get_atom_mut(ii).set_xyz(
                            bret.betas[ii*3]+xyz.0
                            ,bret.betas[ii*3+1]+xyz.1
                            ,bret.betas[ii*3+2]+xyz.2
                        );
                    }
                    cb_energyset.update_distance();
                    let currentenergy:f64 = cb_energyset.calc_energy(&mut cbdummy);
                    if prevenergy - currentenergy < 0.1{
                        break;
                    }
                    println!("->: {}",currentenergy);
                }
                            
                let mut lines:Vec<String> = vec![];
                for aa in cb_energyset.evoef2_env.md_envset.atoms.iter(){
                    let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                    lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                }
                write_to_file((format!("test/6F3H_B_d1_built_cbs{}_{}_a.pdb",tt,kkk)).as_str(),lines);
                
                let ssvec = cb_energyset.calc_energy_sep(&mut cbdummy);
                for sss in ssvec.1.iter().enumerate(){
                    eprintln!("energy:{} {}",sss.0,sss.1);
                }

                let mut cbfix:HashSet<usize> = HashSet::new();
                for (cii,cbm) in cbmapper.iter().enumerate(){
                    if *cbm > -1{
                        cbfix.insert(cii);
                    }
                }
                charmm_based_energy::CharmmEnv::accept_subenv_atom(&mut sub_energyset.evoef2_env.md_envset,&cb_energyset.evoef2_env.md_envset,&cbmapper);

                //sub_energyset.atom_distance_energy.1 = 1.0/(1.0+kkk as f64/2.0);
                mc_iter_array(&mut sub_energyset
                    ,5000
                    ,&cbfix//CB の位置は固定
                    ,Some(123)
                    ,10
                    ,true
                    ,5.0
                    ,(0.3,0.5,200)
                    ,&sub_bond_restrictor
                    ,&vec![]
                    ,(1.0,1.0)
                    ,(0.0,0)
                    ,0
                    ,0
                    ,None
                    );
            }
            //for xx in vec![(100000,5.0,5000,false),(100000,1.0,5000,false),(100,1.0,5000,false),(10,1.0,5000,false),(10,0.5,5000,false)]{
            for xx in vec![(100000,3,true,0.5),
            ]{
                mc_iter_array(&mut sub_energyset
                    ,xx.0
                    ,&fixedatoms
                    ,Some(123)
                    ,xx.1
                    ,xx.2
                    ,xx.3
                    ,(0.3,0.5,200)
                    ,&sub_bond_restrictor,&vec![]
                    ,(1.0,1.0)
                    ,(0.0,0)
                    ,0
                    ,0
                    ,None
                );
                
                let ssvec = sub_energyset.calc_energy_sep(&mut subdummy);
                for sss in ssvec.1.iter().enumerate(){
                    eprintln!("energy:{} {}",sss.0,sss.1);
                }
                
                let mut lines:Vec<String> = vec![];
                for aa in sub_energyset.evoef2_env.md_envset.atoms.iter(){
                    let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                    lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                }
                write_to_file((format!("test/6F3H_B_d1_built_cb_mc_{}_{}_a{}_{}.pdb",tt,kkk,xx.0,xx.1)).as_str(),lines);
            }
            sub_energyset.update_distance();
            let prevenergy:f64 = sub_energyset.calc_energy(&mut subdummy);
            println!("current: {}",prevenergy);
            
            sub_energyset.weights.atom_binned_distance_energy = 0.0;
            for pp in 0..5{
                let prevenergy:f64 = sub_energyset.calc_energy(&mut subdummy);
                println!("current: {}",prevenergy);
                let bret = energy_lbfgs::run_lbfgs(
                    &mut sub_energyset
                    ,(0..subatomnum).collect()
                    ,100
                    ,0.001
                    ,2
                    ,&None
                );
                for ii in 0..subatomnum{
                    let xyz:(f64,f64,f64) = sub_energyset.get_atom(ii).get_xyz();
                    sub_energyset.get_atom_mut(ii).set_xyz(
                        bret.betas[ii*3]+xyz.0
                        ,bret.betas[ii*3+1]+xyz.1
                        ,bret.betas[ii*3+2]+xyz.2
                    );
                }
                sub_energyset.update_distance();
                let currentenergy:f64 = sub_energyset.calc_energy(&mut subdummy);
                if prevenergy - currentenergy < 0.1{
                    break;
                }
                println!("->: {}",currentenergy);
                let mut lines:Vec<String> = vec![];
                for aa in sub_energyset.evoef2_env.md_envset.atoms.iter(){
                    let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                    lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                }
                write_to_file((format!("test/6F3H_B_d1_built_cb_mc_{}_{}_{}.pdb",tt,kkk,pp)).as_str(),lines);
                let ssvec = sub_energyset.calc_energy_sep(&mut subdummy);
                for sss in ssvec.1.iter().enumerate(){
                    eprintln!("energy:{} {}",sss.0,sss.1);
                }
                
            }
            let mut lines:Vec<String> = vec![];
            for aa in sub_energyset.evoef2_env.md_envset.atoms.iter(){
                let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
            }
            write_to_file((format!("test/6F3H_B_d1_built_cb_mc_{}_{}.pdb",tt,kkk)).as_str(),lines);
        }
    }
}




#[test]
fn group_fit_test(){//エラーが起きないか見るだけ
    for tt in 100..101{
        let mut atoms:Vec<charmm_based_energy::MDAtom> = vec![];
        let atom1 = charmm_based_energy::MDAtom::dummy();
        atoms.push(atom1);
        let mut rgen:StdRng = SeedableRng::seed_from_u64(123);
        let randx:(f64,f64,f64) = (rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        for ii in 0..10{
            atoms.push(charmm_based_energy::MDAtom{
                atom_type:"C".to_string(),
                atom_name:"C".to_string(),
                chain_name:"A".to_string(),
                residue_name:"UNK".to_string(),
                residue_number:100+ii,
                residue_ins_code:"".to_string(),
                
                residue_index_in_chain:1,
                atom_index:0,
                x:randx.0+ii as f64*1.5,
                y:randx.1+(ii%2) as f64 /2.0,
                z:randx.2,
                mass:1.0,
                charge:1.0,
                created:true,
                unplaced:false,
                nb_epsilon:1.0,
                nb_r1_2:1.0,//Rmin1/2
                nb_14flag:false,
                nb_14_epsilon:1.0,
                nb_14_r1_2:1.0,
                hybrid_orbit:charmm_based_energy::HybridOrbit::SP2,
                nterminal:false,
                cterminal:true,
            });
        }
        
        let randx:(f64,f64,f64) = (rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        for ii in 0..10{
            atoms.push(charmm_based_energy::MDAtom{
                atom_type:"C".to_string(),
                atom_name:"C".to_string(),
                chain_name:"A".to_string(),
                residue_name:"UNK".to_string(),
                residue_number:300+ii,
                residue_ins_code:"".to_string(),
                
                residue_index_in_chain:1,
                atom_index:0,
                x:randx.0+ii as f64*1.5,
                y:randx.1+(ii%2) as f64 /2.0,
                z:randx.2,
                mass:1.0,
                charge:1.0,
                created:true,
                unplaced:false,
                nb_epsilon:1.0,
                nb_r1_2:1.0,//Rmin1/2
                nb_14flag:false,
                nb_14_epsilon:1.0,
                nb_14_r1_2:1.0,
                hybrid_orbit:charmm_based_energy::HybridOrbit::SP2,
                nterminal:false,
                cterminal:true,
            });
        }

        
        let randx:(f64,f64,f64) = (rgen.gen_range(-10.0,10.0),rgen.gen_range(-10.0,10.0),rgen.gen_range(-10.0,10.0));
        for ii in 0..1{
            atoms.push(charmm_based_energy::MDAtom{
                atom_type:"C".to_string(),
                atom_name:"C".to_string(),
                chain_name:"A".to_string(),
                residue_name:"UNK".to_string(),
                residue_number:400+ii,
                residue_ins_code:"".to_string(),
                
                residue_index_in_chain:1,
                atom_index:0,
                x:randx.0+ii as f64*1.5,
                y:randx.1+(ii%2) as f64 /2.0,
                z:randx.2,
                mass:1.0,
                charge:1.0,
                created:true,
                unplaced:false,
                nb_epsilon:1.0,
                nb_r1_2:1.0,//Rmin1/2
                nb_14flag:false,
                nb_14_epsilon:1.0,
                nb_14_r1_2:1.0,
                hybrid_orbit:charmm_based_energy::HybridOrbit::SP2,
                nterminal:false,
                cterminal:true,
            });
        }

        let atom_last = charmm_based_energy::MDAtom{
            atom_type:"C".to_string(),
            atom_name:"C".to_string(),
            chain_name:"A".to_string(),
            residue_name:"UNK".to_string(),
            residue_number:999,
            residue_ins_code:"".to_string(),
            
            residue_index_in_chain:1,
            atom_index:0,
            
            x:0.0,
            y:30.0,
            z:30.0,
            mass:1.0,
            charge:1.0,
            created:true,
            unplaced:false,
        
            nb_epsilon:1.0,
            nb_r1_2:1.0,//Rmin1/2
        
            nb_14flag:false,
            nb_14_epsilon:1.0,
            nb_14_r1_2:1.0,
            hybrid_orbit:charmm_based_energy::HybridOrbit::SP2,
        
            nterminal:false,
            cterminal:true,
        };
        atoms.push(atom_last);

        let group1:Vec<usize> = (1..11).into_iter().collect();
        let group2:Vec<usize> = (11..21).into_iter().collect();
        let group3:Vec<usize> = (21..22).into_iter().collect();
        
        let mut lines:Vec<String> = vec![];
        for aa in atoms.iter(){
            let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
            lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
        }
        write_to_file(format!("test/before_fix.pdb").as_str(),lines);

        fix_group_bond(
        &mut atoms
        ,&(vec![group1,group2,group3])
        ,
        &(vec![
            charmm_based_energy::BondVars{atoms:(0,1),b0:2.0,kb:0.0}
            ,charmm_based_energy::BondVars{atoms:(10,11),b0:2.0,kb:0.0}
            ,charmm_based_energy::BondVars{atoms:(20,21),b0:2.0,kb:0.0}
            ,charmm_based_energy::BondVars{atoms:(21,22),b0:2.0,kb:0.0}
            ])
        ,&(vec![
            charmm_based_energy::BondVars{atoms:(11,12),b0:0.1,kb:0.0}
        ])
        ,&((1..22).into_iter().collect())
        ,&vec![]
        ,0.2
        ,tt);

        let mut lines:Vec<String> = vec![];
        for aa in atoms.iter(){
            let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
            lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
        }
        write_to_file(format!("test/fixed_c{}.pdb",tt).as_str(),lines);
    }
}



#[test]
fn group_connection_test(){
    let atomnum = 100;
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(10);
    for pp in 0..100{
        let mut groups:Vec<Vec<usize>> = vec![];
        let mut hss:HashSet<usize> = (0..atomnum).into_iter().collect();
        let mut bonds:Vec<charmm_based_energy::BondVars> = vec![];
        for _ in 0..5{
            let mut group1:Vec<usize> = vec![];
            while group1.len() == 0{
                for _ in 0..15{
                    let uc:usize = rgen.gen_range(0,atomnum);
                    if hss.contains(&uc){
                        hss.remove(&uc);
                        group1.push(uc);
                    }
                }
            }
            groups.push(group1);
        }
        groups.shuffle(&mut rgen);
        for ii in 0..4{//グループ間の結合を指定する
            let r1 = rgen.gen_range(0,groups[ii].len());
            let t1 = groups[ii][r1];
            let r2 = rgen.gen_range(0,groups[ii+1].len());
            let t2 = groups[ii+1][r2];
            let b1 = charmm_based_energy::BondVars{
                atoms:(t1,t2),
                b0:0.0,
                kb:0.0
            };
            bonds.push(b1);
        }

        if pp%2 == 0{//グループメンバーの重なりを加える
            hss.insert(groups[rgen.gen_range(0,5)][0]);
        }
        groups.push(hss.into_iter().collect());
        
        let res = get_connected_group(
            atomnum
            ,&groups
            ,&bonds
        );
        if pp%2 == 0{
            assert_eq!(res.len(),1);
        }else{
            assert_eq!(res.len(),2);
        }
    }
}


#[test]
fn peptide_group_test(){
    let mut bonds_all:Vec<charmm_based_energy::BondVars> = vec![];
    let mut atoms_all:Vec<charmm_based_energy::MDAtom> = vec![];
    let mut lastatomnum:usize = 0;
    for ll in 0..3{
        let mut bonds:Vec<charmm_based_energy::BondVars> = vec![];
        let mut atoms:Vec<charmm_based_energy::MDAtom> = vec![];
        for aa in 0..5 as usize{
            let mut n = charmm_based_energy::MDAtom::dummy();
            n.atom_name = "N".to_owned();

            let mut ca = charmm_based_energy::MDAtom::dummy();
            ca.atom_name = "CA".to_owned();

            let mut c = charmm_based_energy::MDAtom::dummy();
            c.atom_name = "C".to_owned();
            
            let mut o = charmm_based_energy::MDAtom::dummy();
            o.atom_name = "O".to_owned();
            
            let mut h = charmm_based_energy::MDAtom::dummy();
            h.atom_name = "H".to_owned();

            let mut ha = charmm_based_energy::MDAtom::dummy();
            ha.atom_name = "HA".to_owned();
            
            n.residue_index_in_chain = aa as i64;
            ca.residue_index_in_chain = aa as i64;
            c.residue_index_in_chain = aa as i64;
            o.residue_index_in_chain = aa as i64;
            h.residue_index_in_chain = aa as i64;
            ha.residue_index_in_chain = aa as i64;

            atoms.push(n);
            atoms.push(ca);
            atoms.push(c);
            atoms.push(o);
            atoms.push(h);
            atoms.push(ha);

            bonds.push(
                charmm_based_energy::BondVars{
                    atoms:(aa*6,aa*6+1),kb:0.0,b0:3.4
                }
            );
            
            bonds.push(
                charmm_based_energy::BondVars{
                    atoms:(aa*6+1,aa*6+2),kb:0.0,b0:3.4
                }
            );
            
            bonds.push(
                charmm_based_energy::BondVars{
                    atoms:(aa*6+2,aa*6+3),kb:0.0,b0:3.4
                }
            );
            
            bonds.push(
                charmm_based_energy::BondVars{
                    atoms:(aa*6+0,aa*6+4),kb:0.0,b0:3.4
                }
            );
            
            bonds.push(
                charmm_based_energy::BondVars{
                    atoms:(aa*6+1,aa*6+5),kb:0.0,b0:3.4
                }
            );
                
            if aa > 0{
                bonds.push(
                    charmm_based_energy::BondVars{
                        atoms:((aa-1)*6+2,aa*6),kb:0.0,b0:3.4
                    }
                );
            }
        }
        if ll == 0{
            let res:Vec<Vec<usize>> = make_peptide_atom_groups(&atoms,&(bonds.iter().collect()));
            let check_hs:HashSet<Vec<usize>> = res.into_iter().collect();

            assert!(check_hs.contains(&vec![0,4]));
            assert!(check_hs.contains(&vec![0,1,4,5]));
            assert!(check_hs.contains(&vec![0,1,2,3,4,5]));
            assert!(check_hs.contains(&vec![0,1,2,3,4,5,6,10]));
            assert!(check_hs.contains(&vec![20,21,24,25,26,27,28,29]));
            assert!(check_hs.contains(&vec![26,27]));
            assert!(!check_hs.contains(&vec![24,25]));//中途半端だが。。。
        }
        //residue number について考えてないので必要になったら追加
        for aa in 30..=51{
            let mut a = charmm_based_energy::MDAtom::dummy();
            a.residue_index_in_chain = aa as i64;
            a.atom_name = "X".to_owned();
            atoms.push(a);
        }

        let aa = 0;//H が重複してしまうがテストなので気にしない
        bonds.push(
            charmm_based_energy::BondVars{
                atoms:(aa*6,31),kb:0.0,b0:3.4
            }
        );
        bonds.push(
            charmm_based_energy::BondVars{
                atoms:(aa*6,32),kb:0.0,b0:3.4
            }
        );
        bonds.push(
            charmm_based_energy::BondVars{
                atoms:(aa*6,33),kb:0.0,b0:3.4
            }
        );
        
        let aa = 1;
        bonds.push(charmm_based_energy::BondVars{atoms:(aa*6+1,35),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(35,36),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(36,37),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(37,38),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(38,39),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(39,40),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(40,35),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(40,41),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(39,42),kb:0.0,b0:3.4});
        
        let aa = 2;
        bonds.push(charmm_based_energy::BondVars{atoms:(aa*6+1,43),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(43,44),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(44,45),kb:0.0,b0:3.4});

        let aa = 3;
        bonds.push(charmm_based_energy::BondVars{atoms:(aa*6+1,46),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(46,47),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(47,48),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(48,49),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(49,50),kb:0.0,b0:3.4});
        bonds.push(charmm_based_energy::BondVars{atoms:(48,51),kb:0.0,b0:3.4});
        for aa in atoms.iter_mut(){
            if ll == 0{
                aa.chain_name = "A".to_owned();
            }else if ll == 1{
                aa.chain_name = "B".to_owned();
            }else if ll == 2{
                aa.chain_name = "C".to_owned();
            }else{
                panic!();
            }
        }
        for bb in bonds.iter_mut(){
            bb.atoms.0 += lastatomnum;
            bb.atoms.1 += lastatomnum;
        }
        let anum = atoms.len();
        atoms_all.append(&mut atoms);
        bonds_all.append(&mut bonds);
        let res:Vec<Vec<usize>> = make_peptide_atom_groups(&atoms_all,&(bonds_all.iter().collect()));
        //println!("{:?}",res);
        let check_hs:HashSet<Vec<usize>> = res.into_iter().collect();
        for kk in 0..=ll{
            let ashif = anum*kk;
            assert!(check_hs.contains(&vec![0+ashif,31+ashif]));
            assert!(check_hs.contains(&vec![0+ashif,32+ashif]));
            assert!(check_hs.contains(&vec![0+ashif,33+ashif]));
            assert!(!check_hs.contains(&vec![0+ashif,31+ashif,32+ashif,33+ashif]));

            assert!(check_hs.contains(&vec![0+ashif,4+ashif,31+ashif,32+ashif,33+ashif]));
            assert!(check_hs.contains(&vec![0+ashif,1+ashif,4+ashif,5+ashif,31+ashif,32+ashif,33+ashif]));
            assert!(check_hs.contains(&vec![0+ashif,1+ashif,2+ashif,3+ashif,4+ashif,5+ashif,31+ashif,32+ashif,33+ashif]));
            assert!(check_hs.contains(&vec![0+ashif,1+ashif,2+ashif,3+ashif,4+ashif,5+ashif,6+ashif,10+ashif,31+ashif,32+ashif,33+ashif]));
            assert!(check_hs.contains(&vec![39+ashif,42+ashif]));
            assert!(check_hs.contains(&vec![40+ashif,41+ashif]));
            assert!(!check_hs.contains(&vec![39+ashif,40+ashif,42+ashif]));
            assert!(!check_hs.contains(&vec![38+ashif,39+ashif,40+ashif]));
            assert!(check_hs.contains(&vec![35+ashif,36+ashif,37+ashif,38+ashif,39+ashif,40+ashif,41+ashif,42+ashif]));
            assert!(check_hs.contains(&vec![7+ashif,35+ashif,36+ashif,37+ashif,38+ashif,39+ashif,40+ashif,41+ashif,42+ashif]));

            assert!(check_hs.contains(&vec![13+ashif,43+ashif,44+ashif,45+ashif]));
            assert!(check_hs.contains(&vec![43+ashif,44+ashif,45+ashif]));
            assert!(!check_hs.contains(&vec![13+ashif,43+ashif]));
            assert!(check_hs.contains(&vec![44+ashif,45+ashif]));
            assert!(!check_hs.contains(&vec![43+ashif,44+ashif]));
            
            assert!(check_hs.contains(&vec![19+ashif,46+ashif,47+ashif,48+ashif,49+ashif,50+ashif,51+ashif]));
            assert!(!check_hs.contains(&vec![19+ashif,46+ashif,47+ashif,48+ashif,49+ashif,50+ashif]));
            assert!(!check_hs.contains(&vec![19+ashif,46+ashif,47+ashif,48+ashif,49+ashif,51+ashif]));
            assert!(check_hs.contains(&vec![46+ashif,47+ashif,48+ashif,49+ashif,50+ashif,51+ashif]));
            assert!(check_hs.contains(&vec![47+ashif,48+ashif,49+ashif,50+ashif,51+ashif]));
            assert!(check_hs.contains(&vec![48+ashif,49+ashif,50+ashif,51+ashif]));
            assert!(check_hs.contains(&vec![48+ashif,51+ashif]));
            assert!(check_hs.contains(&vec![48+ashif,49+ashif,50+ashif]));
            assert!(check_hs.contains(&vec![49+ashif,50+ashif]));
        }
        lastatomnum = atoms_all.len();
    }
}



/*

pub fn tbm_and_loop_test(){
    chain_builder::prepare_static();

    let allaa:Vec<String> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].iter().map(|m|m.to_string()).collect();
    let pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/1a4w_part.pdb");

    let residues:Vec<&pdbdata::PDBResidue> = pdbb.chains[0].residues.iter().collect();
    let mut dummystring:Vec<String> = vec![];
    for rr in residues.iter(){
        dummystring.push(rr.get_name().to_string());
    }
    let mut dummystring = chain_builder::convert_aa_3_to_1(&dummystring);
    let mut dummyquery:Vec<String> = vec!["K".to_string();dummystring.len()];
    dummyquery.push("T".to_string());
    dummyquery.push("T".to_string());
    dummyquery.push("T".to_string());
    dummyquery.push("T".to_string());
    dummyquery.push("T".to_string());
    dummyquery.push("T".to_string());

    dummystring.insert(10,"-".to_string());
    dummystring.insert(10,"-".to_string());
    dummystring.insert(10,"-".to_string());
    dummystring.insert(10,"-".to_string());
    dummystring.insert(10,"-".to_string());
    dummystring.insert(10,"-".to_string());
    eprintln!("Loading backbone samples...");
    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    for aa in allaa.iter(){
        eprintln!("{} {}",aa,bset.get_num(aa));
    }

    eprintln!("Loading sidechain samples...");
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);
    for aa in allaa.iter(){
        eprintln!("{} {}",aa,sset.get_num(aa));
    }


    let ress_and_flag:Vec<(pdbdata::PDBResidue,bool)> = chain_builder::build_from_alignment(&dummyquery,&dummystring,&residues,&bset,&sset);
    let mut ress:Vec<pdbdata::PDBResidue> = vec![];
    let mut flags:Vec<bool> = vec![];
    for (rr,ff) in ress_and_flag.into_iter(){
        ress.push(rr);
        flags.push(ff);
    }
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    

    charmm_based_energy::MDAtom::change_to_charmmnames(&mut ress);
    
    let mut chain:pdbdata::PDBChain = pdbdata::PDBChain::new("A".to_string());
    for rr in ress.into_iter(){
        chain.add_residue(rr,true);
    }

    let (mut md_envset,mut md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&mut chain,&parr,true);
    
    let mut aligned_atoms:HashSet<usize> = HashSet::new();
    for (aii,aa) in md_envset.atoms.iter().enumerate(){
        if flags[aa.residue_index_in_chain as usize]{
            aligned_atoms.insert(aii);
        }
    }

    let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped(&(debug_env::RESOURCE_DIR.to_string()+"/"+"angle_distribution_energy.dat"));
    let (torsion_general,omegas_general) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,true,false,false));
    //ここから
    //拘束がある場合は Vec に入れてある場合ない場合両方対応できるように

    //let pvec:HashMap<String,peptide_backbone_dihedral_energy::PlainDistribution> = peptide_backbone_dihedral_energy::PlainDistribution::load_name_mapped("test/6F3H_B.dihed.dat");
    //let (torsion,_) = peptide_backbone_dihedral_energy::PlainDistribution::create_energy_instance(&pvec,&md_envset,(false,false,true,false));
    
    //let dist_vvec:Vec<(String,usize,String,usize,distance_energy::AtomDistanceEnergy)> = distance_energy::AtomDistanceEnergy::load_name_mapped("test/6F3H_B.cbdist.dat");
    //let diste:Vec<distance_energy::AtomDistanceEnergy> = distance_energy::AtomDistanceEnergy::assign_atom_ids(&md_envset,dist_vvec);

    //binned dist
    //let dist_bvvec:Vec<(String,usize,String,usize,distance_energy::AtomBinnedDistanceEnergy)> = distance_energy::AtomBinnedDistanceEnergy::load_name_mapped("test/6F3H_B.cbdistbin.dat");
    //let distbe:Vec<distance_energy::AtomBinnedDistanceEnergy> = distance_energy::AtomBinnedDistanceEnergy::assign_atom_ids(&md_envset,dist_bvvec);

    let masked:Vec<usize> = peptide_backbone_dihedral_energy::get_overwrapping_dihed(&mut md_envset,&mut md_varset.dihedvec,&torsion_general,&omegas_general);

    let mut dummy:Vec<f64> = vec![0.0;md_envset.atoms.len()];
    let scoresum = charmm_based_energy::calc_energy(&mut md_envset,&md_varset,&mut dummy);
        
    println!("energy:{:?}",scoresum);
   
    let mut weight_dihed = vec![1.0;md_varset.dihedvec.len()];
    for ii in masked.into_iter(){
        weight_dihed[ii] = 0.2;
    }
    let ppenergyset = PPEnergySet{ 
        charmm_bond:(md_varset.bondvec,1.0),
        charmm_angle:(md_varset.anglevec,1.0),
        charmm_ub:(md_varset.ubvec,1.0),
        charmm_dihed:(md_varset.dihedvec,1.0),
        charmm_impr:(md_varset.imprvec,1.0),
        charmm_lj:(md_varset.lj_energy_calculator,1.0),
        charmm_electro:(md_varset.electrostatic_energy_calculator,0.0),
        dihed_weight_charmm:weight_dihed,//backbone と重複するやつにウエイトを掛けるための係数
        backbone_energy_omega:(omegas_general,1.0),
        //backbone_energy_phi_psi:(torsion,1.0),
        backbone_energy_phi_psi:(torsion_general,1.0),
        atom_distance_energy:(vec![],1.0),
        //atom_binned_distance_energy:(distbe,100.0),
        atom_binned_distance_energy:(vec![],100.0),
    };
    
    eprintln!("{}",ppenergyset.calc_energy(&md_envset,&mut dummy));
    
    let mut sparse_atoms:Vec<usize> = vec![];
    
    for (aii,aa) in md_envset.atoms.iter().enumerate(){
        if aa.atom_name == "CA"
        || aa.atom_name == "CB" 
        || aa.atom_name == "C" 
        || aa.atom_name == "N"{
           sparse_atoms.push(aii);
        }
    }
    //ToDo
    //N 末だけでなく中間からもできるようにする
    //作った後にクラスタリングとかするか。
    //動きについて Bond の回転等だけで実現するものを作る

    let (mut subenv,mapper):(charmm_based_energy::CharmmEnv,Vec<i64>) = md_envset.make_sub_env(&sparse_atoms);
    let mut aligned_atoms_sub:HashSet<usize> = HashSet::new();
    for (mii,mm) in mapper.iter().enumerate(){
        if aligned_atoms.contains(&mii){
            aligned_atoms_sub.insert(*mm as usize);
        }
    }

    let subatomnum:usize = subenv.atoms.len();
    //ditribute_atoms_globular(&mut subenv,50.0,Some(1234_u64));
    let sub_energyset:PPEnergySet = ppenergyset.make_set_for_subenv(&mapper);
    subenv.update_distance();
    subenv.update_edges(&sub_energyset.charmm_bond.0,5);
    let ignores:Vec<usize> = vec![];
    let mut lines:Vec<String> = vec![];
    for aa in subenv.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file((format!("test/tbmtest_start1.pdb")).as_str(),lines);
    assign_floating_atoms(&mut subenv,&sub_energyset,&((1..subatomnum).into_iter().collect()),&ignores,0.5);
    
    let mut lines:Vec<String> = vec![];
    for aa in subenv.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file((format!("test/tbmtest_start2.pdb")).as_str(),lines);
    let mut subdummy:Vec<f64> = vec![0.0;subenv.atoms.len()];

    eprintln!("!!!{}",sub_energyset.calc_energy(&subenv,&mut subdummy));
    for tt in 0..30{
        for kkk in 0..100{
            //sub_energyset.atom_distance_energy.1 = 1.0/(1.0+kkk as f64/2.0);
            mc_iter_array(&mut subenv
                ,&sub_energyset
                ,100
                ,&aligned_atoms_sub
                ,Some(123)
                ,100000
                ,false
                ,3.0
            ,(0.3,0.5,200));
            
            let ssvec = sub_energyset.calc_energy_sep(&subenv,&mut subdummy);
            for sss in ssvec.1.iter().enumerate(){
                eprintln!("energy:{} {}",sss.0,sss.1);
            }
            
            let mut lines:Vec<String> = vec![];
            for aa in subenv.atoms.iter(){
                let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
            }
            write_to_file((format!("test/tbmtest_{}_{}_a.pdb",tt,kkk)).as_str(),lines);
            
            subenv.update_distance();
            let prevenergy:f64 = sub_energyset.calc_energy(&subenv, &mut subdummy);
            println!("current: {}",prevenergy);
            
            for pp in 0..3{
                let prevenergy:f64 = sub_energyset.calc_energy(&subenv, &mut subdummy);
                println!("current: {}",prevenergy);
                //全 ATOM aligned については BACKBONE のみ FIX にするか。
                let allsub:Vec<usize> = (0..subatomnum).collect();
                let mut movatom:Vec<usize> = vec![];
                for aa in allsub.into_iter(){
                    if !aligned_atoms_sub.contains(&aa){
                        movatom.push(aa);
                    }
                }
                let bret = energy_lbfgs::run_lbfgs(
                    &mut subenv
                    ,&sub_energyset
                    ,movatom.clone()
                    ,100
                    ,0.001
                    ,2
                    ,&None
                );
                for ii_ in 0..movatom.len(){
                    let ii = movatom[ii_];
                    let xyz:(f64,f64,f64) = subenv.atoms[ii].get_xyz();
                    subenv.atoms[ii].set_xyz(
                        bret.betas[ii_*3]+xyz.0
                        ,bret.betas[ii_*3+1]+xyz.1
                        ,bret.betas[ii_*3+2]+xyz.2
                    );
                }
                subenv.update_distance();
                let currentenergy:f64 = sub_energyset.calc_energy(&subenv, &mut subdummy);
                if prevenergy - currentenergy < 0.1{
                    break;
                }
                println!("->: {}",currentenergy);
                let mut lines:Vec<String> = vec![];
                for aa in subenv.atoms.iter(){
                    let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                    lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
                }
                write_to_file((format!("test/tbm_test_{}_{}_{}.pdb",tt,kkk,pp)).as_str(),lines);
                let ssvec = sub_energyset.calc_energy_sep(&subenv,&mut subdummy);
                for sss in ssvec.1.iter().enumerate(){
                    eprintln!("energy:{} {}",sss.0,sss.1);
                }
                
            }
            let mut lines:Vec<String> = vec![];
            for aa in subenv.atoms.iter(){
                let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
                lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
            }
            write_to_file((format!("test/tbm_test_{}_{}_fin.pdb",tt,kkk)).as_str(),lines);
        }
    }
}

*/
extern crate rand;
use super::charmm_based_energy;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
use self::rand::prelude::*;
use super::geometry::Vector3D;
#[allow(unused_imports)]
use super::misc_util::*;

#[allow(unused_imports)]
use super::pdbdata;
#[allow(unused_imports)]
use super::charmm_param;
#[allow(unused_imports)]
use super::process_3d;
#[allow(unused_imports)]
use super::debug_env;
#[allow(unused_imports)]
use super::mmcif_process;




pub fn drmsd(atoms:&Vec<(f64,f64,f64)>,reference_dist:&Vec<Vec<f64>>)->f64{
    let alen:usize = atoms.len();
    let mut ret:f64 = 0.0;
    let mut lcou:usize = 0;
    for ii in 0..(alen-1){
        for jj in (ii+1)..alen{
            if reference_dist[ii][jj] > 0.0{
                lcou += 1;
                let d:f64 = process_3d::distance(&atoms[ii],&atoms[jj]);
                ret += (d-reference_dist[ii][jj]).powf(2.0);
                //println!("{} {}",d,reference_dist[ii][jj]);
            }
        }
    }
    if ret > 0.0{
        return (ret/(lcou as f64)).sqrt();
    }
    return 0.0;
}


//referencedist[0][1] = distance(atoms[atom_indices[0]],atoms[atom_indices[1]])
pub fn drmsd_sparse(envv:&charmm_based_energy::CharmmEnv,reference_dist:&Vec<Vec<f64>>
,atom_indices:&Vec<usize>)->f64{
    let alen:usize = atom_indices.len();
    let mut ret:f64 = 0.0;
    let mut lcou:usize = 0;
    for ii in 0..(alen-1){
        for jj in (ii+1)..alen{
            if reference_dist[ii][jj] > 0.0{
                lcou += 1;
                let d:f64 = envv.atoms[atom_indices[0]].distance(&envv.atoms[atom_indices[1]]);
                ret += (d-reference_dist[ii][jj]).powf(2.0);
                //println!("{} {}",d,reference_dist[ii][jj]);
            }
        }
    }
    if ret > 0.0{
        return (ret/(lcou as f64)).sqrt();
    }
    return 0.0;
}


pub fn drmsd_array_sparse(envv:&mut charmm_based_energy::CharmmEnv
,reference_dist:&Vec<Vec<f64>>
,atom_indices:&Vec<usize>
,drmsd_buff:&mut Vec<f64>){
    let alen:usize = atom_indices.len();
    for ii in 0..alen{
        drmsd_buff[ii] = 0.0;
    }
    for ii in 0..(alen-1){
        let sii = atom_indices[ii];
        for jj in (ii+1)..alen{
            let sjj = atom_indices[jj];
            if reference_dist[ii][jj] > 0.0{
                let d:f64 = envv.atoms[sii].distance(&envv.atoms[sjj]);
                let dd:f64 = (d-reference_dist[ii][jj]).powf(2.0);
                if dd > 0.0{
                    let dd = dd.sqrt();
                    drmsd_buff[ii] += dd;
                    drmsd_buff[jj] += dd;
                }
            }
        }
    }
    
    for ii in 0..alen{
        drmsd_buff[ii] /= alen as f64;
    }
}


pub fn random_movement(envv:&mut charmm_based_energy::CharmmEnv,maxval:f64,targetnum:usize
    ,rgen:&mut StdRng
    ,mov_weight:&Vec<f64>){
    let alen = envv.atoms.len();
    let mut indices:Vec<usize> = (0..alen).collect();
    indices.shuffle(rgen);
    for ii in 0..targetnum{
        if mov_weight[ii] <= 0.0{
            continue;
        }
        let mw = mov_weight[ii];
        let xyz = envv.atoms[indices[ii]].get_xyz();
        envv.atoms[indices[ii]].set_xyz(
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.0,
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.1,
            rgen.gen_range(-0.5*maxval,maxval*0.5)*mw+xyz.2
        );
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



//原子ごとに評価して、原子ごとに MH 基準で採択する
pub fn mc_iter_array(md_envset:&mut charmm_based_energy::CharmmEnv
    ,md_varset:&mut charmm_based_energy::CharmmVars
    ,reference_dist:&Vec<Vec<f64>>
    ,distatom_indices:&Vec<usize>
    ,iternum:usize
    ,random_seed:Option<u64>
    ,num_moveatom:usize){
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
    let mut tmpatoms:Vec<(f64,f64,f64)> =vec![];
    let mut tmpatoms_prevpos:Vec<(f64,f64,f64)> =vec![];

    for aa in md_envset.atoms.iter(){
        tmpatoms.push((aa.get_x(),aa.get_y(),aa.get_z()));
        tmpatoms_prevpos.push((aa.get_x(),aa.get_y(),aa.get_z()));
    }



    let mut current_drmsds:Vec<f64> = vec![std::f64::NEG_INFINITY;md_envset.atoms.len()];
    let mut prev_scores:Vec<f64> = vec![std::f64::NEG_INFINITY;md_envset.atoms.len()];
    let mut current_scores:Vec<f64> = vec![std::f64::NEG_INFINITY;md_envset.atoms.len()];

    let mut mov_weight:Vec<f64> = vec![1.0;md_envset.atoms.len()];
    let mut unplacedflag:bool = false;
    for ii in 0.. md_envset.atoms.len(){
        if md_envset.atoms[ii].unplaced{
            unplacedflag = true;
            break;
        }
    }
    for _ii in 0..iternum{
        let mov:f64 = 0.1;
        //動かす最大距離。ちゃんとオプションとして渡すように
        if _ii as f64 > iternum as f64 * 0.5{
        //    mov = 0.025;
        }


        let mut sorter:Vec<(usize,f64)> = prev_scores.iter().enumerate().map(|m|(m.0,*m.1)).collect();
        sorter.sort_by(|a,b|a.1.partial_cmp(&b.1).unwrap_or_else(||panic!("error in sort!")));
        //スコアが低い奴は大きく動かす
        for ii in 0.. sorter.len(){
            if ii < sorter.len()/20{
        //        mov_weight[sorter[ii].0] = 1.0;
            }else{
        //        mov_weight[sorter[ii].0] = 0.01;
            }
        }
        
        if unplacedflag && (_ii as f64) < iternum as f64*0.25{ 
            for ii in 0.. md_envset.atoms.len(){
                if md_envset.atoms[ii].unplaced{
                    mov_weight[ii] = 1.0;
                }else{
                    mov_weight[ii] = 0.0;
                }
            }
        }else{
            for ii in 0.. md_envset.atoms.len(){
                mov_weight[ii] = 1.0;
               
            }
        }
        for ii in 0.. md_envset.atoms.len(){
            tmpatoms_prevpos[ii].0 = md_envset.atoms[ii].get_x();
            tmpatoms_prevpos[ii].1 = md_envset.atoms[ii].get_y();
            tmpatoms_prevpos[ii].2 = md_envset.atoms[ii].get_z();
        }
        random_movement(md_envset,mov,num_moveatom,& mut rgen,&mov_weight);

        drmsd_array_sparse(md_envset, reference_dist,&distatom_indices,&mut current_drmsds);
        
        //println!("Calculated drmsd.");
        let current_drmsd = drmsd_sparse(md_envset, reference_dist,&distatom_indices);
        
        //この段階では高い方が良くない値が current_scores に入る
        let scoresum = charmm_based_energy::calc_energy(md_envset,md_varset,&mut current_scores);
        
        //println!("Calculated energy.");
        println!("iter: {} drmsd: {} energy:{:?}",_ii,current_drmsd,scoresum);
        
        //高い方が良い値に変換する
        for uu in 0..current_scores.len(){
            current_scores[uu] = 1.0/(current_scores[uu] +1000.0*md_envset.atoms.len() as f64);
        }
        for ii in 0..distatom_indices.len(){
            let sparseidx:usize = distatom_indices[ii];
            //1000.0 とかいまのところ適当
            current_scores[sparseidx] += 1000.0*1.0/(current_drmsds[sparseidx]+0.00001);
        }
        if false{//全体でスコアが上がらないと採択しない
            let csum:f64 = current_scores.iter().fold(0.0,|s,m| s+m);
            let psum:f64 = prev_scores.iter().fold(0.0,|s,m| s+m);
            let mut accept = csum > psum;
            if !accept{
                let dd = rgen.gen_range(0.0,1.0);
                //10.0 とかいまのところ適当
                if 1.0/(1.0+(psum/csum).exp()*10.0) > dd{
                    accept = true;
                }
            }
            if accept{
                for ii in 0..current_scores.len(){
                    prev_scores[ii] = current_scores[ii];
                }
            }else{
                for ii in 0..current_scores.len(){
                        md_envset.atoms[ii].set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                }
            }
        }else{
            for ii in 0..current_scores.len(){
                let current_score = current_scores[ii];
                let mut accept = current_score > prev_scores[ii];
                if !accept{
                    let dd = rgen.gen_range(0.0,1.0);
                    //10.0 とかいまのところ適当
                    if 1.0/(1.0+(prev_scores[ii]/current_score).exp()*10.0) > dd{
                        accept = true;
                    }
                }
                if accept{
                    prev_scores[ii] = current_score;
                }else{
                    //前の状態に戻す
                    md_envset.atoms[ii].set_xyz(tmpatoms_prevpos[ii].0,tmpatoms_prevpos[ii].1,tmpatoms_prevpos[ii].2);
                }
            }
        }
        if current_drmsd < 0.01{
            break;
        }
    }
    let current_drmsd = drmsd_sparse(&md_envset, reference_dist,&distatom_indices);
    println!("result: {}",current_drmsd);
}


#[test]

fn ca_mc_md_test(){
    let mut pdbb:pdbdata::PDBEntry = mmcif_process::load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/1a4w_part.pdb");
    let mut topp:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_top_all22_inp((debug_env::CHARMM_DIR.to_string()+"\\top_all22_prot.rtf").as_str());
    let mut parr:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_param((debug_env::CHARMM_DIR.to_string()+"\\par_all22_prot.prm").as_str());
    parr.resi.append(&mut topp.resi);
    

    let mut ca_atoms_pos:Vec<(usize,(f64,f64,f64))> = vec![];
    
    charmm_based_energy::MDAtom::change_to_charmmnames(&mut pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_comps().map(|m|*m).collect());

    let (mut md_envset,mut md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars)
     = charmm_based_energy::MDAtom::chain_to_atoms(&vec![pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0)],&parr,true);

    let mut ca_md:Vec<&charmm_based_energy::MDAtom> = vec![];
    for (_rii,aa) in md_envset.atoms.iter().enumerate(){
        if aa.atom_name == "CA"{
        //後のステップでCA が全ての AA について 1 つずつ存在し、Residue 順に Vector に加えられているとみなしている。
            ca_md.push(aa);
        }
    }


    for (_rii,rr) in pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_comps().enumerate(){
        if rr.get_comp_id() == "HOH"{
            continue;
        }

        //CA が全ての AA について 1 つずつ存在し、Residue 順に Vector に加えられているとみなしている。
        //名前だけチェック。ずれていたら panic する。その場合もっと保証できる方法に変更。
        //CA distance の時のマップに必要
        assert_eq!(ca_md[_rii].residue_name,rr.get_comp_id());
        for (_aii,aa) in rr.iter_atoms().enumerate(){
            if aa.atom_code == "CA"{
                ca_atoms_pos.push((_rii,(aa.get_x(),aa.get_y(),aa.get_z())));
            }
        }
    }


    let mut ref_dist:Vec<Vec<f64>> = vec![vec![-1.0;ca_atoms_pos.len()];ca_atoms_pos.len()];
    let alen:usize = ca_atoms_pos.len();
    let mut ca_indices:Vec<usize> = vec![];
    for ii in 0..alen{
        ca_indices.push(ca_atoms_pos[ii].0);
        for jj in 0..alen{
            let ddist = process_3d::distance(&ca_atoms_pos[ii].1,&ca_atoms_pos[jj].1);
            if  ddist < 12.0{
                ref_dist[ii][jj] = ddist;
            }
        }
    }

    let num_atoms = md_envset.atoms.len();
    println!("Params were loaded.");
    mc_iter_array(&mut md_envset,&mut md_varset,&ref_dist,&ca_indices,10,Some(123),num_atoms);
    let mut lines:Vec<String> = vec![];
    for aa in md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file("test/mdmcout.pdb",lines);
}
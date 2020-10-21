extern crate rand;
use super::pdbdata::*;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
use self::rand::prelude::*;
use super::geometry::Vector3D;


pub fn drmsd(atoms:&Vec<PDBAtom>,reference_dist:&Vec<Vec<f64>>)->f64{
    let alen:usize = atoms.len();
    let mut ret:f64 = 0.0;
    let mut lcou:usize = 0;
    for ii in 0..(alen-1){
        for jj in (ii+1)..alen{
            if reference_dist[ii][jj] > 0.0{
                lcou += 1;
                let d:f64 = atoms[ii].distance(&atoms[jj]);
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


pub fn drmsd_array(atoms:&Vec<PDBAtom>,reference_dist:&Vec<Vec<f64>>,drmsd_buff:&mut Vec<f64>){
    let alen:usize = atoms.len();
    for ii in 0..alen{
        drmsd_buff[ii] = 0.0;
    }
    for ii in 0..(alen-1){
        for jj in (ii+1)..alen{
            if reference_dist[ii][jj] > 0.0{
                let d:f64 = atoms[ii].distance(&atoms[jj]);
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


pub fn random_movement(atoms:&mut Vec<PDBAtom>,maxval:f64,targetnum:usize,rgen:&mut StdRng){
    let alen = atoms.len();
    let mut indices:Vec<usize> = (0..alen).collect();
    indices.shuffle(rgen);
    for ii in 0..targetnum{
        let aa:&mut PDBAtom = &mut atoms[indices[ii]];
        let x = aa.get_x()+rgen.gen_range(-0.5*maxval,maxval*0.5);
        let y = aa.get_y()+rgen.gen_range(-0.5*maxval,maxval*0.5);
        let z = aa.get_z()+rgen.gen_range(-0.5*maxval,maxval*0.5);
        aa.set_xyz(x,y,z);
    }
}
pub fn copy_xyz(src:&Vec<PDBAtom>,dest:&mut Vec<PDBAtom>){

    for (ii,aa) in src.iter().enumerate(){
        dest[ii].set_x(aa.get_x());
        dest[ii].set_y(aa.get_y());
        dest[ii].set_z(aa.get_z());
    }
}

//評価は全体の構造
pub fn mc_iter(atoms:&mut Vec<PDBAtom>
    ,reference_dist:&Vec<Vec<f64>>
    ,iternum:usize,random_seed:Option<u64>
    ,num_moveatom:usize//random movement で一回に動かす数
){
    let mut rgen:StdRng;//乱数生成器
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }
    let mut tmpatoms:Vec<PDBAtom> =vec![];
    for aa in atoms.iter(){
        let mut natom = PDBAtom::new();
        natom.set_x(aa.get_x());
        natom.set_y(aa.get_y());
        natom.set_z(aa.get_z());
        tmpatoms.push(natom);
    }
    let mut prev_drmsd = drmsd(&tmpatoms, reference_dist);
    //let mut prev_score:f64 = 1.0/(drmsd(atoms, reference_dist)+0.001); 
    let mut prev_score:f64 = -1.0*prev_drmsd; 
    for _ii in 0..iternum{
        let mut mov:f64 = 10.0;
        if _ii > iternum/2{
            mov = 1.0;
        }
        if _ii as f64 > iternum as f64 * 0.8{
            mov = 0.25;
        }
        
        random_movement(& mut tmpatoms,mov,num_moveatom,& mut rgen);
        let current_drmsd = drmsd(&tmpatoms, reference_dist);
        //let current_score:f64  = 1.0/(current_drmsd+0.001); 
        let current_score:f64  = 1.0/(current_drmsd+0.0001); 
        println!("{} {} scores: current= {}, prev= {}",_ii,prev_drmsd,current_score,prev_score);
        if current_score > prev_score{
            copy_xyz(&tmpatoms,atoms);
            prev_score = current_score;
            prev_drmsd = current_drmsd;
            //println!("swap!");
        }else{
            let dd = rgen.gen_range(0.0,1.0);
            if 1.0/(1.0+(prev_score/current_score).exp()*10.0) > dd{
                copy_xyz(&tmpatoms,atoms);
                prev_score = current_score;
                prev_drmsd = current_drmsd;
                //println!("swap!");
            }else{
                copy_xyz(atoms,&mut tmpatoms);
            }
        }
        if prev_drmsd < 1.0{
            break;
        }
    }
}


//原子ごとに評価して、原子ごとに MH 基準で採択する
pub fn mc_iter_array(atoms:&mut Vec<PDBAtom>,reference_dist:&Vec<Vec<f64>>,iternum:usize,random_seed:Option<u64>,num_moveatom:usize){
    let mut rgen:StdRng;//乱数生成器
    match random_seed{
        Some(x)=>{
            rgen = SeedableRng::seed_from_u64(x);
        },
        None => {
            rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
        }
    }
    let mut tmpatoms:Vec<PDBAtom> =vec![];
    for aa in atoms.iter(){
        let mut natom = PDBAtom::new();
        natom.set_x(aa.get_x());
        natom.set_y(aa.get_y());
        natom.set_z(aa.get_z());
        tmpatoms.push(natom);
    }
    
    let mut current_drmsds:Vec<f64> = vec![std::f64::NEG_INFINITY;atoms.len()];
    let mut prev_scores:Vec<f64> = vec![std::f64::NEG_INFINITY;atoms.len()];

    
    for _ii in 0..iternum{
        let mut mov:f64 = 10.0;
        if _ii > iternum/2{
            mov = 1.0;
        }
        if _ii as f64 > iternum as f64 * 0.8{
            mov = 0.25;
        }
        random_movement(& mut tmpatoms,mov,num_moveatom,& mut rgen);

        drmsd_array(&tmpatoms, reference_dist,&mut current_drmsds);
        let current_drmsd = drmsd(&tmpatoms, reference_dist);
        println!("{} {}",_ii,current_drmsd);
        for ii in 0..atoms.len(){
        
            let current_score = 1.0/(current_drmsds[ii]+0.00001);
            let mut accept = current_score > prev_scores[ii];
            if !accept{
                let dd = rgen.gen_range(0.0,1.0);
                if 1.0/(1.0+(prev_scores[ii]/current_score).exp()*10.0) > dd{
                    accept = true;
                }
            }
            if accept{
                prev_scores[ii] = current_score;
                atoms[ii].set_xyz(tmpatoms[ii].get_x(),tmpatoms[ii].get_y(),tmpatoms[ii].get_z());
            }else{
                tmpatoms[ii].set_xyz(atoms[ii].get_x(),atoms[ii].get_y(),atoms[ii].get_z());
            }
        }
        if current_drmsd < 0.01{
            break;
        }
    }
    let current_drmsd = drmsd(&tmpatoms, reference_dist);
    println!("result: {}",current_drmsd);
}


#[test]

fn ca_mc_test(){
    let pdbb:PDBEntry = load_pdb("D:/dummy/vbox_share/bioo/tbmtest/2r75.pdb");
    let mut ca_atoms:Vec<PDBAtom> = vec![];
    for (_cii,cc) in pdbb.chains.iter().enumerate(){
        for (_rii,rr) in cc.residues.iter().enumerate(){
            assert_eq!(_cii as i64,rr.parent_chain.unwrap());
            if rr.residue_name == "HOH"{
                continue;
            }
            for (_aii,aa) in rr.iter_atoms().enumerate(){
                if aa.atom_code == "CA"{
                    let mut tmpp:PDBAtom = PDBAtom::new();
                    tmpp.set_xyz(aa.get_x(),aa.get_y(),aa.get_z());
                    ca_atoms.push(tmpp);
                }
            }
        }
    }
    
    let alen:usize = ca_atoms.len();
    let mut ref_dist:Vec<Vec<f64>> = vec![vec![-1.0;alen];alen];
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    let mut atoms_unfolded:Vec<PDBAtom> = vec![];
    for ii in 0..alen{
        let mut tmpp:PDBAtom = PDBAtom::new();
        tmpp.set_xyz(rgen.gen_range(0.0,1.0), rgen.gen_range(0.0,1.0), rgen.gen_range(0.0,1.0));
        atoms_unfolded.push(tmpp);
        for jj in 0..alen{
            if ca_atoms[ii].distance(&ca_atoms[jj]) < 12.0{
                ref_dist[ii][jj] = ca_atoms[ii].distance(&ca_atoms[jj]);
            }
        }
    }
//mc_iter_array(&mut atoms_unfolded,&ref_dist,10000,Some(123),1);
    let num_atoms = atoms_unfolded.len();
    mc_iter_array(&mut atoms_unfolded,&ref_dist,10,Some(123),num_atoms);
    
}
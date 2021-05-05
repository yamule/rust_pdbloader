use core::f64;
use std::fs;
use std::cmp::Ordering;
use std::io::{BufWriter,Write};
use regex::Regex;
use std::collections::{HashMap,HashSet};
use crate::{geometry::{Point3D, Vector3D}, mmcif_process, pdbdata::{PDBAsym, PDBAtom,PDBComp, PDBEntry}, process_3d, structural_alignment::align};
use super::structural_alignment;
use super::matrix_process;
use chrono::{Local};

lazy_static! {
    static ref RESIDUES_DEFAULT:Vec<String> =  vec![
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
        ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].into_iter().map(|m|m.to_owned()).collect();
}


pub fn align_residue(t1:&dyn Vector3D,t2:&dyn Vector3D,t3:&dyn Vector3D
    ,s1:&dyn Vector3D,s2:&dyn Vector3D,s3:&dyn Vector3D,
    vset:&mut Vec<&mut dyn Vector3D>)->bool{
    let mut movmatbuff:Vec<Vec<f64>> = vec![vec![0.0;6];3];
    let mut tempmatbuff:Vec<Vec<f64>> = vec![vec![0.0;6];3];
    let mut movcenter:(f64,f64,f64) = (0.0,0.0,0.0);
    let mut tempcenter:(f64,f64,f64) = (0.0,0.0,0.0);

    for (ii,tt) in vec![t1,t2,t3].into_iter().enumerate(){
        movmatbuff[0][ii] = tt.get_x();
        movmatbuff[1][ii] = tt.get_y();
        movmatbuff[2][ii] = tt.get_z();
        movcenter.0 += tt.get_x()/3.0;
        movcenter.1 += tt.get_y()/3.0;
        movcenter.2 += tt.get_z()/3.0;
    }
    for (ii,tt) in vec![s1,s2,s3].into_iter().enumerate(){
        tempmatbuff[0][ii] = tt.get_x();
        tempmatbuff[1][ii] = tt.get_y();
        tempmatbuff[2][ii] = tt.get_z();
        tempcenter.0 += tt.get_x()/3.0;
        tempcenter.1 += tt.get_y()/3.0;
        tempcenter.2 += tt.get_z()/3.0;
    }
    
    let t4 = process_3d::calc_norm_t(&t1.get_xyz(),&t2.get_xyz(),&t3.get_xyz());
    let s4 = process_3d::calc_norm_t(&s1.get_xyz(),&s2.get_xyz(),&s3.get_xyz());
    //3 点だと QR 分解がうまくいかないので追加する
    for i in 3..6{
        movmatbuff[0][i] = t4.0+movmatbuff[0][i-3];
        movmatbuff[1][i] = t4.1+movmatbuff[1][i-3];
        movmatbuff[2][i] = t4.2+movmatbuff[2][i-3];

        tempmatbuff[0][i] = s4.0+tempmatbuff[0][i-3];
        tempmatbuff[1][i] = s4.1+tempmatbuff[1][i-3];
        tempmatbuff[2][i] = s4.2+tempmatbuff[2][i-3];
    }



    for ii in 0..6{
        movmatbuff[0][ii] -= movcenter.0;
        movmatbuff[1][ii] -= movcenter.1;
        movmatbuff[2][ii] -= movcenter.2;
        
        tempmatbuff[0][ii] -= tempcenter.0;
        tempmatbuff[1][ii] -= tempcenter.1;
        tempmatbuff[2][ii] -= tempcenter.2;
    }

    let rotmat_ = structural_alignment::get_rotation_matrix(&movmatbuff
        ,&tempmatbuff
        ,6);
    if let None = rotmat_{
        return false;
    }else{
        let rotmat = rotmat_.unwrap();
        let checker:Vec<&dyn Vector3D> = vec![s1,s2,s3];
        for (ii,tt) in vec![t1,t2,t3].into_iter().enumerate(){
    
            let mres = matrix_process::matrix_multi(&rotmat
            ,&vec![
            vec![tt.get_x()-movcenter.0]
            ,vec![tt.get_y()-movcenter.1]
            ,vec![tt.get_z()-movcenter.2]
            ]);
            
            let xx = mres[0][0]+tempcenter.0;
            let yy = mres[1][0]+tempcenter.1;
            let zz = mres[2][0]+tempcenter.2;
            let chkdist = process_3d::distance(&checker[ii].get_xyz(), &(xx,yy,zz));
            //println!("{}",chkdist);
           // println!("{},{},{}",xx,yy,zz);
        }    
        
        //println!("{:?}",s1.get_xyz());
        //println!("{:?}",s2.get_xyz());
        //println!("{:?}",s3.get_xyz());
        for vv in vset.iter_mut(){
            let mres = matrix_process::matrix_multi(&rotmat
                ,&vec![
                vec![vv.get_x()-movcenter.0]
                ,vec![vv.get_y()-movcenter.1]
                ,vec![vv.get_z()-movcenter.2]
                ]);
            let xx = mres[0][0]+tempcenter.0;
            let yy = mres[1][0]+tempcenter.1;
            let zz = mres[2][0]+tempcenter.2;
            vv.set_xyz(xx, yy, zz);
        }
    }
    return true;
}

pub fn generate_intermediate_files(inputdirname:&str
    ,outputdirname:&str
    ,targets:&Vec<String>
    ,cluster_diff_threshold:f64
    ,cover_ratio:f64//threshold でクラスタリングして大きい順から全体の ratio を満たすまで取る
    ){
    let paths = fs::read_dir(inputdirname).unwrap();
    let exx =  Regex::new(r"(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    let mut entries_:Vec<String> = vec![];
    for path in paths {
        if let Ok(a) = path{
            if let Some(b) = a.path().to_str(){
                if let Some(x) = exx.find(b){
                    entries_.push(b.to_string());
                }
            }
        }
    }
    entries_.sort();
    let mut entries:Vec<(String,PDBEntry)> = vec![];//filename, pdbentry
    for ee in entries_.into_iter(){
        if let Some(x) = exx.captures(&ee){
            let ext1:String = x.get(1).unwrap().as_str().to_string();
            let is_gzip = 
            if let Some(_) = x.get(2){
                true
            }else{
                false
            };
            println!("Loading {}.",ee);
            if &ext1 == ".ent" || &ext1 == ".pdb"{
                entries.push((ee.clone(),mmcif_process::load_pdb(&ee,is_gzip)));
            }else if &ext1 == ".cif"{
                entries.push((ee.clone(),mmcif_process::MMCIFEntry::load_mmcif(&ee,is_gzip)));

            }
        }
    }
    let mut targets_hm:HashMap<String,Vec<PDBComp>> = HashMap::new();
    
    for tt in targets.iter(){
        targets_hm.insert(tt.clone(),vec![]);
    }
    let dd = std::path::Path::new(outputdirname);
    if dd.exists(){
        if  !dd.is_dir(){
            panic!("{} is not a directory.",outputdirname);
        }
    }else{
        if let Err(e) = std::fs::create_dir(outputdirname){
            panic!("Cannot make {}. {:?}",outputdirname,e);
        }
    }
    let mut compcount:i64 = 0;
    for (ii,(ff,ee)) in entries.iter_mut().enumerate(){
        let asyms:Vec<&PDBAsym> = ee.get_all_asyms();
        for (jj,aa) in asyms.iter().enumerate(){
            for (_kk,cc) in aa.iter_comps().enumerate(){
                 if targets_hm.contains_key(cc.get_comp_id()){
                    let mut compp = cc.get_copy_wo_parents();
                    //set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,asym_id:Option<i64>,index:i64){
                    compp.set_parent(Some(ii as i64),None,None,compcount as i64);//sort のために使う
                    compcount += 1;
                    targets_hm.get_mut(cc.get_comp_id()).unwrap().push(compp);
                 }
            } 
        }
    }
    //CHARMM19
    let angle_radian:f64 = 117.5/180.0*std::f64::consts::PI;
    let c_bond_length:f64 = 1.52;
    let base_n:Point3D = Point3D::new(1.45,0.0,0.0);
    let base_ca:Point3D = Point3D::new(0.0,0.0,0.0);
    let base_c:Point3D = Point3D::new(
        c_bond_length*(angle_radian.cos()),
        c_bond_length*(angle_radian.sin()),
       0.0);
    let checker_threshold = 0.15;//N CA C のいずれかのずれがこれより大きい場合使用しない
    for tt in targets.iter(){
        let mut atomcounter:HashMap<String,usize> = HashMap::new();
        let numcomps:f64 = targets_hm.get(tt).unwrap().len() as f64;
        for cc in targets_hm.get(tt).unwrap().iter(){
            for aa in cc.iter_atoms(){
                let code:String = aa.get_atom_code().to_string();
                if !atomcounter.contains_key(&code){
                    atomcounter.insert(code.clone(),0);
                }
                *atomcounter.get_mut(&code).unwrap() += 1;
            }
        }

        let tthreshold = 0.9;//これより大きい割合のメンバーが持っている原子がない場合、使用しない
        //面倒なので OXT についても特別な処理はしていない
        let mut canonical_atoms:HashSet<String> = HashSet::new();
        for aa in atomcounter.iter(){
            if *aa.1 as f64/numcomps > tthreshold{
                canonical_atoms.insert(aa.0.to_string());
            }
        }
        let mut vcc:Vec<String> = canonical_atoms.iter().map(|m|m.to_string()).collect();
        vcc.sort();
        println!("{}\t{:?}",tt,vcc);
        let re_avoid = Regex::new("[^a-zA-Z0-9\\.\\-]").unwrap();
        let mut validcomps_ :Vec<(f64,i64,PDBComp)> = vec![];
        let mut validcomps:Vec<PDBComp> = vec![];
        let canonical_atoms_v:Vec<String> = canonical_atoms.iter().map(|m|m.to_string()).collect();
        let mut compss:Vec<PDBComp> = vec![];
        compss.append(targets_hm.get_mut(tt).unwrap());
        'outc:for mut cc in compss.into_iter(){
            let mut remover:Vec<String> = vec![];
            for aa in canonical_atoms_v.iter(){
                if let None = cc.get_first_atom_by_name(aa){
                    continue 'outc;
                }                
            }
            let mut bfactor_sum:f64 = 0.0;
            for aa in cc.iter_atoms(){
                if !canonical_atoms.contains(aa.get_atom_code()){
                    remover.push(aa.get_atom_code().to_string());
                }else{
                    bfactor_sum += aa.temp_factor;
                }
            }
            for rr in remover.iter(){
                cc.remove_atom_by_name(rr);
            }
            
            validcomps_.push((bfactor_sum,cc.get_index(),cc));
        } 
        validcomps_.sort_by(|a,b|{
            match a.0.partial_cmp(&b.0).unwrap() {
                Ordering::Equal => a.1.cmp(&b.1),
                other => other,
            }   
        });
        if validcomps_.len() > 0{
            let mut qr_failed:usize = 0;
            let mut checker_failed:usize = 0;
            for mut coo in validcomps_.into_iter(){
                
                let n:Point3D = Point3D::from_vector3d(coo.2.get_N().unwrap_or_else(||panic!("Can not get N for {}.",tt)));
                let ca:Point3D = Point3D::from_vector3d(coo.2.get_CA().unwrap_or_else(||panic!("Can not get CA for {}.",tt)));
                let c:Point3D = Point3D::from_vector3d(coo.2.get_C().unwrap_or_else(||panic!("Can not get C for {}.",tt)));
                let mut vvec:Vec<&mut dyn Vector3D> = vec![];
                for vv in coo.2.get_all_atoms().into_iter(){
                    vvec.push(vv);
                }

                //qr 分解で揃えるブロックもあるが今は使ってない
                let mut qralign = true;
                if true{
                    process_3d::fit_to_vector(&n,&ca,&c,&base_n, &base_ca,&base_c
                        ,&mut vvec);
                    let mut c:Point3D = Point3D::from_vector3d(coo.2.get_C().unwrap_or_else(||panic!("Can not get C for {}.",tt)));
                    let ddis = process_3d::distance(&(0.0,0.0,0.0),&c.get_xyz());
                    if ddis > 0.0{
                        c.set_xyz(c.get_x()/ddis,c.get_y()/ddis ,c.get_z()/ddis );
                    }
                    let aac = c.get_x().acos();
                    if aac.is_nan(){
                    }else{
                        let rr = (angle_radian-aac)/2.0;
                        for vv in coo.2.get_all_atoms().into_iter(){
                            let xx = rr.cos()*vv.get_x()-rr.sin()*vv.get_y();
                            let yy = rr.sin()*vv.get_x()+rr.cos()*vv.get_y();
                            vv.set_x(xx);
                            vv.set_y(yy);
                        }
                    }
                }else{
                    qralign = align_residue(&n,&ca,&c
                    ,&base_n, &base_ca,&base_c
                    ,&mut vvec);
                }
                if qralign{
                    let n:Point3D = Point3D::from_vector3d(coo.2.get_N().unwrap_or_else(||panic!("Can not get N for {}.",tt)));
                    let ca:Point3D = Point3D::from_vector3d(coo.2.get_CA().unwrap_or_else(||panic!("Can not get CA for {}.",tt)));
                    let c:Point3D = Point3D::from_vector3d(coo.2.get_C().unwrap_or_else(||panic!("Can not get C for {}.",tt)));
                    if n.distance(&base_n) > checker_threshold
                    || ca.distance(&base_ca) > checker_threshold
                    || c.distance(&base_c) > checker_threshold
                    {  
                        checker_failed += 1;
                    }else{
                        validcomps.push(coo.2);
                    }
                    
                }else{
                    qr_failed += 1;
                }
                
            }
            let num_validcomps:usize = validcomps.len();
            
            let mut cluster:Vec<(PDBComp,usize)> = vec![];
            for vv in validcomps.into_iter(){
                let mut merge_to:i64 = -1;
                for (ii,ccc) in cluster.iter().enumerate(){
                    let mut mergeflag = true;
                    for an in canonical_atoms_v.iter(){
                        if ccc.0.get_first_atom_by_name(an).unwrap().distance(vv.get_first_atom_by_name(an).unwrap()) >= cluster_diff_threshold{
                            mergeflag = false;
                            break;
                        }
                    }
                    if mergeflag{
                        merge_to = ii as i64;
                        break;
                    }
                }
                if merge_to == -1{
                    cluster.push((vv,1));
                }else{
                    cluster[merge_to as usize].1 += 1;
                }
                /*
                
                */
            }
            cluster.sort_by(|a,b|b.1.cmp(&a.1));
            let filename = outputdirname.to_string()+"/"+re_avoid.replace_all(tt,"_").to_string().as_str()+".rotamer.dat";
            let mut f = BufWriter::new(fs::File::create(filename).unwrap());
            f.write_all(format!("#compound: {}\n",tt).as_bytes()).unwrap();
            f.write_all(format!("#date: {}\n",Local::now()).as_bytes()).unwrap();
            f.write_all(format!("#num_validcomps: {}\n",num_validcomps).as_bytes()).unwrap();
            f.write_all(format!("#covered: {}\n",cover_ratio).as_bytes()).unwrap();
            f.write_all(format!("#threshold: {}\n",cluster_diff_threshold).as_bytes()).unwrap();
            f.write_all(format!("#results: accepted: {}, failed: {} (qr: {}, checker: {})\n",num_validcomps,qr_failed+checker_failed,qr_failed,checker_failed).as_bytes()).unwrap();
            f.write_all("#note: Asyms are changed from the original.\n".as_bytes()).unwrap();
            
            let mut ccount:usize = 0;
            for mut cc in cluster.into_iter(){
                f.write_all(format!("//==\n").as_bytes()).unwrap();
                f.write_all(format!("#file: {}\n",entries[cc.0.get_parent_entry().unwrap() as usize].0).as_bytes()).unwrap();
                f.write_all(format!("#ratio: {}\n",(cc.1 as f64)/(num_validcomps as f64)).as_bytes()).unwrap();
                for aaa in cc.0.iter_mut_atoms(){
                    if aaa.get_serial_number() > 99999{
                        aaa.set_serial_number(99999);
                    }
                    f.write_all(format!("{}",aaa.get_pdb_atom_line_string("A",tt, 1,"")+"\n").as_bytes()).unwrap();
                }
                f.write_all("//\n".as_bytes()).unwrap();
                ccount += cc.1;
                if ccount as f64 >= (num_validcomps as f64)*cover_ratio{
                    break;
                }
            }

            
        }
    }

    
    /*
    for tt in targets.iter(){
        

    }
    */
}

#[test]
fn dirloadtest(){
    generate_intermediate_files("example_files","example_files/example_output",&RESIDUES_DEFAULT,0.5,0.8);
    let re_avoid = Regex::new("[^a-zA-Z0-9\\.\\-]").unwrap();
}

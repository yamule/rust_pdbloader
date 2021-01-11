
#[allow(unused_imports)]
use super::debug_env;
use super::pdbdata::*;

#[allow(unused_imports)]
use super::mmcif_process;
use super::geometry::Vector3D;
use super::geometry::Point3D;
use std::collections::HashSet;
#[allow(unused_imports)]
use std::collections::HashMap;

#[derive(Clone,Copy)]
pub struct BackboneAtoms{
    missing:bool,
    n:Point3D,
    h:Option<Point3D>,
    ca:Point3D,
    c:Point3D,
    o:Option<Point3D>,
}

impl BackboneAtoms{
    pub fn new()->BackboneAtoms{
        return BackboneAtoms{
            missing:false,
            n:Point3D::new(0.0,0.0,0.0),
            h:None,
            ca:Point3D::new(0.0,0.0,0.0),
            c:Point3D::new(0.0,0.0,0.0),
            o:None,
        };
    }
    pub fn set_h(&mut self,p:Point3D){
        self.h = Some(p);
    }
    pub fn set_o(&mut self,p:Point3D){
        self.o = Some(p);
    }


}
//Threshold は 3.0 くらいかな
pub fn assing_secstr(chains:&Vec<Vec<&PDBComp>>,hbond_threshold:f64)->Vec<(String,Vec<Vec<usize>>)>{
    let mut ret:Vec<(String,Vec<Vec<usize>>)> = vec![];
    let region_length_threshold:usize = 3;//アサインされる最小の長さ
    //O と H は再構築する
    let estimate_oh =|a1:&dyn Vector3D,a2:&dyn Vector3D,a3:&dyn Vector3D,len:f64|->Point3D{
        let mut p1 = Point3D::new(
            a1.get_x() - a2.get_x(),
            a1.get_y() - a2.get_y(),
            a1.get_z() - a2.get_z());
        let mut p3 = Point3D::new(
            a3.get_x() - a2.get_x(),
            a3.get_y() - a2.get_y(),
            a3.get_z() - a2.get_z()
        );
        p1.standardize();
        p3.standardize();

        let mut p0 =Point3D::new(
            - p1.get_x()/2.0 - p3.get_x()/2.0,
            - p1.get_y()/2.0 - p3.get_y()/2.0,
            - p1.get_z()/2.0 - p3.get_z()/2.0
        );
        p0.standardize();
        p0.set_x(p0.get_x()*len+a2.get_x());
        p0.set_y(p0.get_y()*len+a2.get_y());
        p0.set_z(p0.get_z()*len+a2.get_z());
        return p0;
    };
    for residues in chains.iter(){
        let rlen:usize = residues.len();
        let mut atoms:Vec<BackboneAtoms> = vec![BackboneAtoms::new();residues.len()];
        let mut betasheet_segments:Vec<Vec<usize>> = vec![];
        for ri in 0..rlen{

            let prev_c:Option<&PDBAtom> = if ri > 0{
                residues[ri-1].get_C()
            }else{
                None
            };
            
            let next_n:Option<&PDBAtom> = if ri < rlen-1{
                residues[ri+1].get_N()
            }else{
                None
            };

            let rr:&PDBComp = &residues[ri];
            let aa:&mut BackboneAtoms = &mut atoms[ri];
            
            if let Some(nn) = rr.get_N(){
                aa.n.set_vector(nn);
                if let Some(caa) = rr.get_CA(){
                    if let Some(cc) = prev_c{
                        aa.set_h(estimate_oh(cc,nn,caa,0.98));
                        //if let Some(x) = rr.get_first_atom_by_name("H"){
                        //    println!("{}",aa.h.unwrap().distance(x));
                        //}
                    }
                }
            }else{
                aa.missing = true;
            }
            
            if let Some(x) = rr.get_CA(){
                aa.ca.set_vector(x);
            }else{
                aa.missing = true;
            }
            
            if let Some(cc) = rr.get_C(){
                aa.c.set_vector(cc);
                if let Some(caa) = rr.get_CA(){
                    if let Some(nn) = next_n{
                        if let Some(x) = rr.get_O(){
                            aa.set_o(Point3D::from_vector3d(x));
                        }else{
                            aa.set_o(estimate_oh(caa,cc,nn,1.23));//ちょっとズレるがどっちに動かしていいのか分からない・・・
                        }
                    }
                }
            }else{
                aa.missing = true;
            }
            
            //charmm19
            //H    NH1    405.0       0.98!  GELIN AND IR STRETCH 3200 CM 1
            //C    O      580.0       1.23
            //ANGLE
            //NH1  C    O       65.0     121.0
            //CH1E C    O       85.0     121.5
            //CH1E NH1  H       35.0     120.0
            //C    NH1  H       30.0     120.0
            //X    C    CH1E X        0.0       3       0.0! FROM GELIN THESIS AMIDES

        }
        let mut acc_donn:Vec<HashSet<usize>> = vec![HashSet::new();rlen];
        let mut donn_acc:Vec<HashSet<usize>> = vec![HashSet::new();rlen];
        for r1 in 0..rlen{
            for r2 in 0..rlen{
                if (r1 as i64 - r2 as i64).abs() < 3{
                    continue;
                }
                if let Some(x) = atoms[r1].h{
                    if let Some(y) = atoms[r2].o.as_ref(){
                        if x.distance(y) < hbond_threshold{
                            acc_donn[r2].insert(r1);
                            donn_acc[r1].insert(r2);
                        }
                    }
                }
            }
        }


        let mut helix_a:Vec<Vec<usize>> = vec![];
        let mut helix_3:Vec<Vec<usize>> = vec![];
        let mut helix_p:Vec<Vec<usize>> = vec![];
        for xskip in vec![3,4,5]{
            let mut counted:HashSet<usize> = HashSet::new();
            for r1 in 0..(rlen-1){
                if counted.contains(&r1){
                    continue;
                }
                if !acc_donn[r1].contains(&(r1+xskip)){
                    continue;
                }
                let mut hend:usize = r1;
                for r2 in (r1+1)..rlen{
                    counted.insert(r2);
                    if !acc_donn[r2].contains(&(r2+xskip)){
                        hend = r2-1;
                        break;
                    }
                    hend = r2;
                }
                if hend-r1+1 >= region_length_threshold{
                    let mut vt:Vec<usize> = vec![];
                    for ii in (r1+1)..=(hend+xskip-1){
                        vt.push(ii);
                    }
                    if xskip == 4{
                        helix_a.push(vt);
                    }else if xskip == 3{
                        helix_3.push(vt);
                    }else if xskip == 5{
                        helix_p.push(vt);
                    }
                }
            }
        }
        
        let mut beta_sheet_flag:Vec<usize> = vec![0;rlen];
        for rk in 1..rlen{//beta sheets
            let r1:usize = rlen-rk;
            let mut lasthit:usize = r1;
            for vv in donn_acc[r1].iter(){//anti-parallel
                let mut vreg:HashSet<usize> = HashSet::new();
                vreg.insert(r1);
                let mut vreg_rev:HashSet<usize> = HashSet::new();
                vreg_rev.insert(*vv);
                if !donn_acc[*vv].contains(&(r1 as usize)){
                    continue;
                }

                for s in (2..(r1+1)).step_by(2){
                    let ps = *vv as i64 + s as i64;
                    if ps > -1  && ps < rlen as i64{
                        let r2:usize = r1-s;
                        let mut bflag:bool = false;
                        if donn_acc[r2].contains(&(ps as usize)) && donn_acc[ps as usize].contains(&(r2 as usize)){
                            bflag = true;
                        }
                        if bflag{
                            vreg.insert(r2);
                            vreg_rev.insert(ps as usize);
                            lasthit = r2;
                        }else{
                            break;
                        }
                    }
                }
                if lasthit != r1{
                    let maxv:usize = vreg.iter().fold(0,|s,m|s.max(*m));
                    let maxr:usize = vreg_rev.iter().fold(0,|s,m|s.max(*m));
                    let minv:usize = vreg.iter().fold(rlen,|s,m|s.min(*m));
                    let minr:usize = vreg_rev.iter().fold(rlen,|s,m|s.min(*m));
                    betasheet_segments.push((minv..(maxv+1)).into_iter().collect());
                    for i in minv..(maxv+1){
                        beta_sheet_flag[i] = 1;
                    }

                    betasheet_segments.push((minr..(maxr+1)).into_iter().collect());
                    for i in minr..(maxr+1){
                        beta_sheet_flag[i] = 1;
                    }
                }
            }

            for vv in acc_donn[r1].iter(){//parallel
                let mut vreg:HashSet<usize> = HashSet::new();
                vreg.insert(r1);
                let mut vreg_rev:HashSet<usize> = HashSet::new();
                if *vv as usize + 2 < rlen{
                    if donn_acc[r1].contains(&(*vv as usize +2)){
                        vreg_rev.insert(*vv+2);
                    }
                }else{
                    continue;
                }
                vreg_rev.insert(*vv);
                for s in (2..(r1+1)).step_by(2){
                    let ps = *vv as i64 - s as i64;
                    let ps2 = *vv as i64 - s as i64 -2;
                    if ps > -1  && ps2 > -1{
                        let r2:usize = r1-s;
                        let mut bflag:bool = false;
                        if acc_donn[r2].contains(&(ps as usize)) && donn_acc[r2].contains(&(ps2 as usize)){
                            bflag = true;
                        }
                        if bflag{
                            vreg.insert(r2);
                            vreg_rev.insert(ps as usize);
                            vreg_rev.insert(ps2 as usize);
                            lasthit = r2;
                        }else{
                            if acc_donn[r2].contains(&(ps as usize)){
                                vreg.insert(r2);
                                vreg_rev.insert(ps as usize);
                            }
                            break;
                        }
                    }
                }
                if lasthit != r1{
                    let maxv:usize = vreg.iter().fold(0,|s,m|s.max(*m));
                    let maxr:usize = vreg_rev.iter().fold(0,|s,m|s.max(*m));
                    let minv:usize = vreg.iter().fold(rlen,|s,m|s.min(*m));
                    let minr:usize = vreg_rev.iter().fold(rlen,|s,m|s.min(*m));
                    betasheet_segments.push((minv..(maxv+1)).into_iter().collect());
                    for i in minv..(maxv+1){
                        beta_sheet_flag[i] = 1;
                    }
                    betasheet_segments.push((minr..(maxr+1)).into_iter().collect());
                    for i in minr..(maxr+1){
                        beta_sheet_flag[i] = 1;
                    }
                }
            }
        }
        
        let mut marks:Vec<String> = vec!["C".to_owned();rlen];
        for (_vii,v) in helix_3.iter().enumerate(){
            for vv in v.iter(){
                //marks[*vv] = "[".to_owned()+&vii.to_string()+"]";
                marks[*vv] = "3".to_owned();
            }
        }
        
        for (_vii,v) in helix_a.iter().enumerate(){
            for vv in v.iter(){
                //marks[*vv] = "[".to_owned()+&vii.to_string()+"]";
                marks[*vv] = "H".to_owned();
            }
        }
        
        for (_vii,v) in helix_p.iter().enumerate(){
            for vv in v.iter(){
                //marks[*vv] = "[".to_owned()+&vii.to_string()+"]";
                marks[*vv] = "P".to_owned();
            }
        }
        for (vii,vv) in beta_sheet_flag.iter().enumerate(){
            if *vv != 0{
                marks[vii] = "E".to_owned();
            }
        }
        let mut segments:Vec<Vec<usize>> = vec![];
        segments.append(&mut helix_3);
        segments.append(&mut helix_a);
        segments.append(&mut helix_p);

        for ii in 0..betasheet_segments.len(){
            let mut bflag = true;
            for jj in 0..betasheet_segments.len(){//冗長だが眠い
                if ii == jj{
                    continue;
                }
                if betasheet_segments[ii][0] >= betasheet_segments[jj][0] 
                    && betasheet_segments[ii].last().unwrap() <= betasheet_segments[jj].last().unwrap(){
                    if betasheet_segments[ii][0] == betasheet_segments[jj][0] 
                        && betasheet_segments[ii].last().unwrap() == betasheet_segments[jj].last().unwrap(){
                            if ii < jj{
                                break;
                            }
                    }
                    bflag = false;
                    break;
                }
            }
            if bflag{
                segments.push(betasheet_segments[ii].clone());
            }
        }
        ret.push((marks.iter().fold("".to_owned(),|s,m|s+m),segments));
    }
    return ret;
}

#[test]
fn secstr_test(){
    //let pdb = load_pdb((debug_env::EXAMPLE_DIR.to_string()+"3rgk_A.pdb").as_str());
    //let pdb = load_pdb((debug_env::EXAMPLE_DIR.to_string()+"2gx4_A.pdb").as_str());
    let pdb = mmcif_process::load_pdb((debug_env::EXAMPLE_DIR.to_string()+"1EFH_A.pdb").as_str());
    //let pdb = load_pdb((debug_env::EXAMPLE_DIR.to_string()+"6iws_model1.pdb").as_str());
    
    let mut ress:Vec<Vec<&PDBComp>> = vec![];
    for cc in pdb.get_model_at(0).get_entity_at(0).iter_asyms(){
        let mut rss:Vec<&PDBComp> = vec![];
        for rr in cc.iter_comps(){
            if !rr.get_atom_at(0).is_ligand{
                rss.push(rr);
            }
        }
        ress.push(rss);
    }
    let res = assing_secstr(&ress,2.8);
    for rr in res.iter(){
        println!("{}",rr.0);
    }
}
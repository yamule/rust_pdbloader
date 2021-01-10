extern crate regex;
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use std::fs;
#[allow(unused_imports)]
use std::slice::IterMut;
#[allow(unused_imports)]
use std::slice::Iter;
#[allow(unused_imports)]
use std::collections::HashSet;
#[allow(unused_imports)]
use std::collections::HashMap;
#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};

#[allow(unused_imports)]
use rand::SeedableRng;
#[allow(unused_imports)]
use rand::rngs::StdRng;
#[allow(unused_imports)]
use rand::prelude::*;


#[allow(unused_imports)]
use self::regex::Regex;//module ファイル内で extern crate した場合 self が必要。lib.rs で extern crate すると必要ない。

#[allow(unused_imports)]
use super::geometry;
use super::geometry::Vector3D;
use super::geometry::Point3D;
use super::pdbdata;
use super::charmm_param;
use super::energy_function;

#[allow(unused_imports)]
use super::debug_env;

#[allow(unused_imports)]
use super::backbone_sample;

#[allow(unused_imports)]
use super::side_chain_sample;

#[allow(unused_imports)]
use super::chain_builder;

#[allow(unused_imports)]
use super::process_3d;

#[allow(unused_imports)]
use super::misc_util::*;

use super::energy_function::*;

use std::f64::consts::PI;
const DEGREE_TO_RADIAN2:f64 = PI/180.0*PI/180.0;


pub const EPSILON:f64 = 1.0e-20;


#[derive(Debug,Clone)]
pub enum HybridOrbit{
    SP3,
    SP2,
    Other
}


impl Vector3D for MDAtom{

    fn get_xyz(&self)-> (f64,f64,f64){
        return (self.get_x(),self.get_y(),self.get_z());
    }

    fn set_vector(&mut self,v:&dyn Vector3D){
        self.set_x(v.get_x());
        self.set_y(v.get_y());
        self.set_z(v.get_z());
    }

    fn set_x(&mut self,xx:f64){
        self.x = xx;
    }
    fn set_y(&mut self,yy:f64){
        self.y = yy;
    }
    fn set_z(&mut self,zz:f64){
        self.z = zz;
    }
    fn get_x(&self)->f64{
        return self.x;
    }
    fn get_y(&self)->f64{
        return self.y;
    }
    fn get_z(&self)->f64{
        return self.z;
    }
    
    fn set_xyz(&mut self,xx:f64,yy:f64,zz:f64){
        self.set_x(xx);
        self.set_y(yy);
        self.set_z(zz);
    }
    fn distance(&self,point:&dyn Vector3D)->f64{
        let dd = (self.get_x() - point.get_x()).powf(2.0)+
        (self.get_y() - point.get_y()).powf(2.0)+
        (self.get_z() - point.get_z()).powf(2.0);
        if dd == 0.0{
            return 0.0;
        }
        return dd.sqrt();
    }
}

#[derive(Debug,Clone)]
pub struct BondVars{
    pub atoms:(usize,usize),
    pub kb:f64,
    pub b0:f64
}
impl energy_function::EnergyFunction for BondVars{
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        let dss = mdenv.dist[self.atoms.0][self.atoms.1]-self.b0;
        let dsc = self.kb*dss*dss*weight;
        atom_level_energy[self.atoms.0] += dsc/2.0;
        atom_level_energy[self.atoms.1] += dsc/2.0;
        return dsc;
    }
}

#[derive(Debug,Clone)]
pub struct AngleVars{
    pub atoms:(usize,usize,usize),
    ktheta:f64,
    theta0:f64,
}
impl energy_function::EnergyFunction for AngleVars{
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        let dss = calc_angle(&mdenv.atoms[self.atoms.0]
            , &mdenv.atoms[self.atoms.1]
            ,&mdenv.atoms[self.atoms.2])-self.theta0;
        let dsc = self.ktheta*dss*dss *DEGREE_TO_RADIAN2*weight;
        //let dsc = self.ktheta*dss*dss*weight;
        atom_level_energy[self.atoms.0] += dsc/3.0;
        atom_level_energy[self.atoms.1] += dsc/3.0;
        atom_level_energy[self.atoms.2] += dsc/3.0;
        return dsc;
    }
}

#[derive(Debug,Clone)]
pub struct UBVars{
    pub atoms:(usize,usize,usize),
    kub:f64,
    s0:f64,
}
impl energy_function::EnergyFunction for UBVars{
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        let dsc = self.kub*(mdenv.dist[self.atoms.0][self.atoms.2] - self.s0)*(mdenv.dist[self.atoms.0][self.atoms.2] - self.s0)*weight;
        atom_level_energy[self.atoms.0] += dsc/3.0;
        atom_level_energy[self.atoms.1] += dsc/3.0;
        atom_level_energy[self.atoms.2] += dsc/3.0;
        
        return dsc;
    }
}

#[derive(Debug,Clone)]
pub struct DihedralVars{
    pub atoms:(usize,usize,usize,usize),
    
    kchi:f64,
    n:f64,//usize だが変換が面倒なので
    delta:f64,
    debug_string:String
}
impl energy_function::EnergyFunction for DihedralVars{
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{

        let phi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.0], &mdenv.atoms[self.atoms.1],&mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3]);
        let dsc:f64 = self.kchi*(1.0+((self.n*phi-self.delta)/180.0*PI).cos())*weight;
        
        atom_level_energy[self.atoms.0] += dsc/4.0;
        atom_level_energy[self.atoms.1] += dsc/4.0;
        atom_level_energy[self.atoms.2] += dsc/4.0;
        atom_level_energy[self.atoms.3] += dsc/4.0;
        
        return dsc;
    }
}

#[derive(Debug,Clone)]
pub struct IMPRVars{
    pub atoms:(usize,usize,usize,usize),
    kpsi:f64,
    psi0:f64
}

impl energy_function::EnergyFunction for IMPRVars{
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        let psi = calc_dihedral_angle(
            &mdenv.atoms[self.atoms.0]
            ,&mdenv.atoms[self.atoms.1]
            ,&mdenv.atoms[self.atoms.2]
            ,&mdenv.atoms[self.atoms.3])*-1.0;
        let dsc = self.kpsi*(psi-self.psi0)*(psi-self.psi0) *DEGREE_TO_RADIAN2*weight;
        atom_level_energy[self.atoms.0] += dsc/4.0;
        atom_level_energy[self.atoms.1] += dsc/4.0;
        atom_level_energy[self.atoms.2] += dsc/4.0;
        atom_level_energy[self.atoms.3] += dsc/4.0;
        return dsc;
    }
}


#[derive(Debug,Clone)]
#[allow(dead_code)]
pub struct CMAPVars{
    atoms:(usize,usize,usize,usize),
    matrixindex:usize,//Matrix 全部コピーするのは冗長なので
}

impl energy_function::EnergyFunction for CMAPVars{
    
#[allow(unused_variables)]
    fn calc_energy(&self,mdenv:&CharmmEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        //CMAP は別の方法使った方が良いと思う
        panic!("not implemented yet");
    }
}



//Unplaced な原子の場所を Param に登録された Torsion angle、angle、Bond Length から見積もる
pub fn estimate_positions_unplaced(mdenv:&mut CharmmEnv,cvars:&CharmmVars){
    let mut unplaced:Vec<usize> = vec![];
    let mut all_placed_points:Vec<usize> = vec![];
    for ii in 0..mdenv.atoms.len(){
        if mdenv.atoms[ii].unplaced{
            unplaced.push(ii);
        }else{
            all_placed_points.push(ii);
        }
    }
    if unplaced.len() == 0{
        return;
    }
    let anum = mdenv.atoms.len();
    let mut dihedd:Vec<Vec<&DihedralVars>> = vec![vec![];anum];
    let mut imprr:Vec<Vec<&IMPRVars>> = vec![vec![];anum];
    let mut anglee:Vec<Vec<&AngleVars>> = vec![vec![];anum];
    let mut bondd:Vec<Vec<&BondVars>> = vec![vec![];anum];

    for (_eii,ee) in cvars.dihedvec.iter().enumerate(){
        dihedd[ee.atoms.0].push(ee);
        dihedd[ee.atoms.1].push(ee);
        dihedd[ee.atoms.2].push(ee);
        dihedd[ee.atoms.3].push(ee);
    }




    for (_eii,ee) in cvars.imprvec.iter().enumerate(){
        imprr[ee.atoms.0].push(ee);
        imprr[ee.atoms.1].push(ee);
        imprr[ee.atoms.2].push(ee);
        imprr[ee.atoms.3].push(ee);
    }
    
    for (_eii,ee) in cvars.anglevec.iter().enumerate(){
        anglee[ee.atoms.0].push(ee);
        anglee[ee.atoms.1].push(ee);
        anglee[ee.atoms.2].push(ee);
    }
    for (_eii,ee) in cvars.bondvec.iter().enumerate(){
        bondd[ee.atoms.0].push(ee);
        bondd[ee.atoms.1].push(ee);
    }
    loop{
        let mut assigned:i64 = -1;
        'outer:for (uii,uu) in unplaced.iter().enumerate(){
            let dihed:&Vec<&DihedralVars> = &dihedd[*uu];
            for dd in dihed.iter(){
                let multiplicity = (dd.n+0.00001) as usize;
                let muldeg:f64 = 360.0/dd.n;

                if dd.atoms.3 == *uu{
                    
                    if mdenv.atoms[dd.atoms.0].unplaced 
                    ||  mdenv.atoms[dd.atoms.1].unplaced 
                    ||  mdenv.atoms[dd.atoms.2].unplaced 
                    {
                        continue;
                    }
                    let blen_ = get_bond_length(dd.atoms.2,dd.atoms.3,&bondd[*uu]);
                    let ang_ = get_angle(dd.atoms.1,dd.atoms.2,dd.atoms.3,&anglee[*uu]);

                    if let None = blen_{
                        continue;
                    }
                    if let None = ang_{
                        continue;
                    }
                    let blen = blen_.unwrap();
                    let mut ang = ang_.unwrap();
                    if ang.1{
                        ang.0 *= -1.0;
                    }
                    let mut maxdist:f64 = -9999.0;
                    let mut maxpos:(f64,f64,f64) = (0.0,0.0,0.0);
                    let mut negcheck_debug:f64 = -9999.0;
                    //デバッグ中
                    for mull in 1..=multiplicity{
                        for negg in vec![1.0,-1.0].iter(){
                                let poss = estimate_3_dihed(
                            &mdenv.atoms[dd.atoms.0].get_xyz()
                            , &mdenv.atoms[dd.atoms.1].get_xyz()
                            , &mdenv.atoms[dd.atoms.2].get_xyz()
                            , 1.0
                            , 120.0
                            , dd.delta+muldeg*(mull as f64 -1.0)
                            , ang.0*negg 
                            , blen.0);
                            let vpp:Point3D = Point3D::new(poss.0,poss.1,poss.2);
                            let mut mclose:f64 = std::f64::INFINITY;
                            for ll in all_placed_points.iter(){

                                //結合している原子は取らない
                                if *ll == dd.atoms.2{
                                    continue;
                                }
                                let ddis = mdenv.atoms[*ll].distance(&vpp);
                                if ddis < mclose{
                                    mclose = ddis;
                                }
                            }
                            if maxdist < 0.0{
                                maxdist = mclose;
                                maxpos = poss;
                                negcheck_debug = *negg;
                            }else if maxdist < mclose{
                                maxdist = mclose;
                                maxpos = poss;
                                negcheck_debug = *negg;
                            }
                        }

                    }
                    mdenv.atoms[*uu].set_xyz(maxpos.0,maxpos.1,maxpos.2);
                    mdenv.atoms[*uu].unplaced = false;
                    all_placed_points.push(*uu);
                    assigned = uii as i64;
                    
                    println!("***{} {} {} {} {}",negcheck_debug,multiplicity,mdenv.atoms[*uu].residue_name,mdenv.atoms[*uu].atom_name,dd.debug_string);
                }else if dd.atoms.0 == *uu{
                    if mdenv.atoms[dd.atoms.3].unplaced 
                    ||  mdenv.atoms[dd.atoms.2].unplaced 
                    ||  mdenv.atoms[dd.atoms.1].unplaced 
                    {
                        continue;
                    }
                    let blen_ = get_bond_length(dd.atoms.0,dd.atoms.1,&bondd[*uu]);
                    let ang_ = get_angle(dd.atoms.2,dd.atoms.1,dd.atoms.0,&anglee[*uu]);
                    if let None = blen_{
                        continue;
                    }
                    if let None = ang_{
                        continue;
                    }
                    let blen = blen_.unwrap();
                    let mut ang = ang_.unwrap();
                    if ang.1{
                        ang.0 *= -1.0;
                    }
                    
                    let multiplicity = (dd.n+0.00001) as usize;
                    let mut maxdist:f64 = -999.0;
                    let mut maxpos:(f64,f64,f64) = (0.0,0.0,0.0);
                    let mut negcheck_debug:f64 = -9999.0;
                    for mull in 1..=multiplicity{
                        
                        for negg in vec![1.0,-1.0].iter(){
                            let poss = estimate_3_dihed(
                            &mdenv.atoms[dd.atoms.3].get_xyz()
                            , &mdenv.atoms[dd.atoms.2].get_xyz()
                            , &mdenv.atoms[dd.atoms.1].get_xyz()
                            , 1.0
                            , 120.0
                            , dd.delta+muldeg*(mull as f64 -1.0)
                            , ang.0*negg
                            , blen.0);
                            let vpp:Point3D = Point3D::new(poss.0,poss.1,poss.2);
                            let mut mclose:f64 = std::f64::INFINITY;
                            for ll in all_placed_points.iter(){

                                //結合している原子は取らない
                                if *ll == dd.atoms.1{
                                    continue;
                                }
                                let ddis = mdenv.atoms[*ll].distance(&vpp);
                                if ddis < mclose{
                                    mclose = ddis;
                                }
                            }

                            if maxdist < 0.0{
                                maxdist = mclose;
                                maxpos = poss;
                                negcheck_debug = *negg;
                            }else if maxdist < mclose{
                                maxdist = mclose;
                                maxpos = poss;
                                negcheck_debug = *negg;
                            }
                        }
                    }
                    mdenv.atoms[*uu].set_xyz(maxpos.0,maxpos.1,maxpos.2);
                    all_placed_points.push(*uu);
                    mdenv.atoms[*uu].unplaced = false;
                    assigned = uii as i64;
                    println!("****{} {} {} {}",negcheck_debug,mdenv.atoms[*uu].residue_name,mdenv.atoms[*uu].atom_name,dd.debug_string);
                }
                if assigned > -1{
                    break 'outer;
                }
            }

            //https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node25.html
            //によれば ic セクションとは計算の仕方が別?
            //The final bond-like terms in the parameter file are impropers, which are used exclusively and explicitly in the molecular topology to maintain planarity. As such, the harmonic form $K_\psi (\psi - \psi_0)^2$ with a large spring constant and $\psi_0$ typically zero is used to restrain deformations among an atom and three atoms bonded to it. As with dihedrals, $\psi $ is angle between the plane containing the first three atoms and the plane containing the last three. Notice below that wildcard atom types occur in the second and third positions, rather than the first and fourth as in dihedrals. い
            // 0-1-2 と 1-2-3 のなす角
            //ちょっとこれは計算が面倒くさすぎるので保留。
            /*
            let impp:&Vec<&IMPRVars> = &imprr[*uu];
            //println!("{}",impp.len());
            for dd in impp.iter(){
                
                if dd.atoms.3 == *uu{
                    
                    //par の場合は first atom が central atom とのこと
                    //http://www.charmm-gui.org/charmmdoc/parmfile.html
                    if mdenv.atoms[dd.atoms.0].unplaced 
                    ||  mdenv.atoms[dd.atoms.1].unplaced 
                    ||  mdenv.atoms[dd.atoms.2].unplaced 
                    {
                        continue;
                    }
                    let blen_ = get_bond_length(dd.atoms.0,dd.atoms.3,&bondd[*uu]);
                    let ang_ = get_angle(dd.atoms.2,dd.atoms.0,dd.atoms.3,&anglee[*uu]);
                    if let None = blen_{
                        continue;
                    }
                    if let None = ang_{
                        continue;
                    }
                    let blen = blen_.unwrap();
                    let mut ang = ang_.unwrap();
                    if ang.1{
                        ang.0 *= -1.0;
                    }
                    let mut maxdist:f64 = -9999.0;
                    let mut maxpos:(f64,f64,f64) = (0.0,0.0,0.0);
                    for negg in vec![1.0,-1.0].iter(){
                        let poss = estimate_3_impr(
                        &mdenv.atoms[dd.atoms.1].get_xyz()
                        , &mdenv.atoms[dd.atoms.2].get_xyz()
                        , &mdenv.atoms[dd.atoms.0].get_xyz()
                        , 1.0
                        , 120.0
                        , dd.psi0
                        , ang.0*negg 
                        , blen.0);
                        let vpp:Point3D = Point3D::new(poss.0,poss.1,poss.2);
                        let mut mclose:f64 = std::f64::INFINITY;
                        for ll in all_placed_points.iter(){

                            //結合している原子は取らない
                            if *ll == dd.atoms.0{
                                continue;
                            }
                            let ddis = mdenv.atoms[*ll].distance(&vpp);
                            if ddis < mclose{
                                mclose = ddis;
                            }
                        }
                        if maxdist < 0.0{
                            maxdist = mclose;
                            maxpos = poss;
                        }else if maxdist < mclose{
                            maxdist = mclose;
                            maxpos = poss;
                        }
                    }

                
                    mdenv.atoms[*uu].set_xyz(maxpos.0,maxpos.1,maxpos.2);
                    mdenv.atoms[*uu].unplaced = false;
                    all_placed_points.push(*uu);
                    assigned = uii as i64;
                    println!("!!{} {} {}",mdenv.atoms[*uu].residue_name,mdenv.atoms[*uu].atom_name,"**");

                }else if dd.atoms.1 == *uu{
                    if mdenv.atoms[dd.atoms.0].unplaced 
                    ||  mdenv.atoms[dd.atoms.2].unplaced 
                    ||  mdenv.atoms[dd.atoms.3].unplaced 
                    {
                        continue;
                    }
                    let blen_ = get_bond_length(dd.atoms.0,dd.atoms.1,&bondd[*uu]);
                    let ang_ = get_angle(dd.atoms.2,dd.atoms.0,dd.atoms.1,&anglee[*uu]);
                    if let None = blen_{
                        continue;
                    }
                    if let None = ang_{
                        continue;
                    }
                    let blen = blen_.unwrap();
                    let mut ang = ang_.unwrap();
                    if ang.1{
                        ang.0 *= -1.0;
                    }
                    let mut maxdist:f64 = -9999.0;
                    let mut maxpos:(f64,f64,f64) = (0.0,0.0,0.0);
                    for negg in vec![1.0,-1.0].iter(){
                        let poss = estimate_3_impr(
                        &mdenv.atoms[dd.atoms.3].get_xyz()
                        , &mdenv.atoms[dd.atoms.2].get_xyz()
                        , &mdenv.atoms[dd.atoms.0].get_xyz()
                        , 1.0
                        , 120.0
                        , dd.psi0*-1.0
                        , ang.0*negg 
                        , blen.0);
                        let vpp:Point3D = Point3D::new(poss.0,poss.1,poss.2);
                        let mut mclose:f64 = std::f64::INFINITY;
                        for ll in all_placed_points.iter(){

                            //結合している原子は取らない
                            if *ll == dd.atoms.0{
                                continue;
                            }
                            let ddis = mdenv.atoms[*ll].distance(&vpp);
                            if ddis < mclose{
                                mclose = ddis;
                            }
                        }
                        if maxdist < 0.0{
                            maxdist = mclose;
                            maxpos = poss;
                        }else if maxdist < mclose{
                            maxdist = mclose;
                            maxpos = poss;
                        }
                    }

                
                    mdenv.atoms[*uu].set_xyz(maxpos.0,maxpos.1,maxpos.2);
                    mdenv.atoms[*uu].unplaced = false;
                    all_placed_points.push(*uu);
                    assigned = uii as i64;
                    println!("!!{} {} {}",mdenv.atoms[*uu].residue_name,mdenv.atoms[*uu].atom_name,"**");

                }
                if assigned > -1{
                    break 'outer;
                }
            }
            */

        }
        
        if assigned > -1{
            unplaced.remove(assigned as usize);
            if unplaced.len() == 0{
                break;
            }
        }else{
            break;
        }
    }
}
//BondVars に結合の長さが設定されている場合 Some を返す。
pub fn get_bond_length(p1:usize,p2:usize,bondvec:&Vec<&BondVars>)->Option<(f64,bool)>{
    for bb in bondvec.iter(){
        if bb.atoms.0 == p1 && bb.atoms.1 == p2{
            return Some((bb.b0,false));
        }
        if bb.atoms.1 == p1 && bb.atoms.0 == p2{
            return Some((bb.b0,true));
        }
    }
    return None;
}


//Constraint として登録されている AngleVars Vector から、指定の Angle があるかどうか走査し、あった場合 Some を返す。
//与えられた Angle の順と逆の場合 True が帰る
pub fn get_angle(p0:usize,p1:usize,p2:usize,anglevec:&Vec<&AngleVars>)->Option<(f64,bool)>{
    for bb in anglevec.iter(){
        if bb.atoms.0 == p0 && bb.atoms.1 == p1 && bb.atoms.2 == p2{
            return Some((bb.theta0,false));
        }
        if bb.atoms.0 == p2 && bb.atoms.1 == p1 && bb.atoms.2 == p0{
            return Some((bb.theta0,true));
        }
    }
    return None;
}


pub fn get_dihedral_angle(p0:usize,p1:usize,p2:usize,p3:usize,dihedvec:&Vec<&DihedralVars>)->Option<(f64,bool)>{
    for bb in dihedvec.iter(){
        if bb.atoms.0 == p0 && bb.atoms.1 == p1 && bb.atoms.2 == p2&& bb.atoms.3 == p3{
            return Some((bb.delta,false));
        }
        if bb.atoms.0 == p3 && bb.atoms.1 == p2 && bb.atoms.2 == p1 && bb.atoms.3 == p0{
            return Some((bb.delta,true));
        }
    }
    return None;
}


//残基と N 末 C 末処理のキャップのエントリへの参照を渡すと
// DELETE 文に基づいて諸情報を削除してキャップエントリの情報を結合する。
//capentry 内での重複は見ていないので注意
pub fn create_capped<'a>(target:&'a charmm_param::Param_RESI,capentries:Vec<&'a charmm_param::Param_RESI>)
->charmm_param::Param_RESI_merged<'a>{
    let mut deleten:HashSet<String> = HashSet::new();
    let mut atom_params:Vec<&'a  charmm_param::AtomResiParam>= vec![];
    
    //cap に既にある Atom は挿入しない
    
    let mut added:HashSet<String> = HashSet::new();
    let mut bond:Vec<&'a (charmm_param::PosAtom,charmm_param::PosAtom,usize)>= vec![];
    let mut donor:Vec<&'a (String,String)>= vec![];
    let mut acceptor:Vec<&'a (String,String)>= vec![];
    let mut ic_params:Vec<&'a charmm_param::ICResiParam> = vec![];
    let mut impr:Vec<&'a Vec<charmm_param::PosAtom>> = vec![];
    let mut cmap:Vec<&'a Vec<charmm_param::PosAtom>> = vec![];
    let mut dihe:Vec<&'a Vec<charmm_param::PosAtom>> = vec![];

    for capentry in capentries.into_iter(){
        for aa in capentry.delete.iter(){
            for bb in aa.iter(){
                deleten.insert(bb.to_string());
            }
        }
        for aa in capentry.atom_params.iter(){
            atom_params.push(aa);
            added.insert(aa.atom_name.clone());
        }

        for aa in capentry.bond.iter(){
            bond.push(aa);
        }
        for aa in capentry.donor.iter(){
            donor.push(aa);
        }
        
        for aa in capentry.acceptor.iter(){
            acceptor.push(aa);
        }

        for aa in capentry.ic_params.iter(){
            ic_params.push(aa);
        }
        
        for aa in capentry.impr.iter(){
            impr.push(aa);
        }

        for aa in capentry.cmap.iter(){
            cmap.push(aa);
        }

        for aa in capentry.dihe.iter(){
            dihe.push(aa);
        }
    }


    for aa in target.atom_params.iter(){
        if added.contains(&aa.atom_name){
            continue;
        }
        if !deleten.contains(&aa.atom_name){
            atom_params.push(aa);
        }
    }
    
    for aa in target.bond.iter(){
        if !deleten.contains(&aa.0.atom_name)
        && !deleten.contains(&aa.1.atom_name){
            bond.push(aa);
        }
    }
    for aa in target.cmap.iter(){
        let mut flag:bool = true;
        for i in aa.iter(){
            if deleten.contains(&i.atom_name){
                flag = false;
                break;
            }
        }
        if flag{
            cmap.push(aa);
        }
    }
    for aa in target.donor.iter(){
        if !deleten.contains(&aa.0)
        && !deleten.contains(&aa.1){
            donor.push(aa);
        }
    }
    for aa in target.acceptor.iter(){
        if !deleten.contains(&aa.0)
        && !deleten.contains(&aa.1){
            acceptor.push(aa);
        }
    }

    for aa in target.ic_params.iter(){
        if !deleten.contains(&aa.atoms.0.atom_name)
        && !deleten.contains(&aa.atoms.1.atom_name)
        && !deleten.contains(&aa.atoms.2.atom_name)
        && !deleten.contains(&aa.atoms.3.atom_name){
            ic_params.push(aa);
        }
    }
    for aa in target.impr.iter(){
        let mut flag:bool = true;
        for i in aa.iter(){
            if deleten.contains(&i.atom_name){
                flag = false;
                break;
            }
        }
        if flag{
            impr.push(aa);
        }
    }
    for aa in target.dihe.iter(){
        let mut flag:bool = true;
        for i in aa.iter(){
            if deleten.contains(&i.atom_name){
                flag = false;
                break;
            }
        }
        if flag{
            dihe.push(aa);
        }
    }
    
    let residue_name = target.residue_name.clone();
    let charge = target.charge.clone();
    return charmm_param::Param_RESI_merged{
        residue_name,
        atom_params,
        charge,
        ic_params,
        bond,
        impr,
        cmap,
        dihe,
        donor,
        acceptor
    };
}


pub struct CharmmVars{
    pub bondvec:Vec<BondVars>,
    pub anglevec:Vec<AngleVars>,
    pub ubvec:Vec<UBVars>,
    pub dihedvec:Vec<DihedralVars>,
    pub imprvec:Vec<IMPRVars>,
    pub hb_donors:Vec<(usize,usize)>,
    pub hb_acceptors:Vec<(usize,usize)>,
    pub lj_energy_calculator:LJEnergyCalculator,
    pub electrostatic_energy_calculator:ElectrostaticEnergyCalculator,
}

impl CharmmVars{
    pub fn make_sub_env_vars(&self,mapper:&Vec<i64>)->CharmmVars{
        let mut ret:CharmmVars = CharmmVars{
            bondvec:vec![],
            anglevec:vec![],
            ubvec:vec![],
            dihedvec:vec![],
            imprvec:vec![],
            lj_energy_calculator:self.lj_energy_calculator.clone(),
            electrostatic_energy_calculator:self.electrostatic_energy_calculator.clone(),
            hb_donors:vec![],
            hb_acceptors:vec![],
        };

        for cc in self.bondvec.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            {
                continue;
            }
            let mut pcc:BondVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize
            );
            ret.bondvec.push(pcc);
        }
        
        for cc in self.anglevec.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            {
                continue;
            }
            let mut pcc:AngleVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize
            );
            ret.anglevec.push(pcc);
        }
        

        for cc in self.ubvec.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            {
                continue;
            }
            let mut pcc:UBVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize
            );
            ret.ubvec.push(pcc);
        }
        
        let mut dihedmap:Vec<usize> = vec![];
        for (cii,cc) in self.dihedvec.iter().enumerate(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            || mapper[cc.atoms.3] < 0
            {
                continue;
            }
            let mut pcc:DihedralVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize,
                mapper[cc.atoms.3] as usize
            );
            dihedmap.push(cii);
            ret.dihedvec.push(pcc);
        }
        
        /*ToDo: 呼び出し側で再計算させる
        if self.dihed_weight_charmm.len() > 0{
            for dd in dihedmap.into_iter(){
                ret.dihed_weight_charmm.push(self.dihed_weight_charmm[dd]);
            }
        }
        */

        for cc in self.imprvec.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            || mapper[cc.atoms.3] < 0
            {
                continue;
            }
            let mut pcc:IMPRVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize,
                mapper[cc.atoms.3] as usize
            );
            ret.imprvec.push(pcc);
        }
        
        for cc in self.hb_acceptors.iter(){
            if mapper[cc.0] < 0
            || mapper[cc.1] < 0
            {
                continue;
            }
            ret.hb_acceptors.push((mapper[cc.0] as usize,mapper[cc.1] as usize));
        }
        
        for cc in self.hb_donors.iter(){
            if mapper[cc.0] < 0
            || mapper[cc.1] < 0
            {
                continue;
            }
            ret.hb_donors.push((mapper[cc.0] as usize,mapper[cc.1] as usize));
        }
        
        return ret;
    }
}


//うーん
pub struct CharmmEnv{
    pub atoms:Vec<MDAtom>,
    pub dist:Vec<Vec<f64>>,
    pub num_edges:Vec<Vec<u64>>,//ある原子とある原子の最短距離（存在するエッジの数）。一定数以上についてはカウントしない。今はその一定数は 5。
}
impl CharmmEnv{
    

    pub fn accept_subenv_atom(source:&mut CharmmEnv,subenv:&CharmmEnv,mapper:&Vec<i64>){
        if mapper.len() != source.atoms.len(){
            panic!("Mapper must have the same length with the source atom array. mapper[atom index in source] -> atom index in sub.");
        }
        for (ii,ss_) in mapper.iter().enumerate(){
            if *ss_ < 0{
                continue;
            }
            let ss:usize = *ss_ as usize;
            source.atoms[ii].set_xyz(subenv.atoms[ss].get_x(),subenv.atoms[ss].get_y(),subenv.atoms[ss].get_z());
        } 
    }
    
    //retain だけの原子を取り出した Env を作成する。
    //Tuple の第二要素として元の Atom と同一の長さの配列を返し、使用された Atom のインデクス番目には新しい Env におけるコピーされた原子の Index 番号が入っている。
    pub fn make_sub_env(&self,retain:&Vec<usize>)->(CharmmEnv,Vec<i64>){
        let mut ret:CharmmEnv = CharmmEnv{
            atoms:vec![],
            dist:vec![vec![0.0;retain.len()];retain.len()],
            num_edges:vec![vec![0;retain.len()];retain.len()]
        };
        let mut mapper:Vec<i64> = vec![-1;self.atoms.len()];
        let mut lastval:i64 = -1;
        for rr in retain.iter(){
            if lastval > *rr as i64{
                panic!("The array of atom index must have been sorted! {} {}\n{:?}",lastval,rr,retain);
            }
            lastval = *rr as i64;
            mapper[*rr] = ret.atoms.len() as i64;
            ret.atoms.push(self.atoms[*rr].clone());
        }
        for ii in 0..ret.atoms.len(){
            ret.atoms[ii].atom_index = ii;
        }
        return (ret,mapper);
    }

    pub fn update_distance(&mut self){
        let num_atoms = self.atoms.len();
        for aa in 0..num_atoms-1{
            for bb in (aa+1)..num_atoms{
                let dd = self.atoms[aa].distance(&self.atoms[bb]);
                self.dist[aa][bb] = dd;
                self.dist[bb][aa] = dd;
            }
        }
    }

    pub fn update_edges(&mut self,bondvec:&Vec<BondVars>,maxcon:u64){
        for ii in 0..self.atoms.len(){
            for jj in 0..self.atoms.len(){
                self.num_edges[ii][jj] = maxcon;
            }
        }
        
        MDAtom::mask_connection(bondvec,&mut self.num_edges,maxcon);
    }
    
    pub fn gen_pseudo_edges(&mut self,conn:u64){
        for ii in 0..self.atoms.len(){
            for jj in 0..self.atoms.len(){
                self.num_edges[ii][jj] = conn;
            }
        }
    }

    //原子一個しか動かない場合
    pub fn update_distance_one(&mut self,bb:usize){
        let num_atoms = self.atoms.len();
        for aa in 0..num_atoms{
            let dd = self.atoms[aa].distance(&self.atoms[bb]);
            self.dist[aa][bb] = dd;
            self.dist[bb][aa] = dd;
        }
    }
}



#[derive(Debug,Clone)]
pub struct MDAtom{
    pub atom_type:String,
    pub atom_name:String,
    pub chain_name:String,
    pub residue_name:String,
    pub residue_number:i64,
    pub residue_ins_code:String,
    
    pub residue_index_in_chain:i64,
    pub atom_index:usize,
    
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub mass:f64,
    pub charge:f64,
    pub created:bool,
    pub unplaced:bool,

    pub nb_epsilon:f64,
    pub nb_r1_2:f64,//Rmin1/2

    pub nb_14flag:bool,//無いの多い
    pub nb_14_epsilon:f64,
    pub nb_14_r1_2:f64,
    pub hybrid_orbit:HybridOrbit,

    pub nterminal:bool,
    pub cterminal:bool,
}


impl MDAtom{
    pub fn dummy()->MDAtom{
        return MDAtom{
            atom_type:"C".to_string(),
            atom_name:"C".to_string(),
            chain_name:"A".to_string(),
            residue_name:"UNK".to_string(),
            residue_number:1,
            residue_ins_code:"".to_string(),
            
            residue_index_in_chain:1,
            atom_index:0,
            
            x:1.0,
            y:1.0,
            z:1.0,
            mass:1.0,
            charge:1.0,
            created:true,
            unplaced:false,
        
            nb_epsilon:1.0,
            nb_r1_2:1.0,//Rmin1/2
        
            nb_14flag:false,
            nb_14_epsilon:1.0,
            nb_14_r1_2:1.0,
            hybrid_orbit:HybridOrbit::SP2,
            nterminal:true,
            cterminal:false,
        };
    }
    pub fn to_pdbatom(&self)->(String,(String,i64,String),pdbdata::PDBAtom){
        let ret = pdbdata::PDBAtom{
            parent_entry:None,
            parent_entity:None,
            parent_asym:None,
            parent_comp:None,
            index:-100,
            serial_number:self.atom_index as i64,
            x:self.x,
            y:self.y,
            z:self.z,
            atom_symbol:"".to_string(),
            atom_code:self.atom_name.clone(),
            alt_code:"".to_string(),
            occupancy:1.0,
            temp_factor:0.0,
            charge:None,
            dummy:false,
            het:false,
            alt:false,
            is_ligand:false,
            atom_site_key:-1
        };
        return (self.chain_name.clone(),(self.residue_name.clone(),self.residue_number,self.residue_ins_code.clone()),ret);
    }
    pub fn get_line_representation(&self)->String{
        return format!("{}.{}.{}.{}.{}.{}",
        self.atom_type,
        self.atom_name,
        self.chain_name,
        self.residue_name,
        self.residue_number,
        self.residue_ins_code);
    }

    //残基や原子の名前を CHARMM に合うように変更する
    pub fn change_to_charmmnames(residues:&mut Vec<pdbdata::PDBComp>){
        for rr in residues.iter_mut(){
            if rr.get_name() == "HIS"{
                eprintln!("HIS was changed to HSD.");
                rr.set_name("HSD");
            }
            if rr.get_name() == "ILE"{
                for aa in rr.iter_mut_atoms(){
                    if aa.atom_code == "CD1"{
                        eprintln!("ILE CD1 was changed to CD.");
                        aa.atom_code = "CD".to_string();
                    }
                }
            }
        }
        

        let mut hflag:bool = false;
        for aa in residues[0].iter_atoms(){
            if aa.atom_code == "H"{
                hflag = true;
                eprintln!("The first H {} is removed.",aa.serial_number);
            }
        }
        if hflag{
            //H があるのはおそらく別のエントリかアミノ酸モノマーから作られた Chain
            residues[0].remove_atom_by_name("H");
        }

        let rnum = residues.len();
        let mut o1flag:i64 = -1;
        let mut o2flag:i64 = -1;
        let mut oxtflag:i64 = -1;
        for (ii,aa) in residues[rnum-1].iter_atoms().enumerate(){
            if aa.atom_code == "O" {
                o1flag = ii as i64;
            }
            if aa.atom_code == "OXT"|| aa.atom_code == "OT" {
                o2flag = ii as i64;
                if aa.atom_code == "OXT"{
                    oxtflag = ii as i64;
                }
            }
        }
        if o1flag > -1 && o2flag > -1{
            eprintln!("The last O {} is changed to OT1.",residues[rnum-1].get_atom_at(o1flag as usize).serial_number);
            eprintln!("The last OXT {} is changed to OT2.",residues[rnum-1].get_atom_at(o2flag as usize).serial_number);
            residues[rnum-1].get_mut_atom_at(o1flag as usize).atom_code = "OT1".to_string();
            residues[rnum-1].get_mut_atom_at(o2flag as usize).atom_code = "OT2".to_string();
        }else if o1flag > -1{
            eprintln!("The last O {} is removed",residues[rnum-1].get_atom_at(o1flag as usize).serial_number);
            residues[rnum-1].remove_atom_by_name("O");
            //residues[rnum-1].get_mut_atom_at(o1flag as usize).atom_code = "OT1".to_string();
            
        }else if o2flag > -1{
            eprintln!("The last O(X)T {} is removed",residues[rnum-1].get_atom_at(o2flag as usize).serial_number);
            if oxtflag > -1{
                residues[rnum-1].remove_atom_by_name("OXT");
            }
            residues[rnum-1].remove_atom_by_name("OT");
            //residues[rnum-1].get_mut_atom_at(o1flag as usize).atom_code = "OT1".to_string();
            
        }
    }

    

    pub fn chain_to_atoms(chains:&Vec<&pdbdata::PDBAsym>
    ,parr:&charmm_param::CHARMMParam
    ,add_nc_cap:bool)->(CharmmEnv,CharmmVars){



        //MD に使うのはおそらく MDAtom.atom_type
        let mut mdatoms_all:Vec<MDAtom> = vec![];
        let mut bondvec:Vec<BondVars> = vec![];
        let mut anglevec:Vec<AngleVars> = vec![];
        let mut ubvec:Vec<UBVars> = vec![];
        let mut dihedvec:Vec<DihedralVars> = vec![];
        let mut imprvec:Vec<IMPRVars> = vec![];
        
        let mut masmap:HashMap<String,f64> = HashMap::new();
        let mut hb_donors:Vec<(usize,usize)> = vec![];
        let mut hb_acceptors:Vec<(usize,usize)> = vec![];


        //残基名→RESI エントリへのマップ
        let mut resi_map:HashMap<String,&charmm_param::Param_RESI> = HashMap::new();
        for rss in parr.resi.iter(){
            resi_map.insert(rss.residue_name.clone(),rss);
            //println!("{}",rss.residue_name);
        }

        for chain in chains.iter(){
            //残基順に RESI_MERGED を作成する
            let mut resi_merged_vec:Vec<charmm_param::Param_RESI_merged> = vec![];
            let num_residues:usize = chain.num_comps();
            for (rii,rr) in chain.iter_comps().enumerate(){
                let resi_:Option<&&charmm_param::Param_RESI> = resi_map.get(rr.get_comp_id());
                if let None = resi_{
                    panic!("{} is not found in param.",rr.get_comp_id());
                }
                let mut capentries:Vec<&charmm_param::Param_RESI> = vec![];
                if rii == 0 && add_nc_cap{
                    if rr.get_comp_id() == "GLY"{
                        capentries.push(resi_map.get("GLYP").unwrap_or_else(||panic!("GLYP is not in param!")));
                    }else if rr.get_comp_id() == "PRO"{
                        capentries.push(resi_map.get("PROP").unwrap_or_else(||panic!("PROP is not in param!")));
                    }else{
                        capentries.push(resi_map.get("NTER").unwrap_or_else(||panic!("NTER is not in param!")));
                    }
                }
                if rii == num_residues-1 && add_nc_cap{
                    capentries.push(resi_map.get("CTER").unwrap_or_else(||panic!("CTER is not in param!")));
                }
                if capentries.len() > 0{   
                    resi_merged_vec.push(create_capped(resi_.unwrap(),capentries));
                }else{
                    if !add_nc_cap || (rii != num_residues-1 && rii != 0){
                        resi_merged_vec.push(charmm_param::Param_RESI_merged::new(resi_.unwrap()));        
                    }
                }
            }

            let mut atomtype_ljparam_map:HashMap<String,&charmm_param::Param_NONBONDED_ATOM> = HashMap::new();
            for nn in parr.nonbonded.atom_param.iter(){
                atomtype_ljparam_map.insert(nn.atom_type.clone(),nn);
            }

            for rr in parr.mass.iter(){
                masmap.insert(rr.atom_name.clone(),rr.mass);
            }
            //Residue 順の、atom_name->MDAtom のハッシュ
            let mut residue_atomname_to_indexinmdatoms_all:Vec<HashMap<String,usize>> = vec![];
            let resnum = chain.num_comps();
            for (rii,rr) in chain.iter_comps().enumerate(){

                let mut ratoms:HashMap<String,MDAtom> = HashMap::new();
                let resi:&charmm_param::Param_RESI_merged = &resi_merged_vec[rii];

                //atom_name->AtomResiParam のマップ
                let mut amap:HashMap<String,&charmm_param::AtomResiParam> = HashMap::new();
                for aa in resi.atom_params.iter(){
                    amap.insert(aa.atom_name.clone(),aa);
                }


                //与えられた Chain に含まれている原子に基本的な属性をアサインする。
                for (_,aa) in rr.iter_atoms().enumerate(){

                    let atype = amap.get(&aa.atom_code).unwrap_or_else(||panic!("atom {} in residue {} {} not found in map",aa.atom_code,rr.get_comp_id(),rr.get_seq_id())).atom_type.clone();
                    let ljparam:&charmm_param::Param_NONBONDED_ATOM = atomtype_ljparam_map.get(&atype).unwrap_or_else(||panic!("nonbonded param for {} is not defined.",atype));
                    if ljparam.tanford_kirkwood{
                        panic!("tanford_kirkwood is not supported.");
                    }
                    let mut att = MDAtom{
                        atom_type:atype.clone(),
                        atom_name:aa.atom_code.clone(),
                        chain_name:chain.chain_name.clone(),
                        residue_name:rr.get_comp_id().to_string(),
                        residue_ins_code:rr.get_ins_code().to_string(),
                        atom_index:9999999,
                        residue_index_in_chain:rii as i64,
                        residue_number:rr.get_seq_id(),
                        x:aa.get_x(),
                        y:aa.get_y(),
                        z:aa.get_z(),
                        mass:0.0,
                        charge:0.0,
                        created:false,
                        unplaced:false,
                        nb_epsilon:ljparam.eps,
                        nb_r1_2:ljparam.rmin_05,
                    
                        nb_14flag:ljparam.flag_1_4,
                        nb_14_epsilon:ljparam.eps_1_4,
                        nb_14_r1_2:ljparam.rmin_05_1_4,
                        hybrid_orbit:HybridOrbit::Other,
                        nterminal:rii == 0,
                        cterminal:rii == resnum-1,
                    };
                    let atype = amap.get(&aa.atom_code).unwrap_or_else(||panic!("atom not found in map {}",&aa.atom_code));
                    att.charge = atype.partial_charge;
                    att.mass = *masmap.get(&atype.atom_type).unwrap();
                    ratoms.insert(att.atom_name.clone(),att);
                }

                //本来存在しているはずであるが与えられた Chain に Missing な原子を作成する
                
                
                let mut dummypos = (0.0,0.0,0.0);
                if ratoms.len() > 0{
                    for ppp in ratoms.iter(){
                        dummypos.0 = ppp.1.x;
                        dummypos.1 = ppp.1.y;
                        dummypos.2 = ppp.1.z;
                        break;
                    }
                }
                for (kk,vv) in amap.iter(){
                    if !ratoms.contains_key(kk){
                        let ljparam:&charmm_param::Param_NONBONDED_ATOM = atomtype_ljparam_map.get(&vv.atom_type).unwrap_or_else(||panic!("nonbonded param for {} is not defined.",&vv.atom_type));
                        if ljparam.tanford_kirkwood{
                            panic!("tanford_kirkwood is not supported.");
                        }
                        let mut att = MDAtom{
                            atom_type:vv.atom_type.clone(),
                            atom_name:kk.clone(),
                            chain_name:chain.chain_name.clone(),
                            residue_name:rr.get_comp_id().to_string(),
                            residue_ins_code:rr.get_ins_code().to_string(),
                            atom_index:9999999,
                            residue_index_in_chain:rii as i64,
                            residue_number:rr.get_seq_id(),
                            x:dummypos.0,
                            y:dummypos.1,
                            z:dummypos.2,
                            mass:0.0,
                            charge:0.0,
                            created:true,
                            unplaced:true,
                            nb_epsilon:ljparam.eps,
                            nb_r1_2:ljparam.rmin_05,
                        
                            nb_14flag:ljparam.flag_1_4,
                            nb_14_epsilon:ljparam.eps_1_4,
                            nb_14_r1_2:ljparam.rmin_05_1_4,
                            hybrid_orbit:HybridOrbit::Other,
                            nterminal:rii == 0,
                            cterminal:rii == resnum-1,
                        };

                        att.charge = vv.partial_charge;
                        att.mass = *masmap.get(&vv.atom_type).expect("atom not found in masmap");
                        
                        //println!("{} {} {:?}",rr.get_comp_id(),kk,att);
                        ratoms.insert(kk.clone(),att);
                    }
                }
                //println!("{:?}",ratoms);

                //作成した MDAtom エントリを、一つのベクトルにまとめ、
                //Residue の情報は residue_atomname_to_indexinmdatoms_all として保持する
                //Residue-> Atom へのアクセスは residue_atomname_to_indexinmdatoms_all からインデックスを抽出する
                let mut rvatoms:HashMap<String,usize> = HashMap::new();
                let mut rratoms:Vec<(String,MDAtom)> = ratoms.into_iter().collect();
                rratoms.sort_by(|a,b|a.0.cmp(&b.0));//ToDo もっと綺麗に並べる
                rratoms.reverse();
                for (kk,vv) in rratoms.into_iter(){
                    rvatoms.insert(kk.clone(),mdatoms_all.len());
                    mdatoms_all.push(vv);
                }
                for (aii,mm) in mdatoms_all.iter_mut().enumerate(){
                    mm.atom_index = aii;
                }
                for aa in resi.donor.iter(){
                    hb_donors.push((*rvatoms.get(&aa.0).unwrap_or_else(||panic!("{} is not defined.",&aa.0))
                    ,*rvatoms.get(&aa.1).unwrap_or_else(||panic!("{} is not defined.",&aa.1))));
                }
                
                for aa in resi.acceptor.iter(){
                    hb_acceptors.push((*rvatoms.get(&aa.0).unwrap_or_else(||panic!("{} is not defined.",&aa.0))
                    ,*rvatoms.get(&aa.1).unwrap_or_else(||panic!("{} is not defined.",&aa.1))));
                }
                residue_atomname_to_indexinmdatoms_all.push(rvatoms);
                
            }
            
            //atom_type+"\t"+atom_type->Param_BOND エントリのマップを作成する
            let mut bmap:HashMap<String,&charmm_param::Param_BOND> = HashMap::new();
            for bb in parr.bonds.iter(){
                bmap.insert(bb.atoms.0.clone()+"\t"+bb.atoms.1.as_str(),bb);
                bmap.insert(bb.atoms.1.clone()+"\t"+bb.atoms.0.as_str(),bb);
            }

            //angle,dihed,improper については、文字列化してハッシュでアクセスさせる。
            //wildcard がある場合は Regex で認識させる。
            //angles に wildcard はないかも
            let mut code_angle:HashMap<String,Vec<&charmm_param::Param_ANGLE>> = HashMap::new();
            //ハッシュじゃないので、複数登録ある場合は複数回ヒットする
            let mut code_angle_regex:Vec<(Regex,&charmm_param::Param_ANGLE)> = vec![];
            let mut code_dihed:HashMap<String,Vec<&charmm_param::Param_DIHEDRAL>> = HashMap::new();
            let mut code_dihed_regex:Vec<(Regex,&charmm_param::Param_DIHEDRAL)> = vec![];
            let mut code_impr:HashMap<String,Vec<&charmm_param::Param_IMPROPER>> = HashMap::new();
            let mut code_impr_regex:Vec<(Regex,&charmm_param::Param_IMPROPER)> = vec![];
            
            for pp in parr.angles.iter(){
                let mut att:Vec<String> = vec![pp.atoms.0.clone(),pp.atoms.1.clone(),pp.atoms.2.clone()];
                let mut rflag = false;
                for tii in 0..3{
                    if att[tii] == "X"{
                        rflag = true;
                        att[tii] = ".+".to_string();
                    }
                }
                let codestring:String = format!("#{}\t{}\t{}#",att[0],att[1],att[2]);
                if rflag{
                    code_angle_regex.push((Regex::new(codestring.as_str()).unwrap(),pp));
                }else{
                    if !code_angle.contains_key(&codestring){
                        code_angle.insert(codestring.clone(),vec![]);
                    }
                    code_angle.get_mut(&codestring).unwrap().push(pp);
                }
            }
            
            for pp in parr.dihedrals.iter(){
                let mut att:Vec<String> = vec![pp.atoms.0.clone(),pp.atoms.1.clone(),pp.atoms.2.clone(),pp.atoms.3.clone()];
                let mut rflag = false;
                for tii in 0..4{
                    if att[tii] == "X"{
                        rflag = true;
                        att[tii] = ".*".to_string();
                    }
                }
                let codestring:String = format!("#{}\t{}\t{}\t{}#",att[0],att[1],att[2],att[3]);
                if rflag{
                    code_dihed_regex.push((Regex::new(codestring.as_str()).unwrap(),pp));
                }else{
                    if !code_dihed.contains_key(&codestring){
                        code_dihed.insert(codestring.clone(),vec![]);
                    }
                    code_dihed.get_mut(&codestring).unwrap().push(pp);
                }
            }

            for pp in parr.impropers.iter(){
                let mut att:Vec<String> = vec![pp.atoms.0.clone(),pp.atoms.1.clone(),pp.atoms.2.clone(),pp.atoms.3.clone()];
                let mut rflag = false;
                for tii in 0..4{
                    if att[tii] == "X"{
                        rflag = true;
                        att[tii] = ".*".to_string();
                    }
                }
                let codestring:String = format!("#{}\t{}\t{}\t{}#",att[0],att[1],att[2],att[3]);
                if rflag{
                    code_impr_regex.push((Regex::new(codestring.as_str()).unwrap(),pp));
                }else{
                    if !code_impr.contains_key(&codestring){
                        code_impr.insert(codestring.clone(),vec![]);
                    }
                    code_impr.get_mut(&codestring).unwrap().push(pp);
                }
            }
            let atom_relativepos = |excode:&str|{
                if excode == "-"{
                    return -1;
                }else if excode == "+"{
                    return 1;
                }else if excode == ""{
                    return 0;
                }else{
                    panic!("{} was not defined",excode);
                }
            };
            


            let resnum = residue_atomname_to_indexinmdatoms_all.len();
            
            //unplaced な Atom の場所を見積もる
            
            //bond length が 0.0 の場合 Bond からもってくる
            let mut bond_hm:HashMap<String,f64> = HashMap::new();
            for bb in parr.bonds.iter(){
                bond_hm.insert("".to_string()+&bb.atoms.0+"_"+&bb.atoms.1,bb.b0);
                bond_hm.insert("".to_string()+&bb.atoms.1+"_"+&bb.atoms.0,bb.b0);
            }
            
            let mut angle_hm:HashMap<String,f64> = HashMap::new();
            for aa in parr.angles.iter(){
                
                angle_hm.insert("".to_string()+&aa.atoms.0+"_"+&aa.atoms.1+"_"+&aa.atoms.2,aa.theta0);
                angle_hm.insert("".to_string()+&aa.atoms.2+"_"+&aa.atoms.1+"_"+&aa.atoms.0,aa.theta0);
            }


            for ii in 0..resnum{
                
                let _rname:&str = &chain.get_comp_at(ii).get_comp_id();
                let _hh = &residue_atomname_to_indexinmdatoms_all[ii];
                
                let resi:&charmm_param::Param_RESI_merged = &resi_merged_vec[ii];
                
                for icc in resi.ic_params.iter(){

                    let res0:i64 = atom_relativepos(&icc.atoms.0.ex_code)+ii as i64;
                    let res1:i64 = atom_relativepos(&icc.atoms.1.ex_code)+ii as i64;
                    let res2:i64 = atom_relativepos(&icc.atoms.2.ex_code)+ii as i64;
                    let res3:i64 = atom_relativepos(&icc.atoms.3.ex_code)+ii as i64;
                    if res0 < 0 
                    || res1 < 0 
                    || res2 < 0 
                    || res3 < 0 
                    || res0 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 
                    || res1 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 
                    || res2 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 
                    || res3 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 {
                        continue;
                    }
                    //前後の残基を含むのでこういうことをする必要がある
                    let atom0:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&icc.atoms.0.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res0 as usize).get_comp_id()
                    ,residue_atomname_to_indexinmdatoms_all[res0 as usize],icc.atoms.0).as_str())];
                    
                    let atom1:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res1 as usize].get(&icc.atoms.1.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res1 as usize).get_comp_id()
                    ,residue_atomname_to_indexinmdatoms_all[res1 as usize],icc.atoms.1).as_str())];
                    
                    let atom2:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res2 as usize].get(&icc.atoms.2.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res2 as usize).get_comp_id()
                    ,residue_atomname_to_indexinmdatoms_all[res2 as usize],icc.atoms.2).as_str())];
                    
                    let atom3:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res3 as usize).get_comp_id()
                    ,residue_atomname_to_indexinmdatoms_all[res3 as usize],icc.atoms.3).as_str())];
                    
                    let mut l01_or_02:f64 = icc.length01_or_02;
                    let mut a012_or_021:f64 = icc.angle012_or_021;
                    let mut l23:f64 = icc.length23;
                    let mut a123:f64 = icc.angle123;
                    
                    
                    if l01_or_02 == 0.0{//なんか Dihed しか入ってない奴があるのでそういうのは無視する
                        if icc.improper{
                            let acode = atom0.atom_type.clone()+"_"+&atom2.atom_type;
                            if bond_hm.contains_key(&acode){
                                l01_or_02 = *bond_hm.get(&acode).unwrap();
                            }else{
                                println!("{} was not found.",acode);
                            }
                            let acode = atom0.atom_type.clone()+"_"+&atom2.atom_type+"_"+&atom1.atom_type;
                            if angle_hm.contains_key(&acode){
                                a012_or_021 = *angle_hm.get(&acode).unwrap();
                            }else{
                                println!("{} was not found.",acode);
                            }

                        }else{
                            let acode = atom0.atom_type.clone()+"_"+&atom1.atom_type;
                            if bond_hm.contains_key(&acode){
                                l01_or_02 = *bond_hm.get(&acode).unwrap();
                            }else{
                                println!("{} was not found.",acode);
                            }
                            let acode = atom0.atom_type.clone()+"_"+&atom1.atom_type+"_"+&atom2.atom_type;
                            if angle_hm.contains_key(&acode){
                                a012_or_021 = *angle_hm.get(&acode).unwrap();
                            }else{
                                println!("{} was not found.",acode);
                            }
                        }
                    }
                    
                    if l23 == 0.0{//なんか Dihed しか入ってない奴があるのでそういうのは無視する
                        let acode = atom2.atom_type.clone()+"_"+&atom3.atom_type;
                        if bond_hm.contains_key(&acode){
                            l23 = *bond_hm.get(&acode).unwrap();
                        }else{
                            println!("{} was not found.",acode);
                        }

                        let acode = atom1.atom_type.clone()+"_"+&atom2.atom_type+"_"+&atom3.atom_type;
                        if angle_hm.contains_key(&acode){
                            a123 = *angle_hm.get(&acode).unwrap();
                        }else{
                            println!("{} was not found.",acode);
                        }
                    }
                    if l23 == 0.0 || l01_or_02 == 0.0{
                        continue;
                    }

                    /*
                    //デバッグ用コード                            
                    let mcode = *residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).unwrap();
                    if atom3.unplaced || ( mdatoms_all[mcode].residue_name == "PRO" && mdatoms_all[mcode].atom_name == "C"){
                    */
                    if atom3.unplaced{
                        if !atom0.unplaced 
                        && !atom1.unplaced  
                        && !atom2.unplaced {
                            
                            let poss = if icc.improper{
                                //println!("!!{} {} {}",mdenv.atoms[*uu].residue_name,mdenv.atoms[*uu].atom_name,"**");
                                estimate_3_impr(
                                &atom0.get_xyz()
                                , &atom1.get_xyz()
                                , &atom2.get_xyz()
                                , l01_or_02
                                , a012_or_021
                                , icc.dehedral0123
                                , a123
                                ,l23)
                            }else{
                                estimate_3_dihed(
                                    &atom0.get_xyz()
                                , &atom1.get_xyz()
                                , &atom2.get_xyz()
                                , l01_or_02
                                , a012_or_021
                                , icc.dehedral0123
                                , a123
                                ,l23)
                            };
                            if icc.improper{
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).unwrap()].set_xyz(poss.0,poss.1,poss.2);
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).unwrap()].unplaced = false;
                            }else{
                                /*
                                //デバッグ用コード
                                if mdatoms_all[mcode].residue_name == "PRO" && mdatoms_all[mcode].atom_name == "C"{
                                    println!("PPPP{}",process_3d::distance(&(mdatoms_all[mcode].x,mdatoms_all[mcode].y,mdatoms_all[mcode].z),&poss));
                                }
                                */
                                
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).unwrap()].set_xyz(poss.0,poss.1,poss.2);
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res3 as usize].get(&icc.atoms.3.atom_name).unwrap()].unplaced = false;
                            }
                            
                        }
                    }else if atom0.unplaced{
                        if !atom3.unplaced 
                        && !atom2.unplaced  
                        && !atom1.unplaced {
                            let poss = if icc.improper{
                                estimate_3_impr(
                                &atom3.get_xyz()
                                , &atom1.get_xyz()
                                , &atom2.get_xyz()
                                , l23
                                , a123
                                , icc.dehedral0123
                                , a012_or_021
                                , l01_or_02)
                            }else{
                                estimate_3_dihed(
                                &atom3.get_xyz()
                                , &atom2.get_xyz()
                                , &atom1.get_xyz()
                                , l23
                                , a123
                                , icc.dehedral0123
                                , a012_or_021
                                , l01_or_02)

                            };
                            if icc.improper{
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&icc.atoms.0.atom_name).unwrap()].set_xyz(poss.0,poss.1,poss.2);
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&icc.atoms.0.atom_name).unwrap()].unplaced = false;
                            }else{
                                
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&icc.atoms.0.atom_name).unwrap()].set_xyz(poss.0,poss.1,poss.2);
                                mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&icc.atoms.0.atom_name).unwrap()].unplaced = false;
                            }
                        }
                    }
                    
                }
                
            }


            let mut max_bondorder:Vec<usize> =vec![1;mdatoms_all.len()];
            //bond の作成
            for ii in 0..resnum{
                let _rname:&str = &chain.get_comp_at(ii).get_comp_id();
                let _hh = &residue_atomname_to_indexinmdatoms_all[ii];

                let resi:&charmm_param::Param_RESI_merged = &resi_merged_vec[ii];
                let mut bbonds:Vec<BondVars> = vec![];
                

                for bb in resi.bond.iter(){
                    let res0:i64 = atom_relativepos(&bb.0.ex_code)+ii as i64;
                    let res1:i64 = atom_relativepos(&bb.1.ex_code)+ii as i64;
                    if res0 < 0 || res1 < 0 || res0 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 || res1 > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 {
                        continue;
                    }
                    //前後の残基を含むのでこういうことをする必要がある
                    let atom0:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res0 as usize].get(&bb.0.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res0 as usize).get_comp_id(),residue_atomname_to_indexinmdatoms_all[res0 as usize],bb.0).as_str())];
                    let atom1:&MDAtom = &mdatoms_all[*residue_atomname_to_indexinmdatoms_all[res1 as usize].get(&bb.1.atom_name).expect(format!("Can not find atom. {} {:?} {:?} ",chain.get_comp_at(res1 as usize).get_comp_id(),residue_atomname_to_indexinmdatoms_all[res1 as usize],bb.1).as_str())];
                    max_bondorder[atom0.atom_index] = max_bondorder[atom0.atom_index].max(bb.2);
                    max_bondorder[atom1.atom_index] = max_bondorder[atom1.atom_index].max(bb.2);
                    let mut okflag:bool = false;
                    for bz in parr.bonds.iter(){
                        if (bz.atoms.0 == atom0.atom_type && bz.atoms.1 == atom1.atom_type)
                        ||(bz.atoms.1 == atom0.atom_type && bz.atoms.0 == atom1.atom_type)
                        {
                            
                            okflag = true;
                            bbonds.push(
                                BondVars{
                                    atoms:(atom0.atom_index,atom1.atom_index),
                                    b0:bz.b0,
                                    kb:bz.kb
                                }
                            );
                            break;
                        }
                    }

                    if !okflag{
                        eprintln!("The bond was not found!");
                        eprintln!("===\n{:?}",bb);
                        eprintln!("===\n{:?}\n{:?}\n====",atom0,atom1);
                        panic!();
                    }
                    //println!("{}",okflag);
                }

                for bb in resi.impr.iter(){
                    let rz:Vec<i64> = vec![
                    atom_relativepos(&bb[0].ex_code)+ii as i64
                    ,atom_relativepos(&bb[1].ex_code)+ii as i64
                    ,atom_relativepos(&bb[2].ex_code)+ii as i64
                    ,atom_relativepos(&bb[3].ex_code)+ii as i64
                    ];
                    let mut mdavec:Vec<&MDAtom> = vec![];
                    let mut obflag:bool = false;
                    for (zii,zz) in rz.into_iter().enumerate(){
                        if zz < 0 || zz > residue_atomname_to_indexinmdatoms_all.len() as i64 -1 {
                            obflag = true;
                            break;
                        }
                        mdavec.push(&mdatoms_all[*residue_atomname_to_indexinmdatoms_all[zz as usize].get(
                            &bb[zii].atom_name).expect(
                            format!("Atom mapping error {:?} {:?} "
                            ,residue_atomname_to_indexinmdatoms_all[zz as usize],bb[zii]).as_str())]);
                    }
                    if obflag{//終端残基
                        continue;
                    }
                    
                    //前後の残基を含むのでこういうことをする必要がある
                    let mut hitflag:bool = false;
                    
                    let iicode:String = format!("#{}\t{}\t{}\t{}#"
                    ,mdavec[0].atom_type
                    ,mdavec[1].atom_type
                    ,mdavec[2].atom_type
                    ,mdavec[3].atom_type);
                    
                    let iicode_rev:String = format!("#{}\t{}\t{}\t{}#"
                    ,mdavec[3].atom_type
                    ,mdavec[2].atom_type
                    ,mdavec[1].atom_type
                    ,mdavec[0].atom_type);

                    if code_impr.contains_key(&iicode){
                        imprvec.push(
                            IMPRVars{
                                atoms:(mdavec[0].atom_index
                                    ,mdavec[1].atom_index
                                    ,mdavec[2].atom_index
                                    ,mdavec[3].atom_index),
                                kpsi:code_impr.get(&iicode).unwrap()[0].kpsi,
                                psi0:code_impr.get(&iicode).unwrap()[0].psi0
                            }
                        );
                        hitflag = true;
                    }else if code_impr.contains_key(&iicode_rev){
                        imprvec.push(
                            IMPRVars{
                                atoms:(mdavec[3].atom_index
                                    ,mdavec[2].atom_index
                                    ,mdavec[1].atom_index
                                    ,mdavec[0].atom_index),
                                kpsi:code_impr.get(&iicode_rev).unwrap()[0].kpsi,
                                psi0:code_impr.get(&iicode_rev).unwrap()[0].psi0
                            }
                        );
                        hitflag = true;
                    }
                    if !hitflag{
                        for vv in code_impr_regex.iter(){
                            if let Some(_x) = vv.0.captures(&iicode){
                                imprvec.push(IMPRVars{
                                    atoms:(mdavec[0].atom_index
                                        ,mdavec[1].atom_index
                                        ,mdavec[2].atom_index
                                        ,mdavec[3].atom_index),
                                    kpsi:vv.1.kpsi,
                                    psi0:vv.1.psi0
                                });
                                
                                hitflag = true;
                                break;
                            }
                            if let Some(_x) = vv.0.captures(&iicode_rev){
                                imprvec.push(IMPRVars{
                                    atoms:(mdavec[3].atom_index
                                        ,mdavec[2].atom_index
                                        ,mdavec[1].atom_index
                                        ,mdavec[0].atom_index),
                                    kpsi:vv.1.kpsi,
                                    psi0:vv.1.psi0
                                });
                                hitflag = true;
                                break;
                            }
                        }
                    }

                    if !hitflag{
                        eprintln!("The impr was not found!");
                        eprintln!("===\n{}\n{:?}",iicode,resi);
                        panic!();
                    }
                    //println!("{}",okflag);
                }

                bondvec.append(&mut bbonds);
            }

            let all3 = autogenerate_connection(&bondvec,3);
            for acc in all3.iter(){
                //原子番号が rev の場合 0 は true
                let mut zvec:Vec<(bool,&charmm_param::Param_ANGLE)> = vec![];
                let a0:&MDAtom = &mdatoms_all[acc[0]];
                let a1:&MDAtom = &mdatoms_all[acc[1]];
                let a2:&MDAtom = &mdatoms_all[acc[2]];
                let ccode = format!("#{}\t{}\t{}#",a0.atom_type,a1.atom_type,a2.atom_type);
                let ccode_rev = format!("#{}\t{}\t{}#",a2.atom_type,a1.atom_type,a0.atom_type);
                let mut hitflag:bool = false;

                if code_angle.contains_key(&ccode){
                    for vv in code_angle.get(&ccode).unwrap().iter(){
                        zvec.push((false,vv));
                    }
                    hitflag = true;
                }
                
                if code_angle.contains_key(&ccode_rev){
                    for vv in code_angle.get(&ccode_rev).unwrap().iter(){
                        zvec.push((true,vv));
                    }
                    hitflag = true;
                }

                if !hitflag{
                    for vv in code_angle_regex.iter(){
                        if let Some(_x) = vv.0.captures(&ccode){
                            zvec.push((false,vv.1.clone()));
                        }
                        if let Some(_x) = vv.0.captures(&ccode_rev){
                            zvec.push((true,vv.1.clone()));
                        }
                    }
                }

                //for (bbb,vvv) in zvec.iter(){
                if zvec.len() > 0{
                    let (bbb,vvv):(bool,&charmm_param::Param_ANGLE) = zvec[0];
                    let mut atom3 = (a0.atom_index,a1.atom_index,a2.atom_index);
                    if bbb{
                        atom3 = (a2.atom_index,a1.atom_index,a0.atom_index);
                    }

                    if vvv.kubflag{
                        ubvec.push(
                            UBVars{
                                atoms:atom3,
                                kub:vvv.kub,
                                s0:vvv.s0
                            }
                        );
                    }
                    anglevec.push(
                        AngleVars{
                            atoms:atom3,
                            ktheta:vvv.ktheta,
                            theta0:vvv.theta0
                        }
                    );    
                }else{
                    eprintln!("Can not find the angle in param: #{}\t{}\t{}#",a0.atom_type,a1.atom_type,a2.atom_type);
                }
            }

            let all4 = autogenerate_connection(&bondvec,4);
            //ToDo
            //charmm19 の場合 DIHE にある
            for acc in all4.iter(){
                //原子番号が rev の場合 0 は true
                let mut zvec:Vec<(bool,&charmm_param::Param_DIHEDRAL)> = vec![];
                let a0:&MDAtom = &mdatoms_all[acc[0]];
                let a1:&MDAtom = &mdatoms_all[acc[1]];
                let a2:&MDAtom = &mdatoms_all[acc[2]];
                let a3:&MDAtom = &mdatoms_all[acc[3]];
                let ccode = format!("#{}\t{}\t{}\t{}#",a0.atom_type,a1.atom_type,a2.atom_type,a3.atom_type);
                let ccode_rev = format!("#{}\t{}\t{}\t{}#",a3.atom_type,a2.atom_type,a1.atom_type,a0.atom_type);
                let mut hitflag:bool = false;
                let mut debug_string:String = "".to_string();
                if code_dihed.contains_key(&ccode){
                    for vv in code_dihed.get(&ccode).unwrap().iter(){
                        zvec.push((false,vv));
                    }
                    hitflag = true;
                    debug_string = ccode.clone();
                }

                if code_dihed.contains_key(&ccode_rev){
                    for vv in code_dihed.get(&ccode_rev).unwrap().iter(){
                        zvec.push((true,vv));
                    }
                    hitflag = true;
                    debug_string = ccode_rev.clone();
                }

                if !hitflag{
                    for vv in code_dihed_regex.iter(){
                        if let Some(_x) = vv.0.captures(&ccode){
                            zvec.push((false,vv.1.clone()));
                            debug_string = vv.0.to_string();
                        }
                        if let Some(_x) = vv.0.captures(&ccode_rev){
                            zvec.push((true,vv.1.clone()));
                            debug_string = vv.0.to_string();
                        }
                    }
                }

                //for (bbb,vvv) in zvec.iter(){
                    
                if zvec.len() > 0{
                    let (bbb,vvv):(bool,&charmm_param::Param_DIHEDRAL) = zvec[0];
                    let mut atom4 = (a0.atom_index,a1.atom_index,a2.atom_index,a3.atom_index);
                    if bbb{
                        atom4 = (a3.atom_index,a2.atom_index,a1.atom_index,a0.atom_index);
                    }

                    dihedvec.push(
                        DihedralVars{
                            atoms:atom4,
                            kchi:vvv.kchi,
                            n:vvv.n,
                            delta:vvv.delta,
                            debug_string:debug_string.clone()
                        }
                    );
                    
                }else{
                    
                    /*let chkk = format!("#{}:{}:{}:{}#",a0.atom_type,a1.atom_type,a2.atom_type,a3.atom_type);
                    if chkk == "#CR1E:CR1E:CR1E:CR1E#"{
                    panic!("{}",chkk);  
                    }*/

                    eprintln!("Can not find the dihedral angle in param: #{}\t{}\t{}\t{}#",a0.atom_type,a1.atom_type,a2.atom_type,a3.atom_type);
                }
            }
            for ii in 0..max_bondorder.len(){
                if max_bondorder[ii] == 1{
                    mdatoms_all[ii].hybrid_orbit = HybridOrbit::SP3;
                }else if max_bondorder[ii] == 2{
                    mdatoms_all[ii].hybrid_orbit = HybridOrbit::SP2;
                }
            }
        }

        let mlen:usize = mdatoms_all.len();
        let mxval:u64 = mdatoms_all.len() as u64 +1000;
        let mut num_edges:Vec<Vec<u64>> = vec![vec![mxval;mlen];mlen];
        
        MDAtom::mask_connection(&bondvec,&mut num_edges,5);//間にある結合が少ない近隣原子についてはウエイトをかけることがあるので 5 くらいまで距離を計算する
        let dist = vec![vec![0.0;mdatoms_all.len()];mdatoms_all.len()];
        let electrostatic_energy_calculator:ElectrostaticEnergyCalculator = ElectrostaticEnergyCalculator{
            nboption:parr.nonbonded.option_param.clone()
        };
        let lj_energy_calculator:LJEnergyCalculator = LJEnergyCalculator{
            nboption:parr.nonbonded.option_param.clone()
        };

        return (CharmmEnv{ 
            atoms:mdatoms_all,
            dist:dist,
            num_edges:num_edges
        },CharmmVars{
            bondvec,
            anglevec,
            ubvec,
            hb_donors,
            hb_acceptors,
            dihedvec,
            imprvec,
            electrostatic_energy_calculator,
            lj_energy_calculator
        });
    }
    //bond をたどって maxcon 以下の結合数がある場合はその値を入れる。
    //maxcon 以上の結合数がある場合は maxcon を入れる。
    pub fn mask_connection(bonds:&Vec<BondVars>
        ,distmap:&mut Vec<Vec<u64>>,maxcon:u64){
        let num_atoms = distmap.len();
        for ii in 0..num_atoms{
            for jj in 0..num_atoms{
                if ii == jj {
                    distmap[ii][jj] = 0;
                }else{
                    distmap[ii][jj] = maxcon;
                }
            }
        }
        let mut edges:HashMap<usize,Vec<usize>> = HashMap::new();
        for bb in bonds.iter(){
            if !edges.contains_key(&bb.atoms.0){
                edges.insert(bb.atoms.0.clone(),vec![]);
            }
            if !edges.contains_key(&bb.atoms.1){
                edges.insert(bb.atoms.1.clone(),vec![]);
            }
            edges.get_mut(&bb.atoms.0).unwrap().push(bb.atoms.1.clone());
            edges.get_mut(&bb.atoms.1).unwrap().push(bb.atoms.0.clone());
        }
        for ii in 0..num_atoms{
            let mut checked:HashSet<usize> = HashSet::new();
            
            let mut next:Vec<usize> = vec![ii];
            checked.insert(ii);
            for cc in 0..maxcon{
                let mut nnext:Vec<usize> = vec![];
                for nn in next.iter(){
                    if edges.contains_key(nn){
                        let nv:&Vec<usize> = edges.get(nn).unwrap();
                        for n in nv.iter(){
                            if !checked.contains(n){
                                distmap[ii][*n] = cc+1;
                                distmap[*n][ii] = cc+1;
                                nnext.push(*n);
                            }
                            checked.insert(*n);
                        }
                    }
                }
                next = nnext;
            }
        }
    }
}


//v0-v1-v2 の angle を radian で返す。
//あらかじめ standerdize しておく必要はない。
pub fn calc_angle_radian(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D)->f64{
    let mut vx0 = v0.get_x() - v1.get_x();
    let mut vy0 = v0.get_y() - v1.get_y();
    let mut vz0 = v0.get_z() - v1.get_z();

    let mut vx2 = v2.get_x() - v1.get_x();
    let mut vy2 = v2.get_y() - v1.get_y();
    let mut vz2 = v2.get_z() - v1.get_z();

    let mut dis0 = vx0*vx0+vy0*vy0+vz0*vz0;
    let mut dis2 = vx2*vx2+vy2*vy2+vz2*vz2;
    if dis0 == 0.0 || dis2 == 0.0{
        return 0.0;
    }
    dis0 = dis0.sqrt();
    dis2 = dis2.sqrt();
    vx0 /= dis0;
    vy0 /= dis0;
    vz0 /= dis0;
    
    vx2 /= dis2;
    vy2 /= dis2;
    vz2 /= dis2;

    let dz = process_3d::distance(&(vx0,vy0,vz0),&(vx2,vy2,vz2));
    
    return ((2.0-dz*dz)/2.0).max(-1.0).min(1.0).acos();
}

//v0-v1-v2 の angle を 360 degree で返す。
pub fn calc_angle(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D)->f64{
    return calc_angle_radian(v0,v1,v2)*180.0/PI;
}



pub fn calc_bond_energy(mdenv:&CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    let mut ubond:f64 = 0.0;
    for bb in cvars.bondvec.iter(){
        let dsc = bb.calc_energy(mdenv,atom_level_energy,1.0);
        ubond += dsc;
    }
    return ubond;
}


pub fn calc_angle_energy(mdenv:&mut CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    
    let mut uangle = 0.0;
    for aa in cvars.anglevec.iter(){
        //if dist[aa.atoms.0][aa.atoms.1] > 4.0 || dist[aa.atoms.1][aa.atoms.2] > 4.0{
        //    continue;
        //}
        let dsc = aa.calc_energy(mdenv,atom_level_energy,1.0);
        uangle += dsc;
    }
    return uangle;
}

pub fn calc_urey_bradley_energy(mdenv:&mut CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    let mut u_urey_bradley = 0.0;
    for aa in cvars.ubvec.iter(){
        let dsc = aa.calc_energy(mdenv,atom_level_energy,1.0);
        u_urey_bradley += dsc;
    }
    return u_urey_bradley;
}

pub fn calc_dihedral_energy(mdenv:&mut CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    let mut udihed = 0.0;
    for aa in cvars.dihedvec.iter(){
        let dsc = aa.calc_energy(mdenv,atom_level_energy,1.0);
        udihed += dsc;

    }
    return udihed;
}

pub fn calc_improper_energy(mdenv:&CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    let mut uimproper = 0.0;

    for aa in cvars.imprvec.iter(){
        let dsc = aa.calc_energy(mdenv,atom_level_energy,1.0);
        uimproper += dsc;
    }
    return uimproper;
}

#[derive(Debug,Clone)]
pub struct LJEnergyCalculator{
    pub nboption:charmm_param::Param_NONBONDED_OPTION
}
impl EnergyFunction for LJEnergyCalculator{
    fn calc_energy(&self,mdenv:&CharmmEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{
        let num_atoms:usize = atom_level_energy.len();
        let mut ulj = 0.0;
        for aa in 0..num_atoms-1{
            for bb in (aa+1)..num_atoms{
                if mdenv.num_edges[aa][bb] < 3{
                }else{
                    let ddis = mdenv.dist[aa][bb]+EPSILON;
                    //ToDo ctonnb と Sigmoidal function の導入
                    if ddis > self.nboption.ctofnb{
                        continue;
                    }
                    let epss:f64;
                    let rminij:f64;
                    if mdenv.num_edges[aa][bb] == 3{
                        if mdenv.atoms[aa].nb_14flag  && mdenv.atoms[bb].nb_14flag {
                            epss = (mdenv.atoms[aa].nb_14_epsilon*mdenv.atoms[bb].nb_14_epsilon).sqrt();
                            rminij =  mdenv.atoms[aa].nb_14_r1_2+mdenv.atoms[bb].nb_14_r1_2;
                        }else{
                            epss = (mdenv.atoms[aa].nb_epsilon*mdenv.atoms[bb].nb_epsilon).sqrt()*self.nboption.e14fac;
                            rminij =  mdenv.atoms[aa].nb_r1_2+mdenv.atoms[bb].nb_r1_2;
                        }
                    }else{
                        epss = (mdenv.atoms[aa].nb_epsilon
                            *mdenv.atoms[bb].nb_epsilon).sqrt();
                        rminij =  mdenv.atoms[aa].nb_r1_2+mdenv.atoms[bb].nb_r1_2;
                    }
                    
                    let bzz = (rminij/ddis).powf(6.0);
                    let dsc = epss*(bzz.powf(2.0) - 2.0*bzz)*weight;
                    //let dsc = 0.0;
                    atom_level_energy[aa] += dsc/2.0;
                    atom_level_energy[bb] += dsc/2.0;
                    
                    ulj += dsc;
                }
            }
        }
        return ulj;
    }
}

pub fn calc_lj_energy(mdenv:&CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    return cvars.lj_energy_calculator.calc_energy(mdenv,atom_level_energy,1.0);
}

#[derive(Debug,Clone)]
pub struct ElectrostaticEnergyCalculator{
    pub nboption:charmm_param::Param_NONBONDED_OPTION
}

impl EnergyFunction for ElectrostaticEnergyCalculator{
    fn calc_energy(&self,mdenv:&CharmmEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{
        let num_atoms:usize = atom_level_energy.len();
        //https://www.charmm.org/ubbthreads/ubbthreads.php?ubb=showflat&Number=12227
        //332*qi*qj/rij?
        //なんだか他と比べて異常に絶対値の大きい値が出てしまうが・・・
        //EVOEF は VDW a+b *0.8 未満は 0.8 にしているようだ。
        //ROSETTA は 1.45？
        let mut uelec = 0.0;
        for aa in 0..num_atoms-1{
            for bb in (aa+1)..num_atoms{
                if mdenv.num_edges[aa][bb] < 3{
                }else{
                    let ddis = mdenv.dist[aa][bb]+EPSILON;
                    if ddis > self.nboption.ctofnb{
                        continue;
                    }
                    let mut elec:f64;
                    if mdenv.num_edges[aa][bb] == 3{
                        elec = 332.0*mdenv.atoms[aa].charge*mdenv.atoms[bb].charge*self.nboption.e14fac/ddis;
                    }else{
                        elec = 332.0*mdenv.atoms[aa].charge*mdenv.atoms[bb].charge/ddis;
                    }

                    elec *= weight;
                    atom_level_energy[aa] += elec/2.0;
                    atom_level_energy[bb] += elec/2.0;
                    
                    uelec += elec;
                }
            }
        }
        return uelec;
    }
}
pub fn calc_electrostatic_energy(mdenv:&CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> f64{
    return cvars.electrostatic_energy_calculator.calc_energy(mdenv,atom_level_energy,1.0);
}

pub fn calc_energy(mdenv:&mut CharmmEnv,cvars:&CharmmVars,atom_level_energy:&mut Vec<f64>) -> MDEnergyResult{
    mdenv.update_distance();

    let num_atoms:usize = mdenv.atoms.len();
    assert_eq!(num_atoms , atom_level_energy.len());
    for aa in 0..num_atoms{
        atom_level_energy[aa] = 0.0;
    }
    
    let ubond = calc_bond_energy(mdenv,cvars,atom_level_energy);
    let uangle = calc_angle_energy(mdenv,cvars, atom_level_energy);
    let u_urey_bradley = calc_urey_bradley_energy(mdenv,cvars,atom_level_energy);
    let udihed = calc_dihedral_energy(mdenv,cvars, atom_level_energy);
    let uimproper = calc_improper_energy(mdenv,cvars,atom_level_energy);

    //CMAP は後回し
    let ucmap:f64 = 0.0;
    
    let ulj = calc_lj_energy(mdenv,cvars,atom_level_energy);
    let uelec = calc_electrostatic_energy(mdenv,cvars, atom_level_energy);

    return MDEnergyResult{
        ubond,
        uangle,
        udihed,
        u_urey_bradley,
        uimproper,
        ucmap,
        ulj,
        uelec
    };
    
}

#[derive(Debug)]
pub struct MDEnergyResult{
    ubond:f64,
    uangle:f64,
    udihed:f64,
    u_urey_bradley:f64,
    uimproper:f64,
    ucmap:f64,
    ulj:f64,
    uelec:f64
}

pub fn calc_dihedral_angle_radian(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D)->f64{
    let norm0:(f64,f64,f64) = process_3d::calc_norm(
        v0.get_x() - v1.get_x(),
        v0.get_y() - v1.get_y(),
        v0.get_z() - v1.get_z(),
        v2.get_x() - v1.get_x(),
        v2.get_y() - v1.get_y(),
        v2.get_z() - v1.get_z()
    );
    let norm1:(f64,f64,f64) = process_3d::calc_norm(
        v1.get_x() - v2.get_x(),
        v1.get_y() - v2.get_y(),
        v1.get_z() - v2.get_z(),
        v3.get_x() - v2.get_x(),
        v3.get_y() - v2.get_y(),
        v3.get_z() - v2.get_z()
    );
    let norm2:(f64,f64,f64) = process_3d::calc_norm(
        v2.get_x() - v1.get_x(),
        v2.get_y() - v1.get_y(),
        v2.get_z() - v1.get_z(),
        norm0.0,
        norm0.1,
        norm0.2,
    );
    let mut direc:f64 = 1.0;
    if process_3d::distance(&(norm1.0,norm1.1,norm1.2),&(norm2.0,norm2.1,norm2.2))
    >  (2.0_f64).sqrt(){
        direc = -1.0;
    }
    let drad = calc_angle_radian(
        &Point3D{x:norm0.0,y:norm0.1,z:norm0.2},
        &Point3D{x:0.0,y:0.0,z:0.0},
        &Point3D{x:norm1.0,y:norm1.1,z:norm1.2}
    )*direc;

    return drad;
}

pub fn calc_dihedral_angle(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D)->f64{
    return calc_dihedral_angle_radian(v0, v1, v2, v3)*180.0/PI;
}


//原子と原子を繋ぐ Bonds の配列と原子数を渡すと
//原子数分 Bond で繋がる原子のパスをすべて返す
//逆方向順方向はチェックし、同じパスが既にある場合は含めない
//同じ原子を二回通るパスは含めない
pub fn autogenerate_connection(bonds:&Vec<BondVars>,num_nodes:usize)->Vec<Vec<usize>>{
    let mut edges:HashMap<usize,Vec<usize>> = HashMap::new();
    for bb in bonds.iter(){
        if !edges.contains_key(&bb.atoms.0){
            edges.insert(bb.atoms.0.clone(),vec![]);
        }
        if !edges.contains_key(&bb.atoms.1){
            edges.insert(bb.atoms.1.clone(),vec![]);
        }
        edges.get_mut(&bb.atoms.0).unwrap().push(bb.atoms.1.clone());
        edges.get_mut(&bb.atoms.1).unwrap().push(bb.atoms.0.clone());
    }
    let allstart:Vec<usize> = edges.iter().map(|m|*m.0).collect();
    let mut res:Vec<Vec<usize>> = vec![];
    for a in allstart.iter(){
        let rr = get_all_path(&edges,&vec![*a],num_nodes);
        for r in rr.into_iter(){
            let mut dupcheck:bool = false;
            for rii in 0..r.len()-1{
                for rjj in (rii+1)..r.len(){
                    if r[rii] == r[rjj]{
                        dupcheck = true;
                    }
                }
            }
            if dupcheck{
                continue;
            }
            if r[0] < r[r.len()-1]{//全部チェックするので、必ず順方向逆方向二つある
                res.push(r);
            }
        }
    }
    return res;
}

pub fn get_all_path(next_atoms:&HashMap<usize,Vec<usize>>
    ,path:& Vec<usize>,maxlength:usize)->Vec<Vec<usize>>{
    
    let current = path[path.len()-1].clone();
    let next:&Vec<usize> = next_atoms.get(&current).as_ref().unwrap_or_else(||panic!("Can not find next atom! {} ",current));
    let mut ret:Vec<Vec<usize>> = vec![];

    for nn in next.iter(){
        if path.contains(&nn){
        }else{
            let mut ppath:Vec<usize> = path.iter().map(|m| m.clone()).collect();
            ppath.push(*nn);
            if ppath.len() >= maxlength{
                ret.push(ppath);
            }else{
                ret.append(&mut get_all_path(next_atoms,&ppath, maxlength));
            }
        }
    }
    return ret;
}


//0,1,2 の原子が分かっているときに 3 の原子の位置を見積もる
pub fn estimate_3_dihed(
    a0:&(f64,f64,f64),
    a1:&(f64,f64,f64),
    a2:&(f64,f64,f64),
    _bondlen_01:f64
    ,_angle_012_360:f64
    ,dihed_0123_360:f64,angle_123_360:f64
    ,bondlen_123:f64
)->(f64,f64,f64){
    let norm0 = process_3d::calc_norm_t(a0,a1,a2);
    let std_21:(f64,f64,f64) = process_3d::standardize(a1.0-a2.0,a1.1-a2.1,a1.2-a2.2);
    
    
    let v1 = Point3D::new(norm0.0,norm0.1,norm0.2);
    let std_21_p = Point3D::new(std_21.0,std_21.1,std_21.2);
    let mut v2 = Point3D::new(std_21.0,std_21.1,std_21.2);
    let mut v2r = Point3D::new(std_21.0,std_21.1,std_21.2);

    process_3d::rotate_3d(&mut vec![&mut v2],&v1,angle_123_360/180.0*PI);
    process_3d::rotate_3d(&mut vec![&mut v2r],&v1,angle_123_360/180.0*PI*-1.0);
    v2.add(&a2);
    v2r.add(&a2);
    
    if  process_3d::distance(&v2r.get_xyz(),&a0) > process_3d::distance(&v2.get_xyz(),&a0){
        v2 = v2r;
    }
    
    v2.x -= a2.0;
    v2.y -= a2.1;
    v2.z -= a2.2;

    v2.x *= bondlen_123;
    v2.y *= bondlen_123;
    v2.z *= bondlen_123;
    
    //2->1 方向なのでマイナス
    process_3d::rotate_3d(&mut vec![&mut v2],&std_21_p,(-dihed_0123_360-180.0)/180.0*PI);
    
    v2.x += a2.0;
    v2.y += a2.1;
    v2.z += a2.2;
    return (v2.x,v2.y,v2.z);
}

pub fn estimate_3_impr(
    a0:&(f64,f64,f64),
    a1:&(f64,f64,f64),
    a2:&(f64,f64,f64),
    _bondlen_01:f64
    ,_angle_012_360:f64
    ,dihed_0123_360:f64
    ,angle_203_360:f64
    ,bondlen_23:f64
)->(f64,f64,f64){
    let norm0 = process_3d::calc_norm_t(a0,a2,a1);

    let pnorm0 = Point3D::new(norm0.0,norm0.1,norm0.2);
    let angle_203:f64 = angle_203_360/180.0*PI;
    
    let mut tmp3:Point3D = Point3D::new(a1.0-a2.0,a1.1-a2.1,a1.2-a2.2);
    let mut tmp3b:Point3D = Point3D::new(a1.0-a2.0,a1.1-a2.1,a1.2-a2.2);
    tmp3.standardize();
    tmp3b.standardize();
    tmp3.set_x(tmp3.get_x()*bondlen_23);
    tmp3.set_y(tmp3.get_y()*bondlen_23);
    tmp3.set_z(tmp3.get_z()*bondlen_23);
    
    tmp3b.set_x(tmp3b.get_x()*bondlen_23);
    tmp3b.set_y(tmp3b.get_y()*bondlen_23);
    tmp3b.set_z(tmp3b.get_z()*bondlen_23);


    process_3d::rotate_3d(&mut vec![&mut tmp3],&pnorm0,angle_203);
    process_3d::rotate_3d(&mut vec![&mut tmp3b],&pnorm0,angle_203*-1.0);
    
    tmp3.add(&(a2.0,a2.1,a2.2));
    tmp3b.add(&(a2.0,a2.1,a2.2));
    //どっちか一方しか取れないかも？
    if process_3d::distance(&tmp3b.get_xyz(),&a0) 
    >  process_3d::distance(&tmp3.get_xyz(),&a0){
        tmp3 = tmp3b;
    }

    tmp3.add(&(-a2.0,-a2.1,-a2.2));


    let mut tmp1:Point3D = Point3D::new(a1.0-a2.0,a1.1-a2.1,a1.2-a2.2);
    tmp1.standardize();

    //dihed と同じ 0-1-2-3 の角度ならこれで OK
    process_3d::rotate_3d(&mut vec![&mut tmp3],&tmp1,(-dihed_0123_360-180.0)/180.0*PI);
    
    tmp3.add(&(a2.0,a2.1,a2.2));
    
    return tmp3.get_xyz();
}

#[test]
fn testimpr3d(){
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(10);
    let p1:Point3D = Point3D::new(
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0)
    );
    let p2:Point3D = Point3D::new(
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0)
    );
    let p3:Point3D = Point3D::new(
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0),
        rgen.gen_range(-30.0,30.0)
    );
    estimate_3_impr(
        &p1.get_xyz(),
        &p2.get_xyz(),
        &p3.get_xyz()
        ,1.0
        ,1.0
        ,30.0
        ,60.0
        ,20.0
    );

}



#[test]
fn bondcheck(){
    let bonds:Vec<(usize,usize)> = vec![
        (0,1),
        (0,2),
        (2,3),
        (2,4),
        (2,5),
    ];

    let bbonds:Vec<BondVars> = bonds.iter().map(|m|BondVars{atoms:(m.0,m.1),kb:0.0,b0:0.0}).collect();
    let res = autogenerate_connection(&bbonds,4);
    assert_eq!(res.len(), 3);
    assert!(res.contains(&vec![1,0,2,3]));
    assert!(res.contains(&vec![1,0,2,4]));
    assert!(res.contains(&vec![1,0,2,5]));
    
    let res = autogenerate_connection(&bbonds,3);
    assert_eq!(res.len(),7);
    assert!(res.contains(&vec![1,0,2]));
    assert!(res.contains(&vec![0,2,3]));
    assert!(res.contains(&vec![0,2,4]));
    assert!(res.contains(&vec![0,2,5]));
    assert!(res.contains(&vec![3,2,4]));
    assert!(res.contains(&vec![3,2,5]));
    assert!(res.contains(&vec![4,2,5]));
    
    
    let bonds:Vec<(usize,usize)> = vec![
        (0,1),
        (1,2),
        (2,3),
        (3,6),
        (6,5),
        (5,4),
        (4,1),
    ];
    
    let bbonds:Vec<BondVars> = bonds.iter().map(|m|BondVars{atoms:(m.0,m.1),kb:0.0,b0:0.0}).collect();
    let res = autogenerate_connection(&bbonds,4);
    assert_eq!(res.len(),8);
    assert!(res.contains(&vec![0,1,2,3]));
    assert!(res.contains(&vec![0,1,4,5]));
    assert!(res.contains(&vec![1,2,3,6]));
    assert!(res.contains(&vec![1,4,5,6]));
    assert!(res.contains(&vec![2,1,4,5]));
    assert!(res.contains(&vec![2,3,6,5]));
    assert!(res.contains(&vec![3,2,1,4]));
    assert!(res.contains(&vec![3,6,5,4]));

    let bbonds:Vec<BondVars> = bonds.iter().map(|m|BondVars{atoms:(m.0,m.1),kb:0.0,b0:0.0}).collect();
    let res = autogenerate_connection(&bbonds,3);
    assert_eq!(res.len(),8);
    assert!(res.contains(&vec![0,1,2]));
    assert!(res.contains(&vec![0,1,4]));
    assert!(res.contains(&vec![1,2,3]));
    assert!(res.contains(&vec![1,4,5]));
    assert!(res.contains(&vec![2,1,4]));
    assert!(res.contains(&vec![2,3,6]));
    assert!(res.contains(&vec![3,6,5]));
    assert!(res.contains(&vec![4,5,6]));
    

    let bonds:Vec<(usize,usize)> = vec![
        (0,1),
        (1,2),
        (2,0),
    ];
    
    let bbonds:Vec<BondVars> = bonds.iter().map(|m|BondVars{atoms:(m.0,m.1),kb:0.0,b0:0.0}).collect();
    let res = autogenerate_connection(&bbonds,4);
    assert_eq!(res.len(),0);
    let res = autogenerate_connection(&bbonds,3);
    assert_eq!(res.len(),3);
    assert!(res.contains(&vec![0,1,2]));
    assert!(res.contains(&vec![0,2,1]));
    assert!(res.contains(&vec![1,0,2]));
}


#[test]
fn mdprepare_test(){
    let mut topp:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_top_all22_inp((debug_env::CHARMM_DIR.to_string()+"\\top_all22_prot.rtf").as_str());
    let mut parr:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_param((debug_env::CHARMM_DIR.to_string()+"\\par_all22_prot.prm").as_str());
    parr.resi.append(&mut topp.resi);
    let allaa:Vec<String> = vec![
        "ALA","ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","VAL",
    ].iter().map(|m|m.to_string()).collect();

    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);

    let mut ress:Vec<pdbdata::PDBComp> = chain_builder::build_dirty_chain(&allaa,&bset,&sset);
    let mut chain:pdbdata::PDBAsym = pdbdata::PDBAsym::new("A");
    for rr in ress.into_iter(){
        chain.add_comp(rr);
    }

    MDAtom::change_to_charmmnames(&mut chain.iter_mut_comps().map(|m|*m).collect());
    let (res,_cvars):(CharmmEnv,CharmmVars) = MDAtom::chain_to_atoms(&vec![&chain],&parr,true);
    let mut lines:Vec<String> = vec![];
    let mut sp3_lines:Vec<String> = vec![];
    for aa in res.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
        sp3_lines.push(format!("{}\t{}\t{:?}",resname,aa.atom_name,aa.hybrid_orbit));
    }
    write_to_file("test/fullatomout.pdb",lines);
    //pdbdata::write_to_file("test/atom_hybrid.dat",sp3_lines);
}





#[test]
fn add_h_test(){
    let mut pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/6iws_model1.pdb");
    let mut topp:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_top_all22_inp((debug_env::CHARMM_DIR.to_string()+"\\top_all22_prot.rtf").as_str());
    let mut parr:charmm_param::CHARMMParam = charmm_param::CHARMMParam::load_param((debug_env::CHARMM_DIR.to_string()+"\\par_all22_prot.prm").as_str());
    parr.resi.append(&mut topp.resi);

    
    MDAtom::change_to_charmmnames(&mut pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_mut_comps().map(|m|{*m}).collect());
    let mut resvec:Vec<pdbdata::PDBComp> = vec![];
    for rr in pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_mut_comps(){
        let mut ress:pdbdata::PDBComp = pdbdata::PDBComp::new();
        ress.set_comp_id(rr.get_name());
        ress.set_seq_id(rr.get_seq_id());
        for a in rr.iter_atoms(){
            if start_with(&a.atom_code,"H"){
            }else{
                ress.add_atom(a.clone());
            }
        }
        resvec.push(ress);
    }
    let mut dchain:pdbdata::PDBAsym = pdbdata::PDBAsym::new("A");
    for rr in resvec.into_iter(){
        dchain.add_comp(rr);
    }
    let (md_envset,_md_varset):(CharmmEnv,CharmmVars) = MDAtom::chain_to_atoms(&vec![&dchain],&parr,true);
    
    let mut lines:Vec<String> = vec![];
    for aa in md_envset.atoms.iter(){
        let (chainid,(resname,resnum,altcode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file("test/hadded_b.pdb",lines);
}





#[test]
fn check_backbone_dihed(){
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    let mut pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/6iws_model1_noh.pdb");
    MDAtom::change_to_charmmnames(&mut pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_mut_comps().map(|m|*m).collect());
    let (md_envset,md_varset):(CharmmEnv,CharmmVars) = MDAtom::chain_to_atoms(&vec![pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0)],&parr,true);

    let backbone_dihedrals:Vec<&str> = vec![
        "C#N#CA#C",
        "N#CA#C#N",
        "CA#C#N#CA",
        
        "C#CA#N#C",
        "N#C#CA#N",
        "CA#N#C#CA",
        
    ];    
    
    for aa in md_varset.dihedvec.iter(){

        let acode:String = format!("{}#{}#{}#{}",
        md_envset.atoms[aa.atoms.0].atom_name
        ,md_envset.atoms[aa.atoms.1].atom_name
        ,md_envset.atoms[aa.atoms.2].atom_name
        ,md_envset.atoms[aa.atoms.3].atom_name);
        for vv in backbone_dihedrals.iter(){
            if *vv == acode{
                let rescode:String = format!("{}#{}#{}#{}"
                ,md_envset.atoms[aa.atoms.2].chain_name
                ,md_envset.atoms[aa.atoms.2].residue_name
                ,md_envset.atoms[aa.atoms.2].residue_number
                ,md_envset.atoms[aa.atoms.2].residue_ins_code);
                println!("{}",rescode);
            }
            
        }
    }

    /*
    pub atom_type:String,
    pub atom_name:String,
    pub chain_name:String,
    pub residue_name:String,
    pub residue_number:i64,
    pub residue_ins_code:String,
    pub residue_index_in_chain:i64,
    */

}
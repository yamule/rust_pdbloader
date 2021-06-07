#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(unused_imports)]
use std::collections::{HashMap,HashSet};
#[allow(unused_imports)]
use regex::Regex;
#[allow(unused_imports)]
use super::pdbdata;
#[allow(unused_imports)]
use super::misc_util::*;
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use super::debug_env;
#[allow(unused_imports)]
use super::ff_env;
#[allow(unused_imports)]
use super::geometry::Vector3D;
#[allow(unused_imports)]
use super::geometry::Point3D;
use rand::prelude::*;
use rand_distr::{Normal, Distribution};

use super::energy_function::EnergyFunction;
use super::energy_function::Binning;
use super::energy_function::ROUNDING_EPSILON;

use std::cmp::Ordering;
use std::f64::consts::PI;


//-180~180
//opposite direction with charmm
pub fn calc_dihedral_angle(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D)->f64{
    return ff_env::calc_dihedral_angle(v0, v1, v2, v3)*-1.0;
}

pub fn validate_angle(vv:f64)->f64{
    let mut v = vv;
    while v > 180.0{
        v -= 360.0;
    }
    while v < -180.0{
        v += 360.0;
    }
    return v;
}
//与えられた val が入るインデックスを返す。
pub fn get_vec_index_angle(bbin:&dyn Binning,val:f64)-> usize{
    let mut v = validate_angle(val);
    v += -1.0*bbin.get_start_point();
    v -= bbin.get_div_unit()/2.0;
    v += ROUNDING_EPSILON;

    let ret:i64 = (v/ bbin.get_div_unit()).round() as i64;
    return ret.max(0).min(bbin.get_num_bins() as i64 -1) as usize;
}   


#[derive(Debug,Clone)] 
pub struct AnglePhiPsi{
    divunit:f64,
    start_point:f64,
    //c_n_ca_c_n
    pub atoms:(usize,usize,usize,usize,usize),   
    phi_psi:Vec<Vec<f64>>,
    sorted:Vec<(usize,usize)>,
    fixed_value_180:Option<(f64,f64)>
}

impl AnglePhiPsi{
    pub fn new(start_point_divunit:(f64,f64),
        atoms:(usize,usize,usize,usize,usize),   
        phi_psi:Vec<Vec<(f64,bool)>>,//値が無い場合に true
        fixedvalue:Option<(f64,f64)>)->AnglePhiPsi{
        let rnum = phi_psi.len();
        let mut tops:Vec<(usize,usize)> = vec![];
        if rnum > 0{
            let cnum = phi_psi[0].len();
            let mut sorter:Vec<(usize,usize,f64)> = vec![];
            for xx in 0..rnum{
                for yy in 0..cnum{
                    if phi_psi[xx][yy].1{
                        sorter.push((xx,yy,phi_psi[xx][yy].0));
                    }
                }
            }
            sorter.sort_by(|a,b|{a.2.partial_cmp(&b.2).unwrap()});
            tops = sorter.into_iter().map(|m|(m.0,m.1)).collect();
        }
        return AnglePhiPsi{
            start_point:start_point_divunit.0,
            divunit:start_point_divunit.1,
            atoms:atoms,
            phi_psi:phi_psi.iter().map(|m| m.iter().map(|m|m.0).collect()).collect(),
            sorted:tops,
            fixed_value_180:fixedvalue
        };
    }

    fn get_angle_180(&self,phi_u:usize,psi_u:usize)->(f64,f64){
        let phi:f64 = self.divunit*phi_u as f64 + self.start_point +self.divunit/2.0;
        let psi:f64 = self.divunit*psi_u as f64 + self.start_point +self.divunit/2.0;
        return (phi,psi);
    }

    fn get_random_angle_180(&self,normal:&Normal<f64>,rgen:&mut StdRng,restrict_top:usize)->(f64,f64){
        if let Some(x) = self.fixed_value_180{
            return x;
        }    
        if self.sorted.len() == 0{
            return (0.0,0.0);
        }
        let mut v:usize = ((normal.sample(rgen)/3.0*(self.sorted.len() as f64)).abs() as usize).min(self.sorted.len()-1);
        while v >= restrict_top{
            v = ((normal.sample(rgen)/3.0*(self.sorted.len() as f64)).abs() as usize).min(self.sorted.len()-1);
        }
        return self.get_angle_180(self.sorted[v].0,self.sorted[v].1);
    }
    
    pub fn get_random_angle_radian(&self,normal:&Normal<f64>,rgen:&mut StdRng,restrict_top:usize)->(f64,f64){
        let ret = self.get_random_angle_180(normal,rgen,restrict_top);
        return (ret.0/180.0*PI,ret.1/180.0*PI);
    }
}
impl Binning for AnglePhiPsi{
    fn get_div_unit(&self)->f64{
        return self.divunit;
    }

    fn get_start_point(&self)->f64{
        return self.start_point;
    }
    
    fn get_num_bins(&self)->usize{
        return self.phi_psi.len();
    }
}


pub fn get_linear_interpolated_value_phi_psi(bbin:&dyn Binning,val_phi_psi:(f64,f64),bins:&Vec<Vec<f64>>)-> f64{
    let mut v:(f64,f64) = val_phi_psi;
    v.0 += -1.0*bbin.get_start_point();
    v.1 += -1.0*bbin.get_start_point();
    
    let mut ret:(i64,i64) = ((v.0/ bbin.get_div_unit()).floor() as i64
    ,(v.1/ bbin.get_div_unit()).floor() as i64);

    let lb:(f64,f64) = (bbin.get_start_point()+bbin.get_div_unit()*(ret.0 as f64),bbin.get_start_point()+bbin.get_div_unit()*(ret.1 as f64));
    let ub:(f64,f64) = (bbin.get_start_point()+bbin.get_div_unit()*(ret.0 as f64 +1.0),bbin.get_start_point()+bbin.get_div_unit()*(ret.1 as f64 +1.0));
    
    while ret.0 >= bbin.get_num_bins() as i64 -1{
        ret.0 -= bbin.get_num_bins() as i64;
    }
    while ret.1 >= bbin.get_num_bins() as i64 -1{
        ret.1 -= bbin.get_num_bins() as i64;
    }
    
    while ret.0 < 0{
        ret.0 += bbin.get_num_bins() as i64;
    }
    while ret.1 < 0{
        ret.1 += bbin.get_num_bins() as i64;
    }
    
    let lratio:(f64,f64) = ((val_phi_psi.0-lb.0)/bbin.get_div_unit(),(val_phi_psi.1-lb.1)/bbin.get_div_unit());
    let uratio:(f64,f64) = ((ub.0-val_phi_psi.0)/bbin.get_div_unit(),(ub.1-val_phi_psi.1)/bbin.get_div_unit());

    assert!(lratio.0 >= 0.0);
    assert!(uratio.0 >= 0.0);
    assert!(lratio.1 >= 0.0);
    assert!(uratio.1 >= 0.0);

    let uphi:usize = if ret.0 as usize+1 >= bbin.get_num_bins(){0}else{ret.0 as usize+1};
    let upsi:usize = if ret.1 as usize+1 >= bbin.get_num_bins(){0}else{ret.1 as usize+1};

    return ((lratio.0+lratio.1)*bins[uphi][upsi]
    +(lratio.0+uratio.1)*bins[uphi][ret.1 as usize]
    +(uratio.0+lratio.1)*bins[ret.0 as usize][upsi]
    +(uratio.0+uratio.1)*bins[ret.0 as usize][ret.1 as usize]
    )/(lratio.0+uratio.0+lratio.1+uratio.1);
}   

pub fn get_linear_interpolated_value_dihed(bbin:&dyn Binning,val:f64,bins:&Vec<f64>)-> f64{
    let mut v = val;
    v += -1.0*bbin.get_start_point();
    
    let mut ret:i64 = (v/ bbin.get_div_unit()).floor() as i64;

    let lb:f64 = bbin.get_start_point()+bbin.get_div_unit()*(ret as f64);
    let ub:f64 = bbin.get_start_point()+bbin.get_div_unit()*(ret as f64 +1.0);
    while ret >= bbin.get_num_bins() as i64 -1{
        ret -= bbin.get_num_bins() as i64;
    }
    while ret < 0{
        ret += bbin.get_num_bins() as i64;
    }
    let lratio = (val-lb)/bbin.get_div_unit();
    let uratio = (ub-val)/bbin.get_div_unit();

    assert!(lratio >= 0.0);
    assert!(uratio >= 0.0);
    let upp:usize = if ret as usize +1 >= bbin.get_num_bins(){0}else{ret as usize +1 };
    return (lratio*bins[upp]+uratio*bins[ret as usize])/(lratio+uratio);
} 

impl EnergyFunction for AnglePhiPsi{
    /*fn calc_energy(&self,mdenv:&ff_env::FFEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{
        
        let phi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.0], &mdenv.atoms[self.atoms.1],&mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3]);
        let psi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.1], &mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3],&mdenv.atoms[self.atoms.4]);
        let phiindex:usize = get_vec_index_angle(self,phi);
        let psiindex:usize = get_vec_index_angle(self,psi);
        let mut dsc:f64 = self.phi_psi[phiindex][psiindex];
        dsc *= weight;
        
        atom_level_energy[self.atoms.0] += dsc/5.0;
        atom_level_energy[self.atoms.1] += dsc/5.0;
        atom_level_energy[self.atoms.2] += dsc/5.0;
        atom_level_energy[self.atoms.3] += dsc/5.0;
        atom_level_energy[self.atoms.4] += dsc/5.0;
        return dsc;
    }*/
    fn calc_energy(&self,mdenv:&ff_env::FFEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{
        
        let phi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.0], &mdenv.atoms[self.atoms.1],&mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3]);
        let psi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.1], &mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3],&mdenv.atoms[self.atoms.4]);
        let mut dsc:f64 = get_linear_interpolated_value_phi_psi(self,(phi,psi),&self.phi_psi);
        dsc *= weight;
        
        atom_level_energy[self.atoms.0] += dsc/5.0;
        atom_level_energy[self.atoms.1] += dsc/5.0;
        atom_level_energy[self.atoms.2] += dsc/5.0;
        atom_level_energy[self.atoms.3] += dsc/5.0;
        atom_level_energy[self.atoms.4] += dsc/5.0;
        return dsc;
    }
}



#[derive(Debug,Clone)] 
pub struct AngleDihed{
    start_point:f64,
    divunit:f64,
    //ca_c_n_ca
    pub atoms:(usize,usize,usize,usize),
    dihed:Vec<f64>,
    sorted:Vec<usize>,
    fixed_value_180:Option<f64>
}
impl AngleDihed{
    pub fn new(start_point_divunit:(f64,f64),
        atoms:(usize,usize,usize,usize),   
        dihed:Vec<(f64,bool)>,//value, has value
        fixedvalue:Option<f64>)->AngleDihed{
        let rnum = dihed.len();
        let mut tops:Vec<usize> = vec![];
        if rnum > 0{
            let mut sorter:Vec<(usize,f64)> = vec![];
            for xx in 0..rnum{
                if dihed[xx].1{
                    sorter.push((xx,dihed[xx].0));
                }
            }
            sorter.sort_by(|a,b|{a.1.partial_cmp(&b.1).unwrap()});
            tops = sorter.into_iter().map(|m|m.0).collect();
        }
        return AngleDihed{
            start_point:start_point_divunit.0,
            divunit:start_point_divunit.0,
            atoms:atoms,
            dihed:dihed.iter().map(|m|m.0).collect(),
            fixed_value_180:fixedvalue,
            sorted:tops
        };
    }

    pub fn get_angle_180(&self,phi_u:usize)->f64{
        let phi:f64 = self.divunit*phi_u as f64 + self.start_point +self.divunit/2.0;
        return phi;
    }

    pub fn get_random_angle_180(&self,normal:&Normal<f64>,rgen:&mut StdRng,restrict_top:usize)->f64{
        if let Some(x) = self.fixed_value_180{
            return x;
        }
        if self.sorted.len() == 0{
            return 0.0;
        }
        let mut v:usize = ((normal.sample(rgen)/3.0*(self.sorted.len() as f64)).abs() as usize).min(self.sorted.len()-1);
        while v >= restrict_top{
            v = ((normal.sample(rgen)/3.0*(self.sorted.len() as f64)).abs() as usize).min(self.sorted.len()-1);
        }
        return self.get_angle_180(self.sorted[v]);
    }

    pub fn get_random_angle_radian(&self,normal:&Normal<f64>,rgen:&mut StdRng,restrict_top:usize)->f64{
        let ret = self.get_random_angle_180(normal,rgen,restrict_top);
        return ret/180.0*PI;
    }
}

impl Binning for AngleDihed{
    fn get_div_unit(&self)->f64{
        return self.divunit;
    }
    fn get_num_bins(&self)->usize{
        return self.dihed.len();
    }
    fn get_start_point(&self)->f64{
        return self.start_point;
    }
}

impl EnergyFunction for AngleDihed{
    fn calc_energy(&self,mdenv:&ff_env::FFEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{
        
        let phi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.0], &mdenv.atoms[self.atoms.1],&mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3]);
        let mut dsc:f64 = get_linear_interpolated_value_dihed(self,phi,&self.dihed);
        
        dsc *= weight;

        atom_level_energy[self.atoms.0] += dsc/4.0;
        atom_level_energy[self.atoms.1] += dsc/4.0;
        atom_level_energy[self.atoms.2] += dsc/4.0;
        atom_level_energy[self.atoms.3] += dsc/4.0;

        return dsc;
    }
}




#[derive(Debug,Clone)]
pub struct PlainDistribution{
    divunit:f64,
    undefvalue:f64,
    start_point:f64,
    num_bins:usize,
    residue_name:Option<String>,
    residue_number:Option<i64>,
    residue_ins_code:Option<String>,
    phi_psi:Option<Vec<Vec<(f64,bool)>>>,
    prev_omega:Option<Vec<(f64,bool)>>,
    next_omega:Option<Vec<(f64,bool)>>
}
impl Binning for PlainDistribution{
    fn get_div_unit(&self)->f64{
        return self.divunit;
    }
    fn get_num_bins(&self)->usize{
        return self.num_bins;
    }
    fn get_start_point(&self)->f64{
        return self.start_point;
    }
}

impl PlainDistribution{
    pub fn new(resname_:&str,residue_number_:&str,residue_ins_code_:&str,start_point:f64,divunit:f64,undefval:f64)->PlainDistribution{
        let psiz:usize = (360.0/divunit).round() as usize;
        
        let resname:Option<String> = if resname_ != ""{
            Some(resname_.to_string())
        }else{
            None
        };
        let residue_number:Option<i64> = if residue_number_ != ""{
                Some(residue_number_.parse::<i64>().unwrap_or_else(|_| panic!("Can not parse {}.",residue_number_)))
        }else{
            None
        };
        let residue_ins_code:Option<String> = if residue_ins_code_ != ""{
                Some(residue_ins_code_.to_string())
        }else{
            None
        };
        return PlainDistribution{
            residue_name:resname,
            residue_number:residue_number,
            residue_ins_code:residue_ins_code,
            start_point:start_point,
            divunit:divunit,
            num_bins:psiz,
            undefvalue:undefval,
            phi_psi:None,
            prev_omega:None,
            next_omega:None
        }
    }
    
    pub fn prepare_phi_psi(&mut self){
        if let None = self.phi_psi{
            self.phi_psi = Some(vec![vec![(self.undefvalue,false);self.num_bins];self.num_bins]);
        }else{
            panic!("phi_psi has already been assigned!");
        }
    }
    
    pub fn prepare_prev_omega(&mut self){
        if let None = self.prev_omega{
            self.prev_omega = Some(vec![(self.undefvalue,false);self.num_bins]);
        }else{
            panic!("prev_omega has already been assigned!");
        }
    }

    pub fn prepare_next_omega(&mut self){
        if let None = self.next_omega{
            self.next_omega = Some(vec![(self.undefvalue,false);self.num_bins]);
        }else{
            panic!("next_omega has already been assigned!");
        }
    }

    pub fn set_phi_psi(&mut self,px:usize,py:usize,val:f64,has_value:bool){
        self.phi_psi.as_mut().unwrap()[px][py] = (val,has_value);
    }

    pub fn set_prev_omega(&mut self,px:usize,val:f64,has_value:bool){
        self.prev_omega.as_mut().unwrap()[px] = (val,has_value);
    }
    
    pub fn set_next_omega(&mut self,px:usize,val:f64,has_value:bool){
        self.next_omega.as_mut().unwrap()[px] = (val,has_value);
    }

    pub fn match_with(&self,a:&ff_env::FFAtom)->bool{
        if let Some(x) = self.residue_name.as_ref(){
            if *x != a.residue_name{
                return false;
            }
        }
        if let Some(x) = self.residue_number{
            if x != a.residue_number{
                return false;
            }
        }
        if let Some(x) = self.residue_ins_code.as_ref(){
            if *x != a.residue_ins_code{
                return false;
            }
        }
        return true;
    }
    
    fn create_residue_label(chainname:&str,resname:&str,resnum:&str)->String{
        return format!("#{};{};{}#",chainname,resname,resnum);
    }

    //オレオレフォーマットのファイルを読み込む
    pub fn load_name_mapped(filename:&str)->HashMap<String,PlainDistribution>{
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:HashMap<String,PlainDistribution> = HashMap::new();
        let mut divunit:f64 = -1.0;
        let mut start_point_:Option<f64> = None;
        let mut undefvalue_:Option<f64> = None;//指定のない場合にあてられるエネルギー
        let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
        for (_lcount,line) in lines.iter().enumerate() {
            if start_with(line,"#"){
                continue;
            }
            let hs = line_to_hash(line);
            if hs.contains_key("divunit"){
                divunit = hs.get("divunit").unwrap().parse::<f64>().expect("Can not parse divunit.");
            }
            if hs.contains_key("start_point"){
                start_point_ = Some(hs.get("start_point").unwrap().parse::<f64>().expect("Can not parse start_point."));
            }
            if hs.contains_key("undefval"){
                undefvalue_ = Some(hs.get("undefval").unwrap().parse::<f64>().expect("Can not parse undefval."));
            }
        }
        if divunit < 0.0{
            panic!("Can not find divunit. (binsize)");
        }
        
        let start_point:f64 = start_point_.unwrap_or_else(||panic!("Can not find start_point."));
        let undefvalue:f64 = undefvalue_.unwrap_or_else(||panic!("Can not find undefval."));

        for (_lcount,line) in lines.iter().enumerate() {
            if start_with(line,"#"){
                continue;
            }
            let hs = line_to_hash(line);
            if hs.contains_key("values"){
                let chainname_:String = hs.get("chainname").unwrap_or(&("".to_string())).to_string();
                let resname_:String = hs.get("resname").unwrap_or(&("".to_string())).to_string();
                let resnum_:String = hs.get("residx").unwrap_or(&("".to_string())).to_string();
                let inscode_:String = hs.get("inscode").unwrap_or(&("".to_string())).to_string();
                
                if chainname_.len() + resname_.len() + resnum_.len() + inscode_.len() == 0{
                    panic!("The line does not have valid label; chainname, resname, resnum {:?}",hs);
                }
                if inscode_.len() > 0{
                    panic!("The line indicates index of the residues, not res num of pdb!");
                }

                let reslabel:String = PlainDistribution::create_residue_label(&chainname_,&resname_,&resnum_);
                
                if !ret.contains_key(&reslabel){
                    let rr = PlainDistribution::new(&resname_,&resnum_,&inscode_,start_point,divunit,undefvalue);
                    ret.insert(reslabel.clone(),rr);
                }
                let typ:&str =  hs.get("type").expect("Can not find type section.");
                if typ == "phi_psi"{
                    let val:Vec<&str> = hs.get("values").expect("Can not find values section.").split("#").collect();
                    ret.get_mut(&reslabel).unwrap().prepare_phi_psi();
                    for  vv in val.iter(){
                        if let Some(_) = REGEX_NOLINE.captures(vv){
                            continue;
                        }
                        //x_y,value という形になっていることを想定している。
                        let gval:Vec<&str> = vv.split(",").collect();
                        let pval:Vec<&str> = gval[0].split("_").collect();
                        let px:usize = pval[0].parse::<usize>().unwrap();
                        let py:usize = pval[1].parse::<usize>().unwrap();
                        let val:f64 = gval[1].parse::<f64>().unwrap();
                        ret.get_mut(&reslabel).unwrap().set_phi_psi(px,py,val,true);
                    }
                }else if typ == "prev_omega"{
                    let val:Vec<&str> = hs.get("values").expect("Can not find values section.").split("#").collect();
                    ret.get_mut(&reslabel).unwrap().prepare_prev_omega();
                    for  vv in val.iter(){
                        if let Some(_) = REGEX_NOLINE.captures(vv){
                            continue;
                        }
                        let gval:Vec<&str> = vv.split(",").collect();
                        let px:usize = gval[0].parse::<usize>().unwrap();
                        let val:f64 = gval[1].parse::<f64>().unwrap();
                        ret.get_mut(&reslabel).unwrap().set_prev_omega(px,val,true);
                    }
                }else if typ == "next_omega"{
                    let val:Vec<&str> = hs.get("values").expect("Can not find values section.").split("#").collect();
                    ret.get_mut(&reslabel).unwrap().prepare_next_omega();
                    for  vv in val.iter(){
                        if let Some(_) = REGEX_NOLINE.captures(vv){
                            continue;
                        }
                        let gval:Vec<&str> = vv.split(",").collect();
                        let px:usize = gval[0].parse::<usize>().unwrap();
                        let val:f64 = gval[1].parse::<f64>().unwrap();
                        ret.get_mut(&reslabel).unwrap().set_next_omega(px,val,true);
                    }
                }
            }
        }
        return ret;
    }


    //mdenv に入った atom をソートして、Chain 毎に分割し、タンパク質の Backbone にあたる原子のみ抽出、Phi psi, omega を形成する原子のインデクスを持った
    //それぞれのインスタンスを作成する
    pub fn create_energy_instance(distribution:&HashMap<String,PlainDistribution>
        ,mdenv:& ff_env::FFEnv
        ,_keycheck:(bool,bool,bool)//chain name, residue name, residue_index_in_chain どの値が distribution のキーに使用されているかを指定する。古い
        ,fix_omega:bool)->(Vec<AnglePhiPsi>,Vec<AngleDihed>){
        let mut chains:HashMap<String,Vec<&ff_env::FFAtom>> = HashMap::new();
        let mut ret_phi_psi:Vec<AnglePhiPsi> = vec![];
        let mut ret_omega:Vec<AngleDihed> = vec![];

        for aa in mdenv.atoms.iter(){
            if !chains.contains_key(&aa.chain_name){
                chains.insert(aa.chain_name.clone(),vec![]);
            }
            chains.get_mut(&aa.chain_name).unwrap().push(aa);
        }
        let mut chain_names:Vec<String> = chains.iter().map(|m|m.0.to_string()).collect();
        chain_names.sort();
        for (_cname,avec) in chains.iter(){
            if avec.len() <= 1{
                continue;
            }
            let mut atom_n:Vec<&ff_env::FFAtom> = vec![];
            let mut atom_ca:Vec<&ff_env::FFAtom> = vec![];
            let mut atom_c:Vec<&ff_env::FFAtom> = vec![];
            for aa in avec.iter(){
                if aa.atom_name == "N"{
                    atom_n.push(aa);
                }
                if aa.atom_name == "CA"{
                    atom_ca.push(aa);
                }
                if aa.atom_name == "C"{
                    atom_c.push(aa);
                }
            }
            assert_eq!(atom_n.len(),atom_ca.len());
            assert_eq!(atom_ca.len(),atom_c.len());
            atom_n.sort_by(|a,b|{a.residue_index_in_chain.partial_cmp(&b.residue_index_in_chain).unwrap()});
            atom_ca.sort_by(|a,b|{a.residue_index_in_chain.partial_cmp(&b.residue_index_in_chain).unwrap()});
            atom_c.sort_by(|a,b|{a.residue_index_in_chain.partial_cmp(&b.residue_index_in_chain).unwrap()});

            let alen:usize = atom_n.len();
            for ii in 0..alen{
                assert_eq!(atom_n[ii].residue_index_in_chain,atom_ca[ii].residue_index_in_chain);
                assert_eq!(atom_ca[ii].residue_index_in_chain,atom_c[ii].residue_index_in_chain);
            }
            for ii in 0..alen{
                
                let reslabel0:String =  PlainDistribution::create_residue_label(&atom_ca[ii].chain_name,&atom_ca[ii].residue_name,&atom_ca[ii].residue_index_in_chain.to_string());
                let reslabel1:String =  PlainDistribution::create_residue_label(&atom_ca[ii].chain_name,"",&atom_ca[ii].residue_index_in_chain.to_string());
                let reslabel2:String =  PlainDistribution::create_residue_label("","",&atom_ca[ii].residue_index_in_chain.to_string());
                let reslabel3:String =  PlainDistribution::create_residue_label("",&atom_ca[ii].residue_name,"");


                let diss:&PlainDistribution = distribution.get(&reslabel0).unwrap_or_else(||{
                    distribution.get(&reslabel1).unwrap_or_else(||{
                        distribution.get(&reslabel2).unwrap_or_else(||{
                            distribution.get(&reslabel3).unwrap_or_else(||{
                                panic!("Cannot find value for {} {} {} {}!",&reslabel0,&reslabel1,&reslabel2,&reslabel3);
                            })
                        })
                    })
                });
                if ii != alen-1{
                    if let Some(x) = diss.next_omega.as_ref(){
                        ret_omega.push(AngleDihed::new(
                            (diss.start_point,
                            diss.divunit),
                            (atom_ca[ii].atom_index as usize
                                ,atom_c[ii].atom_index as usize
                                ,atom_n[ii+1].atom_index as usize
                                ,atom_ca[ii+1].atom_index as usize),
                            x.clone(),
                            if fix_omega{Some(PI)}else{None}
                        ));
                    }
                }
                if ii != 0 {
                    if let Some(x) = diss.prev_omega.as_ref(){
                        ret_omega.push(
                            AngleDihed::new(
                                (diss.start_point,
                                diss.divunit),
                               (atom_ca[ii-1].atom_index as usize
                                    ,atom_c[ii-1].atom_index as usize
                                    ,atom_n[ii].atom_index as usize
                                    ,atom_ca[ii].atom_index as usize)
                                    ,x.clone()
                                    ,if fix_omega{Some(PI)}else{None}
                                )
                        );
                    }
                }
                if ii != 0 && ii != alen-1{
                    if let Some(x) = diss.phi_psi.as_ref(){
                        ret_phi_psi.push(
                            AnglePhiPsi::new(
                                (diss.start_point,
                                diss.divunit),
                                (atom_c[ii-1].atom_index as usize
                                    ,atom_n[ii].atom_index as usize
                                    ,atom_ca[ii].atom_index as usize
                                    ,atom_c[ii].atom_index as usize
                                    ,atom_n[ii+1].atom_index as usize
                                ),
                                x.clone(),
                                None
                        ));
                        
                    }
                }
            }
        }
        return (ret_phi_psi,ret_omega);
    }
}


//作成したエネルギーユニットと既にあるユニットで見ている角度が重複している DihedVars を返す
pub fn get_overwrapping_dihed(dihedvec:&Vec<ff_env::DihedralVars>
    ,phi_psi:&Vec<AnglePhiPsi>,omegas:&Vec<AngleDihed>)->Vec<usize>{
    let from_md:i64 = 0;
    let from_phipshi:i64 = 1;
    let from_omega:i64 = 2;
    //0 もしくは 3 の原子のインデクス、どのベクトルから来たか、ベクトル上のインデクス
    let mut pos_firstatom:Vec<(usize,i64,usize)> = vec![];
    for (ii,aa) in dihedvec.iter().enumerate(){
        pos_firstatom.push((aa.atoms.0.min(aa.atoms.3),from_md,ii));
    }
    for (ii,aa) in phi_psi.iter().enumerate(){
        pos_firstatom.push((aa.atoms.0.min(aa.atoms.3),from_phipshi,ii));
        pos_firstatom.push((aa.atoms.1.min(aa.atoms.4),from_phipshi,ii));
    }
    for (ii,aa) in omegas.iter().enumerate(){
        pos_firstatom.push((aa.atoms.0.min(aa.atoms.3),from_omega,ii));
    }

    //全部を一つのベクトルに突っ込んでソートする
    pos_firstatom.sort_by(|a,b|{
        match a.0.cmp(&b.0){
            Ordering::Equal=>{
                a.1.cmp(&b.1)
            },
            _ => {
                a.0.cmp(&b.0)
            }
        }
    });
    
    let plen:usize = pos_firstatom.len();
    let mut pii:usize = 0;
    let mut marked:Vec<usize> = vec![];

    while pii < plen{
        //比較対象が MD 由来のものでない場合スキップ
        if pos_firstatom[pii].1 != from_md{
            pii += 1;
            continue;
        }
        
        let forr:(usize,usize,usize,usize) = dihedvec[pos_firstatom[pii].2].atoms.clone();
        let rev:(usize,usize,usize,usize) = (forr.3,forr.2,forr.1,forr.0);

        for pjj in (pii+1)..plen{
            //ソートしてるので値が変わった場合それ以降で等しくなることはないので Break
            if pos_firstatom[pjj].0 != pos_firstatom[pii].0{
                break;
            }
            //比較対象が MD 由来の場合スキップ
            if pos_firstatom[pjj].1 == from_md{
                continue;
            }
            if pos_firstatom[pjj].1 == from_phipshi{
                let pidd:usize = pos_firstatom[pjj].2;
                let pp0:(usize,usize,usize,usize) = (phi_psi[pidd].atoms.0,phi_psi[pidd].atoms.1,phi_psi[pidd].atoms.2,phi_psi[pidd].atoms.3);
                let pp1:(usize,usize,usize,usize) = (phi_psi[pidd].atoms.1,phi_psi[pidd].atoms.2,phi_psi[pidd].atoms.3,phi_psi[pidd].atoms.4);
                if pp0 == forr
                || pp0 == rev
                || pp1 == forr
                || pp1 == rev{
                    marked.push(pos_firstatom[pii].2);
                }
            }
            if pos_firstatom[pjj].1 == from_omega{
                let pidd:usize = pos_firstatom[pjj].2;
                if omegas[pidd].atoms == forr
                || omegas[pidd].atoms == rev{
                    marked.push(pos_firstatom[pii].2);
                }
            }
        }
        pii += 1;
    }
    return marked;
}


#[test]
fn vectest(){
    let a = vec![1,2,3,4];
    let b = vec![1,2,3,4];
    let c = vec![1,2,3,4,4];
    let mut d = vec![2,3,4,1];
    if a == b{
        println!("1");
    }
    if a == c{
        println!("2");
    }
    if a == d{
        println!("3");
    }
    let z = d.remove(3);
    d.insert(0,z);
    
    if &a == &d{
        println!("4");
    }
}

#[test]
fn divunittest(){
    for ff in vec![20.0,10.0,5.0,1.0,0.5].into_iter(){
        let dihed:Vec<(f64,bool)> = vec![(0.0,true);(360.0_f64/ff).round() as usize];
        let chk =  AngleDihed::new(
            (-180.0,ff),(0,0,0,0),dihed,None
        );
        let mut count:Vec<usize> = vec![0;(360.0/chk.divunit).round() as usize];
        for ii in 0..3600{
            let pval:f64 = (((ii as f64)*0.1 - 180.0)*10.0).round() as i64 as f64 / 10.0;
            //println!("{} {}",ii,chk.get_vec_index(pval));
            count[get_vec_index_angle(&chk,pval)] += 1;
            //他の言語でのビン割り振りが一致するか調べるときにこれを出力する
            //println!("{}\t{}\t{}",ff,pval,chk.get_vec_index(pval));
        }
        let plen:usize = count.len();
        for ii in 1..plen{
            if count[ii] != count[0]{
                println!("{}",ii);
            }
            assert_eq!(count[ii],count[0]);
        }
    }
}

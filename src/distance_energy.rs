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
use super::geometry::Vector3D;
#[allow(unused_imports)]
use super::geometry::Point3D;
#[allow(unused_imports)]
use super::process_3d;
#[allow(unused_imports)]
use super::md_env;

use super::energy_function::EnergyFunction;
use super::energy_function::Binning;
use super::energy_function::ROUNDING_EPSILON;

#[derive(Debug,Clone)]
pub struct AtomDistanceEnergy{
    pub atoms:(usize,usize),
    pub dist0:f64,
    pub score:f64
}


impl EnergyFunction for AtomDistanceEnergy{
    fn calc_energy(&self,mdenv:&md_env::MDEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{  
        let dist:f64 = process_3d::distance(&mdenv.atoms[self.atoms.0].get_xyz(),&mdenv.atoms[self.atoms.1].get_xyz());
        let dsc = self.score*(dist-self.dist0)*(dist-self.dist0)*weight;
        atom_level_energy[self.atoms.0] += dsc/2.0;
        atom_level_energy[self.atoms.1] += dsc/2.0;
        return dsc;
    }
}
impl AtomDistanceEnergy{
    
    pub fn assign_atom_ids(mdenv:&md_env::MDEnv,distdata:Vec<(String,usize,String,usize,AtomDistanceEnergy)>)->Vec<AtomDistanceEnergy>{
        let mut ret:Vec<AtomDistanceEnergy> = vec![];
        let mut cbatoms:HashMap<String,usize> = HashMap::new();//chain;residue_index->index
        let mut cbatoms_nochain:HashMap<String,usize> = HashMap::new();//"";residue_index->index
        for (aii,aa) in mdenv.atoms.iter().enumerate(){
            let cbflag:bool = if aa.residue_name == "GLY"{
                aa.atom_name == "CA"
            }else{
                aa.atom_name == "CB"
            };
            if cbflag {
                let code = aa.chain_name.clone()+";"+aa.residue_index_in_chain.to_string().as_str();
                cbatoms.insert(code,aii);
                cbatoms_nochain.insert(";".to_string()+aa.residue_index_in_chain.to_string().as_str(),aii);
            }
        }
        for mut dd in distdata.into_iter(){
            let code1:String = dd.0.clone()+";"+dd.1.to_string().as_str();
            let code2:String = dd.2.clone()+";"+dd.3.to_string().as_str();
            if dd.0.len() == 0{
                assert!(dd.2.len() == 0);//Chain 名は 2 つとも指定しておくべき
                dd.4.atoms = (*cbatoms_nochain.get(&code1).unwrap_or_else(||panic!("{} was not found.",code1)),
                *cbatoms_nochain.get(&code2).unwrap_or_else(||panic!("{} was not found.",code2)));
            }else{
                dd.4.atoms = (*cbatoms.get(&code1).unwrap_or_else(||panic!("{} was not found.",code1)),
                *cbatoms.get(&code2).unwrap_or_else(||panic!("{} was not found.",code2)));
            }
            ret.push(dd.4);
        }
        return ret;
    }
    pub fn load_name_mapped(filename:&str)->Vec<(String,usize,String,usize,AtomDistanceEnergy)>{
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<(String,usize,String,usize,AtomDistanceEnergy)> = vec![];
        let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
        let mut functioncode:String ="none".to_string();
        for (_lcount,line) in lines.iter().enumerate() {
            if start_with(line,"#"){
                continue;
            }
            let hs = line_to_hash(line);
            if hs.contains_key("function"){
                functioncode = hs.get("function").unwrap().to_string();
            }
        }
        if functioncode != "score*(dist-dist0)**2"{
            panic!("The function must be score*(dist-dist0)**2. But {} was found.",functioncode);
        }
        for (_lcount,line) in lines.iter().enumerate() {
            if start_with(line,"#"){
                continue;
            }
            let hs = line_to_hash(line);
            if hs.contains_key("dist0"){
                let c1:String = hs.get("c1").unwrap_or(&("".to_string())).to_string();
                let c2:String = hs.get("c2").unwrap_or(&("".to_string())).to_string();
                let r1:usize = hs.get("r1").expect("Can not find r1.").parse::<usize>().expect("r1 was not parsed.");
                let r2:usize = hs.get("r2").expect("Can not find r2.").parse::<usize>().expect("r2 was not parsed.");
                let dist0:f64 = hs.get("dist0").expect("Can not find dist0.").parse::<f64>().expect("dist0 was not parsed.");
                let score:f64 = hs.get("score").expect("Can not find score.").parse::<f64>().expect("score was not parsed.");
                //atoms には仮の値が入っているのでアサインできるかどうかは別の関数で調べる事
                ret.push(
                    (c1,r1,c2,r2,
                    AtomDistanceEnergy{
                        atoms:(0,0),
                        dist0:dist0,
                        score:score
                    }
                    )
                );
           }
        }
        return ret;
    }
}


pub fn get_vec_index_dist(bbin:&dyn Binning,val:f64)-> usize{
    let mut v = val;
    v += -1.0*bbin.get_start_point();
    v -= bbin.get_div_unit()/2.0;
    v += ROUNDING_EPSILON;

    let ret:i64 = (v/ bbin.get_div_unit()).round() as i64;
    return ret.max(0).min(bbin.get_num_bins() as i64 -1) as usize;
}   

pub fn get_linear_interpolated_value(bbin:&dyn Binning,val:f64,bins:&Vec<f64>)-> f64{
    let mut v = val;
    v += -1.0*bbin.get_start_point();
    
    let ret:i64 = (v/ bbin.get_div_unit()).floor() as i64;

    let lb:f64 = bbin.get_start_point()+bbin.get_div_unit()*(ret as f64);
    let ub:f64 = bbin.get_start_point()+bbin.get_div_unit()*(ret as f64 +1.0);
    if ret >= bbin.get_num_bins() as i64 -1{
        return bins[bbin.get_num_bins() as usize -1];
    }
    if ret < 0{
        return bins[0];
    }
    let lratio = (val-lb)/bbin.get_div_unit();
    let uratio = (ub-val)/bbin.get_div_unit();

    assert!(lratio >= 0.0);
    assert!(uratio >= 0.0);
    return (lratio*bins[ret as usize+1]+uratio*bins[ret as usize])/(lratio+uratio);
}   


#[derive(Debug,Clone)]
pub struct AtomBinnedDistanceEnergy{
    divunit:f64,
    start_point:f64,
    pub atoms:(usize,usize),   
    energy_bins:Vec<f64>
}

impl Binning for AtomBinnedDistanceEnergy{
    fn get_div_unit(&self)->f64{
        return self.divunit;
    }
    fn get_num_bins(&self)->usize{
        return self.energy_bins.len();
    }
    fn get_start_point(&self)->f64{
        return self.start_point;
    }
}

impl EnergyFunction for AtomBinnedDistanceEnergy{
    fn calc_energy(&self,mdenv:&md_env::MDEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{  
        let val:f64 = mdenv.atoms[self.atoms.0].distance(&mdenv.atoms[self.atoms.1]);
        let dsc = get_linear_interpolated_value(self,val,&self.energy_bins)*weight;
        atom_level_energy[self.atoms.0] += dsc/2.0;
        atom_level_energy[self.atoms.1] += dsc/2.0;
        return dsc;
    }
}
impl AtomBinnedDistanceEnergy{
    pub fn assign_atom_ids(mdenv:&md_env::MDEnv,distdata:Vec<(String,usize,String,usize,AtomBinnedDistanceEnergy)>)->Vec<AtomBinnedDistanceEnergy>{
        let mut ret:Vec<AtomBinnedDistanceEnergy> = vec![];
        let mut cbatoms:HashMap<String,usize> = HashMap::new();//chain;residue_index->index
        let mut cbatoms_nochain:HashMap<String,usize> = HashMap::new();//"";residue_index->index
        for (aii,aa) in mdenv.atoms.iter().enumerate(){
            let cbflag:bool = if aa.residue_name == "GLY"{
                aa.atom_name == "CA"
            }else{
                aa.atom_name == "CB"
            };
            if cbflag {
                let code = aa.chain_name.clone()+";"+aa.residue_index_in_chain.to_string().as_str();
                cbatoms.insert(code,aii);
                cbatoms_nochain.insert(";".to_string()+aa.residue_index_in_chain.to_string().as_str(),aii);
            }
        }
        for mut dd in distdata.into_iter(){
            let code1:String = dd.0.clone()+";"+dd.1.to_string().as_str();
            let code2:String = dd.2.clone()+";"+dd.3.to_string().as_str();
            if dd.0.len() == 0{
                dd.4.atoms = (*cbatoms_nochain.get(&code1).unwrap_or_else(||panic!("{} was not found.",code1)),
                *cbatoms_nochain.get(&code2).unwrap_or_else(||panic!("{} was not found.",code2)));
            }else{
                dd.4.atoms = (*cbatoms.get(&code1).unwrap_or_else(||panic!("{} was not found.",code1)),
                *cbatoms.get(&code2).unwrap_or_else(||panic!("{} was not found.",code2)));
            }
            ret.push(dd.4);
        }
        return ret;
    }
    pub fn load_name_mapped(filename:&str)->Vec<(String,usize,String,usize,AtomBinnedDistanceEnergy)>{
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<(String,usize,String,usize,AtomBinnedDistanceEnergy)> = vec![];
        let lines:Vec<String> = reader.lines().into_iter().map(|m| m.unwrap().to_string()).collect();
        let mut divunit:f64 = -1.0;
        let mut start_point_:Option<f64> = None;
        let mut num_bins_:Option<usize> = None;
        let mut undefvalue_:Option<f64> = None;//指定のない場合にあてられるエネルギー
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
            if hs.contains_key("num_bins"){
                num_bins_ = Some(hs.get("num_bins").unwrap().parse::<usize>().expect("Can not parse num_bins."));
            }
        }
        
        if divunit < 0.0{
            panic!("Can not find divunit. (binsize)");
        }

        let start_point:f64 = start_point_.unwrap_or_else(||panic!("Can not find start_point."));
        let num_bins:usize = num_bins_.unwrap_or_else(||panic!("Can not find num_bins."));
        let undefvalue:f64 = undefvalue_.unwrap_or_else(||panic!("Can not find undefval."));

        for (_lcount,line) in lines.iter().enumerate() {
            if start_with(line,"#"){
                continue;
            }
            let hs = line_to_hash(line);
            if hs.contains_key("r1"){
                let c1:String = hs.get("c1").unwrap_or(&("".to_string())).to_string();
                let c2:String = hs.get("c2").unwrap_or(&("".to_string())).to_string();
                let r1:usize = hs.get("r1").expect("Can not find r1.").parse::<usize>().expect("r1 was not parsed.");
                let r2:usize = hs.get("r2").expect("Can not find r2.").parse::<usize>().expect("r2 was not parsed.");
                
                let val_:Vec<&str> = hs.get("values").expect("Can not find values.").split("#").collect();
                let mut binned = AtomBinnedDistanceEnergy{
                    atoms:(0,0),
                    energy_bins:vec![undefvalue;num_bins],
                    start_point:start_point,
                    divunit:divunit
                };
                for  vv in val_.iter(){
                    if let Some(_) = REGEX_NOLINE.captures(vv){
                        continue;
                    }
                    let gval:Vec<&str> = vv.split(",").collect();
                    let px:usize = gval[0].parse::<usize>().unwrap();
                    let val:f64 = gval[1].parse::<f64>().unwrap();
                    binned.energy_bins[px] = val;
                }
                
                //atoms には仮の値が入っているのでアサインできるかどうかは別の関数で調べる事
                ret.push(
                    (c1,r1,c2,r2,binned
                    )
                );
           }
        }
        return ret;
    }
}

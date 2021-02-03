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
use super::mmcif_process;

#[allow(unused_imports)]
use self::regex::Regex;//module ファイル内で extern crate した場合 self が必要。lib.rs で extern crate すると必要ない。

#[allow(unused_imports)]
use super::geometry;
use super::geometry::Vector3D;
use super::geometry::Point3D;
use super::pdbdata;
use super::energy_function;

#[allow(unused_imports)]
use super::debug_env;

#[allow(unused_imports)]
use super::backbone_sample;

#[allow(unused_imports)]
use super::process_3d;

#[allow(unused_imports)]
use super::misc_util::*;

use super::energy_function::*;

use std::f64::consts::PI;
const DEGREE_TO_RADIAN2:f64 = PI/180.0*PI/180.0;


pub const EPSILON:f64 = 1.0e-20;







/*
pub struct OpenFFEnergy{
    pub bondvec:Vec<BondVars>,
    pub anglevec:Vec<AngleVars>,
    pub dihedvec:Vec<DihedralVars>,
    pub imprvec:Vec<IMPRVars>,
    pub lj_energy_calculator:LJEnergyCalculator,
    pub electrostatic_energy_calculator:ElectrostaticEnergyCalculator,
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
*/
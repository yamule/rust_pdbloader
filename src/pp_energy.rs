extern crate rand;
use super::charmm_based_energy;
use super::geometry::Vector3D;

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
use super::energy_function::EnergyFunction;


use super::evoef2_energy;

//ToDo:new Constructor を作って、set_xxx_energy(obj,weight) みたいな感じでセットさせる。
pub struct PPEnergySet{
    pub evoef2_env:evoef2_energy::EvoEF2Env,
    pub backbone_energy_omega:Vec<peptide_backbone_dihedral_energy::AngleDihed>,
    pub backbone_energy_phi_psi:Vec<peptide_backbone_dihedral_energy::AnglePhiPsi>,
    pub atom_distance_energy:Vec<distance_energy::AtomDistanceEnergy>,
    pub atom_binned_distance_energy:Vec<distance_energy::AtomBinnedDistanceEnergy>,
    pub atom_contact_energy:Vec<AtomContactEnergy>,
    pub weights:PPEnergyWeight
}




#[derive(Debug,Clone)]
pub struct AtomContactEnergy{
    pub atoms_a:Vec<usize>,
    pub atoms_b:Vec<usize>,
    pub energy:f64,
    pub threshold:f64,
    pub linear_penalty:f64,
    //pub square_penalty:f64,
}


impl EnergyFunction for AtomContactEnergy{
    fn calc_energy(&self,mdenv:&charmm_based_energy::CharmmEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64{  
        //ToDo:atom_level_energy 考察
        let mut ret_:f64 = 1000000000.0;//適当
        
        for aa in self.atoms_a.iter(){
            //A group, B group に含まれる Atom が閾値未満の時 energy を与える。
            //より大きいときは self.linear_penalty*(ret_ - self.threshold)+self.energy; を与える。
            for bb in self.atoms_b.iter(){
                ret_  = ret_.min(mdenv.dist[*aa][*bb]);
            }
        }
        let mut ret = self.energy;
        if ret_ > self.threshold{
            ret = self.linear_penalty*(ret_ - self.threshold)+self.energy;
        }
        ret *= weight;

        let anum_:usize = self.atoms_a.len() + self.atoms_b.len();
        if anum_ > 0{
            let aene:f64 = ret/(anum_ as f64);
            for aa in self.atoms_a.iter(){
                atom_level_energy[*aa] += aene;
            }
            for bb in self.atoms_b.iter(){
                atom_level_energy[*bb] += aene;
            }
        }
        return ret;
    }
}


#[derive(Clone)]
pub struct PPEnergyWeight{
    pub charmm_bond:f64,
    pub charmm_angle:f64,
    pub charmm_ub:f64,
    pub charmm_dihed:f64,
    pub charmm_impr:f64,
    pub charmm_lj:f64,
    pub charmm_electro:f64,

    pub dihed_weight_charmm:Vec<f64>,
    pub backbone_energy_omega:f64,
    pub backbone_energy_phi_psi:f64,
    pub atom_distance_energy:f64,
    pub atom_binned_distance_energy:f64,
    pub atom_contact_energy:f64,

    pub evoef2_vdw:Vec<f64>,
    pub evoef2_elec:f64,
    pub evoef2_hb:Vec<f64>,
    pub evoef2_desolv_polar:f64,
    pub evoef2_desolv_nonpolar:f64,
    pub evoef2_ss:f64,
/*
sub_energyset.weights.charmm_bond = 0.0;
sub_energyset.weights.charmm_angle = 0.0;
sub_energyset.weights.charmm_ub = 0.0;
sub_energyset.weights.charmm_dihed = 0.0;
sub_energyset.weights.charmm_impr = 0.0;
sub_energyset.weights.charmm_lj = 0.0;
sub_energyset.weights.charmm_electro = 0.0;

sub_energyset.weights.dihed_weight_charmm=vec![];
sub_energyset.weights.backbone_energy_omega = 0.0;
sub_energyset.weights.backbone_energy_phi_psi = 0.0;
sub_energyset.weights.atom_distance_energy = 0.0;
sub_energyset.weights.atom_binned_distance_energy = 0.0;

sub_energyset.weights.evoef2_vdw=vec![0.0;4];
sub_energyset.weights.evoef2_elec = 0.0;
sub_energyset.weights.evoef2_hb=vec![0.0;9];
sub_energyset.weights.evoef2_desolv_polar = 0.0;
sub_energyset.weights.evoef2_desolv_nonpolar = 0.0;
sub_energyset.weights.evoef2_ss = 0.0;
*/
}
impl PPEnergyWeight{
    pub fn new()->PPEnergyWeight{
        return PPEnergyWeight{
            charmm_bond:1.0,
            charmm_angle:1.0,
            charmm_ub:1.0,
            charmm_dihed:1.0,
            charmm_impr:1.0,
            charmm_lj:1.0,
            charmm_electro:1.0,

            dihed_weight_charmm:vec![],
            backbone_energy_omega:1.0,
            backbone_energy_phi_psi:1.0,
            atom_distance_energy:1.0,
            atom_binned_distance_energy:1.0,
            atom_contact_energy:1.0,

            evoef2_vdw:vec![1.0;4],
            evoef2_elec:1.0,
            evoef2_hb:vec![1.0;9],
            evoef2_desolv_polar:1.0,
            evoef2_desolv_nonpolar:1.0,
            evoef2_ss:1.0,
        };
    }
}
impl PPEnergySet{
    pub fn get_num_atoms(&self)->usize{
        return self.evoef2_env.md_envset.atoms.len();
    }

    pub fn get_atom(&self,ii:usize)->&charmm_based_energy::MDAtom{
        return &self.evoef2_env.md_envset.atoms[ii];
    }
    
    pub fn get_atom_mut(&mut self,ii:usize)->&mut charmm_based_energy::MDAtom{
        return &mut self.evoef2_env.md_envset.atoms[ii];
    }

    pub fn make_set_for_subenv(&self,sparse_atoms:&Vec<usize>)->(PPEnergySet,Vec<i64>){
        let (subenv,mapper):(charmm_based_energy::CharmmEnv,Vec<i64>) = self.evoef2_env.md_envset.make_sub_env(sparse_atoms);
        let charmmvars = self.evoef2_env.charmm_vars.make_sub_env_vars(&mapper);
        let subatomnum:usize = subenv.atoms.len();
        let mut new_evo:evoef2_energy::EvoEF2Env = evoef2_energy::EvoEF2Env{
            md_envset:subenv,
            charmm_vars:charmmvars,
            is_backbone:vec![false;subatomnum],
            solvparam:vec![],
            cys_atoms:vec![]
        };

        for ii in 0..self.evoef2_env.is_backbone.len(){
            if mapper[ii] > -1{
                new_evo.is_backbone[mapper[ii] as usize] = self.evoef2_env.is_backbone[ii];
            }
        }

        let mut sparam:Vec<Option<evoef2_energy::SolvParam>> = vec![None;subatomnum];
        for ii in 0..self.evoef2_env.solvparam.len(){
            if mapper[ii] > -1{
                sparam[mapper[ii] as usize] = Some(self.evoef2_env.solvparam[ii].clone());
            }
        }
        for ss in sparam.into_iter(){
            if let None = ss{
                panic!("???");
            }
            new_evo.solvparam.push(ss.unwrap());
        }

        for cc in self.evoef2_env.cys_atoms.iter(){
            if mapper[cc.ca_index] < 0
            || mapper[cc.cb_index] < 0
            || mapper[cc.sg_index] < 0
            {
                continue;
            }
            new_evo.cys_atoms.push(evoef2_energy::CysAtoms{
                chain_name:cc.chain_name.clone(),
                residue_index_in_chain:cc.residue_index_in_chain,
                ca_index:mapper[cc.ca_index] as usize,
                cb_index:mapper[cc.cb_index] as usize,
                sg_index:mapper[cc.sg_index] as usize
            });

        }
        

        let mut ret = PPEnergySet{
            evoef2_env:new_evo,
            backbone_energy_omega:vec![],
            backbone_energy_phi_psi:vec![],
            atom_distance_energy:vec![],
            atom_binned_distance_energy:vec![],
            atom_contact_energy:vec![],
            weights:self.weights.clone()
        };

        let mut dihedmap:Vec<usize> = vec![];
        for (cii,cc) in self.evoef2_env.charmm_vars.dihedvec.iter().enumerate(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            || mapper[cc.atoms.3] < 0
            {
                continue;
            }
            let mut pcc:charmm_based_energy::DihedralVars = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize,
                mapper[cc.atoms.3] as usize
            );
            dihedmap.push(cii);
        }
        if self.weights.dihed_weight_charmm.len() > 0{
            for dd in dihedmap.into_iter(){
                ret.weights.dihed_weight_charmm.push(self.weights.dihed_weight_charmm[dd]);
            }
        }

        for cc in self.backbone_energy_omega.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            || mapper[cc.atoms.3] < 0
            {
                continue;
            }
            let mut pcc:peptide_backbone_dihedral_energy::AngleDihed = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize,
                mapper[cc.atoms.3] as usize
            );
            ret.backbone_energy_omega.push(pcc);
        }
        
        for cc in self.backbone_energy_phi_psi.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            || mapper[cc.atoms.2] < 0
            || mapper[cc.atoms.3] < 0
            || mapper[cc.atoms.4] < 0
            {
                continue;
            }
            let mut pcc:peptide_backbone_dihedral_energy::AnglePhiPsi = cc.clone();
            pcc.atoms = (
                mapper[cc.atoms.0] as usize,
                mapper[cc.atoms.1] as usize,
                mapper[cc.atoms.2] as usize,
                mapper[cc.atoms.3] as usize,
                mapper[cc.atoms.4] as usize
            );

            ret.backbone_energy_phi_psi.push(pcc);
        }

        for cc in self.atom_distance_energy.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            {
                continue;
            }
            let mut pcc:distance_energy::AtomDistanceEnergy = cc.clone();
            pcc.atoms = (mapper[cc.atoms.0] as usize,mapper[cc.atoms.1] as usize);
            ret.atom_distance_energy.push(pcc);
        }
        for cc in self.atom_binned_distance_energy.iter(){
            if mapper[cc.atoms.0] < 0
            || mapper[cc.atoms.1] < 0
            {
                continue;
            }
            let mut pcc:distance_energy::AtomBinnedDistanceEnergy = cc.clone();
            pcc.atoms = (mapper[cc.atoms.0] as usize,mapper[cc.atoms.1] as usize);
            ret.atom_binned_distance_energy.push(pcc);
        }
        
        for cc in self.atom_contact_energy.iter(){
            let mut avec:Vec<usize> = vec![];
            let mut bvec:Vec<usize> = vec![];
            for caa in cc.atoms_a.iter(){
                if mapper[*caa] < 0{
                    continue;
                }
                avec.push(mapper[*caa] as usize);
            }
            for cbb in cc.atoms_b.iter(){
                if mapper[*cbb] < 0{
                    continue;
                }
                bvec.push(mapper[*cbb] as usize);
            }
            if avec.len() == 0 || bvec.len() == 0{
                continue;
            }
            let mut pcc:AtomContactEnergy = cc.clone();
            pcc.atoms_a = avec;
            pcc.atoms_b = bvec;
            ret.atom_contact_energy.push(pcc);
        }

        ret.update_edges(5);
        return (ret,mapper);
    }

    pub fn update_distance(&mut self){
        self.evoef2_env.md_envset.update_distance();
    }

    pub fn update_distance_one(&mut self,i:usize){
        self.evoef2_env.md_envset.update_distance_one(i);
    }

    pub fn update_edges(&mut self,conn:u64){
        self.evoef2_env.md_envset.update_edges(&self.evoef2_env.charmm_vars.bondvec,conn);
    }

    pub fn gen_pseudo_edges(&mut self,conn:u64){
        self.evoef2_env.md_envset.gen_pseudo_edges(conn);
    }

    pub fn calc_energy(&self,atom_level_energy:&mut Vec<f64>)->f64{
        let mut ret:f64 = 0.0;
        
        if self.weights.charmm_bond > 0.0{
            for bb in self.evoef2_env.charmm_vars.bondvec.iter(){
                ret += bb.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_bond);
            }
        }
        
        if self.weights.charmm_angle > 0.0{
            for aa in self.evoef2_env.charmm_vars.anglevec.iter(){
                ret += aa.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_angle);
            }
        }
        
        if self.weights.charmm_ub > 0.0{
            for uu in self.evoef2_env.charmm_vars.ubvec.iter(){
                ret += uu.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_ub);
            }
        }
        
        if self.weights.charmm_dihed > 0.0{
            if self.weights.dihed_weight_charmm.len() > 0{
                for (dii,dd) in self.evoef2_env.charmm_vars.dihedvec.iter().enumerate(){
                    ret += dd.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_dihed*self.weights.dihed_weight_charmm[dii]);
                }
            }else{
                for dd in self.evoef2_env.charmm_vars.dihedvec.iter(){
                    ret += dd.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_dihed);
                }
            }
        }

        if self.weights.charmm_impr > 0.0{
           for ii in self.evoef2_env.charmm_vars.imprvec.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_impr);
            }
        }

        if self.weights.charmm_lj > 0.0{
            ret += self.evoef2_env.charmm_vars.lj_energy_calculator.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_lj);
        }

        if self.weights.charmm_electro >  0.0{
            ret += self.evoef2_env.charmm_vars.electrostatic_energy_calculator.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_electro);
        }

        if self.weights.backbone_energy_omega > 0.0{
            for ii in self.backbone_energy_omega.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.backbone_energy_omega);
            }
        }

        if self.weights.backbone_energy_phi_psi > 0.0{
            for ii in self.backbone_energy_phi_psi.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.backbone_energy_phi_psi);
            }
        }
        
        if self.weights.atom_distance_energy > 0.0{
            for ii in self.atom_distance_energy.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_distance_energy);
            }
        }

        if self.weights.atom_binned_distance_energy > 0.0{
            for ii in self.atom_binned_distance_energy.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_binned_distance_energy);
            }
        }

        if self.weights.atom_contact_energy > 0.0{
            for ii in self.atom_contact_energy.iter(){
                ret += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_contact_energy);
            }
        }

        if self.weights.evoef2_vdw.iter().fold(0.0,|s,m| s+m) > 0.0{
            ret += evoef2_energy::calcEvdw(&self.evoef2_env,atom_level_energy, &self.weights.evoef2_vdw).iter().fold(0.0,|s,m| s+m);
        }

        if self.weights.evoef2_elec > 0.0{
            ret += evoef2_energy::calcEelec(&self.evoef2_env,atom_level_energy, self.weights.evoef2_elec);
        }

        if self.weights.evoef2_hb.iter().fold(0.0,|s,m| s+m) > 0.0{
            ret += evoef2_energy::calcEhb(&self.evoef2_env,atom_level_energy, &self.weights.evoef2_hb).iter().fold(0.0,|s,m| s+m);
        }

        if self.weights.evoef2_desolv_polar+self.weights.evoef2_desolv_nonpolar > 0.0{
            let deso  =evoef2_energy::calcEdesolv(&self.evoef2_env,atom_level_energy, self.weights.evoef2_desolv_polar, self.weights.evoef2_desolv_nonpolar);
            ret += deso.0;
            ret += deso.1;
        }

        if self.weights.evoef2_ss > 0.0{
            ret +=  evoef2_energy::calcEss(&self.evoef2_env,atom_level_energy, self.weights.evoef2_ss);
        }
        if ret.is_nan(){
            println!("nan found!");
            let mut ene = vec![0.0;atom_level_energy.len()];
            let res = self.calc_energy_sep(&mut ene);
            let mut lastidx:usize = 0;
            for (aii,aa) in self.evoef2_env.md_envset.atoms.iter().enumerate(){
                if ene[aii].is_nan(){
                    lastidx = aii;
                    println!("{} atom:{} {:?}?",atom_level_energy[aii],aa.get_line_representation(),aa.get_xyz());
                }
            }
            println!("{:?}",self.evoef2_env.md_envset.dist[lastidx]);
            println!("res:{:?}",res);
            panic!();
        }

        return ret;
    }
    
    pub fn calc_energy_sep(&self,atom_level_energy:&mut Vec<f64>)->(f64,Vec<f64>){
        let mut ret:Vec<f64> = vec![0.0;17];
        let mut eindex:usize = 0;
        
        if self.weights.charmm_bond > 0.0{
            for bb in self.evoef2_env.charmm_vars.bondvec.iter(){
                ret[eindex] += bb.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_bond);
            }
        }
        eindex += 1;
        
        if self.weights.charmm_angle > 0.0{
            for aa in self.evoef2_env.charmm_vars.anglevec.iter(){
                ret[eindex] += aa.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_angle);
            }
        }
        eindex += 1;
        
        if self.weights.charmm_ub > 0.0{
            for uu in self.evoef2_env.charmm_vars.ubvec.iter(){
                ret[eindex] += uu.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_ub);
            }
        }
        eindex += 1;
        
        if self.weights.charmm_dihed > 0.0{
            if self.weights.dihed_weight_charmm.len() > 0{
                for (dii,dd) in self.evoef2_env.charmm_vars.dihedvec.iter().enumerate(){
                    ret[eindex] += dd.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_dihed*self.weights.dihed_weight_charmm[dii]);
                }
            }else{
                for dd in self.evoef2_env.charmm_vars.dihedvec.iter(){
                    ret[eindex] += dd.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_dihed);
                }
            }
        }
        eindex += 1;

        if self.weights.charmm_impr > 0.0{
           for ii in self.evoef2_env.charmm_vars.imprvec.iter(){
            ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_impr);
            }
        }
        eindex += 1;

        if self.weights.charmm_lj > 0.0{
            ret[eindex] += self.evoef2_env.charmm_vars.lj_energy_calculator.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_lj);
        }
        eindex += 1;

        if self.weights.charmm_electro >  0.0{
            ret[eindex] += self.evoef2_env.charmm_vars.electrostatic_energy_calculator.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.charmm_electro);
        }
        eindex += 1;

        if self.weights.backbone_energy_omega > 0.0{
            for ii in self.backbone_energy_omega.iter(){
                ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.backbone_energy_omega);
            }
        }
        eindex += 1;

        if self.weights.backbone_energy_phi_psi > 0.0{
            for ii in self.backbone_energy_phi_psi.iter(){
                ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.backbone_energy_phi_psi);
            }
        }
        eindex += 1;
        
        if self.weights.atom_distance_energy > 0.0{
            for ii in self.atom_distance_energy.iter(){
                ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_distance_energy);
            }
        }
        eindex += 1;
        

        if self.weights.atom_binned_distance_energy > 0.0{
            for ii in self.atom_binned_distance_energy.iter(){
                ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_binned_distance_energy);
            }
        }
        eindex += 1;

        if self.weights.atom_contact_energy > 0.0{
            for ii in self.atom_contact_energy.iter(){
                ret[eindex] += ii.calc_energy(&self.evoef2_env.md_envset,atom_level_energy,self.weights.atom_contact_energy);
            }
        }
        eindex += 1;

        if self.weights.evoef2_vdw.iter().fold(0.0,|s,m| s+m) > 0.0{
            ret[eindex] += evoef2_energy::calcEvdw(&self.evoef2_env,atom_level_energy, &self.weights.evoef2_vdw).iter().fold(0.0,|s,m| s+m);
        }
        eindex += 1;

        if self.weights.evoef2_elec > 0.0{
            ret[eindex] += evoef2_energy::calcEelec(&self.evoef2_env,atom_level_energy, self.weights.evoef2_elec);
        }
        eindex += 1;

        if self.weights.evoef2_hb.iter().fold(0.0,|s,m| s+m) > 0.0{
            ret[eindex] += evoef2_energy::calcEhb(&self.evoef2_env,atom_level_energy, &self.weights.evoef2_hb).iter().fold(0.0,|s,m| s+m);
        }
        eindex += 1;

        if self.weights.evoef2_desolv_polar+self.weights.evoef2_desolv_nonpolar > 0.0{
            let deso  =evoef2_energy::calcEdesolv(&self.evoef2_env,atom_level_energy, self.weights.evoef2_desolv_polar, self.weights.evoef2_desolv_nonpolar);
            ret[eindex] += deso.0;
            ret[eindex] += deso.1;
        }
        eindex += 1;

        if self.weights.evoef2_ss > 0.0{
            ret[eindex] +=  evoef2_energy::calcEss(&self.evoef2_env,atom_level_energy, self.weights.evoef2_ss);
        }

        return (ret.iter().fold(0.0,|s,m|s+m),ret);
    }

}

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


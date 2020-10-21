
use std::fs::File;
#[allow(unused_imports)]
use std::collections::HashSet;
use std::collections::HashMap;
#[allow(unused_imports)]
use std::fs;
#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use super::charmm_based_energy;
use super::geometry::Vector3D;
use std::f64::consts::PI;
use super::misc_util::*;

const SCALE_FACTOR_1_4:f64 = 0.2;

pub const HB_BB_DIST:usize = 0;
pub const HB_BB_THETA:usize = 1;
pub const HB_BB_PHI:usize = 2;
pub const HB_SB_DIST:usize = 3;
pub const HB_SB_THETA:usize = 4;
pub const HB_SB_PHI:usize = 5;
pub const HB_SS_DIST:usize = 6;
pub const HB_SS_THETA:usize = 7;
pub const HB_SS_PHI:usize = 8;

pub const VDW_INTER_ATT:usize = 0;
pub const VDW_INTRA_ATT:usize = 1;
pub const VDW_INTER_REP:usize = 2;
pub const VDW_INTRA_REP:usize = 3;

enum PairType{
    BB,SB,SS
}

lazy_static!{
    static ref PI_3_2:f64 = PI.powf(1.5);//3/2 乗を使う場合があるのであらかじめ計算しておく
}


/*
https://www.charmm.org/ubbthreads/ubbthreads.php?ubb=showflat&Number=33804
によれば　CHARMM FORCEFIELD は PUBLIC DOMAIN のようだ。
修正があったらフォローできるように CHARMM のフォーマットを読めるようにしておいてほしいとのこと。
*/
/*
Huang, Xiaoqiang, Robin Pearce, and Yang Zhang.
"EvoEF2: accurate and fast energy function for computational protein design."
Bioinformatics (2019).
*/

pub struct EvoEF2Env{
   pub md_envset:charmm_based_energy::CharmmEnv,
   pub charmm_vars:charmm_based_energy::CharmmVars,
   pub is_backbone:Vec<bool>,
   pub solvparam:Vec<SolvParam>,
   pub cys_atoms:Vec<CysAtoms>,
}

#[allow(dead_code)]
pub struct CysAtoms{
    pub chain_name:String,
    pub residue_index_in_chain:usize,
    pub ca_index:usize,
    pub cb_index:usize,
    pub sg_index:usize,
}

#[derive(Clone)]
pub struct SolvParam{
    volume:f64,
    dgfree:f64,
    correllen:f64,
    polar:bool,
}

impl EvoEF2Env{
    pub fn new(
        mdd:charmm_based_energy::CharmmEnv
        ,mdv:charmm_based_energy::CharmmVars
        ,resourcedir:&str
        ,skip_ss:bool//SS BOND の計算はしない。SPARSE の時
    )->EvoEF2Env{
        let sp2_atoms:HashSet<(String,String)> = EvoEF2Env::load_sp2_atoms(&(resourcedir.to_string()+"/sp2_atoms.dat"));
        let mut cys_atoms:Vec<CysAtoms> = vec![];
        if !skip_ss{
            let mut cys_atoms_:HashMap<String,(i64,i64,i64)> = HashMap::new();
            //ca/cb/sg
            for (aii,aa) in mdd.atoms.iter().enumerate(){
                if aa.residue_name == "CYS"{
                    let code:String = format!("{}#{}",aa.chain_name,aa.residue_index_in_chain);
                    if ! cys_atoms_.contains_key(&code){
                        cys_atoms_.insert(code.clone(),(-1,-1,-1));
                    }
                    if aa.atom_name == "SG" {
                        cys_atoms_.get_mut(&code).unwrap().2 = aii as i64;
                    }
                    if aa.atom_name == "CA" {
                        cys_atoms_.get_mut(&code).unwrap().0 = aii as i64;
                    }
                    if aa.atom_name == "CB" {
                        cys_atoms_.get_mut(&code).unwrap().1 = aii as i64;
                    }
                }
            }
            let mut cyss:Vec<String> = cys_atoms_.iter().map(|m|m.0.clone()).collect(); 
            cyss.sort();
            for cc in cyss.iter(){
                let cy:&(i64,i64,i64) = cys_atoms_.get(cc).unwrap();
                let ppt: Vec<&str> = cc.split("#").collect();
                let iid:i64 = ppt[1].parse::<i64>().unwrap_or_else(|_|panic!("{} must be <String>#<i64>,",cc));
                if iid < 0{
                    panic!("Residue index must be positive but {}. Did you forget to assing at first?",iid);
                }
                if cy.0 < 0{
                        panic!("CA index is {}. The atom is missing or using different atom name?",cy.0);
                }
                if cy.1 < 0{
                    panic!("CB index is {}. The atom is missing or using different atom name?",cy.1);
                }
                if cy.2 < 0{
                    panic!("SG index is {}. The atom is missing or using different atom name?",cy.2);
                }
                cys_atoms.push(
                    CysAtoms{
                        chain_name:ppt[0].to_string(),
                        residue_index_in_chain:iid as usize,
                        ca_index:cy.0 as usize,
                        cb_index:cy.1 as usize,
                        sg_index:cy.2 as usize,
                    }

                );
            }
        }

        let mut ret:EvoEF2Env = EvoEF2Env{
            md_envset:mdd,
            charmm_vars:mdv,
            is_backbone:vec![],
            solvparam:vec![],
            cys_atoms:cys_atoms,
        };
        ret.assign_sp2_charmm19(&sp2_atoms);
        ret.flag_backbone();
        ret.assign_deltag_free();
        
        //for aa in ret.mdenv.atoms.iter_mut(){
        //    let rpar = r_mapper.get_rosetta_param(aa);
        //    aa.nb_epsilon = rpar.lj_welldepth;
        //    aa.nb_r1_2 = rpar.lj_radii;
        //}
        return ret;
    }
    pub fn flag_backbone(&mut self){
        //charmm 22
        //let batoms:HashSet<String> = (vec!["N","HN","CA","HA1","HA","C","O"]).iter().map(|m|m.to_string()).collect();
        let batoms:HashSet<String> = (vec!["N","H","CA","C","O"]).iter().map(|m|m.to_string()).collect();
        self.is_backbone = vec![false;self.md_envset.atoms.len()];
        for (aii,aa) in self.md_envset.atoms.iter().enumerate(){
            if batoms.contains(&aa.atom_name){
                self.is_backbone[aii] = true;
            }
        }
    }
    
    pub fn assign_sp2_charmm19(&mut self,sp2atoms:&HashSet<(String,String)>){
        for aa in self.md_envset.atoms.iter_mut(){
            let code = (aa.residue_name.clone(),aa.atom_name.clone());
            if sp2atoms.contains(&code){
                aa.hybrid_orbit = charmm_based_energy::HybridOrbit::SP2;
            }
            if aa.nterminal{
                if sp2atoms.contains(&("NTER".to_string(),aa.atom_name.clone())){
                    aa.hybrid_orbit = charmm_based_energy::HybridOrbit::SP2;
                }   
            }
            
            if aa.cterminal{
                if sp2atoms.contains(&("CTER".to_string(),aa.atom_name.clone())){
                    aa.hybrid_orbit = charmm_based_energy::HybridOrbit::SP2;
                }   
            }
        }
    }

    //residuename<whitespace>atomname<whitespace>SP2
    //という形式になった TSV (っぽいもの。whitespace でもよい)を渡す
    pub fn load_sp2_atoms(filename:&str)->HashSet<(String,String)>{
        let mut ret:HashSet<(String,String)> = HashSet::new();
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        for (_lcount,line_) in reader.lines().enumerate() {
            let line_ = line_.unwrap();
            let line =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();

            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }
            
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt[2] == "SP2"{
                ret.insert((ptt[0].clone(),ptt[1].clone()));
            }
        }
        return ret;
    }
    pub fn assign_deltag_free(&mut self){
        /*
        PROTEINS: Structure, Function, and Genetics 35:133–152 (1999)
        Lazaridis and Karplus - Effective energy function for proteins in solution
        */
        let volume_dgfree_:Vec<(&str,f64,f64)> = vec![
            ("C",14.7,0.00)
            ,("CR",8.3,-1.40)
            ,("CH1E",23.7,-0.25)
            ,("CH2E",22.4,0.52)
            ,("CH3E",30.0,1.50)
            ,("CR1E",18.4,0.08)
            ,("NH1",4.4,-8.90)
            ,("NR",4.4,-4.00)
            ,("NH2",11.2,-7.80)
            ,("NH3",11.2,-20.00)
            ,("NC2",11.2,-10.00)
            ,("N",0.0,-1.55)
            ,("OH1",10.8,-6.70)
            ,("O",10.8,-5.85)
            ,("OC",10.8,-10.00)
            ,("S",14.7,-4.10)
            ,("SH1E",21.4,-2.70)
        ];

        let mut v_dgfree:HashMap<String,(f64,f64)> = HashMap::new();
        for dd in volume_dgfree_.into_iter(){
            v_dgfree.insert(dd.0.to_string(),(dd.1,dd.2));
        }        

        self.solvparam = vec![SolvParam{
            volume:0.0,
            dgfree:0.0,
            correllen:0.0,
            polar:false
        };self.md_envset.atoms.len()];

        for (aii,aa) in self.md_envset.atoms.iter().enumerate(){
            if v_dgfree.contains_key(&aa.atom_type){
                let vdgf:(f64,f64) = *v_dgfree.get(&aa.atom_type).unwrap();
                self.solvparam[aii].volume = vdgf.0;
                self.solvparam[aii].dgfree = vdgf.1;
            }
            let mut correl_len:f64 = 3.5;
            //Lazaridis and Karplus 1999 TableII から
            //Partial Charge 0.0 も 6.0 にするんだろうか？
            //EvoEF2 には S、C、は nonpolar と書かれているが・・・？
            if aa.residue_name == "ARG"{
                if aa.atom_name == "CD" 
                || aa.atom_name == "NE" 
                || aa.atom_name == "HE" 
                || aa.atom_name == "CZ" 
                || aa.atom_name == "NH1"
                || aa.atom_name == "NH2"
                || aa.atom_name == "HH11" 
                || aa.atom_name == "HH12" 
                || aa.atom_name == "HH21" 
                || aa.atom_name == "HH22" 
                {
                    correl_len = 6.0;
                }
            }else if aa.residue_name == "LYS"{
                if aa.atom_name == "CE"//Partial charge 0
                 || aa.atom_name == "NZ"
                 || aa.atom_name == "HZ1" 
                 || aa.atom_name == "HZ2" 
                 || aa.atom_name == "HZ3" {
                    correl_len = 6.0;
                }
            }else if aa.residue_name == "ASP"{
                if aa.atom_name == "CB"
                 || aa.atom_name == "CG"
                 || aa.atom_name == "OD1" 
                 || aa.atom_name == "OD2" 
                 {
                    correl_len = 6.0;
                }
            }else if aa.residue_name == "GLU"{
                if aa.atom_name == "CG"
                 || aa.atom_name == "CD"
                 || aa.atom_name == "OE1" 
                 || aa.atom_name == "OE2" 
                 {
                    correl_len = 6.0;
                }
            }
            if aa.nterminal{
                if aa.atom_name == "N"
                || aa.atom_name == "HT1"
                || aa.atom_name == "HT2" 
                || aa.atom_name == "HT3" 
                {
                    correl_len = 6.0;
                }
            }
            if aa.cterminal{
                if aa.atom_name == "C"
                || aa.atom_name == "OT1"
                || aa.atom_name == "OT2" 
                {
                    correl_len = 6.0;
                }
            }
            self.solvparam[aii].correllen = correl_len;
            if start_with(&aa.atom_name,"C") || start_with(&aa.atom_name,"S"){
                self.solvparam[aii].polar = false;
            }else if start_with(&aa.atom_name,"N") || start_with(&aa.atom_name,"O"){
                self.solvparam[aii].polar = true;
            }else{
                if !start_with(&aa.atom_name,"H") {
                    panic!("Solvation parameter for {} is not defined.",aa.atom_name);
                }
            }
        }
    }
    pub fn is_sp2(atom:&charmm_based_energy::MDAtom)->bool{
        match atom.hybrid_orbit{
            charmm_based_energy::HybridOrbit::SP2=>{return true;},
            _=>{return false;}
        }
    }
}

#[allow(non_snake_case)]
pub fn calcEvdw(evoenv:&EvoEF2Env,atom_level_energy:&mut Vec<f64>,weights:&Vec<f64>)
->Vec<f64>{
    let mdenv = &evoenv.md_envset;
    let num_atoms:usize = mdenv.atoms.len();
    let mut _ulj = 0.0;//善合計

    let mut intra_rep = 0.0;
    let mut intra_att = 0.0;
    let mut inter_rep = 0.0;
    let mut inter_att = 0.0;
    for aa in 0..num_atoms-1{
        for bb in (aa+1)..num_atoms{
            if mdenv.atoms[aa].chain_name == mdenv.atoms[bb].chain_name{
                if evoenv.is_backbone[aa] && evoenv.is_backbone[bb]
                //    && (mdenv.atoms[aa].residue_index_in_chain - mdenv.atoms[bb].residue_index_in_chain).abs() < 2{
                    && mdenv.atoms[aa].residue_index_in_chain == mdenv.atoms[bb].residue_index_in_chain{
                        continue;
                }
            }

            let attractive_weight:f64;
            let repulsive_weight:f64;
            let mut in_same_residue:bool = false;
            //inter residue, intra residue の分離
            if mdenv.atoms[aa].chain_name != mdenv.atoms[bb].chain_name{
                attractive_weight = 1.06;
                repulsive_weight = 0.80;
            }else{
                if mdenv.atoms[bb].residue_index_in_chain == mdenv.atoms[aa].residue_index_in_chain{
                    attractive_weight = 0.43;
                    repulsive_weight = 0.06;
                    in_same_residue = true;
                }else{
                    attractive_weight = 1.21;
                    repulsive_weight = 1.28;
                }
            }

            if mdenv.num_edges[aa][bb] < 3{
                //近くの原子はスキップ
            }else{
                let scale_factor:f64 = if mdenv.num_edges[aa][bb] <= 3{
                    SCALE_FACTOR_1_4
                }else{
                    1.0
                };

                let ddis = mdenv.dist[aa][bb]+charmm_based_energy::EPSILON;
                
                if ddis >= 6.0{
                    continue;
                }

                let epss:f64;
                let rminij:f64;
                epss = (mdenv.atoms[aa].nb_epsilon
                    *mdenv.atoms[bb].nb_epsilon).sqrt();
                rminij =  mdenv.atoms[aa].nb_r1_2+mdenv.atoms[bb].nb_r1_2;
                
                //rminij *= 0.95;

                let mut dsc:f64 = 0.0;

                if ddis < 0.8909*rminij{
                    let bzz = (rminij/ddis).powf(6.0);
                    
                    //ソースには別のブロックがあり、5.0 で切らない方がよかったように書かれている。
                    dsc =epss*(bzz.powf(2.0) - 2.0*bzz);
                    if dsc > 5.0*epss{
                        dsc = 5.0*epss+(0.8909*rminij-ddis)*(0.8909*rminij-ddis);//EvoEF2 では 5.0*epss だが、距離に応じた反発がなくなるので変更。ToDo どこかの文献からそれらしき補正項を見つけてくる
                    }
                    dsc *= repulsive_weight*scale_factor;
                    
                    if in_same_residue{
                        dsc = dsc*weights[VDW_INTRA_REP];
                        intra_rep += dsc;
                    }else{
                        dsc = dsc*weights[VDW_INTER_REP];
                        inter_rep += dsc;
                    }
                }else if ddis < 5.0{
                    let bzz = (rminij/ddis).powf(6.0);
                    dsc = epss*(bzz.powf(2.0) - 2.0*bzz);
                    dsc *= attractive_weight*scale_factor;
                    
                    if in_same_residue{
                        dsc = dsc*weights[VDW_INTRA_ATT];
                        intra_att += dsc;
                    }else{
                        dsc = dsc*weights[VDW_INTER_ATT];
                        inter_att += dsc;
                    }
                }else if ddis < 6.0{
                    let bzz = (rminij/5.0).powf(6.0);
                    let bzz12 = bzz.powf(2.0);
                    let a:f64 =  -0.4*epss*bzz12   -1.6*epss*bzz;
                    let b:f64 =   7.8*epss*bzz12  +25.2*epss*bzz;
                    let c:f64 = -50.4*epss*bzz12 -129.6*epss*bzz;
                    let d:f64 = 108.0*epss*bzz12 +216.0*epss*bzz;
                    dsc = a*ddis*ddis*ddis+b*ddis*ddis+c*ddis+d;
                    dsc *= attractive_weight*scale_factor;
                    
                    if in_same_residue{
                        intra_att += dsc;
                    }else{
                        inter_att += dsc;
                    }
                }

                atom_level_energy[aa] += dsc/2.0;
                atom_level_energy[bb] += dsc/2.0;
                
                _ulj += dsc;
            }
        }
    }
    let mut ret = vec![0.0;4];
    ret[VDW_INTRA_ATT] = intra_att;
    ret[VDW_INTER_ATT] = inter_att;
    ret[VDW_INTRA_REP] = intra_rep;
    ret[VDW_INTER_REP] = inter_rep;
    return ret;
}
#[allow(non_snake_case)]
pub fn calcEelec(evoenv:&EvoEF2Env,atom_level_energy:&mut Vec<f64>,weight:f64)
->f64{
    let mdenv = &evoenv.md_envset; 
    let num_atoms:usize = mdenv.atoms.len();
    let mut uelec = 0.0;
    for aa in 0..num_atoms-1{
        for bb in (aa+1)..num_atoms{
            
            if evoenv.is_backbone[aa] && evoenv.is_backbone[bb]
             && (mdenv.atoms[aa].residue_index_in_chain - mdenv.atoms[aa].residue_index_in_chain).abs() < 2 {
                continue;
            }

            if mdenv.num_edges[aa][bb] < 3{
                //近くの原子はスキップ
            }else{
                
                let scale_factor:f64;
                if mdenv.num_edges[aa][bb] == 3{
                    scale_factor = SCALE_FACTOR_1_4;
                }else{
                    scale_factor = 1.0;
                }

                let rminij =  mdenv.atoms[aa].nb_r1_2+mdenv.atoms[bb].nb_r1_2;
            
                let ddis = mdenv.dist[aa][bb]+charmm_based_energy::EPSILON;
                if ddis >= 6.0{
                    continue;
                }
                let mut elec:f64 = 0.0;
                if ddis < 0.8*rminij{
                   elec = 332.0*mdenv.atoms[aa].charge*mdenv.atoms[bb].charge/(0.8*rminij);
                }else if ddis < 6.0{
                   elec = 332.0*mdenv.atoms[aa].charge*mdenv.atoms[bb].charge/ddis/ddis/40.0;
                }
                elec *= scale_factor;
                elec *= weight;
                atom_level_energy[aa] += elec/2.0;
                atom_level_energy[bb] += elec/2.0;
                uelec += elec;
            }
        }
    }
    return uelec;
}
#[allow(non_snake_case)]
pub fn calcEhb(evoenv:&EvoEF2Env
    ,atom_level_energy:&mut Vec<f64>
    ,weights:&Vec<f64>)
->Vec<f64>{
    let dmin = 1.4;
    let dmax = 3.0;
    let mut ret = vec![0.0;9];

    for aa in evoenv.charmm_vars.hb_acceptors.iter(){
        for dd in evoenv.charmm_vars.hb_donors.iter(){
            if evoenv.md_envset.num_edges[aa.0][dd.0] < 3{
            }else{
                let mut pairtype:PairType =  PairType::SS;
                if evoenv.is_backbone[aa.0]
                && evoenv.is_backbone[dd.1]{
                    pairtype = PairType::BB;
                }else if evoenv.is_backbone[aa.0]
                || evoenv.is_backbone[dd.1]{
                    pairtype = PairType::SB;
                }
                
                let ddis = evoenv.md_envset.dist[aa.0][dd.0];
                if ddis > dmax{
                    continue;
                }
                let mut edha:f64;
                if ddis <= 1.9 &&  ddis >= dmin {
                    edha = (PI/2.0*(ddis-1.9)/(1.9-dmin)).cos()*-1.0;
                }else if ddis < dmax{
                    edha = (PI*(ddis-1.9)/(dmax-1.9)).cos()*-0.5 -0.5;
                }else{
                    edha = 0.0;
                }

                let mut etheta_dha = 0.0;
                let dangle:f64 =  charmm_based_energy::calc_angle(&evoenv.md_envset.atoms[dd.1],
                    &evoenv.md_envset.atoms[dd.0],
                    &evoenv.md_envset.atoms[aa.0]
                );
                if dangle >= 90.0{
                    etheta_dha = (dangle/180.0*PI).cos().powf(4.0)*-1.0;
                }
                
                let mut ephi_hab:f64 = 0.0;
                let bangle:f64 =  charmm_based_energy::calc_angle(
                    &evoenv.md_envset.atoms[dd.0],
                    &evoenv.md_envset.atoms[aa.0],
                    &evoenv.md_envset.atoms[aa.1]
                );
                if bangle >= 80.0{
                    if (evoenv.is_backbone[dd.1] && evoenv.is_backbone[dd.0])
                    || EvoEF2Env::is_sp2(&evoenv.md_envset.atoms[aa.0]){
                        ephi_hab = ((bangle-150.0)/180.0*PI).cos().powf(4.0)*-1.0;
                    }else{
                        ephi_hab = ((bangle-135.0)/180.0*PI).cos().powf(4.0)*-1.0;
                    }
                }

                match pairtype{
                    PairType::BB =>{
                        edha = edha*weights[HB_BB_DIST];
                        etheta_dha = etheta_dha*weights[HB_BB_THETA];
                        ephi_hab = ephi_hab*weights[HB_BB_PHI];

                        ret[HB_BB_DIST] += edha;
                        ret[HB_BB_THETA] += etheta_dha;
                        ret[HB_BB_PHI] += ephi_hab;
                    },
                    PairType::SS =>{
                        edha = edha*weights[HB_SS_DIST];
                        etheta_dha = etheta_dha*weights[HB_SS_THETA];
                        ephi_hab = ephi_hab*weights[HB_SS_PHI];

                        ret[HB_SS_DIST] += edha;
                        ret[HB_SS_THETA] += etheta_dha;
                        ret[HB_SS_PHI] += ephi_hab;
                    },
                    PairType::SB =>{
                        edha = edha*weights[HB_SB_DIST];
                        etheta_dha = etheta_dha*weights[HB_SB_THETA];
                        ephi_hab = ephi_hab*weights[HB_SB_PHI];

                        ret[HB_SB_DIST] += edha;
                        ret[HB_SB_THETA] += etheta_dha;
                        ret[HB_SB_PHI] += ephi_hab;
                    }
                }
                
                atom_level_energy[aa.0] += edha/2.0;
                atom_level_energy[dd.0] += edha/2.0;
                
                atom_level_energy[aa.0] += etheta_dha/3.0;
                atom_level_energy[dd.0] += etheta_dha/3.0;
                atom_level_energy[dd.1] += etheta_dha/3.0;
                
                atom_level_energy[aa.0] += ephi_hab/3.0;
                atom_level_energy[aa.1] += ephi_hab/3.0;
                atom_level_energy[dd.0] += ephi_hab/3.0;
                
            }
        }
    }
    return ret;
}

//polar/nonpolar で返す
#[allow(non_snake_case)]
pub fn calcEdesolv(evoenv:&EvoEF2Env,atom_level_energy:&mut Vec<f64>,weight_polar:f64,weight_nonpolar:f64)
->(f64,f64){
    let mut solvenergy_all_polar:f64 = 0.0;
    let mut solvenergy_all_nonpolar:f64 = 0.0;
    let num_atoms = evoenv.md_envset.atoms.len();
    for aa in 0..num_atoms-1{
        for bb in (aa+1)..num_atoms{
            
            if evoenv.md_envset.num_edges[aa][bb] < 3{
                continue;
            }
            let ddis = evoenv.md_envset.dist[aa][bb].max(evoenv.md_envset.atoms[aa].nb_r1_2+evoenv.md_envset.atoms[bb].nb_r1_2);
            if ddis >= 6.0{
                continue;
            }
            if evoenv.solvparam[bb].volume*evoenv.solvparam[aa].volume == 0.0{
                continue;
            }
            let ddis2 = ddis*ddis;
            let mut dsol1 = -1.0*evoenv.solvparam[bb].volume
            *evoenv.solvparam[aa].dgfree
            /(2.0*(*PI_3_2)*evoenv.solvparam[aa].correllen*ddis2)
            *(-1.0*
                ((ddis-evoenv.md_envset.atoms[aa].nb_r1_2)/evoenv.solvparam[aa].correllen)
                *((ddis-evoenv.md_envset.atoms[aa].nb_r1_2)/evoenv.solvparam[aa].correllen)
            ).exp();
            
            let mut dsol2 = -1.0*evoenv.solvparam[aa].volume
            *evoenv.solvparam[bb].dgfree
            /(2.0*(*PI_3_2)*evoenv.solvparam[bb].correllen*ddis2)
            *(-1.0*
                ((ddis-evoenv.md_envset.atoms[bb].nb_r1_2)/evoenv.solvparam[bb].correllen)
                *((ddis-evoenv.md_envset.atoms[bb].nb_r1_2)/evoenv.solvparam[bb].correllen)
            ).exp()
            ;
            

            if evoenv.solvparam[aa].polar{
                dsol1 *= weight_polar;
                solvenergy_all_polar += dsol1;
            }else{
                dsol1 *= weight_nonpolar;
                solvenergy_all_nonpolar += dsol1;
            }
            if evoenv.solvparam[bb].polar{
                dsol2 *= weight_polar;
                solvenergy_all_polar += dsol2;
            }else{
                dsol2 *= weight_nonpolar;
                solvenergy_all_nonpolar += dsol2;
            }
            
            atom_level_energy[aa] += dsol1/2.0;
            atom_level_energy[bb] += dsol1/2.0;
            atom_level_energy[aa] += dsol2/2.0;
            atom_level_energy[bb] += dsol2/2.0;
        }
    }
    //deslvP, deslvH
    return (solvenergy_all_polar,solvenergy_all_nonpolar);
}


//SS bond
#[allow(non_snake_case)]
pub fn calcEss(evoenv:&EvoEF2Env,atom_level_energy:&mut Vec<f64>,weight:f64)
->f64{
    let snum:usize = evoenv.cys_atoms.len(); 
    let mut esum:f64 = 0.0;
    if snum <= 1{
        return 0.0;
    }
    for ss1 in 0..snum-1{
        let c1:&CysAtoms = &evoenv.cys_atoms[ss1];
        for ss2 in (ss1+1)..snum{
            let mut e_ss:f64 = 0.0;
            let c2:&CysAtoms = &evoenv.cys_atoms[ss2];
            
            let ssdist = evoenv.md_envset.atoms[c1.sg_index]
            .distance(&evoenv.md_envset.atoms[c2.sg_index]);
            if (ssdist -1.95)*(ssdist -2.15) >= 0.0{
                continue;
            }
            let e_dis = (1.0
                -(10.0*(ssdist -2.03)).exp()
                )*0.8;
            e_ss += e_dis;

            let thetaij = charmm_based_energy::calc_angle(&evoenv.md_envset.atoms[c1.cb_index],
                &evoenv.md_envset.atoms[c1.sg_index],
                &evoenv.md_envset.atoms[c2.sg_index]
            );
            let thetaji = charmm_based_energy::calc_angle(&evoenv.md_envset.atoms[c2.cb_index],
                &evoenv.md_envset.atoms[c2.sg_index],
                &evoenv.md_envset.atoms[c1.sg_index]
            );

            //これと後の二つ個別にウエイトを考える必要があると思う
            let e_ang = 0.005*(thetaij-105.0)*(thetaij-105.0)
            +0.005*(thetaji-105.0)*(thetaji-105.0);
            e_ss += e_ang;


            let dihed1 = charmm_based_energy::calc_dihedral_angle_radian(
                &evoenv.md_envset.atoms[c1.cb_index],
                &evoenv.md_envset.atoms[c1.sg_index],
                &evoenv.md_envset.atoms[c2.sg_index],
                &evoenv.md_envset.atoms[c2.cb_index]
            );

            //ウエイト検討の余地あり
            let e_dihed1 = (2.0*dihed1).cos()+1.0;
            e_ss += e_dihed1;

            let dihed2 = charmm_based_energy::calc_dihedral_angle_radian(
                &evoenv.md_envset.atoms[c1.ca_index],
                &evoenv.md_envset.atoms[c1.cb_index],
                &evoenv.md_envset.atoms[c1.sg_index],
                &evoenv.md_envset.atoms[c2.sg_index]
            )+PI/1.5;//abs は必要だろうか
            
            let dihed3 = charmm_based_energy::calc_dihedral_angle_radian(
                &evoenv.md_envset.atoms[c2.ca_index],
                &evoenv.md_envset.atoms[c2.cb_index],
                &evoenv.md_envset.atoms[c2.sg_index],
                &evoenv.md_envset.atoms[c1.sg_index]
            )+PI/1.5;

            let e_dihed23 = 1.25*dihed2.sin()-1.75
            +1.25*dihed3.sin()-1.75;
            e_ss += e_dihed23;

            e_ss *= weight;

            atom_level_energy[c1.ca_index] += e_ss/6.0;
            atom_level_energy[c1.cb_index] += e_ss/6.0;
            atom_level_energy[c1.sg_index] += e_ss/6.0;
            atom_level_energy[c2.ca_index] += e_ss/6.0;
            atom_level_energy[c2.cb_index] += e_ss/6.0;
            atom_level_energy[c2.sg_index] += e_ss/6.0;

            esum += e_ss;
        }    
    }   
    return esum;
}





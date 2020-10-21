#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use std::fs;
/*

use std::collections::{HashMap,HashSet};
use super::misc_util::*;
use super::charmm_based_energy;
use super::backbone_sample;
use super::side_chain_sample;
use super::charmm_param;


fn file_to_hash(filepath:&str,tags:&Vec<&str>)->HashMap<String,Vec<HashMap<String,String>>>{
    let mut tags_hm:HashMap<String,Vec<HashMap<String,String>>>
    = tags.iter().fold(HashMap::new(),|mut s,m|{s.insert(m.to_string(),vec![]); s}); 
    
    let file = File::open(filepath).unwrap();
    let reader = BufReader::new(file);
    for (_lcount,line_) in reader.lines().enumerate() {
        let line:String = line_.unwrap();
        if start_with(&line,"#"){
            continue;
        }
        if let Some(_) = REGEX_NOLINE.captures(&line){
            continue;
        }
        let mapp = line_to_hash(&line);
        if mapp.len() == 0{
            eprintln!("{} is ignored.",line);
            continue;
        }
        if !mapp.contains_key("tag"){
            panic!("The line must have 'tag' section.\n {}",line);
        }
        let tagg:String = mapp.get("tag").unwrap().clone();
        if !tags_hm.contains_key(&tagg){
            panic!("Tag {} is not a valid tag.",tagg);
        }
        tags_hm.get_mut(&tagg).unwrap().push(mapp);
    }
    return tags_hm;
}


enum AtomFiltering{
    OnlyCB,
    BackboneAndCB,
    All
}

pub struct MCOption{
    input_file:String,
    num_iterations:usize,
    atom_filtering:AtomFiltering,
    charmm_19_top:String,
    charmm_19_par:String,
    evoef2_resource_dir:String,

    backbone_dihedral_file:String,
    distance_bin_file:String,
    contact_file:String,
    
    atoms_grouping:bool,
    random_seed:Option<u64>,
    value_mov:f64,
    num_mov_operation:usize,
    mov_sorting:bool,
    acceptance_bound:String,
    bond_restrictor:Vec<String>,
    bond_attractor:Vec<String>,
    atom_repulsive_penalty:Vec<String>,
    freeze:Vec<String>,
}

impl MCOption{
    pub fn load(filepath:&str){
        let tags:Vec<&str> = vec![
            "input_file",

            "num_iterations",
            "atom_filtering",
            "charmm_19_top",
            "charmm_19_par",
            "evoef2_resource_dir",

            "backbone_dihedral_file",
            "distance_bin_file",
            "contact_file",
            
            "atoms_grouping",
            "random_seed",
            "value_mov",
            "num_mov_operation",
            "mov_sorting",
            "acceptance_bound",
            "bond_restrictor",
            "bond_attractor",
            "atom_repulsive_penalty",
            "freeze"
        ];
        let tags_hm:HashMap<String,Vec<HashMap<String,String>>>
        = file_to_hash(&filepath,&tags); 
        
    }
}
*/
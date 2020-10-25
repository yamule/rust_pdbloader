extern crate rand;
extern crate rust_pdbloader;

use std::fs;
use std::collections::HashSet;
use rust_pdbloader::charmm_based_energy;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
#[allow(unused_imports)]
use self::rand::prelude::*;
use rust_pdbloader::geometry::Vector3D;
use rust_pdbloader::pdbdata::*;
use rust_pdbloader::charmm_param;
use rust_pdbloader::process_3d;
use rust_pdbloader::debug_env;
use rust_pdbloader::structural_alignment;
#[allow(unused_imports)]
use rust_pdbloader::mcprocess_md;
use rust_pdbloader::rosetta_param;
use rust_pdbloader::evoef2_energy;
use rust_pdbloader::sequence_alignment;
use rust_pdbloader::max_hit_clust;
use rust_pdbloader::matrix_process;
#[allow(unused_imports)]
use rust_pdbloader::pp_energy;
#[allow(unused_imports)]
use rust_pdbloader::pp_energy_mc;
#[allow(unused_imports)]
use rust_pdbloader::distance_alignment;
use std::collections::HashMap;
use regex::Regex;
use std::env;
use std::f64::consts::PI;
use rust_pdbloader::misc_util::*;
use rust_pdbloader::template_based_modelling;

#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(unused_imports)]
use std::fs::File;


fn arg_to_hash(args:&Vec<String>)->HashMap<String,String>{
    let mut ret:HashMap<String,String> = HashMap::new();

    let argreg = Regex::new(r"^-[a-z]").unwrap();
    for ii in 0..args.len(){
        if let Some(_) = argreg.captures(args[ii].as_str()){
            if ret.contains_key(&args[ii]) {
                eprintln!("{} is already in use. {} -> {:?} is discarded.", args[ii], args[ii],ret.get(&args[ii]));
            }
            if args.len() > ii+1{
                ret.insert(args[ii].clone(),args[ii+1].clone());
            }else{
                ret.insert(args[ii].clone(),"1".to_string());
            }
        }
    }
    return ret;
}


#[allow(dead_code)]
fn main__(){
    //pp_energy::plain_dihed_test();
    //distance_alignment::dist_aligntest();
    //pp_energy_mc::build_test();
    //pp_energy_mc::peptide_group_test();
    //pp_energy_mc::group_connection_test();
    //pp_energy::tbm_and_loop_test();
    
    //let outfilename:String = "test/1a4w_mapped.pdb".to_owned();
    //let fastafile= "test/1a4w_h_2odp_a.fas";
    //let bset = backbone_sample::BackboneSet::new("D:/dummy/vscode_projects/rust/rust_pdbloader/resources/pepbuilderj/resources/sampledresidues/");
    //let sset = side_chain_sample::SideChainSet::new("D:/dummy/vscode_projects/rust/rust_pdbloader/resources/pepbuilderj/resources/sampledresidues/");

    /*template_based_modelling::try_mapping(
        &vec!["test/1a4w_h_2odp_a.fas".to_owned()]
        ,"D:/dummy/vscode_projects/rust/rust_pdbloader/resources/pepbuilderj/resources/sampledresidues/"
        ,"D:/dummy/vscode_projects/rust/rust_pdbloader/resources/pepbuilderj/resources/sampledresidues/"
        , Some(123),"test/1a4w_mapped.pdb");
    */
    
}

//Todo ファイルがない場合のエラーメッセージにパス表示
//

fn main(){
        //コマンド例
    //cargo run --release -- -sequence_clustering -in D:/dummy/vbox_share/casp_decoys/get_ref_structure/casp_fas_test.fas -identity 0.95 -coverage_short 0.95 -kmer_filt 2 -fmt_blastclust
    
    let mut argss: Vec<String> = env::args().collect();
    let _exe:String = argss.remove(0);
    let subcommand:String = argss.remove(0);
    //let argss:Vec<String> = "test test".split_whitespace().into_iter().map(|m|m.to_string()).collect();
    let args:HashMap<String,String> = arg_to_hash(&argss);
    match subcommand.as_str(){
        "sequence_alignment" => {
            main_seq_align(args);
            return;
        },
        "sequence_clustering"=>{
            main_sequence_clustering(args);
            return;
        },
        "structural_alignment"=>{
            main_str_align(args);
            return;
        },
        "comparative_domain_split"=>{
            main_comparative_domain_split(args);
            return;
        },
        "phi_psi_angle"=>{
            calc_phi_psi(args);
            return;
        },
        "prepare_structure"=>{
            main_prepare_structure(args);
            return;
        },
        "residue_mapping"=>{
            residue_mapping(args);
            return;
        },
        "refinement"=>{
            refinement(args);
            return;
        },
        "merge_structures"=>{
            merge_structures(args);
            return;
        },
        "calc_energy"=>{
            calc_energies(args);
            return;
        },
        "make_homo_multimer"=>{
            make_homo_multimer(args);
            return;
        },
        "docking"=>{
            docking(args);
            return;
        },
        _=>{
            panic!("not supported command. {:?}",subcommand);
        }
    }
}
fn residue_mapping(args:HashMap<String,String>) {
    let infiles:Vec<String> = args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in.")).split(",").into_iter().map(|m| m.to_string()).collect();
    template_based_modelling::try_mapping(
    &infiles
    , args.get("-sidechains").unwrap_or_else(|| panic!("Please specify -sidechains for sidechain sample dir."))
    , args.get("-backbones").unwrap_or_else(|| panic!("Please specify  -backbones for backbone sample dir."))
    , match args.get("-seed"){Some(x) => {Some(x.parse::<u64>().unwrap())},_=>{None}}
    ,args.get("-out").unwrap_or_else(|| panic!("Please specify  -out for output file.")));

}
fn calc_energies(args:HashMap<String,String>) {
    let toph19:&str = args.get("-toph19").unwrap_or_else(|| panic!("Please specify charmm19 toph19.inp file with -toph19."));
    let param19:&str = args.get("-param19").unwrap_or_else(|| panic!("Please specify charmm19 param19.inp file with -param19."));
    let evoef2resource:&str = args.get("-resource").unwrap_or_else(|| panic!("Please specify resource dir with -resource"));
    let angle:&str = args.get("-angle").unwrap_or_else(|| panic!("Please specify backbone torsion angle file with -angle"));
    template_based_modelling::calc_energies(
    args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in."))
    ,toph19
    ,param19
    ,evoef2resource
    ,angle
    , match args.get("-cb_dist"){Some(x) => {x},_=>{""}}
    ,args.get("-out").unwrap_or_else(|| panic!("Please specify input file with -out."))
    );
}
fn load_refine_params(paramfile:&str)->Vec<template_based_modelling::RefinementParam>{
    let file = File::open(paramfile).unwrap();
    let reader = BufReader::new(file);
    let mut ret:Vec<template_based_modelling::RefinementParam> = vec![];
    let checker:HashSet<String> = vec![
            "num_iter",
            "wdec_factor",
            "inv_edge",
            "inv_dist",
            "num_target_atoms",
            "mov_dist",
            "lbfgs",
            "sorter",
            "shuffle",
            "backbone_cb_only",
            "accept_lb",
            "accept_ub",
            "steps_bound_check",
            "refine_loop_threshold",
            "use_template",
            "refine_retry",
            "group_rotate",
        ].into_iter().map(|m|m.to_owned()).collect();
        
    for (_lcount,line_) in reader.lines().enumerate() {
        let line:String = line_.unwrap();
        if start_with(&line,"#"){
            continue;
        }
        if let Some(_) = REGEX_NOLINE.captures(&line){
            continue;
        }
        let mapp = line_to_hash(&line);
        for (kk,_vv) in mapp.iter(){
            if !checker.contains(kk){
                panic!("{} is not found in dict!",kk);
            }
        }
        let parr:template_based_modelling::RefinementParam = template_based_modelling::RefinementParam{
            num_iter:mapp.get("num_iter").unwrap().split(",").into_iter().map(|m|m.parse::<usize>().unwrap()).collect(),
            inv_edge:mapp.get("inv_edge").unwrap_or(&"0".to_owned()).split(",").into_iter().map(|m|m.parse::<u64>().unwrap()).collect(),
            inv_dist:mapp.get("inv_dist").unwrap_or(&"0.0".to_owned()).split(",").into_iter().map(|m|m.parse::<f64>().unwrap()).collect(),
            wdec_factor:mapp.get("wdec_factor").unwrap_or(&"0.0".to_owned()).split(",").into_iter().map(|m|m.parse::<f64>().unwrap()).collect(),
            num_target_atoms:mapp.get("num_target_atoms").unwrap_or(&"0".to_owned()).split(",").into_iter().map(|m|m.parse::<usize>().unwrap()).collect(),
            mov_dist:mapp.get("mov_dist").unwrap_or(&"0.0".to_owned()).split(",").into_iter().map(|m|m.parse::<f64>().unwrap()).collect(),
            refine_retry:mapp.get("refine_retry").unwrap_or(&"0".to_owned()).split(",").into_iter().map(|m|m.parse::<usize>().unwrap()).collect(),
            lbfgs:mapp.get("lbfgs").unwrap_or(&"0".to_owned()).split(",").into_iter().map(|m|m.parse::<usize>().unwrap()).collect(),
            refine_loop_threshold:mapp.get("refine_loop_threshold").unwrap_or(&"0.0".to_owned()).split(",").into_iter().map(|m|m.parse::<f64>().unwrap()).collect(),
            group_rotate:mapp.get("group_rotate").unwrap_or(&"-1.0".to_owned()).split(",").into_iter().map(|m|m.parse::<f64>().unwrap()).collect(),
            
            sorter:mapp.get("sorter").unwrap_or(&"false".to_owned()).parse::<bool>().unwrap(),
            shuffle:mapp.get("shuffle").unwrap_or(&"false".to_owned()).parse::<bool>().unwrap(),
            backbone_cb_only:mapp.get("backbone_cb_only").unwrap_or(&"false".to_owned()).parse::<bool>().unwrap(),
            accept_bound:(
                mapp.get("accept_lb").unwrap_or(&"0.0".to_owned()).parse::<f64>().unwrap()
                ,mapp.get("accept_ub").unwrap_or(&"0.0".to_owned()).parse::<f64>().unwrap()
            ),
            steps_bound_check:mapp.get("steps_bound_check").unwrap_or(&"0".to_owned()).parse::<usize>().unwrap(),
            use_template:mapp.get("use_template").unwrap_or(&"false".to_owned()).parse::<bool>().unwrap(),
        };
        ret.push(parr);
    }
    return ret;
}

//refinement template loader
fn load_template_list(templatefile:&str)->Vec<(String,String,String)>{
    let file = File::open(templatefile).unwrap();
    let reader = BufReader::new(file);
    let mut ret:Vec<(String,String,String)> = vec![];
    let checker:HashSet<String> = vec![
            "file",
            "template_chain",
            "query_chain",
        ].into_iter().map(|m|m.to_owned()).collect();
        
    for (_lcount,line_) in reader.lines().enumerate() {
        let line:String = line_.unwrap();
        if start_with(&line,"#"){
            continue;
        }
        if let Some(_) = REGEX_NOLINE.captures(&line){
            continue;
        }
        let mapp = line_to_hash(&line);
        for (kk,_vv) in mapp.iter(){
            if !checker.contains(kk){
                panic!("{} is not found in dict!",kk);
            }
        }
        ret.push((
            mapp.get("query_chain").unwrap_or(&("".to_owned())).to_string()
            ,mapp.get("file").unwrap_or_else(||panic!("file must be specified.")).to_string()
            ,mapp.get("template_chain").unwrap_or_else(||panic!("template_chain must be specified.")).to_string()
        ));
    }
    return ret;
}

fn merge_structures(args:HashMap<String,String>) {
    let toph19:&str = args.get("-toph19").unwrap_or_else(|| panic!("Please specify charmm19 toph19.inp file with -toph19."));
    let param19:&str = args.get("-param19").unwrap_or_else(|| panic!("Please specify charmm19 param19.inp file with -param19."));
    let evoef2resource:&str = args.get("-resource").unwrap_or_else(|| panic!("Please specify resource dir with -resource"));
    let angle:&str = args.get("-angle").unwrap_or_else(|| panic!("Please specify backbone torsion angle file with -angle"));

    let base_:PDBEntry = load_pdb(args.get("-base").unwrap_or_else(|| panic!("Please specify input file with -in.")));
    let sample_:PDBEntry = load_pdb(args.get("-sample").unwrap_or_else(|| panic!("Please specify input file with -in.")));
    if base_.chains.len() > 1 || sample_.chains.len() > 1{
        panic!("Number of chains must be one. {} {}",base_.chains.len(),sample_.chains.len());
    }
    let mut base:Vec<PDBResidue> = vec![];
    let mut sample:Vec<PDBResidue> = vec![];
    let mut chain_name:String = "A".to_owned();
    for mut cc in base_.chains.into_iter(){
        chain_name = cc.chain_name.clone();
        cc.remove_alt(None);
        for rr in cc.residues.into_iter(){
            base.push(rr);
        }
    }

    for mut cc in sample_.chains.into_iter(){
        cc.remove_alt(None);
        for rr in cc.residues.into_iter(){
            sample.push(rr);
        }
    }
    
    let outfile:String =  args.get("-out").unwrap_or_else(|| panic!("Please specify input file with -out.")).to_string();
    template_based_modelling::merge_structure(
    &base
    ,&sample
    ,&chain_name
    ,toph19
    ,param19
    ,evoef2resource
    ,angle
    , match args.get("-cb_dist"){Some(x) => {x},_=>{""}}
    , match args.get("-residue_contact"){Some(x) => {x},_=>{""}}
    , &load_refine_params(args.get("-param_file").unwrap_or_else(|| panic!("Please specify refinement parameter file with -param_file")).as_ref())
    , match args.get("-seed"){Some(x) => {Some(x.parse::<u64>().unwrap())},_=>{None}}
    ,&outfile
    ,args.get("-top_x").unwrap_or(&("0".to_owned())).parse::<usize>().unwrap()
    );
}
fn refinement(args:HashMap<String,String>) {
    let checker:HashSet<String> = vec![
    "-toph19",
    "-param19",
    "-resource",
    "-angle",
    "-build_missing_param1",
    "-build_missing_param2",
    "-num_angles_restrict",
    "-num_structurs_step1",
    "-num_structurs_step2",
    "-change_omega",
    "-in",
    "-out",
    "-template_list",
    "-flag",
    "-group",
    "-num_incremental_build",
    "-cb_dist",
    "-residue_contact",
    "-param_file",
    "-rebuild_iter",
    "-rebuild_region",
    "-rebuild_num",
    "-seed",
    "-steps_checkpoint",
    ].into_iter().map(|m|m.to_owned()).collect();    

    for (kk,_vv) in args.iter(){
        if !checker.contains(kk){
            panic!("{} was not found in dict.",kk);
        }
    }

    let toph19:&str = args.get("-toph19").unwrap_or_else(|| panic!("Please specify charmm19 toph19.inp file with -toph19."));
    let param19:&str = args.get("-param19").unwrap_or_else(|| panic!("Please specify charmm19 param19.inp file with -param19."));
    let evoef2resource:&str = args.get("-resource").unwrap_or_else(|| panic!("Please specify resource dir with -resource"));
    let angle:&str = args.get("-angle").unwrap_or_else(|| panic!("Please specify backbone torsion angle file with -angle"));
    
    let param1:String = args.get("-build_missing_param1").unwrap_or(&"".to_string()).to_string();
    let param2:String = args.get("-build_missing_param2").unwrap_or(&"".to_string()).to_string();
    let phipsi_topx:usize = args.get("-num_angles_restrict").unwrap_or(&("100".to_owned())).parse::<usize>().unwrap();
    let num_structures_step1:usize = args.get("-num_structurs_step1").unwrap_or(&("1000".to_owned())).parse::<usize>().unwrap();
    let num_structures_step2:usize = args.get("-num_structurs_step2").unwrap_or(&("10".to_owned())).parse::<usize>().unwrap();

    let backboneparam:Option<template_based_modelling::MissingBuilderParam>
     = if param1.len() > 0 && param2.len() > 0{
        Some(template_based_modelling::MissingBuilderParam{
            param1:load_refine_params(&param1),
            param2:load_refine_params(&param2),
            phipsi_top_x:phipsi_topx,
            num_structures_step1:num_structures_step1,
            num_structures_step2:num_structures_step2,
            fix_omega:!args.contains_key("-change_omega")//Omega には PI を使うのがデフォルト
        })
    }else if param1.len() > 0 || param2.len() > 0{
        panic!("To construct missing region, -build_missing_param1 & build_missing_param2 are needed. {} {} ",param1.len(), param2.len())
    }else{
        None
    };
    let mut pdbb:PDBEntry = load_pdb(args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in.")));
    for cc in pdbb.chains.iter_mut(){
        cc.remove_alt(None);
        charmm_based_energy::MDAtom::change_to_charmmnames(&mut cc.residues);
    }
    let outfile:String =  args.get("-out").unwrap_or_else(|| panic!("Please specify input file with -out.")).to_string();
    let ppenergyset:pp_energy::PPEnergySet =template_based_modelling::refinement(
    pdbb
    ,match args.get("-template_list"){Some(x) => {Some(load_template_list(x.as_str()))},_=>{None}}
    ,match args.get("-flag"){Some(x) => {Some(x.to_string())},_=>{None}}
    ,match args.get("-group"){Some(x) => {Some(x.to_string())},_=>{None}}
    ,args.get("-num_incremental_build").unwrap_or(&("0".to_owned())).parse::<usize>().unwrap()
    ,toph19
    ,param19
    ,evoef2resource
    ,angle
    , match args.get("-cb_dist"){Some(x) => {x},_=>{""}}
    , match args.get("-residue_contact"){Some(x) => {x},_=>{""}}
    , &load_refine_params(args.get("-param_file").unwrap_or_else(|| panic!("Please specify refinement parameter file with -param_file")).as_ref())
    , backboneparam
    , Some(template_based_modelling::RebuildParam{
        rebuild_iter:match args.get("-rebuild_iter"){Some(x) => {x.parse::<usize>().unwrap()},_=>{0}},
        rebuild_region:match args.get("-rebuild_region"){Some(x) => {x.parse::<usize>().unwrap()},_=>{0}},
        rebuild_num:match args.get("-rebuild_num"){Some(x) => {x.parse::<usize>().unwrap()},_=>{0}},
    })
    , match args.get("-seed"){Some(x) => {Some(x.parse::<u64>().unwrap())},_=>{None}}
    , match args.get("-steps_checkpoint"){Some(x) => {x.parse::<usize>().unwrap()},_=>{0}}
    ,&outfile//checkpoint が出力される
    );
    
    let mut lines:Vec<String> = vec![];
    let mut atom_based_energies:Vec<f64> = vec![0.0;ppenergyset.evoef2_env.md_envset.atoms.len()];
    let _pres = ppenergyset.calc_energy_sep(&mut atom_based_energies);
    for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
        att.temp_factor = atom_based_energies[aii].max(-999.0).min(999.0);
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file(format!("{}",outfile).as_str(),lines);
}

    

fn docking(args:HashMap<String,String>) {
    let toph19:&str = args.get("-toph19").unwrap_or_else(|| panic!("Please specify charmm19 toph19.inp file with -toph19."));
    let param19:&str = args.get("-param19").unwrap_or_else(|| panic!("Please specify charmm19 param19.inp file with -param19."));
    let evoef2resource:&str = args.get("-resource").unwrap_or_else(|| panic!("Please specify resource dir with -resource"));
    let angle:&str = args.get("-angle").unwrap_or_else(|| panic!("Please specify backbone torsion angle file with -angle"));
    let mut pdbb:PDBEntry = load_pdb(args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in.")));
    for cc in pdbb.chains.iter_mut(){
        cc.remove_alt(None);
        charmm_based_energy::MDAtom::change_to_charmmnames(&mut cc.residues);
    }
    let outfile:String =  args.get("-out").unwrap_or_else(|| panic!("Please specify input file with -out.")).to_string();
    let ppenergyset:pp_energy::PPEnergySet =template_based_modelling::docking(
    pdbb
    ,match args.get("-group"){Some(x) => {Some(x.to_string())},_=>{None}}
    ,toph19
    ,param19
    ,evoef2resource
    ,angle
    , match args.get("-cb_dist"){Some(x) => {x},_=>{""}}
    , match args.get("-residue_contact"){Some(x) => {x},_=>{""}}
    , &load_refine_params(args.get("-param_file").unwrap_or_else(|| panic!("Please specify refinement parameter file with -param_file")).as_ref())
    , match args.get("-num_structures_step2"){Some(x) => {x.parse::<usize>().unwrap()},_=>{10}}
    , match args.get("-angle_resolution"){Some(x) => {Some(x.parse::<f64>().unwrap())},_=>{None}}
    , match args.get("-seed"){Some(x) => {Some(x.parse::<u64>().unwrap())},_=>{None}}
    );
    
    let mut lines:Vec<String> = vec![];
    let mut atom_based_energies:Vec<f64> = vec![0.0;ppenergyset.evoef2_env.md_envset.atoms.len()];
    let _pres = ppenergyset.calc_energy_sep(&mut atom_based_energies);
    for (aii,aa) in ppenergyset.evoef2_env.md_envset.atoms.iter().enumerate(){
        let  (chainid,(resname,resnum,altcode),mut att) = aa.to_pdbatom();
        att.temp_factor = atom_based_energies[aii].max(-999.0).min(999.0);
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&altcode));
    }
    write_to_file(format!("{}",outfile).as_str(),lines);
}

    
fn make_homo_multimer(args:HashMap<String,String>) {
    let outfilename:String = args.get("-out").unwrap_or_else(|| panic!("Please specify out file with -out.")).to_string();
    let query_pdb:PDBEntry = load_pdb(args.get("-query").unwrap_or_else(|| panic!("Please specify query file with -query.")));
    let template_pdb:PDBEntry = load_pdb(args.get("-template").unwrap_or_else(|| panic!("Please specify template file with -template.")));
    let mut query_chain:HashSet<String> = args.get("-query_chain").unwrap_or(&("".to_owned())).split(",").map(|m|m.to_string()).collect();
    let mut template_chain:HashSet<String> = args.get("-template_chain").unwrap_or(&("".to_owned())).split(",").map(|m|m.to_string()).collect();
    let mut result_string:Vec<String> = vec![];
    let mut group_string:Vec<String> = vec![];
    let mut chaincount:usize = 0;
    let alignment_type_:String = args.get("-type").unwrap_or(&("sw".to_string())).clone();

    let alignment_type:structural_alignment::AlignmentType = if alignment_type_ == "sw"{
        structural_alignment::AlignmentType::SW
    }else if alignment_type_ == "raw"{
        structural_alignment::AlignmentType::RAW
        
    }else if alignment_type_ == "raw_full"{
        structural_alignment::AlignmentType::RAW_FULLRMSD
        
    }else if alignment_type_ == "max"{
        structural_alignment::AlignmentType::MAXIMUM
    }else{
        panic!("{} is not recognized. sw, raw, max, raw_full",alignment_type_);
    };


    query_chain.remove("");
    template_chain.remove("");
    for qcc in query_pdb.chains.iter(){
        if query_chain.len() > 0{
            if !query_chain.contains(&qcc.chain_name){
                continue;
            }
        }
        let mut qresidues:Vec<&PDBResidue> = vec![];
        for rr in qcc.residues.iter(){
            qresidues.push(rr);
        }
        
        for tcc in template_pdb.chains.iter(){
            if template_chain.len() > 0{
                if !template_chain.contains(&tcc.chain_name){
                    continue;
                }
            }
            let mut tresidues:Vec<&PDBResidue> = vec![];
            for rr in tcc.residues.iter(){
                tresidues.push(rr);
            }
            let res_ = structural_alignment::align_pdb(&qresidues,&tresidues,alignment_type.clone(),0.2);
            let res:structural_alignment::StructuralAlignmentResult = res_.unwrap();
            for rr in qcc.residues.iter(){
                for aa_ in rr.iter_atoms(){
                    let mut aa = aa_.clone();
                    let mres = matrix_process::matrix_multi(&res.transform_matrix,&vec![vec![aa_.get_x()],vec![aa_.get_y()],vec![aa_.get_z()],vec![1.0]]);
                    aa.set_xyz(mres[0][0],mres[1][0],mres[2][0]);
                    result_string.push(aa.get_pdb_atom_line_string(&tcc.chain_name,&rr.get_residue_name(),rr.get_residue_number(),rr.get_ins_code()));
                }
                if rr.get_ins_code().len() > 0{
                    eprintln!("Inscode is not supported for grouping.");
                }
                group_string.push(format!("chain:{}\tresidue_name:{}\tresidue_number:{}\tgroup:{}",tcc.chain_name,&rr.get_residue_name(),rr.get_residue_number(),chaincount));
            }
            chaincount += 1;
        }
    }
    write_to_file(&outfilename,result_string);
    write_to_file(&(outfilename+".chaingroup"),group_string);
}
    
fn calc_phi_psi(args:HashMap<String,String>) {
    
    let mut lines:Vec<String> = vec![];
    let mut pdbb:PDBEntry = load_pdb(args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in.")));
    let outfilename:&str = args.get("-out").unwrap_or_else(|| panic!("Please specify output file with -out."));
    for cc in pdbb.chains.iter_mut(){
        cc.remove_alt(None);
        let rnum =cc.residues.len();
        let dist_threshold:f64 = 2.0;
        for rr in 0..rnum{
            
            let mut prevca_:Option<&PDBAtom> = if rr == 0 {None}else{cc.residues[rr-1].get_CA()};
            let mut prevc_:Option<&PDBAtom> = if rr == 0 {None}else{cc.residues[rr-1].get_C()};
            let currentn_:Option<&PDBAtom> = cc.residues[rr].get_N();
            let currentca_:Option<&PDBAtom> = cc.residues[rr].get_CA();
            let currentc_:Option<&PDBAtom> = cc.residues[rr].get_C();
            let mut nextn_:Option<&PDBAtom> = if rr == rnum-1{None}else{cc.residues[rr+1].get_N()};
            let mut nextca_:Option<&PDBAtom> = if rr == rnum-1{None}else{cc.residues[rr+1].get_CA()};

            
            //前の残基とは結合してないとみなす
            if let None = prevc_{
                prevca_ = None;
                prevc_ = None;
            }else if let None = currentn_{
                prevca_ = None;
                prevc_ = None;
            }else{
                if process_3d::distance(&prevc_.clone().unwrap().get_xyz(),&currentn_.clone().unwrap().get_xyz()) > dist_threshold{
                    prevca_ = None;
                    prevc_ = None;
                }
            }

            
            //後ろの残基とは結合してないとみなす
            if let None = nextn_{
                nextn_ = None;
                nextca_ = None;
            }else if let None = currentc_{
                nextn_ = None;
                nextca_ = None;
            }else{
                if process_3d::distance(&currentc_.clone().unwrap().get_xyz(),&nextn_.clone().unwrap().get_xyz()) > dist_threshold{
                    nextn_ = None;
                    nextca_ = None;
                }
            } 


            let mut phi:String = "-".to_string();
            let mut psi:String = "-".to_string();
            let mut prev_omega:String = "-".to_string();
            let mut next_omega:String = "-".to_string();


            //dihedral angle は charmm とは逆回り（なはず）なので -1.0 を掛ける

            //こんな書き方しかできないだろうか。。。
            let mut phiflag = true;
            let mut psiflag = true;
            let mut pomegaflag = true;
            let mut nomegaflag = true;

            if let None = prevca_{
                pomegaflag = false;
            }
            if let None = prevc_{
                pomegaflag = false;
                phiflag = false;
            }
            if let None = currentn_{
                pomegaflag = false;
                phiflag = false;
                psiflag = false;
            }
            if let None = currentca_{
                pomegaflag = false;
                phiflag = false;
                psiflag = false;
                nomegaflag = false;
            }
            if let None = currentc_{
                phiflag = false;
                psiflag = false;
                nomegaflag = false;
            }
            if let None = nextn_{
                psiflag = false;
                nomegaflag = false;
            }
            if let None = nextca_{
                nomegaflag = false;
            }

            if pomegaflag{
                prev_omega = (charmm_based_energy::calc_dihedral_angle_radian(prevca_.clone().unwrap(),prevc_.clone().unwrap(),currentn_.clone().unwrap(),currentca_.clone().unwrap())*-1.0/PI*180.0).to_string();
            }
            if phiflag{
                phi = (charmm_based_energy::calc_dihedral_angle_radian(prevc_.clone().unwrap(),currentn_.clone().unwrap(),currentca_.clone().unwrap(),currentc_.clone().unwrap())*-1.0/PI*180.0).to_string();
            }
            
            if psiflag {
                psi = (charmm_based_energy::calc_dihedral_angle_radian(currentn_.clone().unwrap(),currentca_.clone().unwrap(),currentc_.clone().unwrap(),nextn_.clone().unwrap())*-1.0/PI*180.0).to_string();
            }
            if nomegaflag {
                next_omega = (charmm_based_energy::calc_dihedral_angle_radian(currentca_.clone().unwrap(),currentc_.clone().unwrap(),nextn_.clone().unwrap(),nextca_.clone().unwrap())*-1.0/PI*180.0).to_string();
            }

            lines.push(format!("chain_name:\t{}\tresidue_name:\t{}\tresidue_number:\t{}\tins_code:\t{}\tomega_prev:\t{}\tphi:\t{}\tpsi:\t{}\tomega_next:\t{}"
            ,cc.chain_name
            ,cc.residues[rr].residue_name
            ,cc.residues[rr].get_residue_number()
            ,cc.residues[rr].ins_code
            ,prev_omega,phi,psi,next_omega));
        }
    }
    write_to_file(outfilename,lines);
}



//長い方から短い方に検索かけて Threshold 以上のヒットがあった場合長い方にまとめる。
//まとめられたものは次以降の検索にかからないしかける側にもならない。
//cd-hit と基本同じアルゴリズムになるはず
fn main_sequence_clustering(args:HashMap<String,String>) {
    let fas = sequence_alignment::SeqData::load_fasta(args.get("-in").unwrap_or_else(|| panic!("Please specify input file with -in.")),false);
    let options_expected:HashSet<String>=vec![
        "-in"
        ,"-identity"
        ,"-coverage_long"
        ,"-coverage_short"
        ,"-kmer_filt"
        ,"-fmt_blastclust"
        ,"-out_fasta"
    ].iter().map(|m|m.to_string()).collect();

    for kk in args.iter(){
        if !options_expected.contains(kk.0){
            eprintln!("{} is not found in option dictionaly.",kk.0);
        }
    }
    let identity_threshold = args.get("-identity").unwrap_or(&("0.9".to_string())).parse::<f64>().unwrap_or(0.9);
    let coverage_threshold_long = args.get("-coverage_long").unwrap_or(&("0.9".to_string())).parse::<f64>().unwrap_or(0.9);
    let coverage_threshold_short = args.get("-coverage_short").unwrap_or(&("0.9".to_string())).parse::<f64>().unwrap_or(0.9);
    let kmer_filter:usize = args.get("-kmer_filt").unwrap_or(&("1".to_string())).parse::<usize>().unwrap_or(1);
    
    println!("#file:{}",args.get("-in").unwrap());
    println!("#identity_threshold:{}",identity_threshold);
    println!("#coverage_long:{}",coverage_threshold_long);
    println!("#coverage_short:{}",coverage_threshold_short);
    println!("#kmer_filt:{}",kmer_filter);
    let outfas = if args.contains_key("-out_fasta"){fas.clone()}else{vec![]};
    let res = max_hit_clust::cluster(fas
        ,identity_threshold
        ,coverage_threshold_long
        ,coverage_threshold_short
        ,kmer_filter
    );
    if  args.contains_key("-out_fasta"){
        let mut f = BufWriter::new(fs::File::create(args.get("-out_fasta").unwrap()).unwrap());
        let hss:HashSet<String> = res.iter().map(|m|m.name_representative.clone()).collect();
        for ss in outfas.into_iter(){
            if hss.contains(&ss.name){
                f.write_all(format!(">{} {}\n",ss.name,ss.desc).as_bytes()).unwrap();
                let sseq:String = ss.seq.into_iter().fold("".to_string(),|s,a|s+&a);
                f.write_all(sseq.as_bytes()).unwrap();
                f.write_all("\n".as_bytes()).unwrap();
            }
        }
    }
    if args.contains_key("-fmt_blastclust"){
        for rr in res.iter(){
            print!("{}",rr.name_representative);
            for mm in rr.members.iter(){
                print!(" {}",mm.0);
            }
            print!("\n");
        }
    }else{
        for rr in res.iter(){
            println!(">{}",rr.name_representative);
            for mm in rr.members.iter(){
                println!("-\t{}\t{}",mm.0,mm.1);
            }
        }
    }
}
fn main_seq_align(args:HashMap<String,String>) {
    let alignment_type:usize;
    let mess:&str;
    if args.contains_key("-global"){
        alignment_type = sequence_alignment::ALIGN_GLOBAL;
        mess = "Global alignment";
    }else if args.contains_key("-glocal"){
        mess = "Glocal alignment";
        alignment_type = sequence_alignment::ALIGN_GLOCAL;
    }else{
        mess = "Local alignment";
        alignment_type = sequence_alignment::ALIGN_LOCAL;
    }
    let seq1 = sequence_alignment::SeqData::load_fasta(args.get("-seq1").unwrap_or_else(|| panic!("-seq1 was not found.")),false);
    let seq2 = sequence_alignment::SeqData::load_fasta(args.get("-seq2").unwrap_or_else(|| panic!("-seq1 was not found.")),false);
    let mut sw = sequence_alignment::SequenceAlignment::new(Box::new(sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,alignment_type);
    for ss1 in seq1.iter(){
        for ss2 in seq2.iter(){
            let res = sw.align(ss1,ss2,args.contains_key("-retain_all"));
            let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
            let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);

            if args.contains_key("-fastaout"){
                println!(">{} {} #{} with {} \n{}\n",ss1.name,ss1.desc,mess,ss2.name,r1);
                println!(">{} {} #{} with {} \n{}\n",ss2.name,ss2.desc,mess,ss1.name,r2);
            }else{
                println!("#score:{}",res.2);
                println!("#type:{}",mess);
                println!("{} {}\n",ss1.name,r1);
                println!("{} {}\n",ss2.name, r2);
            }
            
        }
    }
    return;

}
fn main_comparative_domain_split(args:HashMap<String,String>) {
        if args.contains_key("-h"){
            eprintln!("usage:  rust_pdbloader.exe -comparative_domain_split -query <query_pdb> -template  <template_pdb> [-max_dist <float: default 4.0>] [-min_length <int: default 10>] [-out <output_file>] [-out_pdb <aligned_query_pdb>] [-query_chain <str>] [-template_chain <str>] ");
        }
        //raw_full はカットオフを行わず全 CA についての RMSD を出す。
        
        let qfilename = args.get("-query").unwrap_or_else(|| panic!("-query was not found."));
        let tfilename = args.get("-template").unwrap_or_else(|| panic!("-template was not found."));
        let mut query_pdb:PDBEntry = load_pdb(qfilename);
        let mut template_pdb:PDBEntry = load_pdb(tfilename);
        let result_file:String = args.get("-out").unwrap_or(&("".to_string())).clone(); 
        let outfile_pdb:String = args.get("-out_pdb").unwrap_or(&("".to_string())).clone(); 
        //let tmscore_break:f64 = args.get("-tmscore_break").unwrap_or(&("0.0".to_string())).parse::<f64>().unwrap_or_else(|e| panic!("cannot parse tmscore_break {:?}",e)).clone(); 
        let max_dist:f64 = args.get("-max_dist").unwrap_or(&("4.0".to_string())).parse::<f64>().unwrap_or_else(|e| panic!("cannot parse -max_dist {:?}",e)).clone(); 
        let min_length:usize = args.get("-min_length").unwrap_or(&("10".to_string())).parse::<usize>().unwrap_or_else(|e| panic!("cannot parse -min_length {:?}",e)).clone(); 
        let query_chain:String  = args.get("-query_chain").unwrap_or(&("".to_string())).clone(); 
        let template_chain:String  = args.get("-template_chain").unwrap_or(&("".to_string())).clone(); 
        
        let mut q_chainindex_:i64 = -1;
        let mut t_chainindex_:i64 = -1;
        
        if query_chain.len() != 0{
            for cc in query_pdb.chains.iter().enumerate(){
                if cc.1.chain_name == query_chain{
                    q_chainindex_ = cc.0 as i64;
                }
            }
            if q_chainindex_ < 0{
                panic!("Chain {} was not found.",query_chain);
            }
        }else{
            if query_pdb.chains.len() > 1{
                panic!("Number of chains in query pdb must be one or use -query_chain option.");
            }else{
                q_chainindex_ = 0;
            }
        }
        
        if template_chain.len() != 0{
            for cc in template_pdb.chains.iter().enumerate(){
                if cc.1.chain_name == template_chain{
                    t_chainindex_ = cc.0 as i64;
                }
            }
            if t_chainindex_ < 0{
                panic!("Chain {} was not found.",template_chain);
            }
        }else{
            if template_pdb.chains.len() > 1{
                panic!("Number of chains in template pdb must be one or use -template_chain option.");
            }else{
                t_chainindex_ = 0;
            }
        }

        let q_chainindex:usize = q_chainindex_ as usize;
        let t_chainindex:usize = t_chainindex_ as usize;
        
        let mut query_cas:Vec<Vec<f64>> = vec![];
        let mut query_cas_index:Vec<usize> = vec![];
        let mut template_cas:Vec<Vec<f64>> = vec![];
        let mut template_cas_index:Vec<usize> = vec![];
        query_pdb.chains[q_chainindex].remove_alt(None);
        template_pdb.chains[t_chainindex].remove_alt(None);
        for (rii,rr) in query_pdb.chains[q_chainindex].residues.iter().enumerate(){
            let ca_:Option<&PDBAtom> = rr.get_CA();
            if let None = ca_{
                continue;
            }
            let ca:&PDBAtom = ca_.unwrap();
            query_cas.push(vec![ca.get_x(),ca.get_y(),ca.get_z()]);
            query_cas_index.push(rii);
        }


        for (rii,rr) in template_pdb.chains[t_chainindex].residues.iter().enumerate(){
            let ca_:Option<&PDBAtom> = rr.get_CA();
            if let None = ca_{
                continue;
            }
            let ca:&PDBAtom = ca_.unwrap();
            template_cas.push(vec![ca.get_x(),ca.get_y(),ca.get_z()]);
            template_cas_index.push(rii);
        }

        let res = structural_alignment::comparative_domain_split(
        &query_cas
        ,&template_cas
        ,min_length
        ,max_dist
        ,100
        ,4
        ,3.0
        ,0.1
        ,0.2);

        let mut res_str:Vec<String> = vec![];
        for ss in res.iter(){
            let mut sv:String = query_cas_index[ss[0]].to_string();
            let slen:usize = ss.len();
            for s in 1..slen{ 
                sv += ",";
                sv += &(query_cas_index[ss[s]].to_string());
            }
            res_str.push(sv);
        }
        if result_file.len() == 0{
            for ss in res_str.iter(){
                println!("{}",ss);
            }
        }else{
            write_to_file(&result_file,res_str);
        }
        if outfile_pdb.len() > 0{
            let chainname:&str = &query_pdb.chains[q_chainindex].chain_name;
            for (sii,ss) in res.iter().enumerate(){
                let mut res_str:Vec<String> = vec![];
                let slen:usize = ss.len();
                for s in 0..slen{ 
                    let r:&PDBResidue = &query_pdb.chains[q_chainindex].residues[query_cas_index[ss[s]]];
                    for aa in r.iter_atoms(){
                        res_str.push(aa.get_pdb_atom_line_string(
                            chainname,
                            r.get_residue_name(),
                            r.get_residue_number(),
                            r.get_ins_code()
                            )
                        );
                    }
                }
                
                write_to_file(&(outfile_pdb.clone()+"."+&sii.to_string()+".pdb"),res_str);
            }
        }
        
}


fn main_str_align(args:HashMap<String,String>) {
    
    if args.contains_key("-h"){
        eprintln!("usage:  rust_pdbloader.exe -structural_alignment -query <query_pdb> -template  <template_pdb> -type <max, raw, raw_full, sw (default)> [-out <output_file>] [-out_pdb <aligned_query_pdb>]");
    }
    //raw_full はカットオフを行わず全 CA についての RMSD を出す。
    
    let qfilename = args.get("-query").unwrap_or_else(|| panic!("-query was not found."));
    let tfilename = args.get("-template").unwrap_or_else(|| panic!("-template was not found."));
    let mut query_pdb:PDBEntry = load_pdb(qfilename);
    let template_pdb:PDBEntry = load_pdb(tfilename);
    let mut alignment_type_:String = args.get("-type").unwrap_or(&("".to_string())).clone();
    let result_file:String = args.get("-out").unwrap_or(&("".to_string())).clone(); 
    let outfile_pdb:String = args.get("-out_pdb").unwrap_or(&("".to_string())).clone(); 
    let tmscore_break:f64 = args.get("-tmscore_break").unwrap_or(&("0.0".to_string())).parse::<f64>().unwrap_or_else(|e| panic!("cannot parse tmscore_break {:?}",e)).clone(); 
    let mut query_chain:String  = args.get("-query_chain").unwrap_or(&("".to_string())).clone(); 
    let mut template_chain:String  = args.get("-template_chain").unwrap_or(&("".to_string())).clone(); 

    if alignment_type_ == ""{
        alignment_type_ = "sw".to_string();
    }

    //ToDo: 二次構造で SW して Seed はそれ以外使わないようにするモードの追加
    //seed_length と dist_limit をオプションで指定。
    let alignment_type:structural_alignment::AlignmentType = if alignment_type_ == "sw"{
        structural_alignment::AlignmentType::SW
    }else if alignment_type_ == "raw"{
        structural_alignment::AlignmentType::RAW
        
    }else if alignment_type_ == "raw_full"{
        structural_alignment::AlignmentType::RAW_FULLRMSD
        
    }else if alignment_type_ == "max"{
        structural_alignment::AlignmentType::MAXIMUM
    }else{
        panic!("{} is not recognized. sw, raw, max, raw_full",alignment_type_);
    };

    if query_chain.len() == 0{
        if query_pdb.chains.len() == 1{
            query_chain = query_pdb.chains[0].chain_name.clone();
        }else{
            panic!("Number of chains in query pdb must be one or use -query_chain option.");
        }
    }
    
    if template_chain.len() == 0{
        if template_pdb.chains.len() == 1{
            template_chain = template_pdb.chains[0].chain_name.clone();
        }else{
            panic!("Number of chains in template pdb must be one or use -templete_chain option.");
        }
    }


    let mut qresidues:Vec<&PDBResidue> = vec![];
    let mut tresidues:Vec<&PDBResidue> = vec![];
    for cc in query_pdb.chains.iter(){
        if cc.chain_name == query_chain{
            for rr in cc.residues.iter(){
                qresidues.push(rr);
            }
        }
    }
    
    for cc in template_pdb.chains.iter(){
        if cc.chain_name == template_chain{
            for rr in cc.residues.iter(){
                tresidues.push(rr);
            }
        }
    }
    let res_ = structural_alignment::align_pdb(&qresidues,&tresidues,alignment_type,tmscore_break);
    let res:structural_alignment::StructuralAlignmentResult = res_.unwrap();
    let mut res_str:Vec<String> = vec![];
    res_str.push(format!("query:{}",&qfilename));
    res_str.push(format!("template:{}",&tfilename));
    res_str.push(format!("initial_alignment_type:{}",&alignment_type_));

    res_str.push(format!("tmscore_chain1:{}",res.tmscore_chain1));
    res_str.push(format!("tmscore_chain2:{}",res.tmscore_chain2));
    res_str.push(format!("gdt_ts:{}",res.gdt_ts));
    res_str.push(format!("dist_cutoff:{}",res.dist_cutoff));
    res_str.push(format!("num_aligned:{}",res.num_aligned));
    res_str.push(format!("rmsd_aligned:{}",res.rmsd_aligned));

    res_str.push(format!("transform_matrix:====="));
    for ii in 0..4{
        res_str.push(format!("{} {} {} {}",res.transform_matrix[ii][0],res.transform_matrix[ii][1],res.transform_matrix[ii][2],res.transform_matrix[ii][3]));
    }
    res_str.push(format!("====="));

    let qali:String = res.aligned_chain1.iter().fold("".to_string(),|s,m|s+m);
    let tali:String = res.aligned_chain2.iter().fold("".to_string(),|s,m|s+m);

    res_str.push(format!("alignment:====="));
    res_str.push(qali);
    res_str.push(tali);
    res_str.push(format!("====="));

    if result_file.len() == 0{
        for ss in res_str.iter(){
            println!("{}",ss);
        }
    }else{
        write_to_file(&result_file,res_str);
    }
    if outfile_pdb.len() > 0{
        for cc in query_pdb.chains.iter_mut(){
            for rr in cc.residues.iter_mut(){
                for aa in rr.iter_mut_atoms(){
                    let mres = matrix_process::matrix_multi(&res.transform_matrix,&vec![vec![aa.get_x()],vec![aa.get_y()],vec![aa.get_z()],vec![1.0]]);
                    aa.set_xyz(mres[0][0],mres[1][0],mres[2][0]);
                }
            }
        }
        query_pdb.save(&outfile_pdb);
    }

    return;

}


fn main_prepare_structure(args:HashMap<String,String>) {
    if args.contains_key("-h"){
        eprintln!("usage:  rust_pdbloader.exe -prepare_structure -in <pdbfile> -out <outputfile> -resource_dir <directory contains top19.inp & param19.inp> [-random <float value>] [-random_seed <integer value>]");
    }

    let options_expected:HashSet<String>=vec![
        "-in"
        ,"-out"
        ,"-resource_dir"
        ,"-random"
        ,"-random_seed"
    ].iter().map(|m|m.to_string()).collect();

    for kk in args.iter(){
        if !options_expected.contains(kk.0){
            eprintln!("{} is not found in option dictionaly.",kk.0);
        }
    }


    //raw_full はカットオフを行わず全 CA についての RMSD を出す。
    
    let infilename:&str = args.get("-in").unwrap_or_else(|| panic!("-in was not found."));
    let resource_dir:&str = args.get("-resource_dir").unwrap_or_else(|| panic!("-resource_dir was not found."));
    let mut pdbb:PDBEntry = load_pdb(infilename);
    let parr = charmm_param::CHARMMParam::load_chamm19((resource_dir.to_string()+"\\toph19.inp").as_str(),(resource_dir.to_string()+"\\param19.inp").as_str());
    let random_movement:f64 = args.get("-random").unwrap_or(&("0.0".to_string())).parse::<f64>().unwrap_or_else(|_| panic!("can not parse {:?}",args.get("-random")));
    let outfilename = args.get("-out").unwrap_or_else(|| panic!("-out was not found."));
    
    let random_seed:u64 = args.get("-random_seed").unwrap_or(&("0".to_string())).parse::<u64>().unwrap_or_else(|_| panic!("can not parse {:?}",args.get("-random")));
    
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(random_seed);

    let mut lines_all:Vec<String> = vec![];
    let numchains:usize = pdbb.chains.len();
    for ll in 0..numchains{
        charmm_based_energy::MDAtom::change_to_charmmnames(&mut pdbb.chains[ll].residues);
    }
    let (mut md_envset,md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&pdbb.chains,&parr,true);
    let num_atoms = md_envset.atoms.len();
    charmm_based_energy::estimate_positions_unplaced(&mut md_envset,&md_varset);
    for ii in 0..num_atoms{
        
        if md_envset.atoms[ii].unplaced{
            eprintln!("Can not place {:?}",md_envset.atoms[ii]);
        }
        if random_movement > 0.0{
            md_envset.atoms[ii].x += rgen.gen_range(-random_movement,random_movement);
            md_envset.atoms[ii].y += rgen.gen_range(-random_movement,random_movement);
            md_envset.atoms[ii].z += rgen.gen_range(-random_movement,random_movement);
        }
    }

    for aa in md_envset.atoms.iter(){
        let (chainid,(resname,resnum,inscode),att) = aa.to_pdbatom();
        lines_all.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&inscode));
    }
    
    write_to_file(outfilename,lines_all);
}

fn _main_evoef() {
    let mut pdbb:PDBEntry = load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/6iws_model1_noh.pdb");
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    let outfilename = "test/evoef_hadded.pdb";
    let _rpp:rosetta_param::RosettaParamMapper  = rosetta_param::RosettaParamMapper::construct(debug_env::ROSETTA_DIR);

    let mut ca_atoms_pos:Vec<(usize,(f64,f64,f64))> = vec![];
    
    charmm_based_energy::MDAtom::change_to_charmmnames(&mut pdbb.chains[0].residues);

    let (mut md_envset,md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&vec![pdbb.chains.remove(0)],&parr,true);

    let mut ca_md:Vec<&charmm_based_energy::MDAtom> = vec![];
    for (_rii,aa) in md_envset.atoms.iter().enumerate(){
        if aa.atom_name == "CA"{
        //後のステップでCA が全ての AA について 1 つずつ存在し、Residue 順に Vector に加えられているとみなしている。
            ca_md.push(aa);
            
        }
    }

    for (_rii,rr) in pdbb.chains[0].residues.iter().enumerate(){
        if rr.residue_name == "HOH"{
            continue;
        }

        //CA が全ての AA について 1 つずつ存在し、Residue 順に Vector に加えられているとみなしている。
        //名前だけチェック。ずれていたら panic する。その場合もっと保証できる方法に変更。
        //CA distance の時のマップに必要
        assert_eq!(ca_md[_rii].residue_name,rr.residue_name);
        for (_aii,aa) in rr.iter_atoms().enumerate(){
            if aa.atom_code == "CA"{
                ca_atoms_pos.push((_rii,(aa.get_x(),aa.get_y(),aa.get_z())));
            }
        }
    }


    let mut ref_dist:Vec<Vec<f64>> = vec![vec![-1.0;ca_atoms_pos.len()];ca_atoms_pos.len()];
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    let alen:usize = ca_atoms_pos.len();
    let mut ca_indices:Vec<usize> = vec![];
    for ii in 0..alen{
        ca_indices.push(ca_atoms_pos[ii].0);
        for jj in 0..alen{
            let ddist = process_3d::distance(&ca_atoms_pos[ii].1,&ca_atoms_pos[jj].1);
            if  ddist < 12.0{
                ref_dist[ii][jj] = ddist;
            }
        }
    }

    let num_atoms = md_envset.atoms.len();
    println!("Params were loaded.");
    for ii in 0..num_atoms{
        if md_envset.atoms[ii].unplaced{
            md_envset.atoms[ii].x += rgen.gen_range(-0.5,0.5);
            md_envset.atoms[ii].y += rgen.gen_range(-0.5,0.5);
            md_envset.atoms[ii].z += rgen.gen_range(-0.5,0.5);
        }
    }
    charmm_based_energy::estimate_positions_unplaced(&mut md_envset,&md_varset);
    //mcprocess_md::mc_iter_array(&mut md_envset,&ref_dist,&ca_indices,50,Some(123),num_atoms);
    let mut lines:Vec<String> = vec![];
    for aa in md_envset.atoms.iter(){
        let (chainid,(resname,resnum,inscode),att) = aa.to_pdbatom();
        lines.push(att.get_pdb_atom_line_string(&chainid,&resname,resnum,&inscode));
    }
    write_to_file(outfilename,lines);
    let mut evoenv:evoef2_energy::EvoEF2Env = evoef2_energy::EvoEF2Env::new(md_envset,md_varset,debug_env::RESOURCE_DIR,false);
    let mut atom_ene:Vec<f64> = vec![0.0;evoenv.md_envset.atoms.len()];
    evoenv.md_envset.update_distance();
    println!("vdw: {:?}",evoef2_energy::calcEvdw(&evoenv,&mut atom_ene,&(vec![1.0_f64;4])));
    println!("elec: {}",evoef2_energy::calcEelec(&evoenv,&mut atom_ene,1.0));
    println!("hb: {:?}",evoef2_energy::calcEhb(&evoenv,&mut atom_ene,&(vec![1.0_f64;9])));
    println!("desolv: {:?}",evoef2_energy::calcEdesolv(&evoenv,&mut atom_ene,1.0,1.0));
    println!("disulfide: {}",evoef2_energy::calcEss(&evoenv,&mut atom_ene,1.0));
    let mut atom_level_energy:Vec<f64> = vec![0.0;evoenv.md_envset.atoms.len()];
    let mres = charmm_based_energy::calc_energy(&mut evoenv.md_envset,&mut evoenv.charmm_vars,&mut atom_level_energy);
    println!("{:?}",mres);
}



#[allow(unused_imports)]
use super::backbone_sample;
#[allow(unused_imports)]
use super::side_chain_sample;
#[allow(unused_imports)]
use super::chain_builder;
#[allow(unused_imports)]
use super::pdbdata;
#[allow(unused_imports)]
use super::debug_env;
#[allow(unused_imports)]
use super::charmm_based_energy;
#[allow(unused_imports)]
use super::charmm_param;
#[allow(unused_imports)]
use super::process_3d::*;
#[allow(unused_imports)]
use super::geometry::*;
#[allow(unused_imports)]
use std::collections::HashMap; 
#[allow(unused_imports)]
use regex::Regex;
#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(unused_imports)]
use std::fs::File;


lazy_static! {
    static ref REGEX_WS:Regex = Regex::new(r"[\s]").unwrap();
    static ref REGEX_NOLINE:Regex = Regex::new(r"^[\s]*$").unwrap();
    static ref REGEX_TAILBLANK:Regex = Regex::new(r"[\s]*$").unwrap();
    static ref REGEX_COMMENT_HEAD:Regex = Regex::new(r"[\s]+[!#].*").unwrap();
    static ref REGEX_COMMENT:Regex = Regex::new(r"[!#].*").unwrap();
    static ref REGEX_POSCODE:Regex = Regex::new(r"^([+\-*])(.+)").unwrap();
    static ref CHARMM19_TO_ROSETTA:HashMap<String,String> = HashMap::new();
    static ref PDBATOM_TO_ROSETTA:HashMap<String,String> = HashMap::new();
    
}

#[derive(Debug)]
pub struct RosettaAtomParam{
    pub partial_charge:f64,//使わない
    pub lj_radii:f64,
    pub lj_welldepth:f64,
    pub gfree:f64,//使わない
}

pub fn start_with(target:&str,fragment:&str)-> bool{
    if let Some(x) = target.find(fragment){
        if x == 0{
            return true;
        }
    }
    return false;
}
pub struct RosettaParamMapper{
    atom_pdb_to_rosetta:HashMap<String,String>,
    atom_partial_charge:HashMap<String,f64>,
    atom_radii:HashMap<String,f64>,
    atom_welldepth:HashMap<String,f64>,
    atom_gfree:HashMap<String,f64>,
}
impl RosettaParamMapper{
    pub fn construct(rosetta_param_dir:&str)->RosettaParamMapper{
        let atom_map_file:String = "".to_string()+rosetta_param_dir+"/charmm19_to_rosetta.dat";
        let amfile = File::open(atom_map_file).unwrap();
        let amfile_reader = BufReader::new(amfile);
        let mut atom_pdb_to_rosetta:HashMap<String,String> = HashMap::new();
        
        for (_lcount,line_) in amfile_reader.lines().enumerate(){
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT.replace_all(&line_, "")).to_string();
            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }

            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            atom_pdb_to_rosetta.insert(ptt[0].clone()+"_"+ptt[1].as_str(),ptt[3].clone());
        }
        
        let mut atom_partial_charge:HashMap<String,f64> = HashMap::new();
        
        let mut atom_radii:HashMap<String,f64> = HashMap::new();
        let mut atom_welldepth:HashMap<String,f64> = HashMap::new();
        let mut atom_gfree:HashMap<String,f64> = HashMap::new();

        let radiicol:usize = 1;
        let wdcol:usize = 4;
        let gfreecol:usize = 7;
        
        //Supplementary Information: Simultaneous optimization of biomolecular energy function on features from small molecules and macromolecules
        //Table S3 
        let lj_param_file_:String = "".to_string()+rosetta_param_dir+"/opt_nov15_lj_param.dat";
        let lj_param_file = File::open(lj_param_file_).unwrap();
        let ljfile_reader = BufReader::new(lj_param_file);
        
        for (_lcount,line_) in ljfile_reader.lines().enumerate(){
            if _lcount == 0{
                continue;
            }
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT.replace_all(&line_, "")).to_string();
            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }

            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if let Ok(x) = ptt[radiicol].parse::<f64>(){
                atom_radii.insert(ptt[0].clone(),x);
            }else{
                if ptt[radiicol] == "NODATA"{
                }else{
                    eprintln!("{} was not parsed.",ptt[radiicol]);
                }
            }
            if let Ok(x) = ptt[wdcol].parse::<f64>(){
                atom_welldepth.insert(ptt[0].clone(),x);
            }else{
                if ptt[gfreecol] == "NODATA"{
                }else{
                    eprintln!("{} was not parsed.",ptt[wdcol]);
                }
            }
            if let Ok(x) = ptt[gfreecol].parse::<f64>(){
                atom_gfree.insert(ptt[0].clone(),x);
            }else{
                if ptt[gfreecol] == "NODATA"{
                }else{
                    eprintln!("{} was not parsed.",ptt[gfreecol]);
                }
            }
        }

        //Supplementary Information: Simultaneous optimization of biomolecular energy function on features from small molecules and macromolecules
        //Table S4 
        let pt_file_:String = "".to_string()+rosetta_param_dir+"/opt_nov15_partial_charge.dat";
        let pt_file = File::open(pt_file_).unwrap();
        let ptfile_reader = BufReader::new(pt_file);
        
        for (_lcount,line_) in ptfile_reader.lines().enumerate(){
            if _lcount == 0{
                continue;
            }
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT.replace_all(&line_, "")).to_string();
            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }

            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            atom_partial_charge.insert(
                ptt[0].clone(),
                ptt[1].parse::<f64>().unwrap_or_else(|_|{panic!("File parsing error! {}",line)})
            );
        }

        return RosettaParamMapper{
            atom_pdb_to_rosetta,
            atom_partial_charge,
            atom_radii,
            atom_welldepth,
            atom_gfree,
        };
    }

    pub fn get_rosetta_param(&self,atom:&charmm_based_energy::MDAtom)->RosettaAtomParam{
        let default_value:f64 = 999999.0;
        let mut lj_welldepth = default_value;
        let mut lj_radii = default_value;
        let mut partial_charge = default_value;
        let mut gfree = default_value;

        let mut rescode:String = atom.residue_name.clone().clone()+"_"+atom.atom_name.as_str();
        //アンダーバーが無いだけ
        let mut rescode2:String = atom.residue_name.clone().clone()+atom.atom_name.as_str();
        
        if  atom.residue_name == "HSD"{//CHARMM の HSD は ROSETTA の HIS っぽい。。。うーん？？？
            rescode2 = "HIS".to_string()+atom.atom_name.as_str();
        }

        if atom.nterminal{
            if self.atom_pdb_to_rosetta.contains_key(("%NTER".to_string()+"_"+atom.atom_name.as_str()).as_str()){
                rescode = "%NTER".to_string()+"_"+atom.atom_name.as_str();
            }
        }
        if atom.cterminal{
            if self.atom_pdb_to_rosetta.contains_key(("%CTER".to_string()+"_"+atom.atom_name.as_str()).as_str()){
                rescode = "%CTER".to_string()+"_"+atom.atom_name.as_str();
            }
        }

        if let Some(x) = self.atom_pdb_to_rosetta.get(&rescode){
            lj_radii = *self.atom_radii.get(x).unwrap_or_else(||{eprintln!("Vdw radii for {} is not defined.",x);return &default_value;});
            lj_welldepth = *self.atom_welldepth.get(x).unwrap_or_else(||{eprintln!("Well depth for {} is not defined.",x);return &default_value;});
            gfree = *self.atom_gfree.get(x).unwrap_or_else(||{eprintln!("Gfree for {} is not defined.",x);return &default_value;});
        }else{
            //ToDo なんかここでエラーメッセージが出るが正しかっただろうか
            eprintln!("{} is not in the table S3.",rescode);
        }
        
        if let Some(x) = self.atom_partial_charge.get(&rescode2){
            partial_charge = *x;
        }else{
            eprintln!("{} is not in the table S4.",rescode2);
        }


        return RosettaAtomParam{
        lj_welldepth,
        lj_radii,
        partial_charge,
        gfree,  
        };
    }

}


pub fn get_atom_code(residue_name:&str,atom_name_type:&(&str,&str))-> String{
    let atom_name = atom_name_type.0;
    let atom_type = atom_name_type.1;
    let mut pcode = "NONE";

    if atom_name == "CA"{
        pcode = "CAbb";
    }else if atom_name == "C"{
        pcode = "CObb";
    }else if atom_name == "N"{
        pcode = "Nbb";
    }else if atom_name == "O"{
        pcode = "OCbb";
    }else if atom_name == "H"{
        pcode = "HNbb";
    }else if start_with(atom_type,"C"){
        if residue_name == "GLN" && atom_name == "CD"{
            pcode = "CNH2";//だと思うがはっきりした記述が見つけられなかった
        }else if atom_type == "CH1E"{
            pcode = "CH1";
        }else if atom_type == "CH2E"{
            pcode = "CH2";
        }else if atom_type == "CH3E"{
            pcode = "CH3";
        }else if atom_type == "C"{
            if atom_name == "CG" &&  residue_name == "ASP"{
                pcode = "COO";
            }else if atom_name == "CD" &&  residue_name == "GLU"{
                pcode = "COO";
            }else{
                pcode = "CH0";
            }
        }else if atom_type == "CR1E"{
            pcode = "aroC";
        }else{
            println!("{} {} {} is not defined.",residue_name,atom_name,atom_type);
            pcode = "C%_param19";
        }
    }else if start_with(atom_type,"N"){
        if atom_name == "NE2" &&  residue_name  == "GLN" {
            pcode = "NH2O";
        }else if residue_name == "HIS" || residue_name == "HSD" {
            pcode = "Nhis";
        }else if residue_name == "LYS" {
            pcode = "Nlys";
        }else if residue_name == "PRO" {
            pcode = "Npro";
        }else if residue_name == "TRP"  && atom_name == "NE1" {
            pcode = "Ntrp";
        }else if residue_name == "ARG"  && atom_name == "NE" {
            pcode = "NtrR";
        }else if residue_name == "ARG"  && (atom_name == "NH1" || atom_name == "NH2"){
            pcode = "Narg";
        }else{
            println!("{} {} {} is not defined.",residue_name,atom_name,atom_type);
            pcode = "N_consensus";
        }
    }else if start_with(atom_type,"O"){
        if (atom_name == "OD1" || atom_name == "OD2") &&  residue_name  == "ASP" {
                pcode = "OOC";
        }else if (atom_name == "OE1" || atom_name == "OE2") &&  residue_name  == "GLU" {
                pcode = "OOC";
        }else if atom_name == "OG" && residue_name  == "SER" {
                pcode = "OH";
        }else if atom_name == "OH" && residue_name  == "TYR" {
                pcode = "OH";
        }else if atom_name == "OE1" && residue_name  == "GLN" {
                pcode = "ONH2";
        }else{
            println!("{} {} {} is not defined.",residue_name,atom_name,atom_type);
            pcode = "O*_param19";
        }
    }else if start_with(atom_type,"S"){
        if atom_name == "SD" &&  residue_name  == "MET" {
                pcode = "S";
        }else if atom_name == "SG" &&  residue_name  == "CYS" {
                pcode = "SH";
        }else{
            println!("{} {} {} is not defined.",residue_name,atom_name,atom_type);
            pcode = "S*_param19";
        }
    }else if start_with(atom_type,"H"){
        if atom_type == "HC"{
            pcode = "Hpol";
        }else{
            //他の H は来ることは想定していない
            println!("{} {} {} is not defined.",residue_name,atom_name,atom_type);
            pcode = "H_param19";
        }
    }
    return pcode.to_string();
}



#[test]
fn rosetta_assign_test(){
    
    let allaa:Vec<String> = vec![
        "ALA","ALA","ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","VAL","VAL",
    ].iter().map(|m|m.to_string()).collect();

    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);
    let mut ress:Vec<pdbdata::PDBComp> = chain_builder::build_dirty_chain(&allaa,&bset,&sset);
    let mut acount:i64 = 1;
    //let mut pdbstr:Vec<String> = vec![];
    for (ii,rr) in ress.iter_mut().enumerate(){
        rr.set_residue_number((ii+1) as i64);
        for aa in rr.iter_mut_atoms(){
            aa.set_serial_number(acount);
            acount += 1;
        }
    }
    let parr = charmm_param::CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    let mut chain:pdbdata::PDBAsym = pdbdata::PDBAsym::new("A");
    for rr in ress.into_iter(){
        chain.add_residue(rr,true);
    }

    charmm_based_energy::MDAtom::change_to_charmmnames(&mut chain.residues);
    let (md_envset,_md_varset):(charmm_based_energy::CharmmEnv,charmm_based_energy::CharmmVars) = charmm_based_energy::MDAtom::chain_to_atoms(&vec![chain],&parr,true);
    for aa in md_envset.atoms.iter(){
        println!("{} {} {} {}",aa.residue_name,aa.atom_name,aa.atom_type
        ,get_atom_code(aa.residue_name.as_str(), &(aa.atom_name.as_str(),aa.atom_type.as_str())));
    }
    let rpp:RosettaParamMapper  = RosettaParamMapper::construct(debug_env::ROSETTA_DIR);
    for (_ii,aa) in md_envset.atoms.iter().enumerate(){
        let rpar = rpp.get_rosetta_param(aa);
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
        ,aa.residue_name
        ,aa.atom_name
        ,aa.nb_r1_2
        ,aa.nb_epsilon
        ,rpar.lj_welldepth,rpar.lj_radii
        ,aa.charge
        ,rpar.gfree
        );   
     }
}



use std::fs;
use regex::Regex;

use crate::{mmcif_process, pdbdata::PDBEntry};


pub fn generate_intermediate_files(dirname:&str){
    let paths = fs::read_dir(dirname).unwrap();
    let exx =  Regex::new(r"(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    let mut entries_:Vec<String> = vec![];
    for path in paths {
        if let Ok(a) = path{
            if let Some(b) = a.path().to_str(){
                if let Some(x) = exx.find(b){
                    entries_.push(b.to_string());
                }
            }
        }
    }
    let mut entries:Vec<PDBEntry> = vec![];
    for ee in entries_.into_iter(){
        if let Some(x) = exx.captures(&ee){
            let ext1:String = x.get(1).unwrap().as_str().to_string();
            let is_gzip = 
            if let Some(_) = x.get(2){
                true
            }else{
                false
            };

            if &ext1 == ".ent" || &ext1 == ".pdb"{
                entries.push(mmcif_process::load_pdb(&ee,is_gzip));
            }else if &ext1 == ".cif"{
                entries.push(mmcif_process::MMCIFEntry::load_mmcif(&ee,is_gzip));

            }
        }
    }
    ここから
    Residue を全部抽出する
}

#[test]
fn dirloadtest(){
    generate_intermediate_files("example_files");
}

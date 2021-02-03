
use ff_env::{FFAtom, FFBond};
use roxmltree;
use std::collections::HashSet;
use super::debug_env::*;
use super::openff_energy;
use super::ff_env;
use std::collections::HashMap;

pub struct OpenFFEnergy{
}
impl OpenFFEnergy{

//ToDo: IMPROPER をとる奴は作ってなかった気がする
    pub fn assign_energy_terms(atoms:&Vec<&FFAtom>,bonds:&Vec<&FFBond>){
        let angle_candidate:Vec<Vec<usize>> = autogenerate_connection(&bonds,3);
        for bb in angle_candidate.into_iter(){

        }
    }

    pub fn load(filename:&str)->OpenFFEnergy{
        let  ret:OpenFFEnergy = OpenFFEnergy{

        };

        let text = std::fs::read_to_string(filename).unwrap();
        let doc = match roxmltree::Document::parse(&text) {
            Ok(v) => v,
            Err(e) => {
                println!("Error: {}.", e);
                std::process::exit(1);
            }
        };

        let mut bonds_:Vec<roxmltree::Node> = vec![];
        let mut angles_:Vec<roxmltree::Node> = vec![];
        let mut torsions_:Vec<roxmltree::Node> = vec![];
        let mut imprs_:Vec<roxmltree::Node> = vec![];
        let mut atoms_vdw_:Vec<roxmltree::Node> = vec![];

        let child_elements = doc.root_element().children().filter(|n|n.node_type()==roxmltree::NodeType::Element);
        for n in child_elements.clone() {

            //Python 実装
            //https://github.com/openforcefield/openff-toolkit/blob/de8a4a545351301adfe424dff0d879b2dd13bc0b/openff/toolkit/typing/engines/smirnoff/parameters.py
            //BondHandler.create_forces
            //を見ると
            //三個ある場合は最初の二つの Bond に対するパラメータになるように見える
            if n.tag_name().name() == "Bonds"{
                let b = n.children().filter(|n|n.node_type()==roxmltree::NodeType::Element);;
                for aa in b{
                    bonds_.push(aa);
                }
            }
            if n.tag_name().name() == "Angles"{
                let b = n.children().filter(|n|n.node_type()==roxmltree::NodeType::Element);;
                for aa in b{
                    angles_.push(aa);
                }
            }
            
            if n.tag_name().name() == "ProperTorsions"{
                let b = n.children().filter(|n|n.node_type()==roxmltree::NodeType::Element);;
                for aa in b{
                    torsions_.push(aa);
                }
            }
            
            if n.tag_name().name() == "ImproperTorsions"{
                let b = n.children().filter(|n|n.node_type()==roxmltree::NodeType::Element);;
                for aa in b{
                    imprs_.push(aa);
                }
            }
            
            if n.tag_name().name() == "vdW"{
                let b = n.children().filter(|n|n.node_type()==roxmltree::NodeType::Element);;
                for aa in b{
                    atoms_vdw_.push(aa);
                }
            }
        }

        return ret;
    }
}




//原子と原子を繋ぐ Bonds の配列と原子数を渡すと
//原子数分 Bond で繋がる原子のパスをすべて返す
//逆方向順方向はチェックし、同じパスが既にある場合は含めない
//同じ原子を二回通るパスは含めない
pub fn autogenerate_connection(bonds:&Vec<&FFBond>,num_nodes:usize)->Vec<Vec<usize>>{
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


#[test]
fn openff_loadtest(){
    let r:OpenFFEnergy = OpenFFEnergy::load(&(RESOURCE_DIR.to_string()+"/openff/smirnoff99frosst/smirnoff99Frosst-1.1.0.offxml"));

}
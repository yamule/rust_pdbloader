#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use std::collections::HashMap;

use super::pdbdata::*;
use super::geometry::*;
use super::misc_util::*;

#[allow(unused_imports)]
use super::atom_id_map;

#[allow(dead_code,unused_imports)]
use super::debug_env;

#[allow(dead_code,unused_imports)]
use super::process_3d::*;
use std::fs::File;

#[allow(dead_code,unused_imports)]
use std::f64::consts::PI;


pub const N:i64 = 0;
pub const CA:i64 = 1;
pub const C:i64 = 2;
use std::sync::Mutex;



lazy_static! {
    static ref BACKBONE_ATOM_ID:Mutex<HashMap<String,i64>> = Mutex::new(HashMap::new());
    static ref PREPARED:Mutex<bool> = Mutex::new(false);
}


pub fn prepare(){
    BACKBONE_ATOM_ID.lock().unwrap().insert("N".to_string(), N);
    BACKBONE_ATOM_ID.lock().unwrap().insert("C".to_string(), C);
    BACKBONE_ATOM_ID.lock().unwrap().insert("CA".to_string(), CA);
    *PREPARED.lock().unwrap() = true;
}

pub fn get_backbone_atom_index(s:&str) -> i64 {
    if ! *PREPARED.lock().unwrap(){
        prepare();
    }
    if let Some(x) = BACKBONE_ATOM_ID.lock().unwrap().get(s){
        return *x;
    }else{
        return -1;
    }
}


pub struct SideChainSet {
    sidechains:HashMap<String,Vec<SideChainSample>>,
    conformer_dir:String
}

impl SideChainSet{
    pub fn new(conformer_dir:&str)->SideChainSet{
        let aaname:Vec<&str> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        ];
        let mut sidechains:HashMap<String,Vec<SideChainSample>> = HashMap::new();
        for a in aaname.iter(){
            sidechains.insert(a.to_string(),SideChainSample::load((conformer_dir.to_string()+"/"+a+".rotamers.dat").as_str()));
        }
        return SideChainSet{sidechains:sidechains,conformer_dir:conformer_dir.to_string()};
    }

    
    pub fn get_num(&self,sname:&str)-> usize{
        return self.sidechains.get(sname).expect((sname.to_string()+" not found.").as_str()).len();
    }

    pub fn get_sample_at(&self,sname:&str,i:usize)-> &SideChainSample{
        return &self.sidechains.get(sname).expect((sname.to_string()+" not found.").as_str())[i];
    }
    
    pub fn dir_name(&self)->&str{
        return &self.conformer_dir;
    }
}

pub struct SideChainSample{
	pub residue_name:String,
	pub prob:f64,
	pub count:usize,
    pub atoms:HashMap<String,PDBAtom>,
    pub atoms_backbone:Vec<PDBAtom>
	//atoms_scoring:HashMap<String,PDBAtom>//scoring に使用されるセットがあったが、スコアごとにフィルタリングすることにしよう
}
impl SideChainSample{
    pub fn new()->SideChainSample{
        return SideChainSample{
            residue_name:"UNK".to_string()
            ,prob:0.0
            ,count:0
            ,atoms:HashMap::new()
            ,atoms_backbone:vec![]
        }
    }

    //回転して他の Residue に適用することが前提なのでそのものを変更できない。
    //コピーを作るための関数。
    pub fn get_point_copy(&self)->HashMap<String,Point3D>{
        let mut ret:HashMap<String,Point3D> = HashMap::new();
        for (ss,vv) in self.atoms.iter(){
            let a = vv.get_xyz();
            ret.insert(ss.to_string(),Point3D::new(a.0,a.1,a.2));
        }
        for vv in self.atoms_backbone.iter(){
            let a = vv.get_xyz();
            ret.insert(vv.atom_code.clone(),Point3D::new(a.0,a.1,a.2));
        }
        
        return ret;
    }
	pub fn parse_block(block:&Vec<String>)->SideChainSample{
        let mut ret = SideChainSample::new();
        let mut atoms_buff:HashMap<String,PDBAtom> = HashMap::new();

        for ss in block.iter(){
            let caps = PAT_LABEL.captures(ss.as_str());
            match caps{
                Some(cc)=>{
                    let g = cc.get(1).map_or("", |m| m.as_str()).to_lowercase();
                    let vv = cc.get(2).map_or("", |m| m.as_str());
                    if g.as_str() == "residue" {
                        ret.residue_name = vv.to_string();
                    }else if g.as_str() == "count"{
                        ret.count = vv.parse::<usize>().expect(("Error in parsing count ".to_string()+ss).as_str());
                    }else if start_with(g.as_str(),"atom"){
                        let mapp:HashMap<String,String> = line_to_hash(vv);
                        let atomname:String = mapp.get("name").expect("name not found.").to_string();
                        let mut atom:PDBAtom = PDBAtom::new();
                        atom.set_xyz(
                            mapp.get("x").unwrap().parse::<f64>().expect(("x parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("y").unwrap().parse::<f64>().expect(("y parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("z").unwrap().parse::<f64>().expect(("z parse faild ".to_string()+" "+ss).as_str())					
                        );
                        atom.atom_code = mapp.get("name").unwrap().to_string();
                        if mapp.contains_key("atom_symbol"){
                            atom.atom_symbol = mapp.get("atom_symbol").unwrap().to_string();
                        }else{
                            atom.atom_symbol = mapp.get("name").unwrap()[0..1].to_string();
                        }
                        atoms_buff.insert(atomname,atom);
                    }
                },
                None => {
                    if start_with(ss,"###==="){

                    }else{
                        eprintln!("{} was not parsed.",ss);
                    }
                }
        	}
        }
        let mut vectmp:Vec<Option<PDBAtom>> =vec![None,None,None];
        for (ss,vv) in atoms_buff.into_iter(){
            if ss == "CA"{
                vectmp[CA as usize] = Some(vv);
            }else if ss == "N"{
                vectmp[N as usize] = Some(vv);

            }else if ss == "C"{
                vectmp[C as usize] = Some(vv);
            }else{
                ret.atoms.insert(ss,vv);
            }
        }

        ret.atoms_backbone = vectmp.into_iter().map(|m|m.expect("This side chain lacks backbone atom!")).collect();
		return ret;
	}
	pub fn load(filename:&str)->Vec<SideChainSample>{
        
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<SideChainSample> = Vec::new();
        let mut buff:Vec<String> = vec![];
        for (_lcount,line) in reader.lines().enumerate() {
            if start_with(&line.as_ref().unwrap(),"//"){
                ret.push(SideChainSample::parse_block(&buff));
                buff.clear();
            }else{
                buff.push(line.unwrap());
            }
        }
        let mut count_all:usize = 0;
        for bs in ret.iter_mut(){
            count_all += bs.count;
        }
        if count_all > 0{
            for bs in ret.iter_mut(){
                bs.prob = bs.count as f64/count_all as f64;
            }
        }
        ret.sort_by(|a,b|a.prob.partial_cmp(&b.prob).expect("partialcmp failed"));
        ret.reverse();

        return ret;

	}
    pub fn get_prob(&self)->f64{
		return self.prob;
	}
	
	
}



#[test]
fn load_sidechain_test(){
    let filename = format!("{}{}",debug_env::ROTAMER_DIR,"\\ALA.rotamers.dat");
    let sidechains = SideChainSample::load(&filename);
    assert!(sidechains.len() > 0);
    let sset = SideChainSet::new(&debug_env::ROTAMER_DIR);
    assert!(sset.sidechains.len() > 0);

}
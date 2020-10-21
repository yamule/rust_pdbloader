#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use std::collections::HashMap;
use super::pdbdata::*;
use super::geometry::*;
use super::misc_util::*;


#[allow(dead_code,unused_imports)]
use super::process_3d::*;
use std::fs::File;
use std::f64::consts::PI;



pub const N:i64 = 0;
pub const CA:i64 = 1;
pub const C:i64 = 2;
pub const O:i64 = 3;
pub const PREVC:i64 = 4;
pub const NEXTN:i64 = 5;
use std::sync::Mutex;



lazy_static! {
    static ref ATOM_ID:Mutex<HashMap<String,i64>> = Mutex::new(HashMap::new());
    static ref PREPARED:Mutex<bool> = Mutex::new(false);
}

pub fn prepare(){
    ATOM_ID.lock().unwrap().insert("N".to_string(), N);
    ATOM_ID.lock().unwrap().insert("C".to_string(), C);
    ATOM_ID.lock().unwrap().insert("CA".to_string(), CA);
    ATOM_ID.lock().unwrap().insert("O".to_string(), O);
    ATOM_ID.lock().unwrap().insert("PREVC".to_string(), PREVC);
    ATOM_ID.lock().unwrap().insert("NEXTN".to_string(), NEXTN);
    *PREPARED.lock().unwrap() = true;

}


pub fn get_atom_index(s:&str) -> i64 {
    if ! *PREPARED.lock().unwrap(){
        prepare();
    }
    if let Some(x) = ATOM_ID.lock().unwrap().get(s){
        return *x;
    }else{
        return -1;
    }
}





pub struct BackboneSample{
	pub residue_name:String,
	pub prob:f64,
	pub count:usize,
	pub phi:f64,
    pub psi:f64,
    pub atoms:Vec<PDBAtom>,
    pub prev_omega:Option<OmegaSet>,
    pub next_omega:Option<OmegaSet>,
}

impl BackboneSample{
    pub fn new()->BackboneSample{
        return BackboneSample{
            residue_name:"UNK".to_string(),
            prob:0.0,
            count:0,
            phi:0.0,
            psi:0.0,
            atoms:vec![
            PDBAtom::new(),
            PDBAtom::new(),
            PDBAtom::new(),
            PDBAtom::new(),
            PDBAtom::new(),
            PDBAtom::new()
            ],
            prev_omega:None,
            next_omega:None,
        };
    }
    pub fn parse_block(block:&Vec<String>)->BackboneSample{

        let mut ret:BackboneSample = BackboneSample::new();
        for ss in block.iter(){
            let caps = PAT_LABEL.captures(ss.as_str());
            match caps{
                Some(cc)=>{
                    let g = cc.get(1).map_or("", |m| m.as_str()).to_lowercase();
                    let vv = cc.get(2).map_or("", |m| m.as_str());
                    if g == "residue" {
                        ret.residue_name = cc.get(2).unwrap().as_str().to_string();
                    }else if g == "code"{
                        let mapp = line_to_hash(vv);
                        ret.phi = mapp.get("phi").unwrap().parse::<f64>().expect(("phi parse faild ".to_string()+" "+ss).as_str());
                        ret.psi = mapp.get("psi").unwrap().parse::<f64>().expect(("psi parse faild ".to_string()+" "+ss).as_str());
                        
                    }else if g == "count"{
                        ret.count = cc.get(2).unwrap().as_str().parse::<usize>().expect(("count parse faild ".to_string()+" "+ss).as_str());
                    }else if g == "prevc"{
                        let mapp = line_to_hash(vv);
                        let index:usize = get_atom_index("PREVC") as usize;

                        ret.atoms[index] = PDBAtom::new();
                        ret.atoms[index].set_xyz(
                            mapp.get("x").unwrap().parse::<f64>().expect(("x parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("y").unwrap().parse::<f64>().expect(("y parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("z").unwrap().parse::<f64>().expect(("z parse faild ".to_string()+" "+ss).as_str())
                        );

                        ret.atoms[index].atom_symbol = "C".to_string();
                        ret.atoms[index].atom_code = "C".to_string();
                    }else if g == "nextn"{
                        let mapp = line_to_hash(vv);
                        let index:usize = get_atom_index("NEXTN") as usize;                        
                        ret.atoms[index] =PDBAtom::new();
                        ret.atoms[index].set_xyz(
                            mapp.get("x").unwrap().parse::<f64>().expect(("x parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("y").unwrap().parse::<f64>().expect(("y parse faild ".to_string()+" "+ss).as_str()),
                            mapp.get("z").unwrap().parse::<f64>().expect(("z parse faild ".to_string()+" "+ss).as_str())
                        );

                        ret.atoms[index].atom_symbol = "N".to_string();
                        ret.atoms[index].atom_code = "N".to_string();
                    }else{
                        if let Some(x) = g.find("atom"){
                            if x == 0{
                                let mapp = line_to_hash(vv);
                                let index:i64 = get_atom_index(mapp.get("name").unwrap());
                                if index == -1 {
                                    eprintln!("{}",mapp.get("name").unwrap().to_string()+" was not considered as backbone atom.");
                                }else{

                                    ret.atoms[index as usize] = PDBAtom::new();
                                    ret.atoms[index as usize].set_xyz(            
                                        mapp.get("x").unwrap().parse::<f64>().expect(("x parse faild ".to_string()+" "+ss).as_str()),
                                        mapp.get("y").unwrap().parse::<f64>().expect(("y parse faild ".to_string()+" "+ss).as_str()),
                                        mapp.get("z").unwrap().parse::<f64>().expect(("z parse faild ".to_string()+" "+ss).as_str())
                                    );
                                    ret.atoms[index as usize].atom_code = mapp.get("name").unwrap().to_string();
                                    if mapp.contains_key("atom_symbol"){
                                        ret.atoms[index as usize].atom_symbol = mapp.get("atom_symbol").unwrap().to_string();
                                    }else{
                                        ret.atoms[index as usize].atom_symbol = mapp.get("name").unwrap()[0..1].to_string();
                                    }
                                }
                            }
                        }
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
        return ret;
    }
    pub fn get_prob(&self)-> f64{
        return self.prob;
    }
    pub fn parse_omega_block(block:&Vec<String>)->Vec<OmegaSet>{
        let mut ret:Vec<OmegaSet> = vec![];
        let mut val:Vec<f64> = vec![];
        let mut count:Vec<usize> = vec![];
        let mut val_prev:Vec<f64> = vec![];
        let mut count_prev:Vec<usize> = vec![];
        for ss in block.iter(){
            let caps = PAT_LABEL.captures(ss.as_str());
            match caps{
                Some(cc)=>{
                    let g = cc.get(1).map_or("", |m| m.as_str()).to_lowercase();
                    let vv = cc.get(2).unwrap().as_str();
                    if g  == "prevomega"{
                        let mapp:HashMap<String,String> = line_to_hash(vv);
                        val_prev.push((mapp.get("prevomega").unwrap().parse::<f64>().expect(("prevomega parse faild ".to_string()+" "+ss).as_str()))/180.0*PI);
                        count_prev.push(mapp.get("count").unwrap().parse::<usize>().expect(("count parse faild ".to_string()+" "+ss).as_str()));
                    }else{
                        let mapp:HashMap<String,String> = line_to_hash(vv);
                        val.push((mapp.get("omega").unwrap().parse::<f64>().expect(("omega parse faild ".to_string()+" "+ss).as_str()))/180.0*PI);
                        count.push(mapp.get("count").unwrap().parse::<usize>().expect(("count parse faild ".to_string()+" "+ss).as_str()));
                    }
                },
                _=>{
                  eprintln!("{} was not parsed.",ss);  
                }
            }
        }
        
        if val_prev.len() > 0{
            let os:OmegaSet =OmegaSet::new(val_prev,count_prev,true);
            ret.push(os);
        }
        
        if val.len() > 0{
            let os:OmegaSet =OmegaSet::new(val,count,false);
            ret.push(os);
        }
        return ret;
    }

    
    pub fn load(filename:&str)->Vec<BackboneSample>{
        
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<BackboneSample> = Vec::new();
                let mut buff:Vec<String> = vec![];
        let mut prev_omega:Option<OmegaSet> = None;
        let mut next_omega:Option<OmegaSet> = None;

        for (_lcount,line) in reader.lines().enumerate() {
            if start_with(&line.as_ref().unwrap(),"//"){
                if start_with(&buff[0],"###prevomega") || start_with(&buff[0],"###nextomega"){
                    let al:Vec<OmegaSet> = BackboneSample::parse_omega_block(&buff);
                    
                    for os in al.into_iter(){
                        if os.prev {
                            prev_omega = Some(os);
                        }else{
                            next_omega = Some(os);
                        }
                    }

                }else{
                    ret.push(BackboneSample::parse_block(&buff));
                }
                buff.clear();
            }else{
                buff.push(line.unwrap());
            }
        }
        let mut count_all:usize = 0;
        if let None = prev_omega{
            panic!("can not read omega block! {} ",filename);
        }
        if let None = next_omega{
            panic!("can not read omega block! {} ",filename);
        }

        prev_omega.as_mut().unwrap().sort();
        next_omega.as_mut().unwrap().sort();

        for bs in ret.iter_mut(){
            count_all += bs.count;
            bs.prev_omega = Some(prev_omega.as_ref().unwrap().clone());
            bs.next_omega = Some(next_omega.as_ref().unwrap().clone());
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
}


/*public static void main(String[] args){
    InputStream iss = BackBoneSample.class.getResourceAsStream("resources/sampledresidues/ALA.backbones.dat");
    if(iss == null){
        System.out.println("//");
    }
    load(iss);
}*/


#[derive(Clone)]
pub struct OmegaSample{
    value:f64,
    count:usize,
    prob:f64
}

#[derive(Clone)]
pub struct OmegaSet{
    prev:bool,
    omegas:Vec<OmegaSample>
}
impl OmegaSet{
    //ω を多い順に並べ替える

    pub fn sort(&mut self){
        self.omegas.sort_by(|a,b|a.prob.partial_cmp(&b.prob).expect("omega sort failed."));
        self.omegas.reverse();
    }

    fn new(v:Vec<f64>,c:Vec<usize>,isprev:bool) -> OmegaSet{
        assert_eq!(v.len(),c.len());
        let mut ret:OmegaSet=OmegaSet{
            prev:isprev,
            omegas:vec![]
        };
        for ii in 0..v.len(){
            let mut os:OmegaSample = OmegaSample{value:0.0,count:0,prob:0.0};
            os.value = v[ii];
            os.count = c[ii];
            ret.omegas.push(os);
        }
        ret.calc_prob();
        ret.sort();
        return ret;
    }
    
    pub fn filt_count(&mut self,threshold:usize){
        let mut vve:Vec<OmegaSample> = vec![];
        vve.append(&mut self.omegas);
        for vv in vve.into_iter(){
            if vv.count >= threshold{
                self.omegas.push(vv);
            }
        }
        self.calc_prob();
    }

    pub fn calc_prob(&mut self){
        let mut sum :usize = 0;
        let osiz:usize = self.omegas.len();
        for os in self.omegas.iter(){
            sum += os.count;
        }

        if sum == 0 {
            for os in self.omegas.iter_mut(){
                os.prob = 1.0/osiz as f64;
            }
        }else{
            for os in self.omegas.iter_mut(){
                os.prob = os.count as f64/sum as f64;
            }
        }
    }
}


pub struct BackboneSet {
    backbones:HashMap<String,Vec<BackboneSample>>,
    conformer_dir:String
}


impl BackboneSet{
    pub fn new(conformer_dir:&str)->BackboneSet{
        let aaname:Vec<&str> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        ];
        let mut backbones:HashMap<String,Vec<BackboneSample>> = HashMap::new();
        for a in aaname.iter(){
            backbones.insert(a.to_string(),BackboneSample::load((conformer_dir.to_string()+"/"+a+".backbones.dat").as_str()));
        }
        return BackboneSet{backbones:backbones,conformer_dir:conformer_dir.to_string()};
    }
    
    pub fn has_sample(&self,sname:&str)-> bool{
        return self.backbones.contains_key(sname);
    }
    pub fn get_num(&self,sname:&str)-> usize{
            return self.backbones.get(sname).expect((sname.to_string()+" not found.").as_str()).len();
    }

    pub fn get_sample_at(&self,sname:&str,i:usize)-> &BackboneSample{
        return &self.backbones.get(sname).expect((sname.to_string()+" not found.").as_str())[i];
    }

    pub fn dir_name(&self)->&str{
        return &self.conformer_dir;
    }
}

#[test]
fn load_backbone_test(){
    prepare();
    let filename = "D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\pepbuilderj\\resources\\sampledresidues\\ALA.backbones.dat";
    let backbones = BackboneSample::load(filename);
    assert!(backbones.len() > 0);
    let bset = BackboneSet::new("D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\pepbuilderj\\resources\\sampledresidues\\");
    assert!(bset.backbones.len() > 0);
}

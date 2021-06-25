use core::f64;
use core::num;
use std::fs;
use std::cmp::Ordering;
use std::fs::File;
use std::hash::Hash;
use std::io::{BufWriter,Write,BufReader,BufRead};
use regex::Regex;
use std::collections::{HashMap,HashSet};
use crate::predictor::Predictor;
use crate::process_3d::standardize;
use crate::{geometry::{Point3D, Vector3D}, mmcif_process, pdbdata::{PDBAsym, PDBAtom,PDBComp, PDBEntry}, process_3d, structural_alignment::align};
use super::structural_alignment;
use super::matrix_process;
use super::geometry;
use chrono::{Local};
use std::f64::consts::PI;
use rand::prelude::*;
use super::misc_util::*;
use super::decision_tree;
use std::sync::Mutex;


const NUM_AA_INDEX:usize = 20;
const NUM_SEP:usize = 3;
const NUM_PSEUDOCB:usize = 3;
const NUM_REGION:usize = NUM_SEP*3*2;
const NUM_DIFF_PARTNER:usize = (NUM_PSEUDOCB+3)*NUM_AA_INDEX;

lazy_static! {
    static ref RESIDUES_DEFAULT:Vec<(i8,String)> =  vec![
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
        ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].into_iter().enumerate().map(|m|(m.0 as i8,m.1.to_owned())).collect();
    static ref AA_3_TO_1:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
    static ref AA_1_TO_3:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
    //RESIDUES_DEFAULT の順でインデックスが入る
    static ref AA_3_TO_INDEX:Mutex<HashMap<String,i8>> = Mutex::new(HashMap::new());
    static ref AA_1_TO_INDEX:Mutex<HashMap<String,i8>> = Mutex::new(HashMap::new());
}

pub fn prepare_static(){
    if AA_3_TO_1.lock().unwrap().len() != 0{
        return;
    }
    if AA_1_TO_3.lock().unwrap().len() != 0{
        return;
    }
    AA_3_TO_1.lock().unwrap().insert("ALA".to_string(),"A".to_string());
    AA_3_TO_1.lock().unwrap().insert("ARG".to_string(),"R".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASN".to_string(),"N".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASP".to_string(),"D".to_string());
    AA_3_TO_1.lock().unwrap().insert("CYS".to_string(),"C".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLN".to_string(),"Q".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLU".to_string(),"E".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLY".to_string(),"G".to_string());
    
    AA_3_TO_1.lock().unwrap().insert("HIS".to_string(),"H".to_string());
    //AA_3_TO_1.lock().unwrap().insert("HSD".to_string(),"H".to_string());
    //AA_3_TO_1.lock().unwrap().insert("HSC".to_string(),"H".to_string());

    AA_3_TO_1.lock().unwrap().insert("ILE".to_string(),"I".to_string());
    AA_3_TO_1.lock().unwrap().insert("LEU".to_string(),"L".to_string());
    AA_3_TO_1.lock().unwrap().insert("LYS".to_string(),"K".to_string());
    AA_3_TO_1.lock().unwrap().insert("MET".to_string(),"M".to_string());
    AA_3_TO_1.lock().unwrap().insert("PHE".to_string(),"F".to_string());
    AA_3_TO_1.lock().unwrap().insert("PRO".to_string(),"P".to_string());
    AA_3_TO_1.lock().unwrap().insert("SER".to_string(),"S".to_string());
    AA_3_TO_1.lock().unwrap().insert("THR".to_string(),"T".to_string());
    AA_3_TO_1.lock().unwrap().insert("TRP".to_string(),"W".to_string());
    AA_3_TO_1.lock().unwrap().insert("TYR".to_string(),"Y".to_string());
    AA_3_TO_1.lock().unwrap().insert("VAL".to_string(),"V".to_string());
    AA_3_TO_1.lock().unwrap().insert("UNK".to_string(),"X".to_string());
    for a in AA_3_TO_1.lock().unwrap().iter(){
        if a.0 == "HSD" || a.0 == "HSC"{
            continue;
        }
        AA_1_TO_3.lock().unwrap().insert(a.1.to_string(),a.0.to_string());
    }
    for a in RESIDUES_DEFAULT.iter(){
        AA_1_TO_INDEX.lock().unwrap().insert(AA_3_TO_1.lock().unwrap().get(&a.1).unwrap().clone(),a.0);
        AA_3_TO_INDEX.lock().unwrap().insert(a.1.clone(),a.0);
    }
}

const to_radian:f64 = PI/180.0;
#[derive(Clone,Copy)]
pub struct PseudoAtom{
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub radius:f64
}

impl Vector3D for PseudoAtom{

    fn get_xyz(&self)-> (f64,f64,f64){
        return (self.get_x(),self.get_y(),self.get_z());
    }

    fn set_vector(&mut self,v:&dyn Vector3D){
        self.set_x(v.get_x());
        self.set_y(v.get_y());
        self.set_z(v.get_z());
    }

    fn set_x(&mut self,xx:f64){
        self.x = xx;
    }
    fn set_y(&mut self,yy:f64){
        self.y = yy;
    }
    fn set_z(&mut self,zz:f64){
        self.z = zz;
    }
    fn get_x(&self)->f64{
        return self.x;
    }
    fn get_y(&self)->f64{
        return self.y;
    }
    fn get_z(&self)->f64{
        return self.z;
    }
    
    fn set_xyz(&mut self,xx:f64,yy:f64,zz:f64){
        self.set_x(xx);
        self.set_y(yy);
        self.set_z(zz);
    }
    fn distance(&self,point:&dyn Vector3D)->f64{
        let dd = (self.get_x() - point.get_x()).powf(2.0)+
        (self.get_y() - point.get_y()).powf(2.0)+
        (self.get_z() - point.get_z()).powf(2.0);
        if dd == 0.0{
            return 0.0;
        }
        return dd.sqrt();
    }
}
impl PseudoAtom{
    pub fn new(xx:f64,yy:f64,zz:f64)->PseudoAtom{
        return PseudoAtom{
            x:xx,y:yy,z:zz,radius:1.0
        };
    }
    pub fn get_radius(&self)->f64{
        return self.radius;
    }
    pub fn set_radius(&mut self,r:f64){
        self.radius = r;
    }
}


pub struct DecisionTreePCBModel<'a>{
    pub trees:Vec<decision_tree::SimpleDecisionTree<'a>>,
    pub cylinder_length:f64,
    pub cylinder_num_sep:usize,
    pub atom_radius:f64,
    pub cylinder_radius:f64
}
//ここから
//Load 関数とか作る
//最新版に変える

impl<'a> DecisionTreePCBModel<'a>{
    pub fn new()->DecisionTreePCBModel<'a>{
        return DecisionTreePCBModel{
            trees:vec![],
            cylinder_length:10.0,
            cylinder_num_sep:3,
            atom_radius:8.0,
            cylinder_radius:8.0
        };
    }
    pub fn load(filepath:&str)->DecisionTreePCBModel<'a>{
        let mut ret = DecisionTreePCBModel::new();
        prepare_static();
        
        let file = File::open(filepath).unwrap();
        let reader = BufReader::new(file);
        
        let kv =  Regex::new(r"^([^=]+)=(.+)$").unwrap();
        let lines = reader.lines();
        let mut treeflag = false;
        let mut buff:Vec<String> = vec![];
        for line_ in lines.into_iter(){
            let line_ = line_.unwrap();
            let line =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            
            if line.len() == 0{
                continue;
            }

            if start_with(&line, ">>"){
                treeflag = true;
            }
            if start_with(&line, "<<"){
                ret.trees.push(decision_tree::SimpleDecisionTree::from_string(buff));
                buff = vec![];
                continue;
            }
            if treeflag{
                if line.len() > 0{
                    buff.push(line);
                }
                continue;
            }

            let mut ptt:Vec<String> = line.split_ascii_whitespace().into_iter().map(|m|m.to_owned()).collect();
            if &ptt[0] == "option"{
                ptt.remove(0);
                for pp in ptt.iter(){
                    if let Some(x) = kv.captures(pp){
                        let k:String = x.get(1).unwrap().as_str().to_string();
                        let v:String = x.get(2).unwrap().as_str().to_string();
                        if k == "cylinder_length"{
                            ret.cylinder_length = v.parse::<f64>().unwrap();
                        }else if k == "cylinder_radius"{
                            ret.cylinder_radius = v.parse::<f64>().unwrap();
                        }else if k == "atom_radius"{
                            ret.atom_radius = v.parse::<f64>().unwrap();
                        }else if k == "cylinder_num_sep"{
                            ret.cylinder_num_sep = v.parse::<usize>().unwrap();
                        }else{
                            if pp.len() > 0{
                                panic!("can not parse {}.",pp);
                            }
                        }
                    }
                }
                continue;
            }
        }
        if buff.len() > 0{
            println!("{:?} \nwas not parsed.",buff);
        }
        
        return ret;
    }
    pub fn calc_score(&self,entt:&PDBEntry)->f64{
        //ToDo マップの作製は一回だけにする
        let varmap:HashMap<String,usize> = self.trees[0].var_names.as_ref().unwrap().iter().enumerate().map(|m|(m.1.clone(),m.0)).collect();

        let asyms:Vec<&PDBAsym> = entt.get_all_asyms();
        let mut cylinders:Vec<SideChainCylinder> = vec![];
        let mut score:f64 = 0.0;
        for (_jj,aa) in asyms.iter().enumerate(){
            //aa.remove_alt(None);
            let cnum = aa.num_comps();
            if cnum == 0{
                continue;
            }
            let st = 0;
            let en = cnum-1;
            for cc in st..=en{    
                let compp = aa.get_comp_at(cc);
                if AA_3_TO_INDEX.lock().unwrap().contains_key(compp.get_name()){
                    cylinders.push(SideChainCylinder::create_cylinder(compp, self.cylinder_length, self.cylinder_num_sep, self.atom_radius, self.cylinder_radius));
                }else{
                    if compp.get_name() != "HOH"{
                        eprintln!("{} is not found in dict!",compp.get_name());
                    }
                }
            } 
        }
        let rxcode =  Regex::new(r"^([^_]+)_([0-9]+)$").unwrap();
        for (cii,cc) in cylinders.iter().enumerate(){
            let mut space_count:HashMap<String,HashMap<String,usize>> = HashMap::new();
            let rcode = *AA_3_TO_INDEX.lock().unwrap().get(&cc.comp_name).unwrap() as usize;
            for (cjj,ccc) in cylinders.iter().enumerate(){
                if cii != cjj{
                    space_count = cc.count_atoms_str(ccc,Some(space_count));
                }
            }
            let mut inputdat:Vec<f64> = vec![0.0;varmap.len()];
            for pp in space_count.iter(){
                if let Some(x) = rxcode.captures(pp.0){
                    let _rcode = *AA_3_TO_INDEX.lock().unwrap().get(x.get(1).unwrap().as_str()).unwrap();
                    let pcode = x.get(2).unwrap().as_str().parse::<usize>().unwrap();
                    for qq in pp.1.iter(){
                        let zcode = pcode.to_string()+"_"+qq.0;
                        inputdat[*varmap.get(&zcode).unwrap_or_else(||{panic!("{} was not found in key list.",zcode)})] = *qq.1 as f64;
                    }
                }                    
            }
            for tt in self.trees.iter(){
                score += tt.predict_class(&inputdat)[rcode];
            }
        }
        return score;
    }
}

pub struct VerySimplePCBModel{
    pub exist:Vec<Vec<HashMap<String,f64>>>,//ある場所にある PSEUDOCB が有るときに、その COMP である確率
    pub non_exist:Vec<Vec<HashMap<String,f64>>>,//ある場所にある PSEUDOCB が無いときに、その COMP である確率
    pub r_exist:Vec<Vec<HashMap<String,f64>>>,//その COMP であるときにある場所にある PSEUDOCB が有る確率
    pub r_non_exist:Vec<Vec<HashMap<String,f64>>>,//その COMP であるときにある場所にある PSEUDOCB が無い確率
    pub cylinder_length:f64,
    pub cylinder_num_sep:usize,
    pub atom_radius:f64,
    pub cylinder_radius:f64
}
impl VerySimplePCBModel{
    pub fn new()->VerySimplePCBModel{
        let mut exx:Vec<Vec<HashMap<String,f64>>> = vec![];
        let mut nexx:Vec<Vec<HashMap<String,f64>>> = vec![];
        let mut rexx:Vec<Vec<HashMap<String,f64>>> = vec![];
        let mut rnexx:Vec<Vec<HashMap<String,f64>>> = vec![];
        for _ in 0..NUM_AA_INDEX{
            let mut vvhmm:Vec<HashMap<String,f64>> = vec![];
            for _ in 0..NUM_REGION{
                let mut hmm:HashMap<String,f64> = HashMap::new();
                for aa in RESIDUES_DEFAULT.iter(){
                    for pp in 0..(NUM_PSEUDOCB+3){
                        hmm.insert(aa.1.clone()+"_"+pp.to_string().as_str(),0.0);
                    }
                }
                vvhmm.push(hmm.clone());
            }
            exx.push(vvhmm.clone());
            nexx.push(vvhmm.clone());
            rexx.push(vvhmm.clone());
            rnexx.push(vvhmm);
        }
        return VerySimplePCBModel{
            exist:exx,
            non_exist:nexx,
            r_exist:rexx,
            r_non_exist:rnexx,
            cylinder_length:10.0,
            cylinder_num_sep:3,
            atom_radius:8.0,
            cylinder_radius:8.0
        };
    }
    pub fn calc_score(&self,entt:&PDBEntry)->(f64,f64){
        let asyms:Vec<&PDBAsym> = entt.get_all_asyms();
        let mut cylinders:Vec<SideChainCylinder> = vec![];
        let mut score:f64 = 0.0;
        let mut score2:f64 = 0.0;
        for (_jj,aa) in asyms.iter().enumerate(){
            //aa.remove_alt(None);
            let cnum = aa.num_comps();
            if cnum == 0{
                continue;
            }
            let st = 0;
            let en = cnum-1;
            for cc in st..=en{    
                let compp = aa.get_comp_at(cc);
                if AA_3_TO_INDEX.lock().unwrap().contains_key(compp.get_name()){
                    cylinders.push(SideChainCylinder::create_cylinder(compp, self.cylinder_length, self.cylinder_num_sep, self.atom_radius, self.cylinder_radius));
                }else{
                    if compp.get_name() != "HOH"{
                        eprintln!("{} is not found in dict!",compp.get_name());
                    }
                }
            } 
        }
        let rxcode =  Regex::new(r"^([^_]+)_([0-9]+)$").unwrap();
        for (cii,cc) in cylinders.iter().enumerate(){
            let mut space_count:HashMap<String,HashMap<String,usize>> = HashMap::new();
            let mut counted:Vec<HashSet<String>> = vec![HashSet::new();NUM_REGION];
            let rcode = *AA_3_TO_INDEX.lock().unwrap().get(&cc.comp_name).unwrap() as usize;
            for (cjj,ccc) in cylinders.iter().enumerate(){
                if cii != cjj{
                    space_count = cc.count_atoms_str(ccc,Some(space_count));
                }
            }
            for pp in space_count.iter(){
                if let Some(x) = rxcode.captures(pp.0){
                    let _rcode = *AA_3_TO_INDEX.lock().unwrap().get(x.get(1).unwrap().as_str()).unwrap();
                    let pcode = x.get(2).unwrap().as_str().parse::<usize>().unwrap();
                    for qq in pp.1.iter(){
                        if *qq.1 > 0{
                            counted[pcode].insert(qq.0.clone());
                        }
                    }
                }                    
            }
            for pp in 0..self.exist[rcode].len(){
                for qq in self.exist[rcode][pp].iter(){
                    if counted[pp].contains(qq.0){
                        score += qq.1;
                    }
                }
            }
            for pp in 0..self.non_exist[rcode].len(){
                for qq in self.non_exist[rcode][pp].iter(){
                    if !counted[pp].contains(qq.0){
                        score += qq.1;
                    }
                }
            }
            for pp in 0..self.r_exist[rcode].len(){
                for qq in self.r_exist[rcode][pp].iter(){
                    if counted[pp].contains(qq.0){
                        score2 += qq.1;
                    }
                }
            }
            for pp in 0..self.r_non_exist[rcode].len(){
                for qq in self.r_non_exist[rcode][pp].iter(){
                    if !counted[pp].contains(qq.0){
                        score2 += qq.1;
                    }
                }
            }
        }
        return (score,score2);
    }
    pub fn load(filepath:&str)->VerySimplePCBModel{
        prepare_static();

        let mut ret :VerySimplePCBModel = VerySimplePCBModel::new();
        
        let file = File::open(filepath).unwrap();
        let reader = BufReader::new(file);
        
        let mut exx:Vec<Vec<HashMap<String,f64>>> = vec![vec![HashMap::new();NUM_REGION];NUM_AA_INDEX];
        let mut nexx:Vec<Vec<HashMap<String,f64>>> = vec![vec![HashMap::new();NUM_REGION];NUM_AA_INDEX];

        let rxcode =  Regex::new(r"^([^_]+)_([0-9]+)$").unwrap();
        let kv =  Regex::new(r"^([^=]+)=(.+)$").unwrap();
        let mut count_compound:HashMap<usize,f64> = HashMap::new();//region でなく comound 数
        for line_ in reader.lines() {
            let line_ = line_.unwrap();
            let line =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            if line.len() == 0{
                continue;
            }
            let mut ptt:Vec<String> = line.split_ascii_whitespace().into_iter().map(|m|m.to_owned()).collect();
            let mut selfcode:usize = 0;//円柱のインデクス
            let mut selfaacode:usize = 0;//アミノ酸インデクス
            if &ptt[0] == "option"{
                ptt.remove(0);
                for pp in ptt.iter(){
                    if let Some(x) = kv.captures(pp){
                        let k:String = x.get(1).unwrap().as_str().to_string();
                        let v:String = x.get(2).unwrap().as_str().to_string();
                        if k == "cylinder_length"{
                            ret.cylinder_length = v.parse::<f64>().unwrap();
                        }else if k == "cylinder_radius"{
                            ret.cylinder_radius = v.parse::<f64>().unwrap();
                        }else if k == "atom_radius"{
                            ret.atom_radius = v.parse::<f64>().unwrap();
                        }else if k == "cylinder_num_sep"{
                            ret.cylinder_num_sep = v.parse::<usize>().unwrap();
                        }else{
                            if pp.len() > 0{
                                panic!("can not parse {}.",pp);
                            }
                        }
                    }
                }
                continue;
            }
            for pp in ptt.iter().enumerate(){
                if pp.0 == 0{
                    if let Some(x) =  rxcode.captures(pp.1){
                        let name = x.get(1).unwrap().as_str();
                        let code = x.get(2).unwrap().as_str().parse::<usize>().unwrap();
                        selfcode = code;
                        selfaacode = *AA_3_TO_INDEX.lock().unwrap().get(name).unwrap_or_else(||panic!("{} is not a valid code.",name)) as usize;
                    }else{
                        if pp.1.len() > 0{
                            panic!("??? {} ",pp.1);
                        }
                    }
                }else{
                    if let Some(x) = kv.captures(pp.1){
                        let k = x.get(1).unwrap().as_str().to_string();
                        let v = x.get(2).unwrap().as_str().parse::<f64>().unwrap();
                        if &k == "all"{
                            if count_compound.contains_key(&selfaacode){
                                assert_eq!(*count_compound.get(&selfaacode).unwrap(),v);
                            }else{
                                count_compound.insert(selfaacode,v);
                            }
                        }else{
                            exx[selfaacode][selfcode].insert(k,v);
                        }
                    }else{
                        if pp.1.len() > 0{
                            panic!("{} does not have key value pair.",pp.1);
                        }
                    }
                }
            }
        }
        for aa in 0..NUM_AA_INDEX{
            for bb in 0..NUM_REGION{
                for pp in RESIDUES_DEFAULT.iter(){
                    for cc in 0..(NUM_PSEUDOCB+3){
                        let ccode = pp.1.clone()+"_"+cc.to_string().as_str();
                        if !exx[aa][bb].contains_key(&ccode){
                            //下の nexx で使用するので作成する
                            exx[aa][bb].insert(ccode.clone(),0.0);
                        }
                        nexx[aa][bb].insert(ccode.clone(),*count_compound.get(&aa).unwrap_or(&0.0)-*exx[aa][bb].get(&ccode).unwrap());
                    }
                }
            }
        }
        
        for aa in 0..NUM_AA_INDEX{
            for bb in 0..NUM_REGION{
                for pp in RESIDUES_DEFAULT.iter(){
                    for cc in 0..(NUM_PSEUDOCB+3){
                        let ccode = pp.1.clone()+"_"+cc.to_string().as_str();
                        let cnumall = *count_compound.get(&aa).unwrap();
                        if  cnumall > 0.0{
                            ret.r_exist[aa][bb].insert(ccode.clone(), exx[aa][bb].get(&ccode).unwrap()/cnumall);
                            ret.r_non_exist[aa][bb].insert(ccode.clone(), (cnumall -exx[aa][bb].get(&ccode).unwrap())/cnumall);
                        }else{
                            ret.r_exist[aa][bb].insert(ccode.clone(),0.0);
                            ret.r_non_exist[aa][bb].insert(ccode.clone(),0.0);
                        }
                    }
                }
            }
        }
        //全 compound 数
        let asum = count_compound.iter().fold(0.0,|s,m|s+*m.1) as f64;
        for bb in 0..NUM_REGION{
            for pp in RESIDUES_DEFAULT.iter(){
                for cc in 0..(NUM_PSEUDOCB+3){
                    let ccode = pp.1.clone()+"_"+cc.to_string().as_str();
                    //その pseudo atom が存在する compound 数
                    let bsum = exx.iter().fold(0.0,|s,m|s+*m[bb].get(&ccode).unwrap()) as f64;
                    for aa in 0..NUM_AA_INDEX{
                        if bsum > 0.0{
                            ret.exist[aa][bb].insert(ccode.clone(), *exx[aa][bb].get(&ccode).unwrap() as f64 /bsum);
                        }else{
                            ret.exist[aa][bb].insert(ccode.clone(), 0.0);
                        }
                        if asum -bsum > 0.0{
                            ret.non_exist[aa][bb].insert(ccode.clone(), *nexx[aa][bb].get(&ccode).unwrap() as f64 /(asum-bsum));
                        }else{
                            ret.non_exist[aa][bb].insert(ccode.clone(), 0.0);
                        }
                    }
                }
            }
        }
        
        return ret;
    }
}

pub struct SideChainCylinder<'a>{
    pub n:&'a dyn Vector3D,
    pub ca:&'a dyn Vector3D,
    pub c:&'a dyn Vector3D,
    pub atom_radius:f64,
    pub cylinder_radius:f64,
    pub rx:Point3D,//ry と rz の法線であり、nc center->c のベクトル
    pub ry:Point3D,//n ca c の norm を CB 分回転したベクトル 
    pub rz:Point3D,//CA →見積もった CB へのベクトル
    pub length:f64,//CA から最終点の距離
    pub rz_length:Point3D,
    pub num_sep:usize,//length を何個に分割するか//120 度ごと 3*2 つに分割するので、領域はこれ x6
    pub pseudo_atoms:Vec<PseudoAtom>,// CB 方向 cylinder 内に配置される PseudoAtom 群
    pub comp_name:String,
}
impl<'a> SideChainCylinder<'a>{
    pub fn new(n:&'a dyn Vector3D,ca:&'a dyn Vector3D,c:&'a dyn Vector3D,cname:&str,llen:f64,ssep:usize,aradi:f64,cradi:f64)->SideChainCylinder<'a>{
        prepare_static();
        let mut ret = SideChainCylinder{
        n:n,
        ca:ca,
        c:c,
        atom_radius:aradi,
        cylinder_radius:cradi,
        rx:Point3D::new(0.0,0.0,0.0),//n ca c の norm と rz の norm
        ry:Point3D::new(0.0,0.0,0.0),//n ca c の norm
        rz:Point3D::new(0.0,0.0,0.0),//ca->cb の単位ベクトル
        rz_length:Point3D::new(0.0,0.0,0.0),//ca->cb を length 分伸ばした線分。他原子の評価の際に毎回計算する必要があるのであらかじめ計算しておく
        length:llen,//CA から最終点の距離
        num_sep:ssep,
        pseudo_atoms:vec![PseudoAtom::new(0.0,0.0,0.0);ssep],
        comp_name:cname.to_string()
        };
        ret.update();
        return ret;
    }

    pub fn generate_obj(&self)->Vec<(Vec<Point3D>,Vec<geometry::Face>)>{
        let st:(f64,f64,f64) = self.ca.get_xyz();
        let en:(f64,f64,f64) = (
            st.0+self.rz.get_x()*self.length
            ,st.1+self.rz.get_y()*self.length
            ,st.2+self.rz.get_z()*self.length
        );
        let mut ret:Vec<(Vec<Point3D>,Vec<geometry::Face>)> = vec![];
        for pp in self.pseudo_atoms.iter(){
            let mut rres = geometry::Geometry::generate_sphere(&pp.get_xyz(),0.5,8,8);
            geometry::Face::color_faces(&mut rres.1,&vec![0,0,255]);
            ret.push(rres);
        }
        ret.push(geometry::Geometry::generate_cylinder(&st,&en,self.cylinder_radius,8,false));
        ret.push(geometry::Geometry::generate_sphere(&self.ca.get_xyz(),0.5,8,8));
        ret.push(geometry::Geometry::generate_sphere(&self.n.get_xyz(),0.5,8,8));
        ret.push(geometry::Geometry::generate_sphere(&self.c.get_xyz(),0.5,8,8));
        return ret;        
    }

    pub fn update(&mut self){
        let mut nn = self.n.get_xyz();
        let caa = self.ca.get_xyz();
        let mut cc = self.c.get_xyz();
        nn = process_3d::standardize(nn.0-caa.0,nn.1-caa.1,nn.2-caa.2);
        cc = process_3d::standardize(cc.0-caa.0,cc.1-caa.1,cc.2-caa.2);


        
        //CHARMM の角度と二面角と結合長から計算したら微妙に違ったので
        //もう適当でいいやという気分になった
        //C-CA-N から CB 位置を計算する
        //法線方向の距離
        let pyy = (18.5*to_radian).cos()*((60.0*to_radian).sin());
        //standardize された nc 中心点逆方向ベクトルの距離
        let pzz =  (1.0-pyy*pyy).sqrt();


        let yy:(f64,f64,f64) = process_3d::calc_norm_t(
            &nn,
            &(0.0,0.0,0.0),
            &cc
        );
        let nccenter = (nn.0/2.0+cc.0/2.0,nn.1/2.0+cc.1/2.0,nn.2/2.0+cc.2/2.0);
        self.rx.set_xyz(cc.0-nccenter.0,cc.1-nccenter.1,cc.2-nccenter.2);
        self.rx.standardize();
        

        //n-c の中点
        let mut pcx:(f64,f64,f64) = standardize(nn.0/2.0+cc.0/2.0,nn.1/2.0+cc.1/2.0,nn.2/2.0+cc.2/2.0);
        pcx.0 *= -1.0;
        pcx.1 *= -1.0;
        pcx.2 *= -1.0;

        self.rz.set_xyz(
            yy.0*-pyy+pcx.0*pzz
            ,yy.1*-pyy+pcx.1*pzz
            ,yy.2*-pyy+pcx.2*pzz
        );
        self.rz.standardize();

        let ynorm = process_3d::calc_norm(
            self.rz.get_x(),
            self.rz.get_y(),
            self.rz.get_z(),
            self.rx.get_x(),
            self.rx.get_y(),
            self.rx.get_z()
        );
        self.ry.set_xyz(ynorm.0,ynorm.1,ynorm.2);
        
        let stepp:f64 = self.length/(self.num_sep as f64);
        for (pii,pp) in self.pseudo_atoms.iter_mut().enumerate(){
            let rr = ((pii+1) as f64)*stepp;
            pp.set_xyz(
                self.ca.get_x()+self.rz.get_x()*rr
                ,self.ca.get_y()+self.rz.get_y()*rr
                ,self.ca.get_z()+self.rz.get_z()*rr
            );
        }
        self.rz_length.set_xyz(
            self.rz.get_x()*self.length,
            self.rz.get_y()*self.length,
            self.rz.get_z()*self.length
        );
        if false{//デバッグ用コード
            
            
            let candidate1:(f64,f64,f64)
            = (
                (yy.0*-pyy+pcx.0*pzz)*1.52
                ,
                (yy.1*-pyy+pcx.1*pzz)*1.52
                ,
                (yy.2*-pyy+pcx.2*pzz)*1.52
            );
        
            //CHARMM の Force field param から CB 位置を計算する
            let zz = process_3d::standardize(
                -1.0*(nn.0/2.0+cc.0/2.0)
                ,-1.0*(nn.1/2.0+cc.1/2.0)
                ,-1.0*(nn.2/2.0+cc.2/2.0)
                );
            self.rz.set_xyz(zz.0*self.length,zz.1*self.length,zz.2*self.length);
            //C-CA-CB の ANGLE は charmm19  (param19)で 
            //C    CH1E CH3E    70.0     106.5
            //N-CA-CB の ANGLE は 
            //CH3E CH1E NH1     65.0     108.5
            //bond は
            //CH1E CH3E   225.0       1.52
            //N-CA-C は
            //C    CH1E NH1     45.0     111.6
            //ALA の IMPROPER (toph19)は 
            //IC   N    C    *CA  CB     0.0000    0.00  120.00    0.00   0.0000
            //値を入れて計算してみると少しずれる。。。
            let xx:(f64,f64,f64) = process_3d::calc_norm_t(
                &zz,
                &(0.0,0.0,0.0),
                &yy
            );

            let mut cbb = Point3D::from_tuple(&cc);
            cbb.standardize();
            cbb.set_x(cbb.get_x()* 1.52);
            cbb.set_y(cbb.get_y()* 1.52);
            cbb.set_z(cbb.get_z()* 1.52);

            process_3d::rotate_3d(&mut vec![&mut cbb],&Point3D::from_tuple(&yy),-106.5*to_radian);
            process_3d::rotate_3d(&mut vec![&mut cbb],&Point3D::from_tuple(&cc),60.0*to_radian);
            if true{
                let mut ggeo = geometry::Geometry::new();
                for a in vec![&nn,&caa,&cc]{
                    let mut mm = geometry::Geometry::generate_sphere(a,0.25,8,8);
                    for f in mm.1.iter_mut(){
                        f.set_color(&vec![255,0,255]);
                    }
                    ggeo.add_objects(mm);
                }

                let mut mm = geometry::Geometry::generate_sphere(&cbb.get_xyz(),0.25,8,8);
                println!("{:?}",&cbb);
                for f in mm.1.iter_mut(){
                    f.set_color(&vec![255,0,0]);
                }
                ggeo.add_objects(mm);
                println!("cand1: {} {:?}",process_3d::distance(
                    &(-0.949,-0.005,1.177)
                    ,&candidate1
                ),candidate1);
                println!("cand2: {} {:?}",process_3d::distance(
                    &(-0.949,-0.005,1.177)
                    ,&cbb.get_xyz()
                ),&cbb.get_xyz());

                for a in vec![&(0.822,-1.200,-0.000),&(0.000,0.000,0.000),&(-0.949,-0.005,1.177),&(0.860,1.255,0.000)]{
                    let mut mm = geometry::Geometry::generate_sphere(a,0.25,8,8);
                    for f in mm.1.iter_mut(){
                        f.set_color(&vec![255,255,255]);
                    }
                    ggeo.add_objects(mm);
                }


                ggeo.calc_all_norms();
                ggeo.add_colortile_material();
                ggeo.calc_all_norms();
                let mut gv:Vec<geometry::Geometry> = vec![ggeo];
            
                geometry::Geometry::save("test/anglecheck.obj",&mut gv);
            }
        }
    }

    //
    //方向を 3 分割して 半径/2 半径*1 の 6 個ずつの分割 +
    //分割中の近い位置を 0 として何番目の分割内にあるか *6
    //で返す
    pub fn get_position_of(&self,atom:&dyn Vector3D)->Option<u8>{
        if atom.distance(self.ca).powf(2.0) > self.length*self.length + self.cylinder_radius*self.cylinder_radius{
            return None;
        }
        let mut ppos = atom.get_xyz();
        ppos.0 -= self.ca.get_x();
        ppos.1 -= self.ca.get_y();
        ppos.2 -= self.ca.get_z();
        
        let a:f64 = process_3d::distance(&ppos,&self.rz_length.get_xyz());
        let b:f64 = process_3d::distance(&ppos,&(0.0,0.0,0.0));
        if b*b >= a*a+self.length*self.length || a*a >= b*b+self.length*self.length{            
            return None;
        }

        let d = (b*b-a*a+self.length*self.length)/(2.0*self.length);
        assert!(d > 0.0);
        let ratio = d/self.length;
        let rpos:(f64,f64,f64) = (
            self.rz_length.get_x()*ratio,
            self.rz_length.get_y()*ratio,
            self.rz_length.get_z()*ratio
        );

        //ppos を原点中心の相対的な位置に移動
        ppos.0 -= rpos.0;
        ppos.1 -= rpos.1;
        ppos.2 -= rpos.2;

        let mut ddist = ppos.0*ppos.0+ppos.1*ppos.1+ppos.2*ppos.2;//中心点からの距離
        
        if ddist > 0.0{
            ddist = ddist.sqrt();
            ppos.0 /= ddist;
            ppos.1 /= ddist;
            ppos.2 /= ddist;
        }//ppos は単位ベクトルに直した
        if ddist > self.cylinder_radius {
            return None;
        }
        let bdis = process_3d::distance(
            &ppos
            ,&(self.rx.get_xyz())
        );
        let mut rad_1 = (1.0-bdis*bdis/2.0).acos();
        
        let cdis = process_3d::distance(&ppos,&self.ry.get_xyz());

        if cdis > (2.0_f64).sqrt(){
            rad_1 *= -1.0;
        }
        let direc_code:u8 = 
        if rad_1 < PI*2.0/3.0 && rad_1 >= 0.0{
            0
        }else if rad_1 >= -1.0*PI*2.0/3.0 && rad_1 < 0.0{
            1
        }else{
            2
        };
        let direc_code:u8 = if ddist < self.cylinder_radius/2.0{
            direc_code
        }else{
            direc_code + 3
        };
        let rstep = 1.0/(self.num_sep as f64);
        let mut pcode:u8 = 0;
        for ii in 1..=self.num_sep{
            if rstep*(ii as f64) > ratio{
                pcode = (ii -1) as u8 ;
                break;
            }
        }
        
        //自分側は 近い方から 012345 6789,10,11 12,13,14,15,16,17
        return Some(pcode*6+direc_code);
    }


    //count_atoms_str の結果である HashMap<String,HashMap<String,usize>> を受け取って
    //return される vec の var_mapping から get した index に atomstate の最終的な usize の値を入れて返す
    // var_mapping の key は、selfregionid_selfresiduename_selfatomid であることを想定している
    pub fn create_sample_input(atomstate:HashMap<String,HashMap<String,usize>>,var_mapping:&HashMap<String,usize>)->Vec<f64>{
        let mut ret:Vec<f64> = vec![0.0;var_mapping.len()];
        let exx =  Regex::new(r"_(.+)$").unwrap();
        for aa in atomstate.into_iter(){
            let k = if let Some(x) = exx.captures(&aa.0){
                x.get(1).unwrap().as_str().to_string()
            }else{
                panic!("Can not find region index from {}.",aa.0);
            };
            for pp in aa.1.iter(){
                let code = k.clone()+"_"+pp.0;
                ret[*var_mapping.get(&code).unwrap_or_else(||panic!("{} was not found in mapper.",code))] = (*pp.1) as f64;
            }
        }
        return ret;
    }



    pub fn create_cylinder(comp:&PDBComp,llen:f64,ssep:usize,aradi:f64,cradi:f64)->SideChainCylinder{
        let n = comp.get_N();
        let ca = comp.get_CA();
        let c = comp.get_C();
        if let (Some(_),Some(_),Some(_)) = (n,ca,c){
        }else{
            panic!("This compound lacks backbone atom! {:?}",comp);
        }
        let ret = SideChainCylinder::new(n.unwrap(),ca.unwrap(),c.unwrap(),comp.get_name(),llen,ssep,aradi,cradi);
        return ret;
    }
    pub fn count_atoms_str(&self,target:&SideChainCylinder,ret_:Option<HashMap<String,HashMap<String,usize>>>) -> HashMap<String,HashMap<String,usize>>{
        let mut ret:HashMap<String,HashMap<String,usize>> = match ret_{
            Some(x) =>{x},
            _=>{
                HashMap::new()
            }
        };
        let mut zat:Vec<&dyn Vector3D> = target.pseudo_atoms.iter().map(|m| m as &dyn Vector3D).collect();
        zat.push(target.n);
        zat.push(target.ca);
        zat.push(target.c);
        
        for a in zat.into_iter().enumerate(){
            if let Some(x) = self.get_position_of(a.1){
                let zcode = self.comp_name.clone()+"_"+x.to_string().as_str();
                if !ret.contains_key(&zcode){
                    ret.insert(zcode.clone(),HashMap::new());
                }
                let pcode = target.comp_name.clone()+"_"+a.0.to_string().as_str();
                let co = ret.get(&zcode).unwrap().get(&pcode).unwrap_or(&0)+1;
                ret.get_mut(&zcode).unwrap().insert(pcode.clone(),co);
            }
        }
        return ret;
    }
}


pub fn asyms_to_cylinders<'a>(asyms:Vec<&'a mut PDBAsym>
    , cylinder_length:f64
    , cylinder_num_sep:usize
    , atom_radius:f64
    , cylinder_radius:f64
)->Vec<SideChainCylinder<'a>>{
    let mut cylinders:Vec<SideChainCylinder<'a>> = vec![];
    for (jj,aa) in asyms.into_iter().enumerate(){
        aa.remove_alt(None);
        let cnum = aa.num_comps();
        if cnum == 0{
            continue;
        }
        let mut st = 0;
        let mut en = cnum-1;
        for cc in 0..cnum{
            let comp = aa.get_comp_at(cc);
            let n = comp.get_N();
            let ca = comp.get_CA();
            let c = comp.get_C();
            if let (Some(_),Some(_),Some(_)) = (n,ca,c){
                st = cc;
                break;
            }
        }
        for cc_ in 0..cnum{
            let cc = cnum - cc_-1;
            let comp = aa.get_comp_at(cc);
            let n = comp.get_N();
            let ca = comp.get_CA();
            let c = comp.get_C();
            if let (Some(_),Some(_),Some(_)) = (n,ca,c){
                en = cc;
                break;
            }
        }

        //for (_kk,compp) in aa.iter_comps().enumerate(){
        for cc in st..=en{    
            let compp:&'a PDBComp = aa.get_comp_at(cc);
            if AA_3_TO_INDEX.lock().unwrap().contains_key(compp.get_name()){
                cylinders.push(SideChainCylinder::create_cylinder(compp
                    , cylinder_length, cylinder_num_sep
                    , atom_radius, cylinder_radius));
            }else{
                if compp.get_name() != "HOH"{
                    eprintln!("{} is not found in dict!",compp.get_name());
                }
            }
        } 
    }
    return cylinders;
}

pub fn generate_simplemodel_files(
    structure_files:Vec<(String,Option<HashMap<String,f64>>)>//path, weight(chain->weight)
    , cylinder_length:f64
    , cylinder_num_sep:usize
    , atom_radius:f64
    , cylinder_radius:f64
    ,output_filename:&str
    ){
    //自分側は 近い方から 012345 6789,10,11 12,13,14,15,16,17
    //相手側は
    //ARNDCQEGHILKMFPSTWYVX
    //pseudoatoms[0]...+ncac の順でインデクスをつける
    prepare_static();
    let optionline:String = format!("option\tcylinder_length={}\tcylinder_num_sep={}\tatom_radius={}\tcylinder_radius={}\n"
    , cylinder_length
    , cylinder_num_sep
    , atom_radius
    , cylinder_radius);
    let entries_:Vec<String> = structure_files.iter().map(|m|m.0.clone()).collect();
    
    //COMPNAME_POSITION->COMPNAME_PSEUDOATOMCODE->count
    let mut space_count:HashMap<String,HashMap<String,usize>> = HashMap::new();
    let mut type_count:HashMap<String,usize> = HashMap::new();
    let exx =  Regex::new(r"(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    for ee in entries_.into_iter(){
        if let Some(x) = exx.captures(&ee){
            let ext1:String = x.get(1).unwrap().as_str().to_string();
            let is_gzip = 
            if let Some(_) = x.get(2){
                true
            }else{
                false
            };
            println!("Loading {}.",ee);
            let entt_ = if &ext1 == ".ent" || &ext1 == ".pdb"{
                Some(mmcif_process::load_pdb(&ee,is_gzip))
            }else if &ext1 == ".cif"{
                Some(mmcif_process::MMCIFEntry::load_mmcif(&ee,is_gzip))
            }else{
                println!("Skipped {} (unknown extension)",ee);
                None
            };
            if let Some(mut entt) = entt_{
                let cylinders:Vec<SideChainCylinder> = asyms_to_cylinders(
                    entt.get_mut_all_asyms()
                    , cylinder_length
                    , cylinder_num_sep
                    , atom_radius
                    , cylinder_radius
                );
                for (cii,cc) in cylinders.iter().enumerate(){
                    type_count.insert(cc.comp_name.clone(),type_count.get(&cc.comp_name).unwrap_or(&0)+1);
                    for (cjj,ccc) in cylinders.iter().enumerate(){
                        if cii != cjj{
                            space_count = cc.count_atoms_str(ccc,Some(space_count));
                        }
                    }
                }
            }
        }
    }
    let mut skeys:Vec<&String> = space_count.keys().into_iter().collect();
    skeys.sort();
    let mut f = BufWriter::new(fs::File::create(output_filename).unwrap());
    f.write(optionline.as_bytes()).unwrap();
    let grx =  Regex::new(r"^([^_]+)_").unwrap();
    for ss in skeys.into_iter(){
        
        f.write_all(format!("{}",ss).as_bytes()).unwrap();
        
        if let Some(x) = grx.captures(ss){
            f.write_all(format!("\tall={}",type_count.get(x.get(1).unwrap().as_str()).unwrap()).as_bytes()).unwrap();
        }else{
            panic!("not expected! {}",ss);
        }
        let hm = space_count.get(ss.to_string().as_str()).unwrap();
        let mut keys:Vec<&String> = hm.keys().into_iter().collect();
        keys.sort();
        for kk in keys.iter(){
            f.write_all(format!("\t{}={}",kk,hm.get(kk.as_str()).unwrap()).as_bytes()).unwrap();
        }
        f.write_all("\n".as_bytes()).unwrap();
    }
}



pub fn generate_decision_tree_model_files(
    structure_files:Vec<(String,Option<HashMap<String,f64>>)>//path, weight(chain->weight)
    , cylinder_length:f64
    , cylinder_num_sep:usize
    , atom_radius:f64
    , cylinder_radius:f64
    ,output_filename:&str
    ){
    //自分側は 近い方から 012345 6789,10,11 12,13,14,15,16,17
    //相手側は
    //ARNDCQEGHILKMFPSTWYVX
    //pseudoatoms[0]...+ncac の順でインデクスをつける
    prepare_static();
    let optionline:String = format!("option\tcylinder_length={}\tcylinder_num_sep={}\tatom_radius={}\tcylinder_radius={}\n"
    , cylinder_length
    , cylinder_num_sep
    , atom_radius
    , cylinder_radius);
    let entries_:Vec<String> = structure_files.iter().map(|m|m.0.clone()).collect();
    
    let mut varmap:HashMap<String,usize> = HashMap::new();
    let mut varnames:Vec<String> = vec![];
    for ii in 0..(cylinder_num_sep*3*2){
        let idd:String = ii.to_string();
        for rr in RESIDUES_DEFAULT.iter(){
            let ridd = idd.clone()+"_"+rr.1.as_str();
            for jj in 0..(NUM_PSEUDOCB+3){
                let code = ridd.clone()+"_"+jj.to_string().as_str();
                varnames.push(code);
            }
        }
    }
    for vv in varnames.iter().enumerate(){
        varmap.insert(vv.1.clone(),vv.0);
    }
    let mut samples:Vec<decision_tree::SampleData_usize> = vec![];
    let exx =  Regex::new(r"(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    for ee in entries_.into_iter(){
        //COMPNAME_POSITION->COMPNAME_PSEUDOATOMCODE->count
        if let Some(x) = exx.captures(&ee){
            let ext1:String = x.get(1).unwrap().as_str().to_string();
            let is_gzip = 
            if let Some(_) = x.get(2){
                true
            }else{
                false
            };
            println!("Loading {}.",ee);
            let entt_ = if &ext1 == ".ent" || &ext1 == ".pdb"{
                Some(mmcif_process::load_pdb(&ee,is_gzip))
            }else if &ext1 == ".cif"{
                Some(mmcif_process::MMCIFEntry::load_mmcif(&ee,is_gzip))
            }else{
                println!("Skipped {} (unknown extension)",ee);
                None
            };
            if let Some(mut entt) = entt_{
                let cylinders:Vec<SideChainCylinder> = asyms_to_cylinders(
                    entt.get_mut_all_asyms()
                    , cylinder_length
                    , cylinder_num_sep
                    , atom_radius
                    , cylinder_radius
                );
                for (cii,cc) in cylinders.iter().enumerate(){
                    let mut space_count:HashMap<String,HashMap<String,usize>> = HashMap::new();
                    for (cjj,ccc) in cylinders.iter().enumerate(){
                        if cii != cjj{
                            space_count = cc.count_atoms_str(ccc,Some(space_count));
                        }
                    }
                    
                    let vars = SideChainCylinder::create_sample_input(space_count,&varmap);
                    let sam = decision_tree::SampleData_usize::new(
                        "".to_owned(),*AA_3_TO_INDEX.lock().unwrap().get(&cc.comp_name).unwrap() as usize
                        ,vars
                    );
                    samples.push(sam);
                }
            }
        }
    }
    let mut boxed:Vec<Box<& decision_tree::SampleData_usize>> = vec![];
    for ss in samples.iter(){
        boxed.push(Box::new(ss));
    }
    let varnames_:Vec<&str> = varnames.iter().map(|m| m.as_str()).collect();
    let tree = decision_tree::SimpleDecisionTree::build_classifier(&mut boxed
        ,decision_tree::DecisionTreeOptions{ 
	max_depth:Some(10000),
	max_leaf_nodes:Some(100000),
	min_samples_split:Some(1),
	min_samplefraction_split:None,
	feature_fraction :None,
	split_function_type:decision_tree::SplitFunctionType::Gini,
	random_seed:Some(123)
    },Some(&varnames_) );
    let sv = tree.to_string();
    
    decision_tree::write_to_file(output_filename,vec![optionline,sv]);
}



#[test]
fn pseudo_cb_test(){
    let mut n:PseudoAtom = PseudoAtom::new(0.822,-1.200,-0.000);
    let mut ca:PseudoAtom = PseudoAtom::new(0.000,0.000,0.000);
    let mut cb:PseudoAtom = PseudoAtom::new(-0.949,-0.005,1.177);
    let mut c:PseudoAtom = PseudoAtom::new(0.860,1.255,0.000);

    let mut rgen:StdRng =  SeedableRng::seed_from_u64(10);
    

    let mut v1 = Point3D::new(rgen.gen_range(-10.0..10.0),rgen.gen_range(-10.0..10.0),rgen.gen_range(-10.0..10.0));
    v1.standardize();

    process_3d::rotate_3d(&mut vec![&mut n,&mut ca,&mut cb,&mut c],&v1, rgen.gen_range(-360.0..360.0)/180.0*PI);


    let v2 = Point3D::new(rgen.gen_range(-3.0..3.0),rgen.gen_range(-3.0..3.0),rgen.gen_range(-3.0..3.0));

    for vv in vec![&mut n,&mut ca,&mut cb,&mut c].into_iter(){
        let p = vv.get_xyz();

        vv.set_xyz(p.0+v2.get_x(),p.1+v2.get_y(), p.2+v2.get_z());
    }
    let cyl = SideChainCylinder::new(&n,&ca,&c,"UNK",8.0,3,6.0,6.0);
    let obj = cyl.generate_obj();
    let mut geom:geometry::Geometry = geometry::Geometry::new();
    for bb in obj.into_iter(){
        geom.add_objects(bb);
    }
    for _ in 0..10000{
        let spos:(f64,f64,f64) = (rgen.gen_range(-10.0..10.0),rgen.gen_range(-10.0..10.0),rgen.gen_range(-10.0..10.0));
        let mut spp = 
        geometry::Geometry::generate_sphere(&spos,0.125,8,8);
        let res = cyl.get_position_of(&Point3D::new(spos.0,spos.1,spos.2));
        if let Some(x) = res{
            
            let mut r:u8 = 255;
            let mut g:u8 = 255;
            let mut b:u8 = 255;
            /*
            if x%6/3 == 0{
                r = x%3*127;
                b = x%3*127;
            }else{
                g = x%3*127;
                b = x%3*127;
            }*/
            if x/6 == 0{
                g = 0;
                b = 0;
            }else if x/6 == 1{
                r = 0;
                b = 0;
            }else{
                r = 0;
                g = 0;
            }

            /*
            let r:u8 = x.0*120;
            let g = 0;
            */
            geometry::Face::color_faces(&mut spp.1,&vec![r,g,b]);
        }
        geom.add_objects(spp);
    }


    geom.add_objects(geometry::Geometry::generate_sphere(&cb.get_xyz(),0.5,8,8));

    geom.calc_all_norms();
    geom.add_colortile_material();
    geom.calc_all_norms();

    geometry::Geometry::save("test/cylindercheck.obj",&mut vec![geom]);    
}

#[test]
fn decision_tree_model_test(){
    let filename = "resources/scripts/results/target_path.dat";
    let file = File::open(filename).unwrap_or_else(|e|panic!("{} {:?}",filename,e));
    let reader = BufReader::new(file);

        
    let exx =  Regex::new(r"^.+(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    let mut entries_:Vec<(String,Option<HashMap<String,f64>>)> = vec![];
    for (_lcount,line) in reader.lines().enumerate() {
        let path = line.unwrap();
        if let Some(x) = exx.captures(path.as_str()){
            entries_.push((x.get(0).unwrap().as_str().to_string(),None));
        }
        if _lcount > 2{
            break;
        }
    }
    //println!("{:?}",entries_);
    entries_.sort_by(|a,b|a.0.cmp(&b.0));

    //generate_decision_tree_model_files(entries_,10.0,3,8.0,8.0,"example_files/example_output/testpcb_tree_model.dat");
    let treemodel:DecisionTreePCBModel = DecisionTreePCBModel::load("example_files/example_output/testpcb_tree_model.dat");
    
    let files = vec![
        "D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-geom_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/LAW_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/LAW_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MASS_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MUFOLD_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS4",

    ];
    for ff in files.iter(){
        let e = mmcif_process::load_pdb(*ff,false);
        let sc = treemodel.calc_score(&e);
        println!("{}\t{}",ff,sc);
    }
}


#[test]
fn pseudocb_model_test(){
    let filename = "resources/scripts/results/target_path.dat";
    let file = File::open(filename).unwrap_or_else(|e|panic!("{} {:?}",filename,e));
    let reader = BufReader::new(file);

        
    let exx =  Regex::new(r"^.+(\.ent|\.pdb|\.cif)(\.gz)?").unwrap();
    let mut entries_:Vec<(String,Option<HashMap<String,f64>>)> = vec![];
    for (_lcount,line) in reader.lines().enumerate() {
        let path = line.unwrap();
        if let Some(x) = exx.captures(path.as_str()){
            entries_.push((x.get(0).unwrap().as_str().to_string(),None));
        }
    }
    //println!("{:?}",entries_);
    entries_.sort_by(|a,b|a.0.cmp(&b.0));

    generate_simplemodel_files(entries_,10.0,3,8.0,8.0,"example_files/example_output/testpcbmodel.dat");
    let l = VerySimplePCBModel::load("example_files/example_output/testpcbmodel.dat"); 
    
    let files = vec![
        "D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-CHEN_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/AWSEM-Suite_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROBETTA_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/BAKER-ROSETTASERVER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/CATHER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-DeepFolder_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-geom_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FALCON-TBM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FEIG-S_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/FoldX_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/GAPF_LNCC_SERVER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/IntFOLD6_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Kiharalab_Z_Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/LAW_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/LAW_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MASS_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MUFOLD_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CLUSTER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-CONSTRUCT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DEEP_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-DIST_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/MULTICOM-HYBRID_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/QUARK_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/RaptorX_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Seok-server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-CaT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/tFold-IDT_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/TOWER_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_FM_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang_TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Yang-Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang_Ab_Initio_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-CEthreader_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS4",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-Server_TS5",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS1",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS2",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS3",
"D:/dummy/work/CASP14/server_stage2/T1026/Zhang-TBM_TS4",

    ];
    for ff in files.iter(){
        let e = mmcif_process::load_pdb(*ff,false);
        let sc = l.calc_score(&e);
        println!("{}\t{}\t{}",ff,sc.0,sc.1);
    }
}
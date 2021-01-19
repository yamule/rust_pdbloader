extern crate regex;

#[allow(unused_imports)]
use std::fs::File;
use std::slice::IterMut;
use std::slice::Iter;
use std::collections::HashSet;
use std::collections::HashMap;
#[allow(unused_imports)]
use std::fs;
#[allow(unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use self::regex::Regex;//module ファイル内で extern crate した場合 self が必要。lib.rs で extern crate すると必要ない。

use super::geometry::Vector3D;
use super::mmcif_process::*;

#[derive(Debug)]
pub enum StructureHierarchyLevel{
Entry,
Model,
Entity,
Asym,
Comp,
Atom
}


lazy_static! {
  static ref REGEX_WS:Regex = Regex::new(r"[\s]").unwrap();
}


fn write_to_file(filename:&str,contents:Vec<String>){
    let mut f = BufWriter::new(fs::File::create(filename).unwrap());
     for ll in contents{
        f.write_all(ll.as_bytes()).unwrap();
        f.write_all("\n".as_bytes()).unwrap();
    }
}


pub fn try_best_format(f:f64,fulllen:usize,afterpoint:usize)->String{
    let mut afterpoint:i64 = afterpoint as i64;
    let mut ret:String = format!("{:1$.2$}",f,fulllen,afterpoint as usize);
    while fulllen < ret.len(){
        afterpoint -= 1;
        if afterpoint < 0{
            panic!("{} can not become string with length less than {}.",f,fulllen);
        }
        ret = format!("{:1$.2$}",f,fulllen,afterpoint as usize);
    }
    return ret;
}


#[derive(Debug,Clone)]
pub struct PDBAtom{
    pub parent_entry:Option<i64>,
    pub parent_entity:Option<i64>,
    pub parent_asym:Option<i64>,
    pub parent_comp:Option<i64>,
    pub index:i64,
    pub serial_number:i64,
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub atom_symbol:String,
    pub atom_code:String,
    pub alt_code:String,
    pub occupancy:f64,
    pub temp_factor:f64,
    pub charge:Option<String>,
    pub dummy:bool,
    pub het:bool,
    pub alt:bool,
    pub is_ligand:bool,
    pub atom_site_key:i64,//PDBEntry.atom_records 内のどの要素からとられたか。-1 は後生成とか。
}

impl Vector3D for PDBAtom{

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

impl PDBAtom{
    pub fn set_serial_number(&mut self,i:i64){
        self.serial_number = i;
    }
    pub fn get_serial_number(&self)->i64{
        return self.serial_number;
    }

    pub fn set_atom_site_key(&mut self,i:i64){
        self.atom_site_key = i;
    }
    
    pub fn get_atom_site_key(&self)->i64{
        return self.atom_site_key;
    }
    
        /*
        COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
    1 -  6        Record name   "ATOM  "
    7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    */ 
    pub fn get_pdb_atom_line_string(&self,chain_name:&str,res_name:&str,res_pos:i64,ins_code:&str)->String{
        let mut label = "ATOM";
        if self.het{
            label = "HETATM";
        }
        let xx = try_best_format(self.x,8,3);
        let yy = try_best_format(self.y,8,3);
        let zz = try_best_format(self.z,8,3);

        let mut aname = self.atom_code.clone();
        if !self.is_ligand && aname.len() < 4{
            aname = " ".to_string()+aname.as_str();
        }

        //atomname についてはスペース入れたり色々する必要がある。。。Ca と CA の書き分けとか
        //mmcif は方向性のあるものは seq_id をつけて、ないものは asym でモノマー、Entity でポリマーにするっぽい
        let ret = format!("{lab:<6}{serial:>5} {atomname:<4}{alt:1}{resname:<3} {chain:1}{pos:>4}{insertion:1}   {x:>8}{y:>8}{z:>8}{occ:>6}{temp:>6}          {symbol:>2}"
            ,lab=label,serial=self.serial_number
            ,atomname=aname
            ,alt=if self.alt_code == "." || self.alt_code == "?" {" "}else{&self.alt_code}
            ,resname=res_name
            ,chain=chain_name
            ,pos=res_pos
            ,insertion=if ins_code == "." || ins_code == "?"{" "}else{ins_code}
            ,x=xx,y=yy,z=zz,
            occ=try_best_format(self.occupancy,6,2),
            temp=try_best_format(self.temp_factor,6,2),symbol=self.atom_symbol);
        if let Some(x) = self.charge.as_ref(){
            if x != "." && x != "?"{
                return ret+x.as_str();
            }
        }
        return ret;
    }
    pub fn new()->PDBAtom{
        return PDBAtom{
            parent_entry:None,
            parent_entity:None,
            parent_asym:None,
            parent_comp:None,
            index:-1,
            serial_number:-1,
            x:0.0,
            y:0.0,
            z:0.0,
            atom_symbol:"".to_string(),
            atom_code:"".to_string(),
            alt_code:"".to_string(),
            occupancy:1.0,
            temp_factor:0.0,
            charge:None,
            dummy:true,
            het:false,
            alt:false,
            is_ligand:false,
            atom_site_key:-1,
        };
    }

    pub fn get_coordinates(&self)->(f64,f64,f64){
        return (self.x,self.y,self.z);
    }

    pub fn get_index(&self)-> i64{
        return self.index;
    }
    fn set_index(&mut self,index:i64){
        self.index = index;
    }
    fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,asym_id:Option<i64>,comp_id:Option<i64>,index:i64){
        
        self.set_index(index);
        self.parent_comp = comp_id;
        self.parent_asym = asym_id;
        self.parent_entity = entity_id;
        self.parent_entry = entry_id;
    }
}

#[derive(Debug)]
pub struct PDBComp{
    pub parent_entry:Option<i64>,
    pub parent_entity:Option<i64>,
    pub parent_asym:Option<i64>,
    index:i64,
    seq_id:i64,//sequence number
    pub comp_id:String,
    atoms:Vec<PDBAtom>,
    pub ins_code:String,
    pub no_label_seq_id:bool//label_seq_id is "."
}
impl PDBComp{
    pub fn new()->PDBComp{
        return PDBComp{
            parent_entry:None,
            parent_entity:None,
            parent_asym:None,
            index:-1,
            seq_id:-1,
            comp_id:"UNK".to_string(),
            atoms:vec![],
            ins_code:"".to_string(),
            no_label_seq_id:false
        };
    }
    pub fn num_atoms(&self)->usize{
        return self.atoms.len();
    }
    pub fn get_copy_wo_parents(&self)->PDBComp{
        let mut ret = PDBComp::new();
        ret.seq_id = self.seq_id;
        ret.comp_id = self.comp_id.clone();
        ret.ins_code = self.ins_code.clone();
        for aa in self.atoms.iter(){
            ret.add_atom(aa.clone());
        }
        return ret;   
    }

    pub fn set_name(&mut self,n:&str){
        self.comp_id = n.to_string();
    }
    pub fn get_name(&self)->&str{
        return &self.comp_id;
    }
    pub fn get_label(&self)->String{
        return "".to_string()+&self.comp_id+&self.get_seq_id().to_string()+&self.get_ins_code();
    }
    
    pub fn add_atom(&mut self,aa:PDBAtom){
        self.atoms.push(aa);
    }

    pub fn remove_atom_by_name(&mut self,n:&str){
        let mut bflag:bool = false;
        loop{
            let mut pos:Vec<usize> = vec![];
            for (ii,aa) in self.atoms.iter().enumerate(){
                if aa.atom_code == n{
                    pos.push(ii);
                }
            }
            if pos.len() == 0{
                break;
            }
            self.atoms.remove(pos[0]);
            bflag = true;
        }
        if !bflag{
            eprintln!("{} was not found in {}.",n,self.get_label());
        }
    }
    
    pub fn get_first_atom_by_name(&self,n:&str)->Option<&PDBAtom>{
        for aa in self.atoms.iter(){
            if aa.atom_code == n{
                return Some(aa);
            }
        }
        return None;
    }

    pub fn get_index(&self)-> i64{
        return self.index;
    }
    fn set_index(&mut self,index:i64){
        self.index = index;
    }

    pub fn get_all_atoms(&mut self)->Vec<&PDBAtom>{
        let mut ret:Vec<&PDBAtom> = vec![];
        for i in 0..self.atoms.len(){
            ret.push(&self.atoms[i]);
        }
        return ret;
    }

    #[allow(non_snake_case)]
    pub fn get_CA_mut(&mut self)->Option<&mut PDBAtom>{
        for aa in self.atoms.iter_mut(){
            if aa.atom_code == "CA" && aa.atom_symbol != "CA"{
                return Some(aa);
            }
        }
        return None;
    }
    
    #[allow(non_snake_case)]
    pub fn get_O_mut(&mut self)->Option<&mut PDBAtom>{
        for aa in self.atoms.iter_mut(){
            if aa.atom_code == "O"{
                return Some(aa);
            }
        }
        return None;
    }
    
    #[allow(non_snake_case)]
    pub fn get_C_mut(&mut self)->Option<&mut PDBAtom>{
        for aa in self.atoms.iter_mut(){
            if aa.atom_code == "C"{
                return Some(aa);
            }
        }
        return None;
    }

    #[allow(non_snake_case)]
    pub fn get_N_mut(&mut self)->Option<&mut PDBAtom>{
        for aa in self.atoms.iter_mut(){
            if aa.atom_code == "N"{
                return Some(aa);
            }
        }
        return None;
    }

    #[allow(non_snake_case)]
    pub fn get_CA(&self)->Option<&PDBAtom>{
        for aa in self.atoms.iter(){
            if aa.atom_code == "CA"{
                return Some(aa);
            }
        }
        return None;
    }
    
    #[allow(non_snake_case)]
    pub fn get_O(&self)->Option<&PDBAtom>{
        for aa in self.atoms.iter(){
            if aa.atom_code == "O"{
                return Some(aa);
            }
        }
        return None;
    }
    
    #[allow(non_snake_case)]
    pub fn get_C(&self)->Option<&PDBAtom>{
        for aa in self.atoms.iter(){
            if aa.atom_code == "C"{
                return Some(aa);
            }
        }
        return None;
    }

    #[allow(non_snake_case)]
    pub fn get_N(&self)->Option<&PDBAtom>{
        for aa in self.atoms.iter(){
            if aa.atom_code == "N"{
                return Some(aa);
            }
        }
        return None;
    }

    pub fn get_atom_num(&self)->usize{
        return self.atoms.len();
    }
    pub fn get_seq_id(&self)->i64{
        return self.seq_id;
    }
    pub fn set_seq_id(&mut self,i:i64){
        self.seq_id = i;
    }
    pub fn get_comp_id(&self)->&str{
        return &self.comp_id;
    }
    pub fn set_comp_id(&mut self,n:&str){
        self.comp_id = n.to_string();
    }
    pub fn get_ins_code(&self)->&str{
        return &self.ins_code;
    }
    pub fn set_ins_code(&mut self,n:&str){
        self.ins_code = n.to_string();
    }
    pub fn get_atom_at(&self,ii:usize)->&PDBAtom{
        return &self.atoms[ii];
    }
    pub fn get_mut_atom_at(&mut self,ii:usize)->&mut PDBAtom{
        return &mut self.atoms[ii];
    }
    pub fn iter_atoms(&self)->Iter<'_, PDBAtom>{
        return self.atoms.iter();
    }
    pub fn iter_mut_atoms(&mut self)->IterMut<'_, PDBAtom>{
        return self.atoms.iter_mut();
    }
    fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,asym_id:Option<i64>,index:i64){
        self.set_index(index);
        self.parent_asym = asym_id;
        self.parent_entity = entity_id;
        self.parent_entry = entry_id;
        for (ii,aa) in self.atoms.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry.clone(),self.parent_entity.clone(),self.parent_asym.clone(),Some(index),ii as i64);
        }
    }
}


#[derive(Debug)]
pub struct PDBAsym{
    pub parent_entry:Option<i64>,
    pub parent_entity:Option<i64>,
    index:i64,
    pub chain_name:String,
    comps:Vec<PDBComp>
}

impl PDBAsym{
    pub fn new(name:&str)->PDBAsym{
        return PDBAsym{
            parent_entry:None,
            parent_entity:None,
            index:-1,
            chain_name:name.to_string(),
            comps:vec![]
        }
    }

    pub fn num_comps(&self) -> usize{
        return self.comps.len();
    }

    pub fn set_chain_name(&mut self,name:&str){
        self.chain_name = name.to_string();
    }

    pub fn get_index(&self)-> i64{
        return self.index;
    }

    fn set_index(&mut self,index:i64){
        self.index = index;
    }

    pub fn iter_comps(&self) -> Iter<PDBComp>{
        return self.comps.iter();
    }

    pub fn iter_mut_comps(&mut self) -> IterMut<PDBComp>{
        return self.comps.iter_mut();
    }
    /**
     * 例題ファイル
pdb1a46.ent.gz	insertion
pdb1a48.ent.gz	altloc
pdb1a4f.ent.gz	altloc
pdb1a4p.ent.gz	altloc
pdb1a4w.ent.gz	insertion
pdb1a4x.ent.gz	altloc
pdb2a40.ent.gz	altloc
pdb2a41.ent.gz	altloc
pdb2a42.ent.gz	altloc
pdb2a45.ent.gz	insertion
     */
    //特に保持したい ALT の文字がある場合指定する。None の場合最初に出てきた文字が扱われる
    pub fn remove_alt(&mut self,retain:Option<&Vec<&str>>){
        let mut retain_str:HashSet<String> = HashSet::new();
        if let Some(x) = retain{
            retain_str = x.iter().map(|m|m.to_string()).collect();
        }
        for rr in self.comps.iter_mut(){
            let mut atomid_vec:HashMap<String,Vec<usize>> = HashMap::new();
            for (ii,aa) in rr.atoms.iter().enumerate(){
                if !atomid_vec.contains_key(&aa.atom_code){
                    atomid_vec.insert(aa.atom_code.clone(),vec![]);
                }
                atomid_vec.get_mut(&aa.atom_code).unwrap().push(ii);
            }
            if retain_str.len() == 0{
                'outer:for (_kk,vv) in atomid_vec.iter(){
                    if vv.len() > 1{
                        for vii in vv.iter(){ 
                            retain_str.insert(rr.atoms[*vii].alt_code.clone());
                            break 'outer;
                        }
                    }
                }
            }
            let mut atoms_remove:HashSet<usize> = HashSet::new();
            for (_kk,vv) in atomid_vec.iter(){
                if vv.len() > 1{
                    let mut ccid:i64 = -1;//retain する STRING がない場合 -1 でその場合先頭が取られる
                    for vii in vv.iter(){ 
                        if retain_str.contains(&rr.atoms[*vii].alt_code){
                            ccid = *vii as i64;
                        }
                    }
                    if ccid < 0{
                        ccid = *(vv.iter().next().unwrap()) as i64;
                    }
                    for vii in vv.iter(){
                        if ccid != *vii as i64{
                            atoms_remove.insert(*vii);
                        }
                    }
                }
            }
            if atoms_remove.len() == 0{
                continue;
            }
            let mut patoms:Vec<PDBAtom> = vec![];
            patoms.append(&mut rr.atoms);
            for (vii,a) in patoms.into_iter().enumerate(){
                if !atoms_remove.contains(&vii){
                    rr.atoms.push(a);
                }
            }
        }
    }
    
    pub fn get_comp_at(&self, i:usize)->&PDBComp{
        return &self.comps[i];
    }

    pub fn get_mut_comp_at(&mut self, i:usize)->&mut PDBComp{
        return &mut self.comps[i];
    }

    pub fn add_comp(&mut self,a:PDBComp){
        self.comps.push(a);
    }

    pub fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,index:i64){
        self.set_index(index);
        self.parent_entity = entity_id;
        self.parent_entry = entry_id;
        for (ii,aa) in self.comps.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry.clone(),self.parent_entity.clone(),Some(index),ii as i64);
        }
    }

    
}


pub struct PDBEntity{
    index:i64,
    pub parent_entry:Option<i64>,
    pub parent_model:Option<i64>,
    pub entity_id:String,
    asyms:Vec<PDBAsym>,
}

#[allow(dead_code)]
impl PDBEntity{
    pub fn new()->PDBEntity{
        return PDBEntity{
            index:-1,
            parent_entry:None,
            parent_model:None,
            entity_id:"".to_string(),
            asyms:vec![]
            };
    }
    
    pub fn squeeze_asyms(&mut self)->Vec<PDBAsym>{
        let mut ret:Vec<PDBAsym> = vec![];
        ret.append(&mut self.asyms);
        return ret;
    }

    pub fn num_asyms(&self)->usize{
        return self.asyms.len();
    }

    pub fn get_asym_at(&self, i:usize)->&PDBAsym{
        return &self.asyms[i];
    }

    pub fn get_mut_asym_at(&mut self, i:usize)->&mut PDBAsym{
        return &mut self.asyms[i];
    }
    
    pub fn iter_asyms(&self) -> Iter<PDBAsym>{
        return self.asyms.iter();
    }

    pub fn iter_mut_asyms(&mut self) -> IterMut<PDBAsym>{
        return self.asyms.iter_mut();
    }

    pub fn get_index(&self)-> i64{
        return self.index;
    }

    fn set_index(&mut self,index:i64){
        self.index = index;
    }

    
    pub fn add_asym(&mut self,chain:PDBAsym){
        self.asyms.push(chain);
    } 

    pub fn set_parent(&mut self,entry_id:Option<i64>,model_id:Option<i64>,index:i64){
        self.set_index(index);
        self.parent_model = model_id;
        self.parent_entry = entry_id;
        for (ii,aa) in self.iter_mut_asyms().enumerate(){
            aa.set_parent(entry_id,Some(index),ii as i64);
        }
    }

    pub fn create_chain(&mut self,name:&str)->(usize,&PDBAsym){
        let newchain = PDBAsym{
            chain_name:name.to_string(),
            parent_entry:self.parent_entry.clone(),
            parent_entity:Some(self.index as i64),
            index:self.asyms.len() as i64,
            comps:vec![]
        } ;
        let index:i64 = newchain.index;
        self.asyms.push(newchain);
        return (index as usize,&self.asyms[index as usize])
    }

}

pub struct PDBModel{
    parent_entry:Option<i64>,
    entities:Vec<PDBEntity>,
    index:i64,
    model_id:String,
}
impl PDBModel{
    pub fn new()->PDBModel{
        return PDBModel{
            parent_entry:None,
            entities:vec![],
            index:-1,
            model_id:"".to_owned(),
            
        };
    }
    pub fn get_entity_at(&self, i:usize)->&PDBEntity{
        return &self.entities[i];
    }

    pub fn get_mut_entity_at(&mut self, i:usize)->&mut PDBEntity{
        return &mut self.entities[i];
    }

    pub fn num_entities(&self)->usize{
        return self.entities.len();
    }

    pub fn iter_entities(&self) -> Iter<PDBEntity>{
        return self.entities.iter();
    }

    pub fn iter_mut_entities(&mut self) -> IterMut<PDBEntity>{
        return self.entities.iter_mut();
    }

    pub fn add_entity(&mut self,entt:PDBEntity){
        self.entities.push(entt);
    }

    pub fn set_index(&mut self,i:i64){
        self.index = i;
    }
    pub fn get_index(&self)->i64{
        return self.index;
    }
    pub fn set_parent(&mut self,entry_id:Option<i64>,index:i64){
        self.parent_entry = entry_id.clone();
        self.set_index(index);
        for (ii,aa) in self.iter_mut_entities().enumerate(){
            aa.set_parent(entry_id,Some(index),ii as i64);
        }
    }
}
pub struct PDBEntry{
    index:i64,
    pub entry_id:String,
    pub description:String,
    pub mmcif_data:Option<MMCIFEntry>,//PDBAtom が、record_key > -1 を持っている場合、ここからの参照とさせる
    models:Vec<PDBModel>,
}
impl PDBEntry{
    pub fn new()->PDBEntry{
        return PDBEntry{
        index:-1,
        entry_id:"".to_string(),
        description:"".to_string(),
        mmcif_data:None,
        models:vec![]
        };
    }
    

    pub fn get_all_asyms(&self)->Vec<&PDBAsym>{
        let mut ret:Vec<&PDBAsym> = vec![];
        for mm in self.iter_models(){
            for ee in mm.iter_entities(){
                for cc in ee.iter_asyms(){
                    ret.push(cc);
                }
            }
        }
        return ret;
    }


    pub fn get_mut_all_asyms(&mut self)->Vec<&mut PDBAsym>{
        let mut ret:Vec<&mut PDBAsym> = vec![];
        for mm in self.iter_mut_models(){
            for ee in mm.iter_mut_entities(){
                for cc in ee.iter_mut_asyms(){
                    ret.push(cc);
                }
            }
        }
        return ret;
    }


    pub fn squeeze_asyms(&mut self)->Vec<PDBAsym>{
        let mut ret:Vec<PDBAsym> = vec![];
        for mm in self.iter_mut_models(){
            for ee in mm.iter_mut_entities(){
                ret.append(&mut ee.squeeze_asyms());
            }
        }
        return ret;
    }

    
    pub fn prepare_base()->PDBEntry{
        let mut ret =  PDBEntry{
        index:-1,
        entry_id:"".to_string(),
        description:"".to_string(),
        mmcif_data:None,
        models:vec![]
        };
        ret.add_model(PDBModel::new());
        ret.get_mut_model_at(0).add_entity(PDBEntity::new());
        return ret;
    }
    
    pub fn get_aa_sequences(&self)->Vec<(String,Vec<String>)>{
        let mut ret:Vec<(String,Vec<String>)> = vec![];
        let aaname:Vec<(&str,&str)> = vec![
            ("ALA","A"),
            ("ARG","R"),
            ("ASN","N"),
            ("ASP","D"),
            ("CYS","C"),
            ("GLN","Q"),
            ("GLU","E"),
            ("GLY","G"),
            ("HIS","H"),
            ("ILE","I"),
            ("LEU","L"),
            ("LYS","K"),
            ("MET","M"),
            ("PHE","F"),
            ("PRO","P"),
            ("SER","S"),
            ("THR","T"),
            ("TRP","W"),
            ("TYR","Y"),
            ("VAL","V")];
        let mapper:HashMap<String,String> = aaname.iter().map(|m|(m.0.to_string(),m.1.to_string())).collect();
        for mm in self.iter_models(){
            for ee in mm.iter_entities(){
                for cc in ee.iter_asyms(){
                    let mut ss:Vec<String> = vec![];
                    for rr in cc.comps.iter(){
                        ss.push(mapper.get(rr.get_comp_id()).unwrap_or(&("X".to_string())).clone());
                    }
                    ret.push((cc.chain_name.clone(),ss));
                }
            }
        }
        return ret;
    }

    pub fn num_models(&self)->usize{
        return self.models.len();
    }
    pub fn get_model_at(&self, i:usize)->&PDBModel{
        return &self.models[i];
    }

    pub fn get_mut_model_at(&mut self, i:usize)->&mut PDBModel{
        return &mut self.models[i];
    }

    pub fn iter_models(&self) -> Iter<PDBModel>{
        return self.models.iter();
    }

    pub fn iter_mut_models(&mut self) -> IterMut<PDBModel>{
        return self.models.iter_mut();
    }

    pub fn add_model(&mut self,entt:PDBModel){
        self.models.push(entt);
    }
    
    pub fn get_index(&self)->i64{
        return self.index;
    }

    pub fn update_downstream_index(&mut self){
        let index = self.get_index();
        for (ii,aa) in self.iter_mut_models().enumerate(){
            aa.set_parent(Some(index),ii as i64);
        }
    }

    //mmcif から numeric な値の入ったインスタンスを作成する
    pub fn prepare_entry(mmcifdata:MMCIFEntry)->PDBEntry{
        
        let mut ret:PDBEntry = PDBEntry::new();

        //ToDo ソートするか、順序を持っておく
        let models:Vec<(String,Vec<usize>)> = mmcifdata.get_model_map(&((0..mmcifdata.get_num_atoms()).into_iter().collect()));
        for (_modk,modv) in models.into_iter(){
            let mut modell = PDBModel::new();
            let entities:Vec<(String,Vec<usize>)> = mmcifdata.get_entity_map(&(modv));
            for (_entk,entv) in entities.into_iter(){
                let mut entt = PDBEntity::new();
                let asyms:Vec<(String,Vec<usize>)> = mmcifdata.get_asym_map(&(entv));
                for (asymk,asymv) in asyms.into_iter(){
                    let mut asymm = PDBAsym::new(&asymk);
                    let comps:Vec<(String,Vec<usize>)> = mmcifdata.get_comp_map(&(asymv));
                    for (compk,compv) in comps.into_iter(){
                        if compv.len() == 0{
                            panic!("{} has no atom??",compk);
                        }
                        let mut rr:PDBComp = PDBComp{
                            parent_entry:None,
                            parent_entity:None,
                            parent_asym:None,
                            index:-1,
                            seq_id:mmcifdata.get_atom_site(compv[0]).get_label_seq_id().parse::<i64>().unwrap_or_else(|_e|
                                if mmcifdata.get_atom_site(compv[0]).get_label_seq_id() == "."{
                                    NO_SEQ_ID
                                }else{
                                    panic!("Can not parse {}!",mmcifdata.get_atom_site(compv[0]).get_label_seq_id());
                                }
                            ),
                            //.unwrap_or_else(|e|panic!("{:?} {}",e,mmcifdata.get_atom_site(compv[0]).get_label_seq_id())),//sequence number
                            comp_id:mmcifdata.get_atom_site(compv[0]).get_label_comp_id().to_string(),
                            atoms:vec![],
                            ins_code:mmcifdata.get_atom_site(compv[0]).get_pdbx_PDB_ins_code().to_string(),
                            no_label_seq_id:mmcifdata.get_atom_site(compv[0]).get_label_seq_id() == "."
                        };
                        for aa in compv.into_iter(){
                            let att = mmcifdata.get_atom_site(aa).atomsite_to_atom(false);
                            rr.add_atom(
                                att
                            );
                        }
                        asymm.add_comp(rr);
                    }
                    entt.add_asym(asymm);
                }
                modell.add_entity(entt);
            }
            ret.add_model(modell);
        }

        for (mii,mm) in ret.iter_mut_models().enumerate(){
            mm.set_index(mii as i64);
            for (eii,ee) in mm.iter_mut_entities().enumerate(){
                ee.set_index(eii as i64);
                for (aii,aas) in ee.iter_mut_asyms().enumerate(){
                    aas.set_index(aii as i64);
                    for (cii,cc) in aas.iter_mut_comps().enumerate(){
                        cc.set_index(cii as i64);
                        for (aii,aa) in cc.iter_mut_atoms().enumerate(){
                            aa.set_index(aii as i64);
                        }
                    }
                }
            }
        }
        ret.mmcif_data = Some(mmcifdata);
        if let Some(x) = ret.mmcif_data.as_mut(){
            x.update_atom_site_index();
        }
        return ret;
    }

    
    pub fn get_pdb_atom_line_string(&self)->Vec<String>{
        let mut ret:Vec<String>  = vec![];
        for mm in self.iter_models(){
            for ee in mm.iter_entities(){
                for cc in ee.iter_asyms(){
                    let cid:&str = &cc.chain_name;
                    for rr in cc.iter_comps(){
                        let rname:&str = &rr.comp_id;
                        let rpos:i64 = rr.seq_id;
                        let inscode:&str = &rr.ins_code;
                        for aa in rr.iter_atoms(){
                            ret.push(aa.get_pdb_atom_line_string(cid,rname,rpos,inscode));
                        }
                    }
                }
            }
        }
        return ret;
    }

    pub fn save(&self,filename:&str){
        let res:Vec<String> = self.get_pdb_atom_line_string();
        write_to_file(filename,res);
    }

}


pub fn slice_to_string(chrs:&Vec<char>,start:usize,end_pls_one_:usize)->String{
    let mut end_pls_one = end_pls_one_;
    if start >= chrs.len(){
        return "".to_string();
    }

    if end_pls_one > chrs.len(){
        end_pls_one = chrs.len();
    }
    let ret:String = chrs[start..end_pls_one].iter().fold("".to_string(),|s,m|{s+m.to_string().as_str()});
    return (*REGEX_WS.replace_all(&ret, "")).to_string();
}



#[test]
fn slicetest(){
    let chkk:Vec<char> = vec!['a','b','c','d'];
    let q:String = (&chkk[1..3]).iter().fold("".to_string(),|s,m|{s+m.to_string().as_str()});
    println!("{:?}",&chkk[1..3]);
    println!("{}",q);
    
}

#[test]
fn atomlineloadtest(){
    let mut pdb = load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/alt_example_1a48.pdb");
    for mm in pdb.iter_mut_models(){
        for ee in mm.iter_mut_entities(){
            for cc in ee.iter_mut_asyms(){
                for rr in cc.comps.iter(){
                    if rr.seq_id != 17{
                        continue;
                    }
                    for aa in rr.atoms.iter(){
                        println!("{} {} {} {}",rr.comp_id,aa.atom_code,aa.serial_number,aa.alt_code);
                    }
                }
                cc.remove_alt(Some(&vec!["B"]));
                for rr in cc.comps.iter(){
                    if rr.seq_id != 17{
                        continue;
                    }
                    for aa in rr.atoms.iter(){
                        println!("{} {} {} {}",rr.comp_id,aa.atom_code,aa.serial_number,aa.alt_code);
                    }
                }
            }
        }
    }
    pdb.save("test/testout.pdb");
}

extern crate regex;
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

#[derive(Debug)]
enum StructureHierarchyLevel{
Entry,
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
fn start_with(target:&str,fragment:&str)-> bool{
    if let Some(x) = target.find(fragment){
        if x == 0{
            return true;
        }
    }
    return false;
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
    pub atom_record_key:i64,//PDBEntry.atom_records 内のどの要素からとられたか。-1 は後生成とか。
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

    pub fn set_atom_record_key(&mut self,i:i64){
        self.atom_record_key = i;
    }
    
    pub fn get_atom_record_key(&self)->i64{
        return self.atom_record_key;
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
        
        let ret = format!("{lab:<6}{serial:>5} {atomname:<4}{alt:1}{resname:<3} {chain:1}{pos:>4}{insertion:1}   {x:>8}{y:>8}{z:>8}{occ:>6}{temp:>6}          {symbol:>2}"
            ,lab=label,serial=self.serial_number,atomname=aname,alt=self.alt_code,resname=res_name,chain=chain_name
            ,pos=res_pos,insertion=ins_code,x=xx,y=yy,z=zz,
            occ=try_best_format(self.occupancy,6,2),
            temp=try_best_format(self.temp_factor,6,2),symbol=self.atom_symbol);
        if let Some(x) = self.charge.as_ref(){
            return ret+x.as_str();
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
            atom_record_key:-1,
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
    fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,asym_id:Option<i64>,comp_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_comp{
            if x != comp_id.unwrap(){
                if !force{
                    panic!("{:?} already has parent.",self);
                }
            }
        }
        self.set_index(index);
        self.parent_comp = comp_id.clone();
        self.parent_asym = asym_id.clone();
        self.parent_entity = entity_id.clone();
        self.parent_entry = entry_id.clone();
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
            ins_code:"".to_string()
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
    
    pub fn add_atom(&mut self,mut aa:PDBAtom){
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
    fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,asym_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_entry{
            if let Some(y) = entry_id{
                if x != y && !force{
                    panic!("{:?} already has parent entry_id.",self);
                }
            }
        }
        if let Some(x) = self.parent_entity{
            if let Some(y) = entity_id{
                if x != y && !force{
                    panic!("{:?} already has parent entity_id.",self);
                }
            }
        }
        if let Some(x) = self.parent_asym{
            if let Some(y) = asym_id{
                if x != y && !force{
                    panic!("{:?} already has parent asym_id.",self);
                }
            }
        }
        self.set_index(index);
        self.parent_asym = asym_id;
        self.parent_entity = entity_id;
        self.parent_entry = entry_id;
        let sii = self.get_index();
        for (ii,aa) in self.atoms.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry.clone(),self.parent_entity.clone(),self.parent_asym.clone(),Some(sii),ii as i64,force);
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

    pub fn add_comp(&mut self,mut a:PDBComp){
        self.comps.push(a);
    }

    pub fn set_parent(&mut self,entry_id:Option<i64>,entity_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_entry{
            if let Some(y) = entry_id{
                if x != y && !force{
                    panic!("{:?} already has parent.",self);
                }
            }
        }
        if let Some(x) = self.parent_entity{
            if let Some(y) = entity_id{
                if x != y && !force{
                    panic!("{:?} already has parent.",self);
                }
            }
        }
        
        self.set_index(index);
        self.parent_entity = entity_id;
        self.parent_entry = entry_id;
        let sii = self.get_index();
        for (ii,aa) in self.comps.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry.clone(),self.parent_entity.clone(),Some(sii),ii as i64,force);
        }
    }

    
}


pub struct PDBEntity{
    index:i64,
    pub parent_entry:Option<i64>,
    pub entity_id:String,
    asyms:Vec<PDBAsym>,
}

#[allow(dead_code)]
impl PDBEntity{
    pub fn new()->PDBEntity{
        return PDBEntity{
            index:-1,
            parent_entry:None,
            entity_id:"".to_string(),
            asyms:vec![]
            };
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
        for cc in self.asyms.iter(){
            let mut ss:Vec<String> = vec![];
            for rr in cc.comps.iter(){
                ss.push(mapper.get(rr.get_comp_id()).unwrap_or(&("X".to_string())).clone());
            }
            ret.push((cc.chain_name.clone(),ss));
        }
        return ret;
    }
    
    pub fn add_asym(&mut self,mut chain:PDBAsym){
        self.asyms.push(chain);
    } 

    pub fn set_parent(&mut self,entry_id:Option<i64>,index:i64,force:bool){
        self.set_index(index);
        self.parent_entry = entry_id;
        let sii = self.get_index();
        for (ii,aa) in self.iter_mut_asyms().enumerate(){
            aa.set_parent(self.parent_entry.clone(),Some(sii),ii as i64,force);
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


pub struct PDBEntry{
    index:i64,
    pub entry_id:String,
    pub description:String,
    pub atom_records:Vec<AtomRecord>,//PDBAtom が、record_key > -1 を持っている場合、ここからの参照とさせる
    entities:Vec<PDBEntity>,
}
impl PDBEntry{
    pub fn new()->PDBEntry{
        return PDBEntry{
        index:-1,
        entry_id:"".to_string(),
        description:"".to_string(),
        atom_records:vec![],
        entities:vec![]
        };
    }
    
    pub fn get_entity_at(&self, i:usize)->&PDBEntity{
        return &self.entities[i];
    }

    pub fn get_mut_entity_at(&mut self, i:usize)->&mut PDBEntity{
        return &mut self.entities[i];
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

    pub fn update_downstream_index(&mut self){
        for (ii,aa) in self.iter_mut_entities().enumerate(){
            aa.set_parent(Some(self.index),ii as i64,true);
        }
    }

    pub fn prepare_entry(atom_records:Vec<AtomRecord>)->PDBEntry{
        
        let mut ret:PDBEntry = PDBEntry::new();
        ret.atom_records = atom_records;
        for (ee,aa) in ret.atom_records.iter_mut().enumerate(){
            aa.set_index(ee as i64);
        }

        let entities:HashMap<String,Vec<&AtomRecord>> = PDBEntry::get_entity_map(&(ret.atom_records.iter().collect()));
        for (entk,entv) in entities.into_iter(){
            let mut entt = PDBEntity::new();
            let asyms:HashMap<String,Vec<&AtomRecord>> = PDBEntry::get_asym_map(&(entv.iter().map(|m|m.clone()).collect()));
            for (asymk,asymv) in asyms.into_iter(){
                let mut asymm = PDBAsym::new(&asymk);
                let comps:HashMap<String,Vec<&AtomRecord>> = PDBEntry::get_comp_map(&(asymv.iter().map(|m|m.clone()).collect()));
                for (compk,compv) in comps.into_iter(){
                    if compv.len() == 0{
                        panic!("{} has no atom??",compk);
                    }
                    let mut rr:PDBComp = PDBComp{
                           parent_entry:None,
                            parent_entity:None,
                            parent_asym:None,
                            index:-1,
                            seq_id:compv[0].label_seq_id.parse::<i64>().unwrap(),//sequence number
                            comp_id:compv[0].label_comp_id.clone(),
                            atoms:vec![],
                            ins_code:compv[0].pdbx_PDB_ins_code.clone(),
                    };
                    for aa in compv.into_iter(){
                        let mut att = aa.atomrecord_to_atom();
                        rr.add_atom(
                            aa.atomrecord_to_atom()
                        );
                    }
                    asymm.add_comp(rr);
                }
                entt.add_asym(asymm);
            }
            ret.add_entity(entt);
        }

        
        for (eii,ee) in ret.iter_mut_entities().enumerate(){
            ee.set_index(eii as i64);
            for (aii,aas) in ee.iter_mut_asyms().enumerate(){
                aas.set_index(aii as i64);
                for (cii,cc) in aas.iter_mut_comps().enumerate(){
                    cc.set_index(cii as i64);
                    for (aii,aa) in cc.iter_atoms().enumerate(){
                        aa.set_index(aii as i64);
                    }
                }
            }
        }
        
        return ret;
    }


    pub fn get_entity_map<'a>(atom_records:&Vec<&'a AtomRecord>)->HashMap<String,Vec<&'a AtomRecord>>{
        return PDBEntry::get_map(atom_records,StructureHierarchyLevel::Entity);
    }
    
    pub fn get_asym_map<'a>(atom_records:&Vec<&'a AtomRecord>)->HashMap<String,Vec<&'a AtomRecord>>{
        return PDBEntry::get_map(atom_records,StructureHierarchyLevel::Asym);
    }
    
    pub fn get_comp_map<'a>(atom_records:&Vec<&'a AtomRecord>)->HashMap<String,Vec<&'a AtomRecord>>{
        return PDBEntry::get_map(atom_records,StructureHierarchyLevel::Comp);
    }

    pub fn get_map<'a>(atom_records:&Vec<&'a AtomRecord>,lev:StructureHierarchyLevel)->HashMap<String,Vec<&'a AtomRecord>>{
        let mut ret:HashMap<String,Vec<&AtomRecord>> = HashMap::new();
        for aa in atom_records.iter(){
            let ak:String = match lev{
                StructureHierarchyLevel::Entity => aa.label_entity_id.clone(),
                StructureHierarchyLevel::Asym => aa.label_asym_id.clone(),
                StructureHierarchyLevel::Comp => aa.get_unique_comp_label(),
                _=> panic!("{:?} is not allowed in this function.",lev)
            };
            if !ret.contains_key(&ak){
                ret.insert(ak.clone(),vec![]);
            }
            ret.get_mut(&ak).unwrap().push(aa);
        }
        return ret;
    }
    
    pub fn get_pdb_atom_line_string(&self)->Vec<String>{
        let mut ret:Vec<String>  = vec![];
        for ee in self.iter_entities(){
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
/**
 * 
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.pdbx_formal_charge 
_atom_site.auth_seq_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_atom_id 
_atom_site.pdbx_PDB_model_num 
ATOM   1    N N   . MET A 1 1   ? -17.775 41.398 46.220  1.00 26.26  ? 1   MET A N   1 
**/

#[derive(Debug)]
pub struct AtomRecord{
    index_in_entry:i64,//PDBEntry 内の atom_records 内のインデクス
    pub group_PDB:String,
    pub id:String,
    pub type_symbol:String,
    pub label_atom_id:String,
    pub label_alt_id:String,
    pub label_comp_id:String,
    pub label_asym_id:String,
    pub label_entity_id:String,
    pub label_seq_id:String,
    pub pdbx_PDB_ins_code:String,
    pub Cartn_x:String,
    pub Cartn_y:String,
    pub Cartn_z:String,
    pub occupancy:String,
    pub B_iso_or_equiv:String,
    pub pdbx_formal_charge:String,
    pub auth_seq_id:String,
    pub auth_comp_id:String,
    pub auth_asym_id:String,
    pub auth_atom_id:String,
    pub pdbx_PDB_model_num:String,
}

impl AtomRecord{

    pub fn set_index(&mut self,i:i64){
        self.index_in_entry = i;
    }
    
    pub fn get_index(self)-> i64{
        return self.index_in_entry;
    }
    
    pub fn atomrecord_to_atom(&self)->PDBAtom{
        let ret:PDBAtom = PDBAtom{
            parent_entry:None,
            parent_entity:None,
            parent_asym:None,
            parent_comp:None,
            index:-1,
            serial_number:self.label_atom_id.parse::<i64>().unwrap().clone(),
            x:self.Cartn_x.parse::<f64>().unwrap().clone(),
            y:self.Cartn_y.parse::<f64>().unwrap().clone(),
            z:self.Cartn_z.parse::<f64>().unwrap().clone(),
            charge:if self.pdbx_formal_charge.len() > 0{Some(self.pdbx_formal_charge.clone())}else{None},
            occupancy:if self.occupancy.len() > 0{self.occupancy.parse::<f64>().unwrap()}else{1.0},
            temp_factor:if self.B_iso_or_equiv.len() > 0{self.B_iso_or_equiv.parse::<f64>().unwrap()}else{0.0},
            atom_symbol:self.type_symbol.clone(),
            atom_code:self.label_atom_id,
            alt_code:self.label_alt_id.clone(),
            dummy:true,
            het:self.is_het(),
            alt:self.is_alt(),
            is_ligand:false,
            atom_record_key:self.index_in_entry
        };//alt_loc は他の Atom も見ないと処理できないと思う
            
        return ret;
    }

    pub fn is_het(&self)->bool{
        if &self.group_PDB == "ATOM"{
            return false;
        }
        return true;
    }

    pub fn is_alt(&self)->bool{
        if &self.label_alt_id == ""
        || &self.label_alt_id == "."
        || &self.label_alt_id == "?"
        || &self.label_alt_id == " "
        || &self.label_alt_id == "A"{
            return false;
        }
        return true;
    }

    pub fn new()->AtomRecord{
        return AtomRecord{
            index_in_entry:-1,
            group_PDB:"".to_owned(),
            id:"".to_owned(),
            type_symbol:"".to_owned(),
            label_atom_id:"".to_owned(),
            label_alt_id:"".to_owned(),
            label_comp_id:"".to_owned(),
            label_asym_id:"".to_owned(),
            label_entity_id:"".to_owned(),
            label_seq_id:"".to_owned(),
            pdbx_PDB_ins_code:"".to_owned(),
            Cartn_x:"".to_owned(),
            Cartn_y:"".to_owned(),
            Cartn_z:"".to_owned(),
            occupancy:"".to_owned(),
            B_iso_or_equiv:"".to_owned(),
            pdbx_formal_charge:"".to_owned(),
            auth_seq_id:"".to_owned(),
            auth_comp_id:"".to_owned(),
            auth_asym_id:"".to_owned(),
            auth_atom_id:"".to_owned(),
            pdbx_PDB_model_num:"".to_owned(),
        };
    }    
    
    pub fn get_unique_comp_label(&self)->String{
        return "".to_string()+self.label_comp_id.as_str()
        +"#"+self.label_seq_id.as_str()
        +"#"+self.pdbx_PDB_ins_code.as_str();
    }
    pub fn parse_atom_line_pdb(line:&str,model_code:&str)-> AtomRecord{
        let mut ret = AtomRecord::new();

        let pt:Vec<char> = line.chars().collect();
        let label = slice_to_string(&pt,0,6);
        ret.group_PDB = label;

        let idd = slice_to_string(&pt,6,11);
        ret.id = idd;

        let atomid = slice_to_string(&pt,12,16);
        ret.label_atom_id = atomid.clone();
        ret.auth_atom_id = atomid;

        let altloc = slice_to_string(&pt,16,17);
        ret.label_alt_id = altloc;

        let resname = slice_to_string(&pt,17,20);
        ret.label_comp_id = resname.clone();
        ret.auth_comp_id = resname;

        let chainid = slice_to_string(&pt,21,22);
        ret.label_asym_id = chainid.clone();
        ret.auth_asym_id = chainid;

        let respos = slice_to_string(&pt,22,26);
        ret.label_seq_id = respos.clone();
        ret.auth_seq_id = respos;

        let icode = slice_to_string(&pt,26,27);
        ret.pdbx_PDB_ins_code = icode;

        let x = slice_to_string(&pt,30,38);
        ret.Cartn_x = x;

        let y = slice_to_string(&pt,38,46);
        ret.Cartn_y = y;

        let z = slice_to_string(&pt,46,54);
        ret.Cartn_z = z;

        let occupancy = slice_to_string(&pt,54,60);
        if occupancy.len() > 0{
            ret.occupancy = occupancy;
        }

        let tempfactor = slice_to_string(&pt,60,66);
        if tempfactor.len() > 0{
            ret.B_iso_or_equiv = tempfactor;
        }

        let element = slice_to_string(&pt,76,78);
        if element.len() > 0{
            ret.type_symbol = element;
        }

        let charge = slice_to_string(&pt,78,80);
        if charge.len() > 0{
            ret.pdbx_formal_charge = charge;
        }
        ret.pdbx_PDB_model_num = model_code.to_string();

        return ret;
    }
}




pub fn load_pdb(filename:&str) ->PDBEntry{
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    
    //let mut lcount:i64 = 0;
    let mut _ret:Vec<Vec<String>> = Vec::new();

    let _noline = Regex::new(r"^[\r\n]*$").unwrap();
    let mut records:Vec<AtomRecord> = vec![];
    let mut ligand_records:Vec<AtomRecord> = vec![];
    let mut terflag:bool = false;
    let mut possibly_ligand = false;
    let mut current_model_num:String = "".to_owned();
    for (_lcount,line) in reader.lines().enumerate() {

        let sstr = line.unwrap_or_else(|e|panic!("{:?}",e));
        
        if start_with(&sstr,"MODEL"){
            let ppt: Vec<&str> = sstr.split("[\\s;]+").collect();
            current_model_num = ppt[1].to_owned();
            terflag = false;
        }
        if start_with(&sstr,"ATOM") || start_with(&sstr,"HETATM"){
            let mut arecord = AtomRecord::parse_atom_line_pdb(&sstr,&current_model_num);
            //println!("{:?}",arecord);
            if terflag{
                possibly_ligand = true;
            }
            records.push(arecord);
        }else if start_with(&sstr,"TER"){
            terflag = true;
        }else{

        }
    }
    if possibly_ligand{
        eprintln!("There are possibly ligand records. ");
    }

    let mut ret_:PDBEntry = PDBEntry::prepare_entry(records);
    //?? こうしないとボローチェッカー通らない。。。
    
    let mut ret:PDBEntry = PDBEntry::new();
    ret.entities.append(&mut ret_.entities);
    return ret;
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
    for ee in pdb.iter_mut_entities(){
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
    pdb.save("test/testout.pdb");
}

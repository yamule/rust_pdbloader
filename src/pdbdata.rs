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
    pub parent_chain:Option<i64>,
    pub parent_residue:Option<i64>,
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
    pub external_array_pointer:i64,//mmcif とか外部データから取られた場合マップに使用される
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

    pub fn set_external_array_pointer(&mut self,i:i64){
        self.external_array_pointer = i;
    }
    
    pub fn get_external_array_pointer(&self)->i64{
        return self.external_array_pointer;
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
            parent_chain:None,
            parent_residue:None,
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
            external_array_pointer:-1,
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
    fn set_parent(&mut self,entry_id:Option<i64>,chain_id:Option<i64>,residue_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_residue{
            if x != residue_id.unwrap(){
                if !force{
                    panic!("{:?} already has parent.",self);
                }
            }
        }
        self.set_index(index);
        self.parent_residue = residue_id.clone();
        self.parent_chain = chain_id.clone();
        self.parent_entry = entry_id.clone();
    }
}

#[derive(Debug)]
pub struct PDBResidue{
    pub parent_entry:Option<i64>,
    pub parent_chain:Option<i64>,
    index:i64,
    residue_number:i64,//sequence number
    pub residue_name:String,
    atoms:Vec<PDBAtom>,
    pub ins_code:String,
}
impl PDBResidue{
    pub fn new()->PDBResidue{
        return PDBResidue{
            parent_entry:None,
            parent_chain:None,
            index:-1,
            residue_number:-1,
            residue_name:"UNK".to_string(),
            atoms:vec![],
            ins_code:"".to_string()
        };
    }

    pub fn get_copy_wo_parents(&self)->PDBResidue{
        let mut ret = PDBResidue::new();
        ret.residue_number = self.residue_number;
        ret.residue_name = self.residue_name.clone();
        ret.ins_code = self.ins_code.clone();
        for aa in self.atoms.iter(){
            ret.add_atom(aa.clone(), true);
        }
        return ret;   
    }

    pub fn set_name(&mut self,n:&str){
        self.residue_name = n.to_string();
    }
    pub fn get_name(&self)->&str{
        return &self.residue_name;
    }
    pub fn get_label(&self)->String{
        return "".to_string()+&self.residue_name+&self.get_residue_number().to_string()+&self.get_ins_code();
    }
    
    pub fn add_atom(&mut self,mut aa:PDBAtom,force:bool){
        //let llen = self.atoms.len().clone() as i64;
        let alen = self.atoms.len();
        let ai = self.get_index();
        aa.set_parent(self.parent_entry.clone()
        ,self.parent_chain.clone()
        ,Some(ai)
        ,alen.clone() as i64,force);
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
        for (aii,aa) in self.atoms.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry.clone()
            ,self.parent_chain.clone()
            ,Some(index)
            ,aii as i64,true);
        }
    }
    pub fn num_atoms(&self)->usize{
        return self.atoms.len();
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
    pub fn get_residue_number(&self)->i64{
        return self.residue_number;
    }
    pub fn set_residue_number(&mut self,i:i64){
        self.residue_number = i;
    }
    pub fn get_residue_name(&self)->&str{
        return &self.residue_name;
    }
    pub fn set_residue_name(&mut self,n:&str){
        self.residue_name = n.to_string();
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
    fn set_parent(&mut self,entry_id:Option<i64>,chain_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_chain{
            if !force{
                panic!("{:?} already has parent. {:?}",self,x);
            }
        }
        self.set_index(index);
        self.parent_chain = chain_id;
        self.parent_entry = entry_id;
        let sii = self.get_index();
        for (ii,aa) in self.atoms.iter_mut().enumerate(){
            aa.set_parent(entry_id.clone(),chain_id.clone(),Some(sii),ii as i64,force);
        }
    }
}


#[derive(Debug)]
pub struct PDBChain{
    pub parent_entry:Option<i64>,
    index:i64,
    pub chain_name:String,
    pub residues:Vec<PDBResidue>
}

impl PDBChain{
    pub fn new(name:&str)->PDBChain{
        return PDBChain{
            parent_entry:None,
            index:-1,
            chain_name:name.to_string(),
            residues:vec![]
        }
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
        for rr in self.residues.iter_mut(){
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
    
    pub fn add_residue(&mut self,mut a:PDBResidue,force:bool){
        a.set_parent(self.parent_entry.clone(),Some(self.get_index()),self.residues.len() as i64,force);
        self.residues.push(a);
    }

    pub fn set_parent(&mut self,entry_id:Option<i64>,index:i64,force:bool){
        if let Some(x) = self.parent_entry{
            if let Some(y) = entry_id{

                if x != y && !force{
                    panic!("{:?} already has parent.",self);
                }
            }
        }
        self.set_index(index);
        self.parent_entry = entry_id;
        let sii = self.get_index();
        for (ii,aa) in self.residues.iter_mut().enumerate(){
            aa.set_parent(self.parent_entry,Some(sii),ii as i64,force);
        }
    }

    
}

pub struct PDBEntry{
    index:i64,
    pub entry_id:String,
    pub description:String,
    pub chains:Vec<PDBChain>,
}

#[allow(dead_code)]
impl PDBEntry{
    pub fn new()->PDBEntry{
        return PDBEntry{
            index:-1,
            entry_id:"".to_string(),
            description:"".to_string(),
            chains:vec![]
            };
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
        for cc in self.chains.iter(){
            let mut ss:Vec<String> = vec![];
            for rr in cc.residues.iter(){
                ss.push(mapper.get(rr.get_residue_name()).unwrap_or(&("X".to_string())).clone());
            }
            ret.push((cc.chain_name.clone(),ss));
        }
        return ret;
    }
    pub fn update_downstream_index(&mut self){
        for (ii,aa) in self.chains.iter_mut().enumerate(){
            aa.set_parent(Some(self.index),ii as i64,true);
        }
    }
    

    pub fn add_chain(&mut self,mut chain:PDBChain,force:bool){
        chain.set_parent(Some(self.index),self.chains.len() as i64, force);
        self.chains.push(chain);
    } 

    pub fn create_chain(&mut self,name:&str)->(usize,&PDBChain){
        let newchain = PDBChain{
            chain_name:name.to_string(),
            parent_entry:Some(self.index as i64),
            index:self.chains.len() as i64,
            residues:vec![]
        } ;
        let index:i64 = newchain.index;
        self.chains.push(newchain);
        return (index as usize,&self.chains[index as usize])
    }

    pub fn get_pdb_atom_line_string(&self)->Vec<String>{
        let mut ret:Vec<String>  = vec![];
        for cc in self.chains.iter(){
            let cid:&str = &cc.chain_name;
            for rr in cc.residues.iter(){
                let rname:&str = &rr.residue_name;
                let rpos:i64 = rr.residue_number;
                let inscode:&str = &rr.ins_code;
                for aa in rr.atoms.iter(){
                    ret.push(aa.get_pdb_atom_line_string(cid,rname,rpos,inscode));
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

#[derive(Debug)]
pub struct AtomRecord{
    pub is_het:bool,
    pub atom_id:i64,
    pub atom_name:String,
    pub alt_loc:String,
    pub residue_name:String,
    pub chain_id:String,
    pub residue_pos:i64,
    pub ins_code:String,
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub occupancy:Option<f64>,
    pub temp_factor:Option<f64>,
    pub element:Option<String>,
    pub charge:Option<String>,
    pub is_auth_data:bool,
    pub is_ligand:bool,

}

impl AtomRecord{
    pub fn atomrecord_to_atom(&self)->PDBAtom{
        let ret:PDBAtom = PDBAtom{
            parent_entry:None,
            parent_chain:None,
            parent_residue:None,
            index:-1,
            serial_number:self.atom_id.clone(),
            x:self.x.clone(),
            y:self.y.clone(),
            z:self.z.clone(),
            charge:self.charge.clone(),
            occupancy:self.occupancy.unwrap_or(1.0),
            temp_factor:self.temp_factor.unwrap_or(0.0),
            atom_symbol:self.element.clone().unwrap_or("".to_string()),
            atom_code:self.atom_name.clone(),
            alt_code:self.alt_loc.clone(),
            dummy:true,
            het:self.is_het.clone(),
            alt:false,
            is_ligand:self.is_ligand,
            external_array_pointer:-1//出現場所を入れた方が良いか？
        };//alt_loc は他の Atom も見ないと処理できないと思う
            
            return ret;
    }
    pub fn atomrecords_to_entry(atom_records:&Vec<AtomRecord>)->PDBEntry{
        
        let mut ret:PDBEntry = PDBEntry::new();
        let mut chains:HashMap<String,Vec<Vec<&AtomRecord>>> = HashMap::new();
        let mut chains_rindex:HashMap<String,HashMap<String,usize>> = HashMap::new();
        
        for aa in atom_records.iter(){
            if !chains.contains_key(&aa.chain_id){
                chains.insert(aa.chain_id.clone(),vec![]);
                chains_rindex.insert(aa.chain_id.clone(),HashMap::new());
            }

            
            let r_label:String = aa.get_unique_residue_label();
            if !chains_rindex.get(aa.chain_id.as_str()).unwrap().contains_key(r_label.as_str()){
                chains_rindex.get_mut(aa.chain_id.as_str()).unwrap().insert(r_label.clone(),chains.get(aa.chain_id.as_str()).unwrap().len());
                chains.get_mut(aa.chain_id.as_str()).unwrap().push(vec![]);
            }
            let rpos:&usize = chains_rindex.get_mut(aa.chain_id.as_str()).unwrap().get(r_label.as_str()).unwrap();
            
            chains.get_mut(aa.chain_id.as_str()).unwrap()[*rpos].push(aa);
        }
        for (_cc,rr) in chains.into_iter(){
            let mut cc:PDBChain = PDBChain{
                parent_entry:None,
                index:-1,
                chain_name:_cc.to_string(),
                residues:vec![]
            };
            for vv in rr.into_iter(){
                
                let atoms_v:Vec<&AtomRecord> = vv;
                let mut rr:PDBResidue = PDBResidue{
                        parent_entry:None,
                        parent_chain:None,
                        index:-1,
                        residue_number:atoms_v[0].residue_pos.clone(),//sequence number
                        residue_name:atoms_v[0].residue_name.clone(),
                        atoms:vec![],
                        ins_code:atoms_v[0].ins_code.clone(),
                };
                for aa in atoms_v.iter(){
                    let att:PDBAtom = aa.atomrecord_to_atom();
                    rr.add_atom(att,true);
                }
                cc.add_residue(rr,true);
            }
            ret.add_chain(cc,true);
        }

        for (cii,cc) in ret.chains.iter().enumerate(){
            for (rii,rr) in cc.residues.iter().enumerate(){
                assert_eq!(cii as i64,rr.parent_chain.unwrap());
                for (_aii,aa) in rr.iter_atoms().enumerate(){
                    assert_eq!(cii as i64,aa.parent_chain.unwrap());
                    assert_eq!(rii as i64,aa.parent_residue.unwrap());
                    //println!("{} {} {} {} {} {} ",cc.chain_name,rr.residue_name,aa.atom_code,aa.get_x(),aa.get_y(),aa.get_z());
                }
            }
        }
        return ret;
    }


    pub fn new()->AtomRecord{
        return AtomRecord{
            is_het:false,
            atom_id:-999,
            atom_name:"".to_string(),
            alt_loc:"".to_string(),
            residue_name:"".to_string(),
            chain_id:"".to_string(),
            residue_pos:-999,
            ins_code:"".to_string(),
            x:-999.0,
            y:-999.0,
            z:-999.0,
            occupancy:None,
            temp_factor:None,
            element:None,
            charge:None,
            is_auth_data:true,
            is_ligand:false,
        };
    }    
    
    pub fn get_unique_residue_label(&self)->String{
        return "".to_string()+self.residue_name.as_str()
        +"#"+self.residue_pos.to_string().as_str()
        +"#"+self.ins_code.as_str();
    }
    pub fn parse_atom_line(line:&str)-> AtomRecord{
        let mut ret = AtomRecord::new();

        let pt:Vec<char> = line.chars().collect();
        let label = slice_to_string(&pt,0,6);
        if label == "HETATM"{
            ret.is_het = true;
        }else if label == "ATOM"{
        }else{
            panic!("Label {} is not defined.",label);
        }

        let atomid = slice_to_string(&pt,6,11);
        ret.atom_id = (&atomid).parse::<i64>().unwrap_or_else(|_|{panic!("Can not parse {}.",atomid);});

        let atomname = slice_to_string(&pt,12,16);
        ret.atom_name = atomname;

        let altloc = slice_to_string(&pt,16,17);
        ret.alt_loc = altloc;

        let resname = slice_to_string(&pt,17,20);
        ret.residue_name = resname;

        let chainid = slice_to_string(&pt,21,22);
        ret.chain_id = chainid;

        let respos = slice_to_string(&pt,22,26);
        ret.residue_pos = (&respos).parse::<i64>().unwrap_or_else(|_|{panic!("Can not parse {}.",respos);});

        let icode = slice_to_string(&pt,26,27);
        ret.ins_code = icode;

        let x = slice_to_string(&pt,30,38);
        ret.x = (&x).parse::<f64>().unwrap_or_else(|_|{panic!("Can not parse {}.",x);});

        let y = slice_to_string(&pt,38,46);
        ret.y = (&y).parse::<f64>().unwrap_or_else(|_|{panic!("Can not parse {}.",y);});

        let z = slice_to_string(&pt,46,54);
        ret.z = (&z).parse::<f64>().unwrap_or_else(|_|{panic!("Can not parse {}.",z);});

        let occupancy = slice_to_string(&pt,54,60);
        if occupancy.len() > 0{
            ret.occupancy = Some((&occupancy).parse::<f64>().unwrap_or_else(|_|{panic!("Can not parse {}.",occupancy);}));
        }

        let tempfactor = slice_to_string(&pt,60,66);
        if tempfactor.len() > 0{
            ret.temp_factor = Some((&tempfactor).parse::<f64>().unwrap_or_else(|_|{panic!("Can not parse {}.",tempfactor);}));
        }
        let element = slice_to_string(&pt,76,78);
        if element.len() > 0{
            ret.element = Some(element);
        }

        let charge = slice_to_string(&pt,78,80);
        if charge.len() > 0{
            ret.charge = Some(charge);
        }
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
    let mut terflag:bool = false;
    let mut modelflag:bool = false;
    for (_lcount,line) in reader.lines().enumerate() {

        let sstr = line.unwrap_or_else(|e|panic!("{:?}",e));
        
        if start_with(&sstr,"MODEL"){
            if modelflag{
                eprintln!("Multi model is not supported.");
                break;
            }
            modelflag = true;
        }
        if start_with(&sstr,"ATOM"){
            let mut arecord = AtomRecord::parse_atom_line(&sstr);
            //println!("{:?}",arecord);
            if terflag{
                arecord.is_ligand = true;
            }
            records.push(arecord);
        }else if start_with(&sstr,"HETATM"){
            let mut arecord = AtomRecord::parse_atom_line(&sstr);
            //println!("{:?}",arecord);
            if terflag{
                arecord.is_ligand = true;
            }
            records.push(arecord);
        }else if start_with(&sstr,"TER"){
            terflag = true;
        }else{

        }
    }
    let mut ret_:PDBEntry = AtomRecord::atomrecords_to_entry(&records);
    //?? こうしないとボローチェッカー通らない。。。
    
    let mut ret:PDBEntry = PDBEntry::new();
    ret.chains.append(&mut ret_.chains);
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
    for cc in pdb.chains.iter_mut(){
        for rr in cc.residues.iter(){
            if rr.residue_number != 17{
                continue;
            }
            for aa in rr.atoms.iter(){
                println!("{} {} {} {}",rr.residue_name,aa.atom_code,aa.serial_number,aa.alt_code);
            }
        }
        cc.remove_alt(Some(&vec!["B"]));
        for rr in cc.residues.iter(){
            if rr.residue_number != 17{
                continue;
            }
            for aa in rr.atoms.iter(){
                println!("{} {} {} {}",rr.residue_name,aa.atom_code,aa.serial_number,aa.alt_code);
            }
        }
    }
    pdb.save("test/testout.pdb");
}

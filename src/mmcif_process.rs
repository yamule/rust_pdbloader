#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use std::fs::File;

use std::collections::HashMap;
#[allow(unused_imports)]
use std::collections::HashSet;
use regex::Regex;
use super::pdbdata::*;
use super::geometry::Vector3D;


const IS_MULTILINE:i64 = 0;
const IS_SINGLELINE:i64 = 1;

const NO_SEQ_ID:i64 = -999;


const SECTION_NUMLETTER_MAX:usize = 30;//これより長いと複数行モードにされる


const _ATOM_SITE_GROUP_PDB:&str ="_atom_site.group_PDB";
const _ATOM_SITE_ID:&str ="_atom_site.id";
const _ATOM_SITE_TYPE_SYMBOL:&str ="_atom_site.type_symbol";
const _ATOM_SITE_LABEL_ATOM_ID:&str ="_atom_site.label_atom_id";
const _ATOM_SITE_LABEL_ALT_ID:&str ="_atom_site.label_alt_id";
const _ATOM_SITE_LABEL_COMP_ID:&str ="_atom_site.label_comp_id";
const _ATOM_SITE_LABEL_ASYM_ID:&str ="_atom_site.label_asym_id";
const _ATOM_SITE_LABEL_ENTITY_ID:&str ="_atom_site.label_entity_id";
const _ATOM_SITE_LABEL_SEQ_ID:&str ="_atom_site.label_seq_id";
const _ATOM_SITE_PDBX_PDB_INS_CODE:&str ="_atom_site.pdbx_PDB_ins_code";
const _ATOM_SITE_CARTN_X:&str ="_atom_site.Cartn_x";
const _ATOM_SITE_CARTN_Y:&str ="_atom_site.Cartn_y";
const _ATOM_SITE_CARTN_Z:&str ="_atom_site.Cartn_z";
const _ATOM_SITE_OCCUPANCY:&str ="_atom_site.occupancy";
const _ATOM_SITE_B_ISO_OR_EQUIV:&str ="_atom_site.B_iso_or_equiv";
const _ATOM_SITE_PDBX_FORMAL_CHARGE:&str ="_atom_site.pdbx_formal_charge";
const _ATOM_SITE_AUTH_SEQ_ID:&str ="_atom_site.auth_seq_id";
const _ATOM_SITE_AUTH_COMP_ID:&str ="_atom_site.auth_comp_id";
const _ATOM_SITE_AUTH_ASYM_ID:&str ="_atom_site.auth_asym_id";
const _ATOM_SITE_AUTH_ATOM_ID:&str ="_atom_site.auth_atom_id";
const _ATOM_SITE_PDBX_PDB_MODEL_NUM:&str ="_atom_site.pdbx_PDB_model_num";


/*ToDo
テスト構造で
unknown polymer entity '1' near line 27
Unknown polymer entity '2' near line 230
Unknown polymer entity '3' near line 2196
Missing or incomplete entity_poly_seq table. Inferred polymer connectivity.
と言われた。
*/



lazy_static! {
    static ref REGEX_TAILBLANK:Regex = Regex::new(r"[\s]*$").unwrap();
    static ref REGEX_DECIMAL:Regex = Regex::new(r"^([0-9\-]+)\.([0-9]+)$").unwrap();
    static ref REGEX_QUOTATION_NEEDED:Regex = Regex::new("[\\s\"'\\(\\)]").unwrap(); //"
    static ref REGEX_WS:Regex = Regex::new(r"[\s]").unwrap();
}


pub fn write_to_file(filename:&str,contents:Vec<String>){
    let mut f = BufWriter::new(File::create(filename).unwrap());
     for ll in contents{
        f.write_all(ll.as_bytes()).unwrap();
        f.write_all("\n".as_bytes()).unwrap();
    }
}

  
pub fn start_with(target:&str,fragment:&str)-> bool{
    if let Some(x) = target.find(fragment){
        if x == 0{
            return true;
        }
    }
    return false;
}


//一つの生要素を引用符着け、複数行化等 mmcif の値として互換の形にして MMCIF に書き込めるようにする
//複数行にわたる場合は .1 が true
pub fn get_compatible_values(target:&str)->(String,i64){
    let re = regex::Regex::new(r"\r\n").unwrap();
    let re2 = regex::Regex::new(r"[\r\n]").unwrap();
    let mut lines:Vec<String> = vec![];
    for pp in re.split(target){
        for qq in re2.split(pp){
            lines.push(qq.to_string());
        }
    }

    if lines.len() > 1{
        let mut ret:String=";".to_string();
        for ll in lines.into_iter(){
            ret += &ll;
            ret += "\n";
        }
        return (ret+";",IS_MULTILINE);
    }
    let ret = quote(&lines[0]);
    return ret.unwrap_or((lines.remove(0),IS_SINGLELINE));
}

//二番目は一行か複数行かを示す値が入る。
//閾値より長い、もしくは
//ダブルクォート、シングルクォート両方があると複数行に分けられる。
pub fn quote(target:&str)->Option<(String,i64)>{
    if let None = REGEX_QUOTATION_NEEDED.captures(&target){
        return None;
    }
    let q:bool = match target.find("'"){
        Some(_x)=>{true},
        None=>{false},
    };
    
    let dq:bool = match target.find("\""){
        Some(_x)=>{true},
        None=>{false},
    };
    if q && dq || target.len() > SECTION_NUMLETTER_MAX{
        return  Some((";".to_string()+target+"\n;",IS_MULTILINE));
    }else{
        if q{
            return Some(("\"".to_string()+target+"\"",IS_SINGLELINE));
        }else if dq{
            return Some(("'".to_string()+target+"'",IS_SINGLELINE));
        }else{
            return Some(("'".to_string()+target+"'",IS_SINGLELINE));
        }
    }
}


pub fn parse_block(block:&Vec<String>)
->(Vec<String>,Vec<Vec<String>>){
    let mut keys:Vec<String> = vec![];
    let mut values:Vec<Vec<String>> = vec![];

    //From https://github.com/rcsb/py-mmcif/blob/cf7d45a0c178658c856443c14dd2a2cd0b36a87b/mmcif/io/PdbxReader.py#L115
    //rcsb @jdwestbrook
    //APL v2
    
    let mmcif_re = Regex::new(
        ("(?:".to_string()+
        "(?:_(.+?)[.](\\S+))"+//1:親ラベル、2:子ラベル
        "|"+
        "(?:['](.*?)(?:[']\\s|[']$))"+//3:クオート
        "|"+
        "(?:[\"](.*?)(?:[\"]\\s|[\"]$))"+//4:ダブルクオート
        "|"+
        "(?:\\s*#.*$)"+//コメント
        "|"+
        "(\\S+)"+//5:通常文字列
        ")").as_str()
    ).unwrap();

    let get_values = |l:&str|->Vec<String>{
        let mut ret:Vec<String> = vec![];
        let capp = mmcif_re.captures_iter(l);
        for cc in capp.into_iter(){
            if let Some(x) = cc.get(3){
                ret.push(x.as_str().to_string());
            }
            if let Some(x) = cc.get(4){
                ret.push(x.as_str().to_string());
            }
            if let Some(x) = cc.get(5){
                ret.push(x.as_str().to_string());
            }
        }
        return ret;
    };

    let oneline = Regex::new(r"^(_[^\s]+)[\s]+([^\s].*)[\r\n]*$").unwrap();
    let mut multiline_flag = false;
    let mut buff:Vec<String> = vec![];
    let loopmode:bool =  start_with(&block[0],"loop_");
    values.push(vec![]);

    for (lii,line) in block.iter().enumerate(){
        if loopmode && lii == 0{
            continue;
        }
        let mut currentblock:usize  = values.len()-1;
        if !multiline_flag && start_with(line,"_"){
            if let Some(x) = oneline.captures(line){
                let keyy:String = x.get(1).unwrap().as_str().to_string();
                let vall:String = x.get(2).unwrap().as_str().to_string();
                values.get_mut(currentblock).unwrap().push(get_values(vall.as_str())[0].clone());
                keys.push(keyy);
            }else{
                keys.push(line.to_string());
            }
            continue;
        }
        if start_with(line,";"){
            if multiline_flag{
                values[currentblock].push(
                    buff.iter().fold("".to_string(),|s,m|s+m+"\n")
                );
                multiline_flag = false;
                buff.clear();
            }else{
                assert!(buff.len() == 0);
                multiline_flag = true;
                let slen = line.len();
                buff.push((&line[1..slen]).to_string());
            }
            
            if loopmode && values[currentblock].len() == keys.len(){
                values.push(vec![]);
            }
            continue;
        }else{
            if multiline_flag{
                buff.push(line.to_string());
            }else{
                let parsed = get_values(line.as_str());
                for pp in parsed.into_iter(){
                    values.get_mut(currentblock).unwrap().push(pp);
                    
                    if loopmode && values[currentblock].len() == keys.len(){
                        values.push(vec![]);
                        currentblock += 1;
                    }
                }
            }

            //if loopmode && values[currentblock].len() == keys.len(){
            //    values.push(vec![]);
            //}
        }
    }
    for kk in keys.iter_mut(){
        *kk = (*REGEX_TAILBLANK.replace_all(kk, "")).to_string();
    }
    if values[values.len()-1].len() == 0{
        values.pop();
    }

    return (keys,values)
}

pub fn parse_mmcif(filename:&str)->Vec<(Vec<String>,Vec<Vec<String>>)>{
   
    let mut ret:Vec<(Vec<String>,Vec<Vec<String>>)> = vec![];
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut buff:Vec<String> = vec![];
    for (_lcount,line_) in reader.lines().enumerate() {
        let line = line_.unwrap();
        if start_with(&line,"#"){
            if buff.len() > 0{
                ret.push(parse_block(&buff));
                buff.clear();
            }
            continue;
        }else{
            buff.push(line);
        }
    }

    if buff.len() > 0{
        ret.push(parse_block(&buff));
        buff.clear();
    }
    return ret;
}

#[allow(dead_code)]
pub struct MMCIFEntry{
    //二番目要素はベクトルを持っていて、そのベクトルの値に対応するラベルが一番目の要素に入っている
    atom_site:(Vec<String>,Vec<AtomSite>),
    misc_section:Vec<(Vec<String>,Vec<MiscSection>)>,
    entry_id:String,
    header:String,
}


impl MMCIFEntry{
    
    pub fn set_atom_site_section(&mut self,kk:Vec<String>,vll:Vec<AtomSite>){
        /*何をしようとしていたか忘れてしまった。
        let atom_site_map:HashMap<String,usize> = MMCIFEntry::get_key_index_map(&self.atom_site.0);
        let pxx:(i64,i64,i64)=(
            match atom_site_map.get(_ATOM_SITE_CARTN_X){Some(x)=>{*x as i64},None=>{-1}},
            match atom_site_map.get(_ATOM_SITE_CARTN_Y){Some(x)=>{*x as i64},None=>{-1}},
            match atom_site_map.get(_ATOM_SITE_CARTN_Z){Some(x)=>{*x as i64},None=>{-1}}
        );
        */
        self.atom_site =(kk,vll);
        self.update_atom_site_index();
    }

    pub fn update_atom_site_index(&mut self){
        for (aii,aa) in self.atom_site.1.iter_mut().enumerate(){
            aa.set_atom_index(aii as i64);
        }
    }

    //PDB データに見つからなかった AtomSite のインデクスを返す
    pub fn assign_cartn(&mut self,pdbb:&PDBEntry,add_new_atoms:bool)
    ->Vec<usize>{
        let atom_site_map:HashMap<String,usize> = MMCIFEntry::get_key_index_map(&self.atom_site.0);
        let mut newatoms_all:Vec<AtomSite> = vec![];
        let mut newatoms_notassigned:Vec<AtomSite> = vec![];
        let xyz:(usize,usize,usize) = (
         *atom_site_map.get(_ATOM_SITE_CARTN_X).expect("Cartn_x is not found.")
        ,*atom_site_map.get(_ATOM_SITE_CARTN_Y).expect("Cartn_y is not found.")
        ,*atom_site_map.get(_ATOM_SITE_CARTN_Z).expect("Cartn_z is not found.")
        );
        let occ:i64 = match atom_site_map.get(_ATOM_SITE_OCCUPANCY){Some(x) => {*x as i64},None=>{-1}};
        let biso:i64 = match atom_site_map.get(_ATOM_SITE_B_ISO_OR_EQUIV){Some(x) => {*x as i64},None=>{-1}};
        let f_charge:i64 = match atom_site_map.get(_ATOM_SITE_PDBX_FORMAL_CHARGE){Some(x) => {*x as i64},None=>{-1}};
        let mut updated:Vec<bool> = vec![false;self.atom_site.1.len()];

        let mut haspointer_all:i64 = -1;//最後に Model num とかを付ける場合に使用
        for cc in pdbb.chains.iter(){
            for rr in cc.residues.iter(){
                let mut newatoms:Vec<AtomSite> = vec![];
                for aa in rr.iter_atoms(){
                    let p_ = aa.get_external_array_pointer();
                    let mut targetatom_:Option<&mut AtomSite> = None;
                    if p_ > -1{
                        if haspointer_all < 0{
                            haspointer_all = p_;
                        }
                        let p:usize = p_ as usize;
                        updated[p] = true;
                        targetatom_ = Some(&mut self.atom_site.1[p]);   
                    }else if add_new_atoms{
                        let mut naa:AtomSite = AtomSite::new();
                        naa.set_num_values(self.atom_site.0.len());
                        if aa.het{
                            naa.set_value_of(_ATOM_SITE_GROUP_PDB,"HETATM".to_string(),&atom_site_map,true);
                        }else{
                            naa.set_value_of(_ATOM_SITE_GROUP_PDB,"ATOM".to_string(),&atom_site_map,true);
                        }
                        naa.set_value_of(_ATOM_SITE_ID,aa.index.to_string(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_TYPE_SYMBOL,aa.atom_symbol.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_LABEL_ATOM_ID,aa.atom_code.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_LABEL_ALT_ID,aa.alt_code.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_AUTH_ATOM_ID,aa.atom_code.to_string(),&atom_site_map,true);
                        newatoms.push(naa);
                        targetatom_ = Some(newatoms.last_mut().unwrap());
                    }
                    if let None = targetatom_{
                        continue;
                    }
                    let targetatom:&mut AtomSite = targetatom_.unwrap();
                    targetatom.set_xyz(
                        aa.get_x(),
                        aa.get_y(),
                        aa.get_z(),
                        xyz
                    );
                    if occ > -1{
                        targetatom.set_occupancy(aa.occupancy,(occ as usize,));
                    }
                    if biso > -1{
                        targetatom.set_temperature_factor(aa.temp_factor,(biso as usize,));
                    }
                    if f_charge > -1{
                        if let Some(x) = aa.charge.as_ref(){
                            targetatom.set_formal_charge(x.clone(),(f_charge as usize,));
                        }
                    }
                }
                let mut haspointer:i64 = -1;
                for (aii,aa) in rr.iter_atoms().enumerate(){
                    let p_ = aa.get_external_array_pointer();
                    if p_ > -1{
                        haspointer = aii as i64;
                    }
                }
                if haspointer > -1{
                    let refatom:&AtomSite = &self.atom_site.1[haspointer as usize];
                    let ikeys:Vec<&str> = vec![
                    _ATOM_SITE_LABEL_COMP_ID
                    ,_ATOM_SITE_LABEL_ASYM_ID
                    ,_ATOM_SITE_LABEL_ENTITY_ID 
                    ,_ATOM_SITE_LABEL_SEQ_ID
                    ,_ATOM_SITE_PDBX_PDB_INS_CODE
                    ,_ATOM_SITE_AUTH_SEQ_ID
                    ,_ATOM_SITE_AUTH_COMP_ID
                    ,_ATOM_SITE_AUTH_ASYM_ID
                    ,_ATOM_SITE_PDBX_PDB_MODEL_NUM
                    ];
                    for nn in newatoms.iter_mut(){
                        nn.copy_information_from(refatom,&ikeys,&atom_site_map);
                    }
                    newatoms_all.append(&mut newatoms);
                }else{
                    for naa in newatoms.iter_mut(){
                        //同じ Residue 中にマップされた原子が無い場合 residue とかの情報を使って再構成する
                        naa.set_value_of(_ATOM_SITE_LABEL_COMP_ID,rr.residue_name.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_LABEL_ASYM_ID,cc.chain_name.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_LABEL_SEQ_ID,rr.get_residue_number().to_string(),&atom_site_map,true);
                        
                        naa.set_value_of(_ATOM_SITE_AUTH_COMP_ID,rr.residue_name.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_AUTH_ASYM_ID,cc.chain_name.clone(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_AUTH_SEQ_ID,rr.get_residue_number().to_string(),&atom_site_map,true);
                        
                        naa.set_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE,rr.get_ins_code().to_string(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_LABEL_ENTITY_ID,"?".to_string(),&atom_site_map,true);
                        naa.set_value_of(_ATOM_SITE_PDBX_PDB_MODEL_NUM,"?".to_string(),&atom_site_map,true);
                    }
                    newatoms_notassigned.append(&mut newatoms);
                }
            }
        }

        if haspointer_all > -1 && newatoms_notassigned.len() > 0{
            let refatom:&AtomSite = &self.atom_site.1[haspointer_all as usize];
            let ikeys:Vec<&str> = vec![
            _ATOM_SITE_LABEL_ENTITY_ID 
            ,_ATOM_SITE_PDBX_PDB_MODEL_NUM
            ];
            for nn in newatoms_notassigned.iter_mut(){
                nn.copy_information_from(refatom,&ikeys,&atom_site_map);
            }
            newatoms_all.append(&mut newatoms_notassigned);
        }

        self.atom_site.1.append(&mut newatoms_all);
        let mut ret:Vec<usize> = vec![];
        for mm in 0..updated.len(){
            if !updated[mm]{
                ret.push(mm);
            }
        }
        return ret;
    }    
    
    //ToDo: use auth 途中
    pub fn load_mmcif(filename:&str,use_auth:bool)->(MMCIFEntry,PDBEntry){
        let mut blocks:Vec<(Vec<String>,Vec<Vec<String>>)> = parse_mmcif(filename);
        let mut atom_site_block_:Option<(Vec<String>,Vec<Vec<String>>)> = None;
        
        //let mut atom_sites:Vec<AtomSite> = vec![];
        let header = blocks.remove(0);
        //println!("{}",MMCIFEntry::blocks_to_string(&blocks));
        let mut misc_section:Vec<(Vec<String>,Vec<MiscSection>)> = vec![];
        let mut entry_id:String = "NONE".to_string();
        for bb in blocks.into_iter(){
            if bb.0.len() == 0{//entry
                continue;
            }
            if start_with(&bb.0[0],"_atom_site."){
                atom_site_block_ = Some(bb);
            }else if start_with(&bb.0[0],"_entry."){
                let mut eindex:i64 = -1;
                for eii in 0..bb.0.len(){
                    if bb.0[eii] == "_entry.id"{
                        eindex = eii as i64;
                    }
                }
                if eindex > -1{
                    entry_id = bb.1[0][eindex as usize].clone();   
                }
                misc_section.push(
                    (bb.0,bb.1.into_iter().map(|m|MiscSection{values:m}).collect())
                );
            }else{
                misc_section.push(
                    (bb.0,bb.1.into_iter().map(|m|MiscSection{values:m}).collect())
                );
            }
            
        }
        let mut ret = MMCIFEntry{
            atom_site:(vec![],vec![]),
            entry_id:entry_id,
            header:header.1[0][0].clone(),
            misc_section:misc_section
        };
        if let Some(x) = atom_site_block_{
            let (keys,values) = x;
            let mut atoms:Vec<AtomSite> = vec![];
            for vv in values.into_iter(){
                if vv.len() == 0{
                    continue;
                }
                let mut atom:AtomSite = AtomSite::new();
                atom.values = vv;
                atoms.push(atom);
            }
            ret.set_atom_site_section(keys,atoms);
        }else{
            panic!("can not find _atom_site. block!");
        }


        let ret2:PDBEntry = MMCIFEntry::atomsite_to_entry(&ret,false);





        return (ret,ret2);
    }

    pub fn get_key_index_map(vs:&Vec<String>)->HashMap<String,usize>{
        let mut ret:HashMap<String,usize> = HashMap::new();
        for (vii,vss) in vs.iter().enumerate(){
            ret.insert(vss.to_string(),vii);
        }
        return ret;
    }
    pub fn atomsite_to_entry(mmcif:&MMCIFEntry,use_auth:bool)->PDBEntry{
        
        let mut ret:PDBEntry = PDBEntry::new();
        let mut chains:HashMap<String,Vec<Vec<&AtomSite>>> = HashMap::new();
        let mut chains_rindex:HashMap<String,HashMap<String,usize>> = HashMap::new();
        let atom_site_map:HashMap<String,usize> = MMCIFEntry::get_key_index_map(&mmcif.atom_site.0);
    
    
        let asym_id_label:String;
        let _atom_id_label:String;//CA とか CB とか atom_name
        let comp_id_label:String;
        let seq_id_label:String;
        if use_auth{
            asym_id_label = _ATOM_SITE_AUTH_ASYM_ID.to_string();
            comp_id_label = _ATOM_SITE_AUTH_COMP_ID.to_string();
            _atom_id_label = _ATOM_SITE_AUTH_ATOM_ID.to_string();
            seq_id_label = _ATOM_SITE_AUTH_SEQ_ID.to_string();
        }else{
            asym_id_label = _ATOM_SITE_LABEL_ASYM_ID.to_string();
            comp_id_label = _ATOM_SITE_LABEL_COMP_ID.to_string();
            _atom_id_label = _ATOM_SITE_LABEL_ATOM_ID.to_string();
            seq_id_label = _ATOM_SITE_LABEL_SEQ_ID.to_string();
        }
        
        let asym_id_:usize = *atom_site_map.get(&asym_id_label).unwrap_or_else(||panic!("{} is not defined.",asym_id_label));
        let comp_id_:usize = *atom_site_map.get(&comp_id_label).unwrap_or_else(||panic!("{} is not defined.",comp_id_label));
        let seq_id_:usize = *atom_site_map.get(&seq_id_label).unwrap_or_else(||panic!("{} is not defined.",seq_id_label));
        //let atom_id_:usize = *atom_site_map.get(&atom_id_label).unwrap_or_else(||panic!("{} is not defined.",atom_id_label));
    
        for (_aii,aa) in mmcif.atom_site.1.iter().enumerate(){
    
            let asym_id:String = aa.values[asym_id_].clone();
    
            if !chains.contains_key(&asym_id){
                chains.insert(asym_id.clone(),vec![]);
                chains_rindex.insert(asym_id.clone(),HashMap::new());
            }
            
            let r_label:String = aa.get_unique_residue_label(&atom_site_map);
            if !chains_rindex.get(&asym_id).unwrap().contains_key(r_label.as_str()){
                chains_rindex.get_mut(&asym_id).unwrap().insert(
                    r_label.clone(),chains.get(&asym_id).unwrap().len());
                chains.get_mut(&asym_id).unwrap().push(vec![]);
            }
            //一つの残基につき一意になるような Reside label をハッシュキーにしてが同じものは一つの Vector に入れる。
            let rpos:&usize = chains_rindex.get_mut(asym_id.as_str()).unwrap().get(r_label.as_str()).unwrap();
            
            chains.get_mut(asym_id.as_str()).unwrap()[*rpos].push(aa);
        }
        for (_cc,rr) in chains.into_iter(){
            
            //Chain の実体を作成する
            let mut cc:PDBAsym = PDBAsym::new(&_cc);
            for vv in rr.into_iter(){
                //Residue の実体を作成する
                let atoms_v:Vec<&AtomSite> = vv;
                let aa = &atoms_v[0];
                let comp_id:String = aa.values[comp_id_].clone();
                let seq_id:String = aa.values[seq_id_].clone();
                let ins_code:String = aa.get_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE,&atom_site_map).to_string();
                let mut rr:PDBComp = PDBComp::new();
                if seq_id == "."{
                    //HOH に seq_id が設定されていない。。。
                    rr.set_residue_number(NO_SEQ_ID);
                }else{
                    rr.set_residue_number(seq_id.parse::<i64>().unwrap_or_else(|_|panic!("{} can not parse.",seq_id)));
                    if rr.get_residue_number() == NO_SEQ_ID{
                        panic!("{} is not allowed!",NO_SEQ_ID);
                    }
                }
                rr.residue_name = comp_id.clone();
                rr.ins_code = ins_code;
                for aa in atoms_v.iter(){
                    let att:PDBAtom = aa.atom_site_to_pdbatom(&atom_site_map
                        ,use_auth);
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
    

    //それぞれのセクションが持っているのは
    //(Vec<key>,Vec<Vec<value>>) というタプルであるべきで、それを MMCIF フォーマットに整形する。
    pub fn blocks_to_strings(blocks:&Vec<(&Vec<String>,&Vec<&Vec<String>>)>)->Vec<String>{
        let mut ret:Vec<String> = vec![];
        for vv in blocks.iter(){
            ret.push("# ".to_string());
            let lab = &vv.0;
            let val = &vv.1;
            if val.len() == 1{
                let lines:&Vec<String> = &val[0];
                assert_eq!(lab.len(),lines.len());
                for ii in 0..lines.len(){
                    let cval = get_compatible_values(&lines[ii]);
                    if cval.1 == IS_MULTILINE{
                        ret.push(format!("{}     \n{}",lab[ii],&cval.0));
                    }else{
                        ret.push(format!("{}     {} ",lab[ii],&cval.0));
                    }
                }
            }else{
                ret.push("loop_".to_string());
                for ii in 0..lab.len(){
                    ret.push(lab[ii].clone()+" ");
                }

                //カラムサイズが同じでないと駄目らしい。少なくとも UCSF Chimera では
                let mut maxnumletter:Vec<usize> = vec![0;lab.len()];
                for vv in val.iter(){
                    if vv.len() == 0{
                        continue;
                    }
                    assert_eq!(lab.len(),vv.len());
                    for ii in 0..lab.len(){
                        let cval = get_compatible_values(&vv[ii]);
                        if cval.1 == IS_SINGLELINE{
                            maxnumletter[ii] = maxnumletter[ii].max(cval.0.len());
                        }
                    }
                }
                let mut ws_maxnum:Vec<String> = vec![];
                for sii in maxnumletter.iter(){
                    //文字の位置を他の行と合わせる
                    ws_maxnum.push((0..*sii).into_iter().fold(" ".to_string(),|s,_|s+" ")+" ");
                }

                for vv in val.iter(){
                    if vv.len() == 0{
                        continue;
                    }
                    let mut lline:String = "".to_string();
                    let mut prev_cr:i64 = IS_SINGLELINE;
                    for ii in 0..lab.len(){
                        let cval = get_compatible_values(&vv[ii]);
                        if cval.1 == IS_MULTILINE{
                            if ii != 0 && prev_cr == IS_SINGLELINE{
                                lline += "\n";
                            }
                            lline += cval.0.as_str();
                            if ii != lab.len()-1{
                                lline += "\n";
                            }
                        }else{
                            lline += &((format!("{}",cval.0.as_str())+ws_maxnum[ii].as_str())[0..=maxnumletter[ii]]);
                        }
                        prev_cr = cval.1;
                    }

                    ret.push(lline);
                }
            }
        }
        return ret;
    }

    //めちゃ遅い・・・
    pub fn save(&self,filename:&str){
        let mut lines:Vec<String> = vec![];
        lines.push(self.header.clone());
        
        for ll in self.misc_section.iter(){
            let val:Vec<&Vec<String>> = ll.1.iter().map(|m|{&m.values}).collect();
            lines.append(&mut MMCIFEntry::blocks_to_strings(
                &vec![(&ll.0,&val)]
            ));    
        }
        

        let asite:Vec<&Vec<String>> = self.atom_site.1.iter().map(|m|&m.values).collect();
        lines.append(&mut MMCIFEntry::blocks_to_strings(
            &vec![(&self.atom_site.0,&asite)]
        ));
        lines.push("# ".to_string());
        write_to_file(filename,lines);
        
    }
}

pub struct MiscSection{
    pub values:Vec<String>
}

#[allow(non_snake_case)]
pub struct AtomSite{
    pub values:Vec<String>,
    pub atom_index:i64, //MMCIFEntry 内の Vec にあるこのインスタンスのインデクス。PDBAtom から参照する
    pub used_auth:bool
}

impl AtomSite{
    

    //残基名とかチェーン名とか、上位情報をコピーする
    pub fn copy_information_from(&mut self,src:&AtomSite,keys:&Vec<&str>,vmap:&HashMap<String,usize>){
        assert_eq!(self.values.len(),src.values.len());
        for ii in 0..keys.len(){
            if vmap.contains_key(keys[ii]){
                let uii:usize = *vmap.get(keys[ii]).unwrap();
                self.values[uii]
                 = src.values[uii].clone();
            }
        }
    }

    pub fn set_num_values(&mut self,siz:usize){
        assert!(self.values.len() == 0);
        self.values = vec!["?".to_string();siz];
    }

    pub fn new()->AtomSite{
        return AtomSite{
            values:vec![],
            atom_index:-1,
            used_auth:false
        };
    }
    

    //MMCIFEntry 内の Vec にあるこのインスタンスのインデクスを与えてください
    pub fn set_atom_index(&mut self,i:i64){
        self.atom_index = i;
    }
    
    //PDBAtom 以外から変更されることを今のところ想定していない
    //小数点以下の桁数を合わせるためであり、あまり意味はない。
    //set_value とかでもよいと思う。
    fn set_xyz(&mut self,x:f64,y:f64,z:f64,indices:(usize,usize,usize)){
        self.values[indices.0] = format!("{:.3}",x);
        self.values[indices.1] = format!("{:.3}",y);
        self.values[indices.2] = format!("{:.3}",z);
    }
    fn set_occupancy(&mut self,v:f64,indices:(usize,)){
        self.values[indices.0] = format!("{:.2}",v);
    }
    fn set_temperature_factor(&mut self,v:f64,indices:(usize,)){
        self.values[indices.0] = format!("{:.2}",v);
    }
    fn set_formal_charge(&mut self,v:String,indices:(usize,)){
        self.values[indices.0] = v;
    }
    
    
    pub fn set_value_of(&mut self,k:&str,v:String,vmap:&HashMap<String,usize>,dont_panic:bool){
        if !vmap.contains_key(k){
            if !dont_panic{
                panic!("{} is not found in key list.",k);
            }
            return;
        }
        self.values[*vmap.get(k).unwrap()] = v;
    }
    
    
    pub fn atom_site_to_pdbatom(&self,vmap:&HashMap<String,usize>,use_auth:bool)->PDBAtom{
        let mut ret = PDBAtom::new();
        if self.get_atom_index() < 0{
            panic!("update_atom_site_index must have been performed at first.");
        }

        ret.set_external_array_pointer(self.get_atom_index());
        ret.set_xyz(self.get_value_of(_ATOM_SITE_CARTN_X,&vmap).parse::<f64>().unwrap_or_else(|_|panic!("can not parse x"))
        ,self.get_value_of(_ATOM_SITE_CARTN_Y,&vmap).parse::<f64>().unwrap_or_else(|_|panic!("can not parse y"))
        ,self.get_value_of(_ATOM_SITE_CARTN_Z,&vmap).parse::<f64>().unwrap_or_else(|_|panic!("can not parse z"))
        );
       
        ret.serial_number = self.get_value_of(_ATOM_SITE_ID,vmap).parse::<i64>().expect("Cannot parse atom id.");
        ret.atom_symbol = self.get_value_of(_ATOM_SITE_TYPE_SYMBOL,vmap).to_string();
        ret.alt_code = self.get_value_of(_ATOM_SITE_LABEL_ALT_ID,vmap).to_string();
        if ret.alt_code == "?" || ret.alt_code == "."{
            ret.alt_code = "".to_string();
        }

        ret.dummy = false;
        if self.get_value_of(_ATOM_SITE_GROUP_PDB,vmap) == "HETATM"{
            ret.het = true;
        }else{
            ret.het = false;
        }
        
        if vmap.contains_key(_ATOM_SITE_PDBX_FORMAL_CHARGE){
            let fcc:String = self.get_value_of(_ATOM_SITE_PDBX_FORMAL_CHARGE,vmap).to_string();
            if fcc == "?" || fcc == "."{
                ret.charge = None;
            }else{
                ret.charge = Some(fcc);
            }
        }
        if vmap.contains_key(_ATOM_SITE_OCCUPANCY){
            ret.occupancy = self.get_value_of(_ATOM_SITE_OCCUPANCY,vmap).parse::<f64>().unwrap();
        }
        if vmap.contains_key(_ATOM_SITE_B_ISO_OR_EQUIV){
            ret.temp_factor = self.get_value_of(_ATOM_SITE_B_ISO_OR_EQUIV,vmap).parse::<f64>().unwrap();
        }

        if use_auth{
            ret.atom_code = self.get_value_of(_ATOM_SITE_AUTH_ATOM_ID,vmap).to_string();
            
        }else{
            ret.atom_code = self.get_value_of(_ATOM_SITE_LABEL_ATOM_ID,vmap).to_string();
        }
        
        return ret;
    }
    pub fn get_atom_index(&self)->i64{
        return self.atom_index;
    }


    pub fn get_value_of(&self,key:&str,vmap:&HashMap<String,usize>)->&str{
        return match vmap.get(key){
            Some(x)=>self.values[*x].as_str(),
            None=>"?"};
    }
    pub fn get_unique_residue_label(&self,vmap:&HashMap<String,usize>)->String{
        return "".to_string()
        +self.get_value_of(_ATOM_SITE_LABEL_COMP_ID,vmap)
        +"#"
        +self.get_value_of(_ATOM_SITE_LABEL_SEQ_ID,vmap)
        +"#"
        +self.get_value_of(_ATOM_SITE_AUTH_COMP_ID,vmap)
        +"#"
        +self.get_value_of(_ATOM_SITE_AUTH_SEQ_ID,vmap)
        +"#"
        +self.get_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE,vmap)
        ;
    }
}


#[test]
fn regextest(){
    let mmcif_re = Regex::new(
        ("(?:".to_string()+
        "(?:_(.+?)[.](\\S+))"+//1:親ラベル、2:子ラベル
        "|"+
        "(?:['](.*?)(?:[']\\s|[']$))"+//3:クオート
        "|"+
        "(?:[\"](.*?)(?:[\"]\\s|[\"]$))"+//4:ダブルクオート
        "|"+
        "(?:\\s*#.*$)"+//コメント
        "|"+
        "(\\S+)"+//5:通常文字列
        ")").as_str()
    ).unwrap();
    let line = "HETATM 2274 C  \"C1'\"  . 'QWE' F 5 .   ? 16.593 -13.014 18.550 1.00 33.08 ? 373 QWE H \"C1'\"  1 ";
    let capp = mmcif_re.captures_iter(line);
    for cc in capp.into_iter(){
        println!("{:?}",cc);
    }
}

#[test]
fn mmcifloadtest(){
    let pdbentry = MMCIFEntry::load_mmcif("example_files/ins_example_1a4w.cif",false);
    for mm in pdbentry.0.misc_section.iter(){//atom_site を間違うと他でエラーが出ると思う
        for ii in 0..(mm.1).len(){
            assert_eq!(mm.0.len(),mm.1[ii].values.len());
        }
    }
    pdbentry.1.save("test/mmcifout_1a4w.pdb");
    pdbentry.0.save("test/mmcifout_1a4w.cif");
}

#[test]
fn quote_check(){
    //こういう仕様だと思うが間違っているかもしれない。
    assert_eq!(get_compatible_values("test"),("test".to_string(),IS_SINGLELINE));
    assert_eq!(get_compatible_values("(test)"),("'(test)'".to_string(),IS_SINGLELINE));
    assert_eq!(get_compatible_values("t\"est"),("'t\"est'".to_string(),IS_SINGLELINE));
    //ダブルクォートとシングルクォートが混ざっている場合が特に不明
    assert_eq!(get_compatible_values("t\"'est"),(";t\"'est\n;".to_string(),IS_MULTILINE));
    assert_eq!(get_compatible_values("t'est"),("\"t'est\"".to_string(),IS_SINGLELINE));
    assert_eq!(get_compatible_values("te\nst"),(";te\nst\n;".to_string(),IS_MULTILINE));
    assert_eq!(get_compatible_values("te\rst"),(";te\nst\n;".to_string(),IS_MULTILINE));
    assert_eq!(get_compatible_values("te\r\nst"),(";te\nst\n;".to_string(),IS_MULTILINE));
    assert_eq!(get_compatible_values("t'e\r\ns't"),(";t'e\ns't\n;".to_string(),IS_MULTILINE));
}

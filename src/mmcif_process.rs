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


const __ELEMENT_INDEX:&str = "__element_index";
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
    //三番目要素はベクトルを持っていて、そのベクトルの値に対応するラベルが一番目の要素に入っている
    //二番目の要素は逆引き
    pub atom_site:(Vec<String>,HashMap<String,usize>,Vec<Vec<Box<String>>>),
    pub misc_section:Vec<(Vec<String>,HashMap<String,usize>,Vec<Vec<Box<String>>>)>,
    pub entry_id:String,
    pub header:String,
}


impl MMCIFEntry{
    pub fn new()->MMCIFEntry{
        return MMCIFEntry{
            atom_site:(vec![],HashMap::new(),vec![]),
            misc_section:vec![(vec![],HashMap::new(),vec![])],
            entry_id:"".to_owned(),
            header:"".to_owned()
        };
    }

    pub fn get_num_atoms(&self)->usize{
        return self.atom_site.2.len();
    }

    //与えられたインデックスのリストを Entity ごとに分ける
    pub fn get_entity_map(&self,atom_records:&Vec<usize>)->HashMap<String,Vec<usize>>{
        return self.get_map(atom_records,StructureHierarchyLevel::Entity);
    }
    
    //与えられたインデックスのリストを Model ごとに分ける
    pub fn get_model_map(&self,atom_records:&Vec<usize>)->HashMap<String,Vec<usize>>{
        return self.get_map(atom_records,StructureHierarchyLevel::Model);
    }
    
    //与えられたインデックスのリストを asym ごとに分ける
    pub fn get_asym_map(&self,atom_records:&Vec<usize>)->HashMap<String,Vec<usize>>{
        return self.get_map(atom_records,StructureHierarchyLevel::Asym);
    }
    
    
    //与えられたインデックスのリストを comp ごとに分ける
    pub fn get_comp_map(&self,atom_records:&Vec<usize>)->HashMap<String,Vec<usize>>{
        return self.get_map(atom_records,StructureHierarchyLevel::Comp);
    }

    //インデックスのリストを分割する関数の本体
    pub fn get_map(&self,atom_records:&Vec<usize>,lev:StructureHierarchyLevel)->HashMap<String,Vec<usize>>{
        let mut ret:HashMap<String,Vec<usize>> = HashMap::new();
        for aa in atom_records.iter(){
            let ak:String = match lev{
                StructureHierarchyLevel::Model => self.get_atom_site(*aa).get_pdbx_PDB_model_num().to_string(),
                StructureHierarchyLevel::Entity => self.get_atom_site(*aa).get_label_entity_id().to_string(),
                StructureHierarchyLevel::Asym => self.get_atom_site(*aa).get_label_asym_id().to_string(),
                StructureHierarchyLevel::Comp => self.get_atom_site(*aa).get_unique_residue_label(),
                _=> panic!("{:?} is not allowed in this function.",lev)
            };
            if !ret.contains_key(&ak){
                ret.insert(ak.clone(),vec![]);
            }
            ret.get_mut(&ak).unwrap().push(*aa);
        }
        return ret;
    }
    
    //atom_site 情報を設定する
    pub fn set_atom_site_section(&mut self,kk:Vec<String>,vll:Vec<Vec<Box<String>>>){
        self.atom_site =(kk.clone(),MMCIFEntry::get_key_index_map(&kk),vll);
        self.update_atom_site_index();
    }

    pub fn update_atom_site_index(&mut self){
        for (aii,aa) in self.atom_site.2.iter_mut().enumerate(){
            AtomSiteMut::new(&mut self.atom_site.1,aa).set_index(aii as i64);
        }
    }

    pub fn get_atom_site(&self,i:usize)->AtomSite{
        return AtomSite::new(&self.atom_site.1,&self.atom_site.2[i]);
    }
    
    pub fn get_mut_atom_site(&mut self,i:usize)->AtomSiteMut{
        return AtomSiteMut::new(&mut self.atom_site.1,&mut self.atom_site.2[i]);
    }

    //keymap のインデクスにそれぞれの要素が入った String のベクトルを返す。
    pub fn parse_atom_line_pdb(line:&str,model_code:&str,keymap:&HashMap<String,usize>)-> Vec<String>{
        let mut ret = vec!["".to_owned();keymap.len()];

        let pt:Vec<char> = line.chars().collect();
        let label = slice_to_string(&pt,0,6);
        ret[*keymap.get(_ATOM_SITE_GROUP_PDB).unwrap()] = label;
        
        let idd = slice_to_string(&pt,6,11);
        ret[*keymap.get(_ATOM_SITE_ID).unwrap()] = idd;

        let atomid = slice_to_string(&pt,12,16);
        ret[*keymap.get(_ATOM_SITE_LABEL_ATOM_ID).unwrap()] = atomid.clone();
        ret[*keymap.get(_ATOM_SITE_AUTH_ATOM_ID).unwrap()] = atomid;

        let altloc = slice_to_string(&pt,16,17);
        ret[*keymap.get(_ATOM_SITE_LABEL_ALT_ID).unwrap()] = altloc;

        let resname = slice_to_string(&pt,17,20);
        ret[*keymap.get(_ATOM_SITE_LABEL_COMP_ID).unwrap()] = resname.clone();
        ret[*keymap.get(_ATOM_SITE_AUTH_COMP_ID).unwrap()] = resname;

        let chainid = slice_to_string(&pt,21,22);
        ret[*keymap.get(_ATOM_SITE_LABEL_ASYM_ID).unwrap()] = chainid.clone();
        ret[*keymap.get(_ATOM_SITE_AUTH_ASYM_ID).unwrap()] = chainid;

        let respos = slice_to_string(&pt,22,26);
        ret[*keymap.get(_ATOM_SITE_LABEL_SEQ_ID).unwrap()] = respos.clone();
        ret[*keymap.get(_ATOM_SITE_AUTH_SEQ_ID).unwrap()] = respos;

        let icode = slice_to_string(&pt,26,27);
        ret[*keymap.get(_ATOM_SITE_PDBX_PDB_INS_CODE).unwrap()] = icode;

        let x = slice_to_string(&pt,30,38);
        ret[*keymap.get(_ATOM_SITE_CARTN_X).unwrap()] = x;

        let y = slice_to_string(&pt,38,46);
        ret[*keymap.get(_ATOM_SITE_CARTN_Y).unwrap()] = y;

        let z = slice_to_string(&pt,46,54);
        ret[*keymap.get(_ATOM_SITE_CARTN_Z).unwrap()] = z;

        let occupancy = slice_to_string(&pt,54,60);
        if occupancy.len() > 0{
            ret[*keymap.get(_ATOM_SITE_OCCUPANCY).unwrap()] = occupancy;
        }

        let tempfactor = slice_to_string(&pt,60,66);
        if tempfactor.len() > 0{
            ret[*keymap.get(_ATOM_SITE_B_ISO_OR_EQUIV).unwrap()] = tempfactor;
        }

        let element = slice_to_string(&pt,76,78);
        if element.len() > 0{
            ret[*keymap.get(_ATOM_SITE_TYPE_SYMBOL).unwrap()] = element;
        }

        let charge = slice_to_string(&pt,78,80);
        if charge.len() > 0{
            ret[*keymap.get(_ATOM_SITE_PDBX_FORMAL_CHARGE).unwrap()] = charge;
        }
        ret[*keymap.get(_ATOM_SITE_PDBX_PDB_MODEL_NUM).unwrap()] = model_code.to_string();
        
        return ret;
    }

    //ToDo: use auth 途中
    pub fn load_mmcif(filename:&str,use_auth:bool)->(MMCIFEntry,PDBEntry){
        let mut blocks:Vec<(Vec<String>,Vec<Vec<String>>)> = parse_mmcif(filename);
        let mut atom_site_block_:Option<(Vec<String>,Vec<Vec<Box<String>>>)> = None;
        
        //let mut atom_sites:Vec<AtomSiteMut> = vec![];
        let header = blocks.remove(0);
        //println!("{}",MMCIFEntry::blocks_to_string(&blocks));
        let mut misc_section:Vec<(Vec<String>,HashMap<String,usize>,Vec<Vec<Box<String>>>)> = vec![];
        let mut entry_id:String = "NONE".to_string();
        let vecvec_string_in_box = |b:Vec<Vec<String>>|->Vec<Vec<Box<String>>> {
            b.into_iter().map(|m|m.into_iter().map(|mm|Box::new(mm)).collect()).collect()
        };


        for bb in blocks.into_iter(){
            if bb.0.len() == 0{//entry
                continue;
            }
            if start_with(&bb.0[0],"_atom_site."){
                atom_site_block_ = Some((bb.0,vecvec_string_in_box(bb.1)));
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
                    (bb.0.clone(),MMCIFEntry::get_key_index_map(&bb.0),vecvec_string_in_box(bb.1))
                );
            }else{
                misc_section.push(
                    (bb.0.clone(),MMCIFEntry::get_key_index_map(&bb.0),vecvec_string_in_box(bb.1))
                );
            }
            
        }
        let mut ret = MMCIFEntry{
            atom_site:(vec![],HashMap::new(),vec![]),
            entry_id:entry_id,
            header:header.1[0][0].clone(),
            misc_section:misc_section
        };
        if let Some(x) = atom_site_block_{
            let (keys,values) = x;
            let mut atoms:Vec<Vec<Box<String>>> = vec![];
            for vv in values.into_iter(){
                if vv.len() == 0{
                    continue;
                }
                atoms.push(vv);
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
        let mut chains:HashMap<String,Vec<Vec<&AtomSiteMut>>> = HashMap::new();
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
        let atom_sites:Vec<AtomSiteMut> = vec![];
        for (_aii,aa) in mmcif.atom_site.2.iter().enumerate(){
    //かきさし
            //let asym_id:String = aa.values[asym_id_].clone();
    
        }
    /*
        for (cii,cc) in ret.get_model_at(0).get_entity_at(0).iter_asyms().enumerate(){
            for (rii,rr) in cc.iter_comps().enumerate(){
                assert_eq!(cii as i64,rr.parent_chain.unwrap());
                for (_aii,aa) in rr.iter_atoms().enumerate(){
                    assert_eq!(cii as i64,aa.parent_chain.unwrap());
                    assert_eq!(rii as i64,aa.parent_residue.unwrap());
                    //println!("{} {} {} {} {} {} ",cc.chain_name,rr.get_comp_id(),aa.atom_code,aa.get_x(),aa.get_y(),aa.get_z());
                }
            }
        }
    */
        return ret;
    }
    

    //それぞれのセクションが持っているのは
    //(Vec<key>,Vec<Vec<value>>) というタプルであるべきで、それを MMCIF フォーマットに整形する。
    pub fn blocks_to_strings(blocks:&Vec<(&Vec<String>,&Vec<Vec<Box<String>>>)>)->Vec<String>{
        let mut ret:Vec<String> = vec![];
        for vv in blocks.iter(){
            ret.push("# ".to_string());
            let lab = &vv.0;
            let val = &vv.1;
            if val.len() == 1{
                let lines:&Vec<Box<String>> = &val[0];
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
            lines.append(&mut MMCIFEntry::blocks_to_strings(
                &vec![(&ll.0,&ll.2)]
            ));    
        }
        

        lines.append(&mut MMCIFEntry::blocks_to_strings(
            &vec![(&self.atom_site.0,&self.atom_site.2)]
        ));
        lines.push("# ".to_string());
        write_to_file(filename,lines);
        
    }
}

pub struct MiscSection{
    pub values:Vec<String>
}

#[allow(non_snake_case)]
pub struct AtomSiteMut<'a>{
    pub keymap:&'a mut HashMap<String,usize>,
    pub values:&'a mut Vec<Box<String>>,
}

impl<'a> AtomSiteMut<'a>{
    
    pub fn set_index(&mut self,s:i64){
        self.set_value_of(__ELEMENT_INDEX,s.to_string(),false);
    }
    

    //残基名とかチェーン名とか、上位情報をコピーする
    pub fn copy_information_from(&mut self,src:&AtomSiteMut,keys:&Vec<&str>,vmap:&HashMap<String,usize>){
        assert_eq!(self.values.len(),src.values.len());
        for ii in 0..keys.len(){
            if vmap.contains_key(keys[ii]){
                let uii:usize = *vmap.get(keys[ii]).unwrap();
                self.set_value_of(keys[ii],(*src.values[uii]).clone(),true);
            }
        }
    }

    pub fn new(km:&'a mut HashMap<String,usize>,val:&'a mut Vec<Box<String>>)->AtomSiteMut<'a>{
        return AtomSiteMut{
            keymap:km,values:val
        };
    }
//ToDo
//Vec しかないので Collect で HM も作る
//get_model_at(0).get_entity_at(0) で Chain を取ろうとしているところは、
//Get All Chains というような関数を作って Entity から何からまたいだ結果を全部返させる

    pub fn set_group_PDB(&mut self,s:String){self.set_value_of(_ATOM_SITE_GROUP_PDB,s,false);} 
    pub fn set_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_ID,s,false);} 
    pub fn set_type_symbol(&mut self,s:String){self.set_value_of(_ATOM_SITE_TYPE_SYMBOL,s,false);} 
    pub fn set_label_atom_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_ATOM_ID,s,false);} 
    pub fn set_label_alt_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_ALT_ID,s,false);} 
    pub fn set_label_comp_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_COMP_ID,s,false);} 
    pub fn set_label_asym_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_ASYM_ID,s,false);} 
    pub fn set_label_entity_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_ENTITY_ID,s,false);} 
    pub fn set_label_seq_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_LABEL_SEQ_ID,s,false);} 
    pub fn set_pdbx_PDB_ins_code(&mut self,s:String){self.set_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE,s,false);} 
    pub fn set_Cartn_x(&mut self,s:String){self.set_value_of(_ATOM_SITE_CARTN_X,s,false);} 
    pub fn set_Cartn_y(&mut self,s:String){self.set_value_of(_ATOM_SITE_CARTN_Y,s,false);} 
    pub fn set_Cartn_z(&mut self,s:String){self.set_value_of(_ATOM_SITE_CARTN_Z,s,false);} 
    pub fn set_occupancy(&mut self,s:String){self.set_value_of(_ATOM_SITE_OCCUPANCY,s,false);} 
    pub fn set_B_iso_or_equiv(&mut self,s:String){self.set_value_of(_ATOM_SITE_B_ISO_OR_EQUIV,s,false);} 
    pub fn set_pdbx_formal_charge(&mut self,s:String){self.set_value_of(_ATOM_SITE_PDBX_FORMAL_CHARGE,s,false);} 
    pub fn set_auth_seq_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_AUTH_SEQ_ID,s,false);} 
    pub fn set_auth_comp_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_AUTH_COMP_ID,s,false);} 
    pub fn set_auth_asym_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_AUTH_ASYM_ID,s,false);} 
    pub fn set_auth_atom_id(&mut self,s:String){self.set_value_of(_ATOM_SITE_AUTH_ATOM_ID,s,false);} 
    pub fn set_pdbx_PDB_model_num(&mut self,s:String){self.set_value_of(_ATOM_SITE_PDBX_PDB_MODEL_NUM,s,false);} 


    //PDBAtom 以外から変更されることを今のところ想定していない
    //小数点以下の桁数を合わせるためであり、あまり意味はない。
    //set_value とかでもよいと思う。
    fn set_xyz(&mut self,x:f64,y:f64,z:f64){
        self.set_Cartn_x(format!("{:.3}",x));
        self.set_Cartn_y(format!("{:.3}",y));
        self.set_Cartn_z(format!("{:.3}",z));
    }
    pub fn set_value_of(&mut self,k:&str,v:String,dont_panic:bool){
        if !self.keymap.contains_key(k){
            if !dont_panic{
                panic!("{} is not found in key list.",k);
            }
            return;
        }
        *self.values[*self.keymap.get(k).unwrap()] = v;
    }
    /*
    
    pub fn atomsite_to_atom(&self)->PDBAtom{
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
            atom_site_key:self.index_in_entry
        };//alt_loc は他の Atom も見ないと処理できないと思う
            
        return ret;
    }

    */
    
    
}
#[allow(non_snake_case)]
pub struct AtomSite<'a>{
    pub keymap:&'a HashMap<String,usize>,
    pub values:&'a Vec<Box<String>>,
}

impl<'a> AtomSite<'a>{
    
    
    pub fn get_index(&self)->i64{
        return self.get_value_of(__ELEMENT_INDEX).parse::<i64>().unwrap_or(-1);
    }
    
    pub fn atomsite_to_atom(&self)->PDBAtom{
        let ret:PDBAtom = PDBAtom{
            parent_entry:None,
            parent_entity:None,
            parent_asym:None,
            parent_comp:None,
            index:-1,
            serial_number:self.get_label_atom_id().parse::<i64>().unwrap().clone(),
            x:self.get_Cartn_x().parse::<f64>().unwrap().clone(),
            y:self.get_Cartn_y().parse::<f64>().unwrap().clone(),
            z:self.get_Cartn_z().parse::<f64>().unwrap().clone(),
            charge:if self.get_pdbx_formal_charge().len() > 0{Some(self.get_pdbx_formal_charge().to_string())}else{None},
            occupancy:if self.get_occupancy().len() > 0{self.get_occupancy().parse::<f64>().unwrap()}else{1.0},
            temp_factor:if self.get_B_iso_or_equiv().len() > 0{self.get_B_iso_or_equiv().parse::<f64>().unwrap()}else{0.0},
            atom_symbol:self.get_type_symbol().to_string(),
            atom_code:self.get_label_atom_id().to_string(),
            alt_code:self.get_label_alt_id().to_string(),
            dummy:true,
            het:self.is_het(),
            alt:self.is_alt(),
            is_ligand:false,
            atom_site_key:self.get_index()
        };//alt_loc は他の Atom も見ないと処理できないと思う
            
        return ret;
    }

    pub fn is_het(&self)->bool{
        if self.get_group_PDB() == "ATOM"{
            return false;
        }
        return true;
    }

    pub fn is_alt(&self)->bool{
        let lab:&str = self.get_label_alt_id();
        if lab == ""
        || lab == "."
        || lab == "?"
        || lab == " "
        || lab == "A"{
            return false;
        }
        return true;
    }
    
    pub fn new(km:&'a HashMap<String,usize>,val:&'a Vec<Box<String>>)->AtomSite<'a>{
        return AtomSite{
            keymap:km,values:val
        };
    }
    pub fn get_group_PDB(&self) -> &str{return self.get_value_of(_ATOM_SITE_GROUP_PDB);} 
    pub fn get_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_ID);} 
    pub fn get_type_symbol(&self) -> &str{return self.get_value_of(_ATOM_SITE_TYPE_SYMBOL);} 
    pub fn get_label_atom_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_ATOM_ID);} 
    pub fn get_label_alt_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_ALT_ID);} 
    pub fn get_label_comp_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_COMP_ID);} 
    pub fn get_label_asym_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_ASYM_ID);} 
    pub fn get_label_entity_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_ENTITY_ID);} 
    pub fn get_label_seq_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_LABEL_SEQ_ID);} 
    pub fn get_pdbx_PDB_ins_code(&self) -> &str{return self.get_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE);} 
    pub fn get_Cartn_x(&self) -> &str{return self.get_value_of(_ATOM_SITE_CARTN_X);} 
    pub fn get_Cartn_y(&self) -> &str{return self.get_value_of(_ATOM_SITE_CARTN_Y);} 
    pub fn get_Cartn_z(&self) -> &str{return self.get_value_of(_ATOM_SITE_CARTN_Z);} 
    pub fn get_occupancy(&self) -> &str{return self.get_value_of(_ATOM_SITE_OCCUPANCY);} 
    pub fn get_B_iso_or_equiv(&self) -> &str{return self.get_value_of(_ATOM_SITE_B_ISO_OR_EQUIV);} 
    pub fn get_pdbx_formal_charge(&self) -> &str{return self.get_value_of(_ATOM_SITE_PDBX_FORMAL_CHARGE);} 
    pub fn get_auth_seq_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_AUTH_SEQ_ID);} 
    pub fn get_auth_comp_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_AUTH_COMP_ID);} 
    pub fn get_auth_asym_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_AUTH_ASYM_ID);} 
    pub fn get_auth_atom_id(&self) -> &str{return self.get_value_of(_ATOM_SITE_AUTH_ATOM_ID);} 
    pub fn get_pdbx_PDB_model_num(&self) -> &str{return self.get_value_of(_ATOM_SITE_PDBX_PDB_MODEL_NUM);} 

    pub fn atom_site_to_pdbatom(&self,use_auth:bool)->PDBAtom{
        let mut ret = PDBAtom::new();
        if self.get_index() < 0{
            panic!("update_atom_site_index must have been performed at first.");
        }

        ret.set_atom_site_key(self.get_index());
        ret.set_xyz(self.get_value_of(_ATOM_SITE_CARTN_X).parse::<f64>().unwrap_or_else(|_|panic!("can not parse x"))
        ,self.get_value_of(_ATOM_SITE_CARTN_Y).parse::<f64>().unwrap_or_else(|_|panic!("can not parse y"))
        ,self.get_value_of(_ATOM_SITE_CARTN_Z).parse::<f64>().unwrap_or_else(|_|panic!("can not parse z"))
        );
       
        ret.serial_number = self.get_value_of(_ATOM_SITE_ID).parse::<i64>().expect("Cannot parse atom id.");
        ret.atom_symbol = self.get_value_of(_ATOM_SITE_TYPE_SYMBOL).to_string();
        ret.alt_code = self.get_value_of(_ATOM_SITE_LABEL_ALT_ID).to_string();
        if !self.is_alt(){
            ret.alt_code = "".to_string();
        }

        ret.dummy = false;
        if self.get_value_of(_ATOM_SITE_GROUP_PDB) == "HETATM"{
            ret.het = true;
        }else{
            ret.het = false;
        }
        
        if self.keymap.contains_key(_ATOM_SITE_PDBX_FORMAL_CHARGE){
            let fcc:String = self.get_value_of(_ATOM_SITE_PDBX_FORMAL_CHARGE).to_string();
            if fcc == "?" || fcc == "." || fcc == "" {
                ret.charge = None;
            }else{
                ret.charge = Some(fcc);
            }
        }
        if self.keymap.contains_key(_ATOM_SITE_OCCUPANCY){
            ret.occupancy = self.get_value_of(_ATOM_SITE_OCCUPANCY).parse::<f64>().unwrap();
        }
        if self.keymap.contains_key(_ATOM_SITE_B_ISO_OR_EQUIV){
            ret.temp_factor = self.get_value_of(_ATOM_SITE_B_ISO_OR_EQUIV).parse::<f64>().unwrap();
        }

        if use_auth{
            ret.atom_code = self.get_value_of(_ATOM_SITE_AUTH_ATOM_ID).to_string();
            
        }else{
            ret.atom_code = self.get_value_of(_ATOM_SITE_LABEL_ATOM_ID).to_string();
        }
        
        return ret;
    }
    
    pub fn get_value_of(&self,key:&str)->&str{
        return match self.keymap.get(key){
            Some(x)=>self.values[*x].as_str(),
            None=>"?"};
    }
    pub fn get_unique_residue_label(&self)->String{
        return "".to_string()
        +self.get_value_of(_ATOM_SITE_LABEL_COMP_ID)
        +"#"
        +self.get_value_of(_ATOM_SITE_LABEL_SEQ_ID)
        +"#"
        +self.get_value_of(_ATOM_SITE_AUTH_COMP_ID)
        +"#"
        +self.get_value_of(_ATOM_SITE_AUTH_SEQ_ID)
        +"#"
        +self.get_value_of(_ATOM_SITE_PDBX_PDB_INS_CODE)
        ;
    }
}

pub fn load_pdb(filename:&str) ->PDBEntry{
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    
    //let mut lcount:i64 = 0;

    let _noline = Regex::new(r"^[\r\n]*$").unwrap();
    let mut records:Vec<Vec<Box<String>>> = vec![];
    let mut ligand_records:Vec<Vec<String>> = vec![];
    let mut terflag:bool = false;
    let mut possibly_ligand = false;
    let mut current_model_num:String = "".to_owned();
    let keyvec:Vec<String> = (vec![
        _ATOM_SITE_GROUP_PDB,
        _ATOM_SITE_ID,
        _ATOM_SITE_LABEL_ATOM_ID,
        _ATOM_SITE_AUTH_ATOM_ID,
        _ATOM_SITE_LABEL_ALT_ID,
        _ATOM_SITE_LABEL_COMP_ID,
        _ATOM_SITE_AUTH_COMP_ID,
        _ATOM_SITE_LABEL_ASYM_ID,
        _ATOM_SITE_AUTH_ASYM_ID,
        _ATOM_SITE_LABEL_SEQ_ID,
        _ATOM_SITE_AUTH_SEQ_ID,
        _ATOM_SITE_PDBX_PDB_INS_CODE,
        _ATOM_SITE_CARTN_X,
        _ATOM_SITE_CARTN_Y,
        _ATOM_SITE_CARTN_Z,
        _ATOM_SITE_OCCUPANCY,
        _ATOM_SITE_B_ISO_OR_EQUIV,
        _ATOM_SITE_TYPE_SYMBOL,
        _ATOM_SITE_PDBX_FORMAL_CHARGE,
        _ATOM_SITE_PDBX_PDB_MODEL_NUM,
        ]).into_iter().map(|m|m.to_string()).collect();
    let keymap:HashMap<String,usize> = keyvec.iter().enumerate().map(|m|(m.1.clone(),m.0)).collect();




    for (_lcount,line) in reader.lines().enumerate() {

        let sstr = line.unwrap_or_else(|e|panic!("{:?}",e));
        
        if start_with(&sstr,"MODEL"){
            let ppt: Vec<&str> = sstr.split("[\\s;]+").collect();
            current_model_num = ppt[1].to_owned();
            terflag = false;
        }
        if start_with(&sstr,"ATOM") || start_with(&sstr,"HETATM"){
            let arecord = MMCIFEntry::parse_atom_line_pdb(&sstr,&current_model_num,&keymap);
            //println!("{:?}",arecord);
            if terflag{
                possibly_ligand = true;
            }
            records.push(arecord.into_iter().map(|m|Box::new(m)).collect());
        }else if start_with(&sstr,"TER"){
            terflag = true;
        }else{

        }
    }
    if possibly_ligand{
        eprintln!("There are possibly ligand records. ");
    }
    let mut mmcifdata = MMCIFEntry::new();
    mmcifdata.set_atom_site_section(keyvec,records);
    let mut ret:PDBEntry = PDBEntry::prepare_entry(mmcifdata);
    return ret;
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
            assert_eq!(mm.0.len(),mm.2[ii].len());
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

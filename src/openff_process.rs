
use roxmltree;
use std::collections::HashSet;
use super::debug_env::*;
use super::openff_energy;
use super::smirks_data;
use super::ff_env::*;
use regex::Regex;
use std::collections::HashMap;


//https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

/*
exclamation 	!e1 	not e1
ampersand 	e1&e2 	a1 and e2 (high precedence（優先）)
comma 	e1,e2 	e1 or e2
semicolon 	e1;e2 	a1 and e2 (low precedence)

[CH2] 	aliphatic carbon with two hydrogens (methylene carbon)
[!C;R] 	( NOT aliphatic carbon ) AND in ring
[!C;!R0] 	same as above ("!R0" means not in zero rings)
[n;H1] 	H-pyrrole nitrogen
[n&H1] 	same as above
[nH1] 	same as above
[c,n&H1] 	any arom carbon OR H-pyrrole nitrogen
[X3&H0] 	atom with 3 total bonds and no H's
[c,n;H1] 	(arom carbon OR arom nitrogen) and exactly one H
[Cl] 	any chlorine atom
[35*] 	any atom of mass 35
[35Cl] 	chlorine atom of mass 35
[F,Cl,Br,I] 	the 1st four halogens.
*/


//ネスト可能な文字列の結合を示す中間データ
#[derive(Debug)]
pub struct StringAtomConnector{
    pub atom:String,
    pub bonds_prev:Vec<String>,
    pub bonds_next:Vec<String>,
    pub connected_prev:Vec<usize>,
    pub connected_next:Vec<usize>,
    pub index_tokenlist:usize,//最初に分解した token 内の index
    pub index_in_vec:usize,//Atom 内の INDEX
    pub stereo_center:String
}
impl StringAtomConnector{
    fn new(s:&str,i:usize,i2:usize)->StringAtomConnector{
        assert!(i>=i2);
        return StringAtomConnector{
            atom:s.to_string(),
            bonds_prev:vec![],
            bonds_next:vec![],
            connected_next:vec![],
            connected_prev:vec![],
            index_tokenlist:i,
            index_in_vec:i2,
            stereo_center:"".to_owned()
            
        };
    }
    pub fn add_next_atom(&mut self,i:usize){
        self.connected_next.push(i);
    }
    pub fn add_prev_atom(&mut self,i:usize){
        self.connected_prev.push(i);
    }
}


pub struct SMIRKAtom{
    //"*", any
    //"a", aromatic
    //"A", aliphatic
    //"a,A","A,a", aromatic and aliphatic
    pub atom_type:String,
    pub atom_index:String,
    pub num_substituents:i64
}


impl SMIRKAtom{
    pub fn match_with(&self,atomcode:&str)->bool{
        return false;
    }
}

//start の括弧に対応する閉じ括弧のインデクスを返す
pub fn get_string_range(v:&Vec<String>,start:usize)->usize{
    assert_eq!(v[start],"(");
    let mut xcount:usize = 0;
    for s in (start+1)..v.len(){
        if v[s] == "("{
            xcount += 1;
        }
        if v[s] == ")"{
            if xcount == 0{
                return s;
            }else{
                xcount -= 1;
            }
        }
    }
    let mut zret:Vec<&str> = vec![];
    for s in start..v.len(){
        zret.push(&v[s]);
    }
    panic!("Couldn't find closing parenthesis! {:?}",zret);
}
pub fn tree_print(tmp_atomstring:&Vec<StringAtomConnector>,ppos:usize,parent_length:usize){
    print!("-{}",tmp_atomstring[ppos].atom);
    let nexx:&Vec<usize> = &tmp_atomstring[ppos].connected_next;
    //println!("{} {:?}",ppos,nexx);
    for (ii,nn) in nexx.iter().enumerate(){
        if ii != 0{
            println!("");
            print!("{}",vec![" ";parent_length+tmp_atomstring[ppos].atom.len()+1].into_iter().fold("".to_owned(),|s,m|s+m));
            tree_print(tmp_atomstring,*nn,parent_length+tmp_atomstring[ppos].atom.len()+1);
        }else{
            tree_print(tmp_atomstring,*nn,parent_length+tmp_atomstring[ppos].atom.len()+1);
        }
    }
}
pub fn smirks_to_molecule(smirks:&str)->FFMolecule{
    let mut ret = FFMolecule{
        atoms:vec![],bonds:vec![]
    };

    let mut cvec:Vec<String> = smirks.chars().map(|m|m.to_string()).collect();
    cvec.reverse();

    let mut token:Vec<String> = vec![];
    let ii:usize = 0;
    loop{
        let lcc:String = cvec.pop().unwrap();
        let mut vcc:Vec<String> = vec![];
        if lcc == "["{
            let mut pcc:String = cvec.pop().unwrap();
            loop{
                let pcc_:String = cvec.pop().unwrap();
                if pcc_ == "]"{
                    break;
                }
                
                pcc += pcc_.as_str();
            }
            vcc.push("[".to_owned()+&pcc+"]");
        }else{
            vcc.push(lcc);
        }
        token.push(vcc.into_iter().fold("".to_owned(),|s,m|s+&m));
        if cvec.len() == 0{
            break;
        }
    }
    let numtoken:usize = token.len();
    let mut tmp_atomstring:Vec<StringAtomConnector> = vec![];
    let regex_nonatom:Regex = Regex::new(r"^(-|=|#|$|\(|\)|/|\\|@)$").unwrap();
    let regex_bond:Regex = Regex::new(r"^(-|=|#|$|/|\\|@)$").unwrap();
    
    let slen = token.len();
    let mut iimap:Vec<i64> = vec![-1;slen];
    for (ii,tt) in token.iter().enumerate(){
        match regex_nonatom.captures(tt){
            Some(_x)=>{
            },
            _=>{
                let tii = tmp_atomstring.len();
                tmp_atomstring.push(StringAtomConnector::new(tt,ii,tii));
                iimap[ii] = tii as i64;
            }
        }
    }
    let tlen = tmp_atomstring.len();
    for kk in 0..tlen{
        let mut start:usize = tmp_atomstring[kk].index_tokenlist+1;
        
        if start >= slen{
            break;
        }
        
        if token[start] == ")"{
            continue;
        }
        loop{
            if token[start] == "("{
                let e = get_string_range(&token, start);
                for pp in start..slen{
                    if iimap[pp] > -1{
                        tmp_atomstring[kk].add_next_atom(iimap[pp] as usize);
                        break;
                    }
                }
                start = e+1;
                if start >= slen{
                    break;
                }
            }else{
                for pp in start..slen{
                    if iimap[pp] > -1{
                        tmp_atomstring[kk].add_next_atom(iimap[pp] as usize);
                        break;
                    }
                }
                break;
            }
            if start >= slen{
                break;
            }
        }
    }
    for tt in 0..tlen{
        let mut st = tmp_atomstring[tt].index_tokenlist+1;
        let mut bvec:Vec<String> = vec![];
        //後方についている bond を判定する
        while st < slen{
            match regex_bond.captures(&token[st]){
                Some(x) =>{
                    let b:&str = x.get(1).map(|m| m.as_str()).unwrap();
                    bvec.push(b.to_owned());
                },
                _ =>{
                    break;
                }
            }
            st += 1;
        }
        if bvec.len() > 0{
            if bvec[0] == "@"{
                if bvec.len() > 1{
                    if bvec[1] == "@"{
                        tmp_atomstring[tt].stereo_center = "@@".to_owned();
                        if bvec.len() > 2{
                            tmp_atomstring[tt].bonds_next.push(bvec[2].clone());
                        }
                    }else{
                        tmp_atomstring[tt].stereo_center = "@".to_owned();
                        tmp_atomstring[tt].bonds_next.push(bvec[1].clone());
                    }
                }
            }else{
                assert!(bvec.len() == 1);
                tmp_atomstring[tt].bonds_next.push(bvec[0].clone());
            }
        }

        let mut st:i64= tmp_atomstring[tt].index_tokenlist as i64 -1;
        let mut bvec:Vec<String> = vec![];
        //前方についている bond を判定する
        while st >= 0{
            match regex_bond.captures(&token[st as usize]){
                Some(x) =>{
                    let b:&str = x.get(1).map(|m| m.as_str()).unwrap();
                    bvec.push(b.to_owned());
                },
                _ =>{
                    break;
                }
            }
            st += 1;
        }
        assert!(bvec.len() <= 1);
        if bvec.len() == 1{
            tmp_atomstring[tt].bonds_next.push(bvec[0].clone());
            ここから
            原子インデクスで指定している場合があるはず
            つまり index_tokenlist は vec でないといけない。
        }
    }


    tree_print(&tmp_atomstring,0,0);
    
    let regex_num:Regex = Regex::new(r"([0-9]+)").unwrap();
    for ii in 0..numtoken{
        let start:usize= ii;
        let mut end:usize = ii;
        let mut buff:String = "".to_owned();
        for jj in ii..numtoken{
            if smirks_data::ELEMENT_NAME_TO_NUM.contains_key(&buff){
                end = jj;
            }
        }
    }

    return ret;
}


pub struct SMIRKBond{
    /*https://ja.wikipedia.org/wiki/SMILES%E8%A8%98%E6%B3%95
    結合は一次から順に-、=、#、$ で表される（ただし一重結合-は通常省略される）。二重結合=につながっている一重結合の向きを/、\で表すことでシス-トランス異性体を区別する。たとえばC/C=C\C、C/C=C/Cはそれぞれシス・トランス2-ブテンである。結合がないことは.で表現される（たとえば過酸化水素OOに対しO.Oは水2分子）。 
    */
    /*
    一重結合の向きを/、\で表すことでシス-トランス異性体を区別する。
    （二重三重結合を挟んだ状態でしか出現しないことが保証されているのだろうか？）
    */

    /*
    二桁のラベルを表すには%を前置する（たとえばC%12はラベル12）。 
    不斉中心には@または@@を後置し、根の方向から見てそれぞれ左回り・右回りに後続の原子団が並んでいることを表す（@が左回りのため）。たとえばS-アラニンのSMILESは、アミノ基を根にするとN[C@@H](C)C(=O)Oである（N[C@@]([H])(C)C(=O)Oのように書いてもよい）。 
    */

    //-,=,#,$
    pub bond_type:String,   
}
impl SMIRKBond{
    pub fn match_with(){

    }
}


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

#[test]
fn smirkstest(){
    smirks_to_molecule("C[C@@H](C(=O)O)N");
    println!("");
    smirks_to_molecule("C[C@H](C(=O)O)N");
    println!("");
}

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
    pub bonds_next:Vec<String>,
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
            bonds_next:vec![],
            connected_next:vec![],
            index_tokenlist:i,
            index_in_vec:i2,
            stereo_center:"".to_owned()
            
        };
    }
    pub fn add_next_atom(&mut self,atomindex:usize,ss:String){
        self.connected_next.push(atomindex);
        self.bonds_next.push(ss);
    }
    pub fn set_bond(&mut self,i:usize,ss:String){
        self.bonds_next[i] = ss;
    }
    pub fn to_molecule2d(atoms:&Vec<StringAtomConnector>)->Molecule2D{
        let mut ret = Molecule2D::new();
        let mut bonds_hs:HashMap<usize,HashSet<(usize,String)>> = HashMap::new();
        
        for (ii,aa) in atoms.iter().enumerate(){
            for (bii,bss) in aa.connected_next.iter().zip(aa.bonds_next.iter()){
                let ls = *(bii.min(&ii));
                let us = ii.min(*bii);
                if !bonds_hs.contains_key(&ls){
                    bonds_hs.insert(ls,HashSet::new());
                }
                let bpair = (us,bss.clone());
                if !bonds_hs.get(&ls).unwrap().contains(&bpair){
                    bonds_hs.get_mut(&ls).unwrap().insert(bpair);
                }else{
                    if &bonds_hs.get(&ls).unwrap().get(&bpair).unwrap().1 != bss{
                        panic!("Bond letter mismatch! {} vs {} ",&bonds_hs.get(&ls).unwrap().get(&bpair).unwrap().1,bss);
                    }
                }
            }
        }
        
        //atomindex1,atomindex2,bondindex の双方向のマップ
        //面倒なので大小は考えない
        let mut atom_to_bonds:HashMap<(usize,usize),usize> = HashMap::new();
        for uu  in bonds_hs.iter(){
            for hh in uu.1.iter(){
                let mut b = Bond2D::new();
                b.atom1 = *uu.0 as i64;
                b.atom2 = hh.0 as i64;
                b.bond_type.push(hh.1.clone());
                
                let a1 = b.atom1 as usize;
                let a2 = b.atom2 as usize;
                let bondindex = ret.bonds.len();
                if !atom_to_bonds.contains_key(&(a1,a2)){
                    atom_to_bonds.insert((a1,a2), bondindex);
                    atom_to_bonds.insert((a2,a1), bondindex);
                }
                ret.bonds.push(b);
            }
        }
        

        for (ii,aa) in atoms.iter().enumerate(){
            let mut atom:Atom2D =Atom2D::new();
            atom.atom_type = aa.atom.clone();
            atom.atom_index = ii as i64;
            atom.stereo_center = aa.stereo_center.clone();
            for bb in aa.connected_next.iter(){
                let code :(usize,usize) = (ii,*bb);
                if !atom_to_bonds.contains_key(&code){
                    panic!("{:?} is not found in the dictionary. {:?}",code,atom_to_bonds);
                }
                atom.bonds.push(*atom_to_bonds.get(&code).unwrap());
            }

            ret.atoms.push(
                atom
            );
        }

        return ret;
    }
}


#[derive(Debug)]
pub struct Bond2D{
    /*
    - 	single bond (aliphatic)
    / 	directional bond "up"1
    \ 	directional bond "down"1
    /? 	directional bond "up or unspecified"
    \? 	directional bond "down or unspecified"
    = 	double bond
    # 	triple bond
    : 	aromatic bond
    ~ 	any bond (wildcard)
    @ 	any ring bond1
    */
    pub bond_type:Vec<String>,
    pub atom1:i64,
    pub atom2:i64,
}
impl Bond2D{
    pub fn new() -> Bond2D{

        return Bond2D{
            bond_type:vec![],
            atom1:-1,atom2:-1
        }
    }
}

#[derive(Debug)]
pub struct Molecule2D{
    atoms:Vec<Atom2D>,
    bonds:Vec<Bond2D>
}
impl Molecule2D{
    pub fn new()->Molecule2D{
        return Molecule2D{atoms:vec![],bonds:vec![]};

    }
}

#[derive(Debug)]
pub struct Atom2D{
    pub atom_type:String,
    pub atom_index:i64,
    pub stereo_center:String,
    pub bonds:Vec<usize>,
}


impl Atom2D{
    pub fn new()->Atom2D{
        return Atom2D{
            atom_type: "".to_owned(),
            atom_index: -1,
            stereo_center: "".to_owned(),
            bonds:vec![]
        };
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
pub fn smirks_to_molecule(smirks:&str)->Molecule2D{

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
    let regex_nonatom:Regex = Regex::new(r"^(-|=|#|$|\(|\)|/|\\|@|[0-9])$").unwrap();
    let regex_bond:Regex = Regex::new(r"^(-|=|#|$|/|\\|@)$").unwrap();
    let regex_distant_connection_index:Regex = Regex::new(r"^[0-9]$").unwrap();
    
    let slen = token.len();
    let mut token_to_atom:Vec<i64> = vec![-1;slen];
    
    for (ii,tt) in token.iter().enumerate(){
        match regex_nonatom.captures(tt){
            Some(_x)=>{
            },
            _=>{
                let tii = tmp_atomstring.len();
                tmp_atomstring.push(StringAtomConnector::new(tt,ii,tii));
                token_to_atom[ii] = tii as i64;
            }
        }
    }
    
    //現在 Connection index （atom の後ろにある数値）を持っている Atom の atom List 上の Index と bondorder を示す文字列
    //前方にある記号を bond string として入れているが・・・
    let mut connectionindex_to_atomindex:HashMap<usize,(i64,Vec<String>)> = HashMap::new();
    let mut distant_connections:Vec<(usize,usize,Vec<String>)> = vec![];
    for ss in 0..slen{
        match regex_distant_connection_index.captures(&token[ss]){
            Some(x)=>{
                let mut st:i64 = ss as i64-1;
                let mut bond_string:Vec<String> = vec![];
                while st >= 0{//0 は atom であるはずだが一応
                    match regex_nonatom.captures(&token[ss]){
                        Some(_x) => {//前方にある Atom でないものをとる
                            match regex_bond.captures(&token[ss]){
                                Some(_y) => {//前方にある bond をとる 
                                        bond_string.push(token[ss].clone());
                                },
                                _ => {}
                            }
                        },
                        _ =>{
                            break;
                        }
                    }
                    st -= 1;
                }

                let mut st:i64 = ss as i64-1;
                while st > 0{
                    if token_to_atom[st as usize] > -1{
                        let cindex = x.get(1).map(|m| m.as_str()).unwrap().parse::<usize>().unwrap();
                        let defval:(i64,Vec<String>) = (-1,vec![]);//その Connection Index が空の場合は -1 
                        let cc = connectionindex_to_atomindex.get(&cindex).unwrap_or(&defval);
                        if cc.0 == -1{
                            //二桁のインデクスを持つものはないようだ。
                            //二桁ある場合、それは一桁＋一桁
                            //よって同じ数字が何度も出てくることがある
                            //文法としては一番最後に出現した同じ数字のインデクスを持つ Atom ということでよいのだろうか。
                            
                            connectionindex_to_atomindex.insert(cindex,(token_to_atom[st as usize],bond_string));
                        }else{
                            if cc.1.len() > 0 && bond_string.len() > 0{
                                if cc.1 != bond_string{
                                    panic!("different bond parameter! {:?} , {:?} ",cc.1,bond_string);
                                }
                            }
                            distant_connections.push((cc.0 as usize,ss,bond_string));
                            connectionindex_to_atomindex.insert(cindex,(-1,vec![]));
                        }
                        break;
                    }
                    st -= 1;
                }
            },
            _=>{
            }
        }
    }

    let mut normal_connection:Vec<(usize,usize)> = vec![];
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
            ここから
            やっぱり Bond も Connection と一緒に設定する
            if token[start] == "("{
                let e = get_string_range(&token, start);
                for pp in start..slen{
                    if token_to_atom[pp] > -1{
                        tmp_atomstring[kk].add_next_atom(token_to_atom[pp] as usize,"".to_owned());
                        normal_connection.push((token_to_atom[pp] as usize,kk));
                        break;
                    }
                }
                start = e+1;
                if start >= slen{
                    break;
                }
            }else{
                for pp in start..slen{
                    if token_to_atom[pp] > -1{
                        tmp_atomstring[kk].add_next_atom(token_to_atom[pp] as usize,"".to_owned());
                        normal_connection.push((token_to_atom[pp] as usize,kk));
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

    //bond は Atom-Atom というより単体の Atom についていると考えた方が楽そうなので、Atom と別に Parse する
    for tt in 0..tlen{
        let mut st = tmp_atomstring[tt].index_tokenlist+1;
        let mut bvec:Vec<String> = vec![];
        //stereo senter を判定する
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
                    }else{
                        tmp_atomstring[tt].stereo_center = "@".to_owned();
                    }
                }
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
        }
    }

    for (cc,ss,mut strr) in distant_connections.into_iter(){
        tmp_atomstring[cc].connected_next.push(ss);
        tmp_atomstring[cc].bonds_next.append(&mut strr);
    }




    tree_print(&tmp_atomstring,0,0);
    

    return StringAtomConnector::to_molecule2d(&tmp_atomstring);
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
    let mol = smirks_to_molecule("C[C@@H](C(=O)O)N");
    println!("{:?}",mol);
    smirks_to_molecule("C[C@H](C(=O)O)N");
    println!("");
}

use roxmltree;
use core::panic;
use std::collections::HashSet;
use std::mem;
use std::os::windows::process;
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
    

    //start の括弧に対応する閉じ括弧のインデクスを返す
    pub fn get_string_range(start:usize,v:&Vec<String>)->usize{
        assert_eq!(v[start],"(");
        let mut xcount:usize = 0;
        for s in start..v.len(){
            if v[s] == "("{
                xcount += 1;
            }
            if v[s] == ")"{
                xcount -= 1;
                if xcount == 0{
                    return s;
                }
            }
        }
        let mut zret:Vec<&str> = vec![];
        for s in start..v.len(){
            zret.push(&v[s]);
        }
        panic!("Couldn't find closing parenthesis! {:?}",zret);
    }

    pub fn create_group(start_pos:usize,end_pos:usize,tokens:&Vec<String>)->(Vec<usize>,Vec<Vec<usize>>){
        let mut ssin:Vec<usize> = vec![start_pos];
        let mut ggro:Vec<Vec<usize>> = vec![];

        let mut current_atom:i64 = start_pos as i64;
        let mut processed:HashSet<usize> = HashSet::new();
        processed.insert(start_pos);
        let mut multi_start:usize = start_pos;
        for ii in start_pos..=end_pos{
            multi_start = ii;
            if tokens[ii] == "("{
                break;
            }
            ssin.push(ii);
        }
        if multi_start <= end_pos{
            loop{
                if tokens[multi_start] != "("{
                    let mut pgro:Vec<usize> = vec![];
                    for ii in multi_start..=end_pos{
                        pgro.push(ii);
                    }
                    ggro.push(pgro);
                    break;
                }
                let inex = StringAtomConnector::get_string_range(multi_start, &tokens);
                assert_eq!(tokens[inex],")");
                let mut pgro:Vec<usize> = vec![];
                for ii in (multi_start+1)..=(inex-1){
                    pgro.push(ii);
                }
                ggro.push(pgro);
                if inex < end_pos{
                    multi_start += 1;
                }else{
                    break;
                }
            }
        }
        return (ssin,ggro);
    }

    pub fn smirks_to_molecule(smirks:&str)->Molecule2D{

        let mut cvec:Vec<String> = smirks.chars().map(|m|m.to_string()).collect();
        cvec.reverse();
    
        let mut token:Vec<String> = vec![];
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
        let mut tmp_atomstring:Vec<StringAtomConnector> = vec![];
        let regex_nonatom:Regex = Regex::new(r"^(-|=|#|$|\(|\)|/|\\|@|[0-9])$").unwrap();
        let regex_special:Regex = Regex::new(r"^(-|=|#|$|/|\\|@)$").unwrap();
        let regex_distant_connection_index:Regex = Regex::new(r"^([0-9])$").unwrap();
        
        let mut token_to_atom:Vec<i64> = vec![-1;token.len()];
        
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

        
        let mut connection_normal_:Vec<(usize,usize)> = vec![];
        let mut updated_groups:Vec<(i64,(usize,usize))> = vec![(-1,(0,token.len()-1))];
        while updated_groups.len() > 0{
            let tt = updated_groups.pop().unwrap();
            let start_pos = (tt.1).0;
            let end_pos = (tt.1).1;
            let (singles,multis) = StringAtomConnector::create_group(start_pos, end_pos, &token);
            let mut parentid = tt.0;
            assert!(singles.len() > 0);
            if parentid > -1{
                connection_normal_.push((parentid as usize,singles[0]));
            }
            parentid = singles[0] as i64;
            for ss in 1..(singles.len()){
                connection_normal_.push((parentid as usize,singles[ss]));
                parentid = singles[ss] as i64;
            }
            if multis.len() > 0{
                for mm in multis.into_iter(){
                    updated_groups.push((parentid,(mm[0],mm[mm.len()-1])));
                }
            }
        }
        //atom 以外の記号が結合パートナーとして入っている可能性があるので注意
        let mut connection_normal_hs:HashSet<(usize,usize)> = HashSet::new();
        for cc in connection_normal_.into_iter(){
            let mut st = cc.0 as i64;
            let mut en = cc.1;
            while st >= 0{
                if token_to_atom[st as usize] > -1{
                    break;
                }
                if st == 0{
                    panic!("??");
                }
                st -= 1;
            }
            while en < token.len(){
                if token_to_atom[en] > -1{
                    break;
                }
                if en == token.len()-1{
                    panic!("????");
                }
                en += 1;
            }
            connection_normal_hs.insert((st as usize,en));
        }
        
        let mut connection_normal:Vec<(usize,usize)> = connection_normal_hs.into_iter().collect();
        connection_normal.sort();

        //例として
        //http://www.chemspider.com/Chemical-Structure.65325987.html

        let mut ex_prev:Vec<Vec<String>> = vec![vec![];token.len()];
        let mut ex_next:Vec<Vec<String>> = vec![vec![];token.len()];
        //どちらに入っているか理解してないので後でテストすること
        //特に離れた原子に対する結合で特殊な結合の場合どちらに入っているのか・・・
        let regex_prev:Regex = Regex::new(r"^(-|=|#|$|/|\\)$").unwrap();
        let regex_next:Regex = Regex::new(r"^(@|[0-9]|-|=|#|$|/|\\)$").unwrap();
        for ii in 0..token.len(){
            if token_to_atom[ii] > -1{
                let mut prev = ii as i64 -1;
                while prev > -1{
                    match regex_prev.captures(&token[prev as usize]){
                        Some(_x)=>{
                            ex_prev[ii].push(token[prev as usize].clone());
                        }
                        ,
                        _=>{break;}
                    }
                    prev -= 1;
                }
                let mut nex:usize = ii +1;
                while nex < token.len(){
                    match regex_next.captures(&token[nex]){
                        Some(_x)=>{
                            ex_next[ii].push(token[nex].clone());
                        }
                        ,
                        _=>{break;}
                    }
                    nex += 1;
                }
            }
        }


        //インデックスで結合パートナーを指定する Bond の処理
        //現在 Connection index （atom の後ろにある数値）を持っている Atom の atom List 上の Index と bondorder を示す文字列
        //数字より前方にある記号を bond string として入れているが・・・
        let mut connectionindex_to_tokenindex:HashMap<usize,(i64,Vec<String>)> = HashMap::new();
        let mut distant_connections:Vec<(usize,usize,Vec<String>)> = vec![];
        let mut bond_string:Vec<String> = vec![];
        for ss in 0..token.len(){
            if ex_next[ss].len() == 0{
                continue;
            }
            for pp in 0..ex_next[ss].len(){
                if let None =  regex_distant_connection_index.find(&ex_next[ss][pp]){
                    bond_string.push(ex_next[ss][pp].clone());
                    continue;
                }
                let cindex:usize = ex_next[ss][pp].parse::<usize>().unwrap();
                let defval:(i64,Vec<String>) = (-1,vec![]);//その Connection Index が空の場合は -1 
                let cc = connectionindex_to_tokenindex.get(&cindex).unwrap_or(&defval);
                if cc.0 == -1{
                    //二桁のインデクスを持つものはないようだ。
                    //二桁ある場合、それは一桁＋一桁
                    //よって同じ数字が何度も出てくることがある
                    //文法としては一番最後に出現した同じ数字のインデクスを持つ Atom ということでよいのだろうか。
                    
                    connectionindex_to_tokenindex.insert(cindex,(ss as i64,bond_string));
                    bond_string = vec![];
                }else{
                    //既にペアになる原子が出現している
                    if cc.1.len() > 0 && bond_string.len() > 0{
                        if cc.1 != bond_string{
                            panic!("different bond parameter! {:?} , {:?} ",cc.1,bond_string);
                        }
                    }
                    //ペアとして登録する
                    if bond_string.len() > 0{
                        distant_connections.push((cc.0 as usize,ss,bond_string));
                    }else{
                        distant_connections.push((cc.0 as usize,ss,cc.1.clone()));
                    }


                    //使用された原子の登録は初期化する
                    connectionindex_to_tokenindex.insert(cindex,(-1,vec![]));
                    bond_string = vec![];
                }
            }
        }
    
        
        let mut connected_atoms:Vec<(usize,usize,Vec<String>)> = vec![];
        let tlen = tmp_atomstring.len();
        for kk in 0..tlen{
            let nxx = &ex_next[tmp_atomstring[kk].index_tokenlist];
            for n in nxx.iter(){
                if n == "@"{
                    tmp_atomstring[kk].stereo_center += n;
                }
            }
        }

        for cc in connection_normal.into_iter(){
            tmp_atomstring[token_to_atom[cc.0] as usize].connected_next.push(token_to_atom[cc.1] as usize);
            for pp in ex_prev[cc.1].iter(){
                tmp_atomstring[token_to_atom[cc.0] as usize].bonds_next.push(pp.clone());
            }
        }
    
        for (cc,ss,mut strr) in distant_connections.into_iter(){
            tmp_atomstring[cc].connected_next.push(ss);
            tmp_atomstring[cc].bonds_next.append(&mut strr);
        }
    
    
    
    
        //tree_print(&tmp_atomstring,0,0);
        
    
        return StringAtomConnector::to_molecule2d(&tmp_atomstring);
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
                let us = ii.max(*bii);
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

    pub fn create_group(start_atom:usize,candidates:&HashSet<usize>,connected:&HashMap<usize,Vec<usize>>)->(Vec<usize>,Vec<Vec<usize>>){
        let mut ssin:Vec<usize> = vec![start_atom];
        let mut ggro:Vec<Vec<usize>> = vec![];

        let mut current_atom:i64 = start_atom as i64;
        let mut processed:HashSet<usize> = HashSet::new();
        processed.insert(start_atom);
        loop{
            let catom:usize = current_atom as usize;
            if !connected.contains_key(&catom){
                break;
            }
            let mut zvec:Vec<usize> =  connected.get(&catom).unwrap().iter().filter(|m| candidates.contains(m)).map(|m| *m).collect();
            if zvec.len() == 0{
                break;
            }
            if zvec.len() == 1{
                let nex = zvec[0];
                if processed.contains(&nex){
                    break;
                }
                current_atom = nex as i64;
                
                ssin.push(nex);
                processed.insert(nex);
            }else{
                while zvec.len() > 0{
                    let mut updated:Vec<usize> = vec![zvec.pop().unwrap()];
                    let mut members:Vec<usize> = vec![];
                    while updated.len() > 0{
                        let uu = updated.pop().unwrap();
                        if !candidates.contains(&uu){
                            continue;
                        }
                        if processed.contains(&uu){
                            continue;
                        }
                        processed.insert(uu);
                        members.push(uu);
                        if !connected.contains_key(&uu){
                            continue;
                        }
                        for pp in connected.get(&uu).unwrap().iter(){
                            updated.push(*pp);
                        }
                    }
                    if members.len() > 0{
                        ggro.push(members);
                    }
                }
                break;
            }
        }
        return (ssin,ggro);
    }

    pub fn print_members(&self,targetid:usize,members:&Vec<Vec<usize>>,children:&Vec<Vec<usize>>)->String{
        let mut ret:String = members[targetid].iter().map(|m|&self.atoms[*m].atom_type).fold("".to_owned()
        ,|s,m|{s+m});
        for ll in 0..children[targetid].len(){
            ret = ret+"("+&self.print_members(children[targetid][ll],members,children)+")";
        }
        return ret;
    }
    pub fn to_smirks(&self)->String{
//キラリティ入ってない
//Circler 入ってない
        let mut members:Vec<Vec<usize>> = vec![];
        let mut children:Vec<Vec<usize>> = vec![];

        let mut bondss:HashMap<usize,Vec<usize>> = HashMap::new();
        for bb in self.bonds.iter(){
            let a1:usize = bb.atom1 as usize;
            let a2:usize = bb.atom2 as usize;
            if !bondss.contains_key(&a1){
                bondss.insert(a1, vec![]);
            }
            if !bondss.contains_key(&a2){
                bondss.insert(a2, vec![]);
            }
            bondss.get_mut(&a1).unwrap().push(a2);
            bondss.get_mut(&a2).unwrap().push(a1);
        }
        println!("\n{:?}",self.bonds);
        let mut start:usize;
        let mut target_groups:Vec<(i64,Vec<usize>)> = vec![(-1,(0..self.atoms.len()).into_iter().collect())];
        while target_groups.len() > 0{
            let tg = target_groups.pop().unwrap();
            println!("\n{:?}",tg);
            if tg.1.len() == 0{
                continue;
            }
            start = tg.1[0];
            let candidates:HashSet<usize> = tg.1.iter().map(|m| *m).collect();
            let (singlebond,groups) = Molecule2D::create_group(start,&candidates,&bondss);
            //println!("{:?}\n{:?}",singlebond,groups);
            let nodeid:usize = members.len();
            if tg.0 > -1{
                children[tg.0 as usize].push(nodeid);
            }
            members.push(singlebond);
            children.push(vec![]);
            if groups.len() > 0{
                for gg in groups.into_iter(){
                    target_groups.push((nodeid as i64,gg));
                }
            }
        }
        println!("{:?}\n{:?}",members,children);
        return self.print_members(0,&members,&children);

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
    let smirkstring = "C[C@@H](C(=O)O)N";
    println!("{}",smirkstring);
    let mol = StringAtomConnector::smirks_to_molecule(smirkstring);
    //ここから
    //色々間違ってそう
    println!("\n+++++++{}+++++",mol.to_smirks());
    //println!("\n{:?}",mol);
    StringAtomConnector::smirks_to_molecule("C[C@H](C(=O)O)N");
    //println!("");

}
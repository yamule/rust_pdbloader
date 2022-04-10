
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

const LEFT_PARENTHESIS:i64 = -1;
const RIGHT_PARENTHESIS:i64 = -2;

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
        let mut ssin:Vec<usize> = vec![];
        let mut ggro:Vec<Vec<usize>> = vec![];

        //let mut processed:HashSet<usize> = HashSet::new();//いらんかな
        //processed.insert(start_pos);
        let mut multi_start:i64 = -1;
        for ii in start_pos..=end_pos{
            if tokens[ii] == "("{
                multi_start = ii as i64;
                break;
            }
            ssin.push(ii);
        }
        if multi_start > -1{
            let mut multi_start :usize = multi_start as usize;
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
                    multi_start = inex+1;
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
            if lcc == "%"{//二桁ラベル
                vcc.push(lcc+&cvec.pop().unwrap()+&cvec.pop().unwrap());
            }else if lcc == "["{
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
        //大括弧だけ分離した

        let mut tmp_atomstring:Vec<StringAtomConnector> = vec![];
        let regex_nonatom:Regex = Regex::new(r"^(-|=|#|$|\(|\)|/|\\|@|[0-9])$").unwrap();
        let _regex_special:Regex = Regex::new(r"^(-|=|#|$|/|\\|@)$").unwrap();
        let regex_distant_connection_index:Regex = Regex::new(r"^([0-9]|%[0-9][0-9])$").unwrap();
        
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
        //atom である token には token_to_atom に 0 以上の数を入れた
        
        let mut connection_normal_:Vec<(usize,usize)> = vec![];
        let mut updated_groups:Vec<(i64,(usize,usize))> = vec![(-1,(0,token.len()-1))];
        let mut bond_listorder:Vec<Vec<usize>> = vec![vec![];token.len()];
        while updated_groups.len() > 0{
            let tt = updated_groups.pop().unwrap();
            let parentid_ = tt.0;
            let start_pos = (tt.1).0;
            let end_pos = (tt.1).1;
            if start_pos == end_pos{
                if parentid_ > -1 && token_to_atom[start_pos] > -1{
                    connection_normal_.push((parentid_ as usize,start_pos));
                    bond_listorder[parentid_ as usize].push(start_pos);
                }
                continue;
            }
            let (singles,multis) = StringAtomConnector::create_group(start_pos, end_pos, &token);
            assert!(singles.len() > 0);
            let mut filtered:Vec<usize> = vec![];
            for ss in 0..singles.len(){
                if token_to_atom[singles[ss]] > -1{
                    filtered.push(singles[ss]);
                }
            }
            if parentid_ > -1{
                connection_normal_.push((parentid_ as usize,filtered[0]));
                //bond_listorder には既に入っているはずなので加えない
            }
            if filtered.len() == 0{
                panic!("????");
            }
            let mut parentid_next = filtered[0];

            for ss in 1..filtered.len(){
                connection_normal_.push((parentid_next,filtered[ss]));
                bond_listorder[parentid_next].push(filtered[ss]);
                parentid_next = filtered[ss] ;
            }
            if multis.len() > 0{
                for mm in multis.into_iter(){
                    let mut natt: i64 = -1;
                    for nn in mm.iter(){
                        if token_to_atom[*nn] > -1{
                            natt = *nn as i64;
                            break;
                        }
                    }
                    if natt > -1{
                        bond_listorder[parentid_next].push(natt as usize);
                    }else{
                        panic!("??!");
                    }
                    updated_groups.push((parentid_next as i64,(mm[0],mm[mm.len()-1])));
                }
            }
        }
        //connection_normal_ には atom 以外の記号が結合パートナーとして入っているので後方前方を探して原子同士にする
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
                en += 1;
                if en == token.len(){
                    break;
                }
            }
            if en == token.len(){//Distant connection に対する記号が入っているだけのはず
                panic!("?????");//入らないようにしたかな？
            }
            connection_normal_hs.insert((st as usize,en));
        }
        
        let mut connection_normal:Vec<(usize,usize)> = connection_normal_hs.into_iter().collect();
        connection_normal.sort();
        
        let mut ex_prev:Vec<Vec<String>> = vec![vec![];token.len()];
        let mut ex_next:Vec<Vec<String>> = vec![vec![];token.len()];
        //どちらに入っているか理解してないので後でテストすること
        //特に離れた原子に対する結合で特殊な結合の場合どちらに入っているのか・・・
        let regex_prev:Regex = Regex::new(r"^(-|=|#|$|/|\\)$").unwrap();
        let regex_next:Regex = Regex::new(r"^(@|[0-9]|-|=|#|$|/|\\|%[0-9][0-9])$").unwrap();
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
        
        let mut distant_connections_atomid_index_nextatomid:Vec<HashMap<usize,usize>> = vec![HashMap::new();token.len()];
        for ss in 0..token.len(){
            if ex_next[ss].len() == 0{
                continue;
            }

            let mut bond_string:Vec<String> = vec![];
            for pp in 0..ex_next[ss].len(){
                if let None =  regex_distant_connection_index.find(&ex_next[ss][pp]){
                    //原子インデクス以外の記号
                    bond_string.push(ex_next[ss][pp].clone());
                    continue;
                }
                let cindex_ = if ex_next[ss][pp].len() == 1{
                    ex_next[ss][pp].clone()
                }else{
                    //二桁の場合
                    let cz:Vec<char> = ex_next[ss][pp].chars().into_iter().collect();
                    cz[1].to_string()+&cz[2].to_string()
                };
                let cindex:usize = cindex_.parse::<usize>().unwrap();//connection のインデクス
                let defval:(i64,Vec<String>) = (-1,vec![]);//その Connection Index が空の場合は -1 
                //connection の相手と、結合修飾文字列
                let cc = connectionindex_to_tokenindex.get(&cindex).unwrap_or(&defval);
                if cc.0 == -1{
                    //二桁のインデクスを持つものは % が接頭辞として付く。
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
                    distant_connections_atomid_index_nextatomid[cc.0 as usize].insert(cindex,ss);
                    distant_connections_atomid_index_nextatomid[ss].insert(cindex,cc.0 as usize);

                    //使用された原子の登録は初期化する
                    connectionindex_to_tokenindex.insert(cindex,(-1,vec![]));
                    bond_string = vec![];
                }
            }
        }
    
        
        let tlen = tmp_atomstring.len();
        for kk in 0..tlen{
            let nxx = &ex_next[tmp_atomstring[kk].index_tokenlist];
            let mut bonder:Vec<usize> = vec![];//distant connection がある場合に順番通りに詰めていく
            for n in nxx.iter(){
                if let Some(x) =  regex_distant_connection_index.captures(n){
                    let s = x.get(1).unwrap().as_str();
                    let cindex_ = if s.len() == 1{
                        s.to_string()
                    }else{
                        //二桁の場合
                        let cz:Vec<char> = s.chars().into_iter().collect();
                        assert!(cz.len()==3,"Index must be made from 2 letters. {} ",s);
                        cz[1].to_string()+&cz[2].to_string()
                    };
                    let cindex = cindex_.parse::<usize>().unwrap();
                    bonder.push(*distant_connections_atomid_index_nextatomid[tmp_atomstring[kk].index_tokenlist].get(&cindex).unwrap_or_else(
                        ||panic!("position: {} id: {} is not found. {:?}",kk,cindex,distant_connections_atomid_index_nextatomid)));
                }
                if n == "@"{
                    tmp_atomstring[kk].stereo_center += n;
                }    
            }
            if bonder.len() > 0{
                bonder.append(&mut bond_listorder[kk]);
                bond_listorder[kk].append(&mut bonder);//前に詰める
            }
            
            ここから
            
            テストの出力を RDKIT とかに読ませて同じものが出るか見る
            水素を加える
            get_clockwise_atoms(axis//こちら側にくる原子)
            みたいな関数を作って、stereocenter の処理をする
            N@ もあるっぽい。ダミーの H を入れるか。-> canonical smiles では取れているようだ
            

        }

        for cc in connection_normal.into_iter(){
            tmp_atomstring[token_to_atom[cc.0] as usize].connected_next.push(token_to_atom[cc.1] as usize);
            let bondstring:String = ex_prev[cc.1].iter().fold("".to_owned(),|s,m|s+m);
            if bondstring.len() == 0{
                tmp_atomstring[token_to_atom[cc.0] as usize].bonds_next.push("-".to_owned());
            }else{
                tmp_atomstring[token_to_atom[cc.0] as usize].bonds_next.push(bondstring);
            }
        }
    
        for (cc,ss,strr) in distant_connections.into_iter(){
            tmp_atomstring[token_to_atom[cc] as usize].connected_next.push(token_to_atom[ss] as usize);
            let bondstring:String = strr.iter().fold("".to_owned(),|s,m|s+m);
            tmp_atomstring[token_to_atom[cc] as usize].bonds_next.push(bondstring);
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
        let mut bonds_v:Vec<(usize,HashSet<(usize,String)>)> = bonds_hs.into_iter().map(|m|(m.0,m.1)).collect();
        bonds_v.sort_by(|a,b|a.0.cmp(&b.0));

        for uu  in bonds_v.into_iter(){
            let target = uu.0;
            let mut partner:Vec<(usize,String)> = uu.1.into_iter().collect();
            partner.sort();
            for hh in partner.iter(){
                let mut b = Bond2D::new();
                b.atom1 = target as i64;
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

    //startatom から初めて、connected を参考に Atom を追加していく
    //結合が二つ以上に分かれた場合、ggro にそれぞれ別れた道の先にある Atom を入れて返す
    //Circler になっている場合についてきちんとまだ確認していない
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
                        //面倒だからこうしたが、遅いだろうか
                        updated.sort();
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

    //members には Atom が分岐なくつながったグループが入っており、これは Molecule 全てを網羅する
    //children には Members の最後の Atom から分岐した場合の分岐先の members のインデクスが入っている
    //これらを SMIRKS 形式に類似した形式で Vec に atoms のインデクスを入れて返す。
    //二番目の要素は、暗黙的に示されるため結合を作成するときに必要でない結合のペアのタプル
    //括弧はそれぞれ LEFT_PARENTHESIS, RIGHT_PARENTHESIS という負の定数
    pub fn array_of_members(&self,targetid:usize,members:&Vec<Vec<usize>>,children:&Vec<Vec<usize>>)
    ->(Vec<i64>,Vec<(usize,usize)>){
        let mut ret:Vec<i64> = vec![];
        let mut implied:Vec<(usize,usize)> = vec![];
        for ii in 0..members[targetid].len(){
            let m = members[targetid][ii];
            if ii < members[targetid].len() -1{
                implied.push((members[targetid][ii],members[targetid][ii+1]));
            }
            ret.push(m as i64);
        }
        for ll in 0..children[targetid].len(){
            implied.push((members[targetid][members[targetid].len()-1],members[children[targetid][ll]][0]));
            if ll != children[targetid].len()-1{
                ret.push(LEFT_PARENTHESIS);
                let (mut arr,mut bb) =  self.array_of_members(children[targetid][ll],members,children);
                ret.append(&mut arr);
                implied.append(&mut bb);
                ret.push(RIGHT_PARENTHESIS);
            }else{
                let (mut arr,mut bb) =  self.array_of_members(children[targetid][ll],members,children);
                ret.append(&mut arr);
                implied.append(&mut bb);
            }
        }
        return (ret,implied);
    }
    pub fn to_smirks(&self)->String{
        //members には Single bond でつながった塊が入っている
        //children は Multi で別れた際にどの members に繋がっているかが入っている
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
        let mut start:usize;
        let mut target_groups:Vec<(i64,Vec<usize>)> = vec![(-1,(0..self.atoms.len()).into_iter().collect())];
        while target_groups.len() > 0{
            let tg = target_groups.pop().unwrap();

            if tg.1.len() == 0{
                continue;
            }
            start = tg.1[0];
            let candidates:HashSet<usize> = tg.1.iter().map(|m| *m).collect();

            let (singlebond,groups) = Molecule2D::create_group(start,&candidates,&bondss);
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
        //return self.print_members(0,&members,&children);
        let arr =self.array_of_members(0, &members, &children);
        let mut bondmap:HashMap<usize,Vec<(usize,String)>> = HashMap::new();
        let mut bondmap_rev:HashMap<usize,Vec<(usize,String)>> = HashMap::new();
        for bb in self.bonds.iter(){
            let a1 = bb.atom1.min(bb.atom2) as usize;
            let a2 = bb.atom1.max(bb.atom2) as usize;
            if !bondmap.contains_key(&a1){
                bondmap.insert(a1,vec![]);
            }
            if !bondmap_rev.contains_key(&a2){
                bondmap_rev.insert(a2,vec![]);
            }
            bondmap.get_mut(&a1).unwrap().push((a2,bb.bond_type.iter().fold("".to_owned(),|s,m|s+m)));
            bondmap_rev.get_mut(&a2).unwrap().push((a1,bb.bond_type.iter().fold("".to_owned(),|s,m|s+m)));
        }
        let mut res:Vec<String> = vec![];
        let mut atom_vecmap:Vec<i64> = vec![-1;self.atoms.len()];
        let implied_bonds:HashSet<(usize,usize)> = arr.1.iter().fold(HashSet::new(),|mut s,m|{s.insert((m.0,m.1)); s.insert((m.1,m.0));s});
        for aa in arr.0.into_iter(){
            if aa == LEFT_PARENTHESIS{
                res.push("(".to_owned());
            }else if aa == RIGHT_PARENTHESIS{
                res.push(")".to_owned());
            }else{
                assert!(aa > -1);
                let aa:usize = aa as usize;
                atom_vecmap[aa] = res.len() as i64;
                res.push(self.atoms[aa].atom_type.clone());
            }
        }
        
        let mut distbond_used:Vec<i64> = vec![-1;res.len()];
        for ii in 0..self.atoms.len(){
            //if bondmap_rev.contains_key(&ii){
            //    let bb = bondmap_rev.get(&ii).unwrap();
            if bondmap.contains_key(&ii){
                let bb = bondmap.get(&ii).unwrap();
                for kk in 0..bb.len(){
                    let ast = (atom_vecmap[bb[kk].0] as usize).min(atom_vecmap[ii] as usize);
                    let aen = (atom_vecmap[bb[kk].0] as usize).max(atom_vecmap[ii] as usize);
                    if implied_bonds.contains(&(ii,bb[kk].0)){
                        if bb[kk].1 != "-"{
                            res[aen] = bb[kk].1.clone()+&res[aen];
                        }
                    }else{
                        let mut distbond_check:HashSet<i64> = HashSet::new();
                        let mut distbond_check_:Vec<i64> = vec![0;100];
                        for pp in 0..ast{
                            if distbond_used[pp] > 0{
                                distbond_check_[distbond_used[pp] as usize] += 1;
                            }
                        }
                        for ii in 0..100{
                            if distbond_check_[ii] %2 == 1{
                                distbond_check.insert(ii as i64);
                            }
                        }

                        for pp in ast..=aen{
                            distbond_check.insert(distbond_used[pp]);
                        }

                        //二つの結合した原子の間にある原子で使われていないインデクスを探す
                        let mut distbond_index = -1;
                        for pp in 1..=99{
                            if !distbond_check.contains(&pp){
                                distbond_index = pp;
                                break;
                            }
                        }
                        if distbond_index < 1{
                            panic!("This molecule cannot express as smirks???");
                        }

                        distbond_used[atom_vecmap[bb[kk].0] as usize] = distbond_index;
                        distbond_used[atom_vecmap[ii] as usize] = distbond_index;

                        let pprefix = if distbond_index >= 10{
                            "%"
                        } else{
                            ""
                        };
                        res[ast] += &(pprefix.to_owned()+distbond_index.to_string().as_str());
                        if bb[kk].1 == "-"{
                            res[aen] += &(pprefix.to_owned()+distbond_index.to_string().as_str());
                        }else{
                            res[aen] += &(bb[kk].1.clone() + pprefix + distbond_index.to_string().as_str());
                        }
                    }
                }
            }
        }
        return res.into_iter().fold("".to_owned(),|s,m|s+&m);
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
    let mut allstart:Vec<usize> = edges.iter().map(|m|*m.0).collect();
    allstart.sort();
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
    /*
    https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00785
    Predictive Multitask Deep Neural Network Models for ADME-Tox Properties: Learning from Large Data Sets
        Jan Wenzel*Hans MatterFriedemann Schmidt
    の Supporting Info
    */
    let examples = vec![
        "CCOc1ccc(Nc2c(C)c(N[C@H]3CCCNC3)nc4ccnn24)cc1",
        "N[C@@H]1CC[C@H](CC1)Nc2cc(Nc3ccc(F)c(Cl)c3)n4nccc4n2",
        "COc1ccc(CC(=O)\\N=C(/N)\\N[C@H](CC2CCCCC2)C(=O)NCc3ccc(cc3)c4nn[nH]n4)cc1OC",
        "COc1ccccc1C(=O)\\C=C\\c2ccccc2C(F)(F)F",
        "COc1cc2ncnc(Nc3c(C)c(O)ccc3F)c2cc1OC",
        "Cn1ncc2cc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC[C@H](O)C5)ccc12",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CCCC5)cc2",
        "COC1CN(C1)C(=O)c2ccc3[nH]nc(Cc4ccc(cc4)c5cnn(C)c5)c3c2",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CCOCC5)cc2",
        "CNC(=O)c1ccc2[nH]nc(Cc3ccc(cc3)c4cnn(C)c4)c2c1",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC[C@@H](O)C5)cc2",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC(O)C5)cc2",
        "CNC(=O)c1ccc2[nH]nc(Cc3ccc4c(cnn4C)c3)c2c1",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC[C@H](O)C5)cc2",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CCNCC5)cc2",
        "Cn1ncc2cc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC(O)C5)ccc12",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CC(F)(F)C5)cc2",
        "Cn1ncc2cc(Cc3[nH]nc4ccc(cc34)C(=O)N5CCC(F)(F)C5)ccc12",
        "COC1CN(C1)C(=O)c2ccc3[nH]nc(Cc4ccc5c(cnn5C)c4)c3c2",
        "Cn1cc(cn1)c2ccc(Cc3[nH]nc4ccc(cc34)C(=O)N5CCC(F)(F)C5)cc2",
        "CCc1cc(nc(N)n1)c2c[nH]c3ncc(cc23)c4cnn(c4)C5CCNCC5",
        "CCc1cc(nc(N)n1)c2c[nH]c3ncc(cc23)c4cnn(CCN(C)C)c4",
        "CCCCc1cc(nc(N)n1)c2c[nH]c3ncc(cc23)c4cnn(c4)C5CCNCC5",
        "C[C@H](CCC(=O)Nc1ccccc1)[C@H]2CC[C@H]3[C@@H]4[C@@H](C[C@@H]5CC6(CC[C@]5(C)[C@H]4C[C@H](OC(=O)C)[C@]23C)OOC7(CCC(CC7)C(=O)Nc8ccccc8)OO6)OC(=O)C",
        "C[C@H](CCC(=O)NCCN(C)C)[C@H]1CC[C@H]2[C@@H]3[C@@H](C[C@@H]4CC5(CC[C@]4(C)[C@H]3C[C@H](OC(=O)C)[C@]12C)OOC6(CCC(CC6)C(=O)NCCN(C)C)OO5)OC(=O)C",
        "C[C@H](CCC(=O)N)[C@H]1CC[C@H]2[C@@H]3[C@@H](C[C@@H]4CC5(CC[C@]4(C)[C@H]3C[C@H](OC(=O)C)[C@]12C)OOC6(CCC(CC6)C(=O)N)OO5)OC(=O)C",
        "O=C(NC1CCN(C\\C\\2=C\\CCCCCC2)CC1)Nc3ccc4ccccc4c3",
        "COc1ccc2c(Oc3ccc(NC(=O)C4=C(C)N(CC(C)(C)O)N(C4=O)c5ccccc5)nc3)ccnc2c1",
        "OC(=O)c1cccc(Cc2cc(Cl)ccc2OCc3ccc(Cl)cc3F)n1",
        "CNS(=O)(=O)c1cc(c2c3C(=O)N(C)C(=O)N(CC4CC4)c3nn2Cc5ccnc6ccc(Cl)cc56)n(C)c1",
        "CNS(=O)(=O)c1cc(c2c3C(=O)N(C)C(=O)N(CC4CC4)c3nn2Cc5c[nH]c6ccc(Cl)cc56)n(C)c1",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4ccnc5ccc(Cl)cc45)c(c6oc(cc6)S(=O)(=O)C)c3C1=O",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4cn(C)c5ccc(Cl)cc45)c(c3C1=O)c6cc(cn6C)S(=O)(=O)C",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4c[nH]c5ccc(Cl)cc45)c(c3C1=O)c6cc(cn6C)S(=O)(=O)C",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4ccnc5ccc(Cl)cc45)c(c6oc(cc6Cl)S(=O)(=O)C)c3C1=O",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4ccnc5ccc(Cl)cc45)c(c3C1=O)c6nc(cn6C)S(=O)(=O)C",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4c[nH]c5ccc(Cl)cc45)c(c3C1=O)c6nc(cn6C)S(=O)(=O)C",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4cn(C)c5ccc(Cl)cc45)c(c6cc(oc6C)S(=O)(=O)C)c3C1=O",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4c[nH]c5ccc(Cl)cc45)c(c3C1=O)c6nc(cn6C)S(=O)(=O)N",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4ccnc5ccc(Cl)cc45)c(c6cc(oc6C)S(=O)(=O)C)c3C1=O",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4c[nH]c5ccc(Cl)cc45)c(c6cc(cn6C)S(=O)C)c3C1=O",
        "CN1C(=O)N(CC2CC2)c3nn(Cc4ccnc5ccc(Cl)cc45)c(c3C1=O)c6cc(cn6C)S(=O)C",
        "Clc1ccc(cc1)C(NC(=O)CNC(=O)CCc2ccccc2)c3ccccc3",
        "CC(C)COc1ccc(Cl)cc1c2ccccc2c3cccc(n3)C(=O)O",
        "CCNc1cc(N2CCCCS2(=O)=O)c(F)c(c1)C(=O)N[C@@H](Cc3ccccc3)[C@H](O)CNCc4cccc(c4)C(F)(F)F",
        "CCNc1cc(N2CCCC2=O)c(F)c(c1)C(=O)N[C@@H](Cc3ccccc3)[C@H](O)CNCc4cccc(c4)C(F)(F)F",
        "CCNc1cc(cc(c1)C(=O)N[C@@H](Cc2ccccc2)[C@H](O)CNCc3cnn(CC)c3)N4CCCCS4(=O)=O",
        "CC(C)(C)NCc1cc(Nc2ccnc3cc(Cl)ccc23)ccc1F",
        "CCNCc1cc(Nc2ccnc3cc(Cl)ccc23)ccc1O",
        "CCc1c(nn(c2ccc(Cl)cc2Cl)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCC5",
        "Cc1c(nn(c2ccc(F)cc2F)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "COc1c(nn(c1c2ccc(cc2)C3CC3)c4ccc(Cl)cc4Cl)C(=O)NN5CCCCC5",
        "CCc1c(nn(c1c2ccc(cc2)C3CC3)c4ccc(Cl)cc4Cl)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccc(Cl)cc2Cl)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccccc2)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccccc2Cl)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccc(Cl)cc2F)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccc(Cl)cc2Cl)c1c3ccc(cc3)C4CCCC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c1c2ccc(cc2)C3CCC3)c4ccc(Cl)cc4Cl)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccc(F)cc2Cl)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Clc1ccc(c(Cl)c1)n2nc(cc2c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "Cc1c(nn(c2ccccc2F)c1c3ccc(cc3)C4CC4)C(=O)NN5CCCCC5",
        "COc1cnc2c(NCc3nnc4ccc(nn34)c5cc(C)ns5)ccnc2c1",
        "COc1cnc2c(NCc3nnc4ccc(nn34)c5cc(F)cc(F)c5)ccnc2c1",
        "NC(=O)c1cccc2cn(nc12)c3ccc4CCNCc4c3",
        "NC(=O)c1cccc2cn(nc12)c3ccc4CNCCc4c3",
        "NC(=O)c1cccc2cn(nc12)c3ccc(CN4CCCCC4)cc3",
        "CN1CCN(Cc2ccc(cc2)n3cc4cccc(C(=O)N)c4n3)CC1",
        "NC(=O)c1cccc2cn(nc12)c3ccc(cc3)C4CCCN4",
        "NC(=O)c1cccc2cn(nc12)c3ccccc3",
        "NC(=O)c1cccc2cn(nc12)c3ccc(cc3)C4CCCNC4",
        "CNC(C)c1ccc(cc1)n2cc3cccc(C(=O)N)c3n2",
        "NC(=O)c1cccc2cn(nc12)c3ccc(cc3)C4CCNCC4",
        "CNCc1cccc(c1)n2cc3cccc(C(=O)N)c3n2",
        "NC(=O)c1cccc2cn(nc12)c3ccc(cc3)[C@H]4CCCNC4",
        "CNC(C)(C)c1ccc(cc1)n2cc3cccc(C(=O)N)c3n2",
        "CNCc1ccc(cc1)n2cc3cccc(C(=O)N)c3n2",
        "NC(=O)c1cccc2cn(nc12)c3ccc(cc3)[C@@H]4CCCNC4",
        "CC(C)NCc1ccc(cc1)n2cc3cccc(C(=O)N)c3n2",
        "CN(C)Cc1ccc(cc1)n2cc3cccc(C(=O)N)c3n2",
        "C[C@H](Oc1ccc(cc1C(=O)N2Cc3ccc(Cl)cc3C2)S(=O)(=O)C)C(F)(F)F",
        "Oc1ccc(cc1)c2sc3cc(O)ccc3c2C(=O)c4ccc(OCCN5CCCCC5)cc4",
        "COc1ccc(OC(F)(F)F)cc1Cn2c(cc3cc(ccc23)C#N)C(=O)NCC(C)(C)CO",
        "Cc1ccc(cc1Cn2c(cc3cc(ccc23)C#N)C(=O)NCCC(C)(C)O)C(F)(F)F",
        "CC(C)N(CCO)Cc1csc(n1)c2cn(CC3CCOCC3)c4c(Cl)cccc24",
        "CC(=O)NCCOc1ccc2C(=O)c3c([nH]c4cc(ccc34)C#N)C(C)(C)c2c1",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)N5CCN(CC5)C6COC6",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)N5CCN(CC5)S(=O)(=O)C",
        "CC1(C)c2cc(OCCN3CCS(=O)(=O)CC3)ccc2C(=O)c4c1[nH]c5cc(ccc45)C#N",
        "CCN(CC)CCOc1ccc2C(=O)c3c([nH]c4cc(ccc34)C#N)C(C)(C)c2c1",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)C5CCN(CC5)C6COC6",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)N5CCCCC5",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)N5CCC(O)CC5",
        "CC(C)N1CCN(CC1)c2ccc3C(=O)c4c([nH]c5cc(ccc45)C#N)C(C)(C)c3c2",
        "CC1(C)c2cc(OCCNC(=O)N)ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N",
        "CC1(C)c2cc(ccc2C(=O)c3c1[nH]c4cc(ccc34)C#N)N5CCOCC5",
        "Cc1cc(Nc2ccc(cc2)C#C)n3ncnc3n1",
        "CN(C)c1ccc(Nc2cc(C)nc3ncnn23)cc1",
        "Cc1cc(Nc2ccc(cc2)C3CC3)n4ncnc4n1",
    ];
    let examples = vec![
        "CCOc1ccc(Nc2c(C)c(N[C@H]3CCCNC3)nc4ccnn24)cc1"];
    for ee in examples.into_iter(){
        let smirkstring = ee;
        println!("\"{}\",",smirkstring);
        let mol = StringAtomConnector::smirks_to_molecule(smirkstring);
        println!("\"{}\",",mol.to_smirks());
    }


}
//例として
//http://www.chemspider.com/Chemical-Structure.65325987.html
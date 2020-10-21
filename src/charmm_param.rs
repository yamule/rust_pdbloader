#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(dead_code,unused_imports)]
use std::collections::HashMap;
#[allow(dead_code,unused_imports)]
use std::collections::HashSet;
use regex::Regex;
#[allow(dead_code,unused_imports)]
use super::pdbdata::*;
#[allow(dead_code,unused_imports)]
use super::geometry::*;
#[allow(dead_code,unused_imports)]
use super::atom_id_map;

use super::misc_util::*;

#[allow(dead_code,unused_imports)]
use super::process_3d::*;
#[allow(dead_code,unused_imports)]
use super::debug_env;
use std::fs::File;
#[allow(dead_code,unused_imports)]
use std::f64::consts::PI;


lazy_static! {
    static ref REGEX_COMMENT:Regex = Regex::new(r"[\s]+!.*").unwrap();
    static ref REGEX_COMMENT_STRICT:Regex = Regex::new(r"!.*").unwrap();
}



#[derive(Debug)]
#[allow(non_camel_case_types)]
pub struct Param_MASS{
    pub atom_name:String,
    pub mass:f64,
    pub index:i64
}


#[derive(Debug)]
pub struct AtomResiParam{
    pub atom_name:String,
    pub atom_type:String,
    pub partial_charge:f64
}


//IC と書かれたセクションだが、データが無かったりするのであまり使われていない
#[derive(Debug)]
pub struct ICResiParam{
    pub atoms:(PosAtom,PosAtom,PosAtom,PosAtom),

    pub length01_or_02:f64,
    pub angle012_or_021:f64,

    pub improper:bool,    

    pub dehedral0123:f64,
    pub angle123:f64,
    pub length23:f64

}

//+ とか - とかついているのでそれに対応する
#[derive(Debug,Clone)]
pub struct PosAtom{
    pub atom_name:String,
    pub ex_code:String,
    pub original:String
}
impl PosAtom{
    pub fn sep_code(code:&str)->PosAtom{
        
        if let Some(x) = REGEX_POSCODE.captures(code){
            let ret =  PosAtom{original:code.to_string(),atom_name:x.get(2).unwrap().as_str().to_string(),ex_code:x.get(1).unwrap().as_str().to_string()};
            
            return ret;
        }
        return PosAtom{original:code.to_string(),atom_name:code.to_string(),ex_code:"".to_string()};
    }
    
}

//NTERM 等の処理がされた RESI
#[derive(Debug)]
#[allow(non_camel_case_types)]
pub struct Param_RESI_merged<'a>{
    pub residue_name:String,
    pub atom_params:Vec<&'a AtomResiParam>,
    pub charge:f64,
    pub ic_params:Vec<&'a ICResiParam>,
    
    pub bond:Vec<&'a (PosAtom,PosAtom,usize)>,//atom0, atom1, order
    pub impr:Vec<&'a Vec<PosAtom>>,
    pub cmap:Vec<&'a Vec<PosAtom>>,
    //pub imph:Vec<Vec<PosAtom>>,
    pub dihe:Vec<&'a Vec<PosAtom>>,//DIHEDRAL//何故あるのか謎。PARAM File から取れるようにか？
    pub donor:Vec<&'a (String,String)>,
    pub acceptor:Vec<&'a (String,String)>,
}
impl<'a> Param_RESI_merged<'a>{
    pub fn new(resi:&'a Param_RESI)->Param_RESI_merged<'a>{
        return Param_RESI_merged{
        residue_name:resi.residue_name.clone(),
        atom_params:resi.atom_params.iter().map(|m|m).collect(),
        charge:resi.charge,
        ic_params:resi.ic_params.iter().map(|m|m).collect(),
        bond:resi.bond.iter().map(|m|m).collect(),
        impr:resi.impr.iter().map(|m|m).collect(),
        cmap:resi.cmap.iter().map(|m|m).collect(),
        //pub imph:Vec<Vec<PosAtom>>,
        dihe:resi.dihe.iter().map(|m|m).collect(),
        donor:resi.donor.iter().map(|m|m).collect(),
        acceptor:resi.acceptor.iter().map(|m|m).collect(),
        };
    }
}


#[derive(Debug)]
#[allow(non_camel_case_types)]
pub struct Param_RESI{
    pub pres:bool,
    pub residue_name:String,
    pub atom_params:Vec<AtomResiParam>,
    pub charge:f64,
    pub ic_params:Vec<ICResiParam>,

    
    pub bond:Vec<(PosAtom,PosAtom,usize)>,//atom0, atom1, order
    pub impr:Vec<Vec<PosAtom>>,
    pub cmap:Vec<Vec<PosAtom>>,
    //pub imph:Vec<Vec<PosAtom>>,
    pub dihe:Vec<Vec<PosAtom>>,//DIHEDRAL//何故あるのか謎。PARAM File から取れるようにか？
    pub delete:Vec<Vec<String>>,//良く分からないので一行分全部入れる。。。
    pub donor:Vec<(String,String)>,
    pub acceptor:Vec<(String,String)>,
}
impl Param_RESI{
    pub fn new()->Param_RESI{
        return Param_RESI{
        residue_name:"UNK".to_string(),
        charge:0.0,
        atom_params:vec![],
        ic_params:vec![],
        bond:vec![],
        impr:vec![],
        //imph:vec![],
        dihe:vec![],
        cmap:vec![],
        delete:vec![],
        donor:vec![],
        acceptor:vec![],
            pres:false
        };
    }
}

#[allow(non_camel_case_types)]
pub struct Param_BOND{
    pub atoms:(String,String),
    pub kb:f64,//force constant
    pub b0:f64//equiblium geometry
}

#[allow(non_camel_case_types)]
pub struct Param_ANGLE{//theta
    pub atoms:(String,String,String),
    pub ktheta:f64,//force constant
    pub theta0:f64,//equiblium geometry
    pub kub:f64,//urey-bradley
    pub s0:f64,
    pub kubflag:bool
}


#[allow(non_camel_case_types)]
pub struct Param_DIHEDRAL{//phi
    pub atoms:(String,String,String,String),
    pub kchi:f64,//force constant
    pub n:f64,//multiplicity
    pub delta:f64//minimum geometry of the dihedral
}

#[allow(non_camel_case_types)]
pub struct Param_IMPROPER{//improper dihedrals?
    pub atoms:(String,String,String,String),
    pub kpsi:f64,//force constant
    pub n:f64,//multiplicity
    pub psi0:f64//minimum geometry of the dihedral
}
#[allow(non_camel_case_types)]
pub struct Param_CMAP{
    pub residue_name:String,
    pub before_proline:bool,
    pub atoms0:(String,String,String,String),
    pub atoms1:(String,String,String,String),
    pub phi_psi_val:Vec<(f64,f64,f64)>
}


#[allow(non_camel_case_types)]
pub struct Param_NONBONDED{
    pub option_param:Param_NONBONDED_OPTION,
    pub atom_param:Vec<Param_NONBONDED_ATOM>
}
#[derive(Clone,Debug)]
#[allow(non_camel_case_types)]
pub struct Param_NONBONDED_OPTION{
    pub cutnb:f64,
    pub ctofnb:f64,
    pub ctonnb:f64,
    pub eps:f64,
    pub e14fac:f64,
    pub wmin:f64,
}


#[derive(Clone,Debug)]
#[allow(non_camel_case_types)]
pub struct Param_NONBONDED_ATOM{
    pub atom_type:String,
    pub polarizability:f64,
    pub number_of_effective_electrons:f64,
    pub eps:f64,
    pub rmin_05:f64,
    pub tanford_kirkwood:bool,

    //polarizability にもあるのだろうか、、、
    pub eps_1_4:f64,
    pub rmin_05_1_4:f64,
    pub flag_1_4:bool
}

pub struct CHARMMParam{
    pub resi:Vec<Param_RESI>,
    pub mass:Vec<Param_MASS>,
    pub bonds:Vec<Param_BOND>,
    pub angles:Vec<Param_ANGLE>,
    pub dihedrals:Vec<Param_DIHEDRAL>,
    pub impropers:Vec<Param_IMPROPER>,
    pub cmaps:Vec<Param_CMAP>,
    pub nonbonded:Param_NONBONDED
}

impl CHARMMParam{

    pub fn load_chamm19(tophfile:&str,param19file:&str)->CHARMMParam{
        let mut topp:CHARMMParam = CHARMMParam::load_toph19_inp(tophfile);
        let mut parr:CHARMMParam = CHARMMParam::load_param19_inp(param19file);
        parr.resi.append(&mut topp.resi);
        parr.mass.append(&mut topp.mass);
        return parr;
    }
    
    pub fn load_chamm22(topall22file:&str,parall22file:&str)->CHARMMParam{
        let mut topp:CHARMMParam = CHARMMParam::load_top_all22_inp(topall22file);
        let mut parr:CHARMMParam = CHARMMParam::load_param(parall22file);
        parr.resi.append(&mut topp.resi);
        return parr;
    }
    pub fn load_param19_inp(filename:&str)->CHARMMParam{
        //param19.inp 以外未確認
        let blocks = CHARMMParam::split_file_param19_inp(filename);
        let mass:Vec<Param_MASS> = vec![];
        let mut bonds:Vec<Param_BOND> = vec![];
        let mut angles:Vec<Param_ANGLE> = vec![];
        let mut dihedrals:Vec<Param_DIHEDRAL> = vec![];
        let mut impropers:Vec<Param_IMPROPER> = vec![];
        let cmaps:Vec<Param_CMAP> = vec![];
        let mut nonbonded:Option<Param_NONBONDED> = None;
        
        for bb in blocks.iter(){
            if start_with(&bb[0],"BOND"){
                bonds.append(&mut CHARMMParam::parse_bond_block(bb));
            }else if start_with(&bb[0],"THETAS"){
                angles.append(&mut CHARMMParam::parse_angle_block(bb));
            }else if start_with(&bb[0],"PHI"){
                dihedrals.append(&mut CHARMMParam::parse_dihedral_block(bb));
            }else if start_with(&bb[0],"IMPHI"){
                impropers.append(&mut CHARMMParam::parse_improper_block(bb));
            }else if start_with(&bb[0],"NONBONDED"){
                //C%       1.65    -0.0262       2.490 1.65 -0.1 1.9 ! includes CT and CM
                //N*       1.1000    -0.2384    1.6000   ! includes N,NC2,NH1,NH2,NH3,NP,and NR
                //O*       0.8400    -0.1591    1.6000   ! includes O, OH1, OM
                //S*       0.3400    -0.0430       1.890 ! includes S and SH1E
                let mut nbent = CHARMMParam::parse_nonbonded_block(bb);
                let mut nvec:Vec<Param_NONBONDED_ATOM> = vec![];
                for nn in nbent.atom_param.iter_mut(){
                    //-----------------------
                    //------------確認の必要あり
                    //nn.rmin_05 /= 2.0;
                    //nn.rmin_05_1_4 /= 2.0;
                    if nn.atom_type == "C%"{
                        for z in vec!["CT","CM"]{
                            let mut cc = nn.clone();
                            cc.atom_type = z.to_string();
                            nvec.push(cc);
                        }
                    }else if nn.atom_type == "N*"{
                        for z in vec!["N","NC2","NH1","NH2","NH3","NP","NR"]{
                            let mut cc = nn.clone();
                            cc.atom_type = z.to_string();
                            nvec.push(cc);
                        }
                    }else if nn.atom_type == "O*"{
                        for z in vec!["O","OH1","OM"]{
                            let mut cc = nn.clone();
                            cc.atom_type = z.to_string();
                            nvec.push(cc);
                        }
                    }else if nn.atom_type == "S*"{
                        for z in vec!["S","SH1E"]{
                            let mut cc = nn.clone();
                            cc.atom_type = z.to_string();
                            nvec.push(cc);
                        }
                    }else{
                        nvec.push(nn.clone());
                    }

                }
                nbent.atom_param = nvec;
                nonbonded = Some(nbent);
            }else{
                eprintln!("Currently, block\n {} \n is not used.",bb[0]);
            }
            //NBFIX とか HBOND 
        }
        let resi:Vec<Param_RESI> = vec![];
        return CHARMMParam{
            resi,
            mass,
            bonds,
            angles,
            dihedrals,
            impropers,
            cmaps,
            nonbonded:nonbonded.unwrap_or_else(||panic!("nonbonded section was not found"))
        };
    }
    pub fn load_param(filename:&str)->CHARMMParam{
    //par_all22_prot.prm 以外未確認
        let blocks = CHARMMParam::split_file(filename);
        let mut mass:Vec<Param_MASS> = vec![];
        let mut bonds:Vec<Param_BOND> = vec![];
        let mut angles:Vec<Param_ANGLE> = vec![];
        let mut dihedrals:Vec<Param_DIHEDRAL> = vec![];
        let mut impropers:Vec<Param_IMPROPER> = vec![];
        let mut cmaps:Vec<Param_CMAP> = vec![];
        let mut nonbonded:Option<Param_NONBONDED> = None;
        
        for bb in blocks.iter(){
            if start_with(&bb[0],"ATOMS"){
                mass.append(&mut CHARMMParam::parse_mass_block(bb));
            }else if start_with(&bb[0],"BONDS"){
                bonds.append(&mut CHARMMParam::parse_bond_block(bb));
            }else if start_with(&bb[0],"ANGLES"){
                angles.append(&mut CHARMMParam::parse_angle_block(bb));
            }else if start_with(&bb[0],"DIHEDRALS"){
                dihedrals.append(&mut CHARMMParam::parse_dihedral_block(bb));
            }else if start_with(&bb[0],"IMPROPER"){
                impropers.append(&mut CHARMMParam::parse_improper_block(bb));
            }else if start_with(&bb[0],"CMAP"){
                cmaps.append(&mut CHARMMParam::parse_cmap_block(bb));
            }else if start_with(&bb[0],"NONBONDED"){
                nonbonded = Some(CHARMMParam::parse_nonbonded_block(bb));
            }else{
                eprintln!("{} is not supported yet.",&bb[0]);
            }
        }
        let resi:Vec<Param_RESI> = vec![];
        return CHARMMParam{
            resi,
            mass,
            bonds,
            angles,
            dihedrals,
            impropers,
            cmaps,
            nonbonded:nonbonded.unwrap_or_else(||panic!("nonbonded section was not found"))
        };
    }

    pub fn parse_bond_block(block:&Vec<String>)->Vec<Param_BOND>{
        let mut ret:Vec<Param_BOND> = vec![];
        for line in block.iter(){
            if start_with(line,"!"){
                continue;
            }
            
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() < 4{
                continue;
            }
            ret.push(
                Param_BOND{
                    atoms:(ptt[0].clone(),ptt[1].clone()),
                    kb:ptt[2].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                    b0:ptt[3].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                }
            );
        }
        return ret;
    }

    
    pub fn parse_angle_block(block:&Vec<String>)->Vec<Param_ANGLE>{
        let mut ret:Vec<Param_ANGLE> = vec![];
        for line in block.iter(){
            if start_with(line,"!"){
                continue;
            }
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() < 5{
                continue;
            }
            let mut ang:Param_ANGLE = Param_ANGLE{
                atoms:(ptt[0].clone(),ptt[1].clone(),ptt[2].clone()),
                ktheta:ptt[3].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                theta0:ptt[4].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                kub:0.0,
                s0:0.0,
                kubflag:false
            };
            if ptt.len() > 5{
                ang.kub = ptt[5].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str());
                ang.s0 = ptt[6].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str());
                ang.kubflag = true;
            }
            ret.push(ang);
            if ptt.len() == 0{
                continue;
            }
        }
        return ret;
    }

    pub fn parse_dihedral_block(block:&Vec<String>)->Vec<Param_DIHEDRAL>{
        let mut ret:Vec<Param_DIHEDRAL> = vec![];
        for line in block.iter(){
            if start_with(line,"!"){
                continue;
            }
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() < 4{//なんか変なデータがあったらエラーになるように
                continue;
            }
            let dihed:Param_DIHEDRAL = Param_DIHEDRAL{
                atoms:(ptt[0].clone(),ptt[1].clone(),ptt[2].clone(),ptt[3].clone()),
                kchi:ptt[4].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                n:ptt[5].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                delta:ptt[6].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
            };
            ret.push(dihed);
        }
        return ret;
    }
    pub fn parse_nonbonded_block(block:&Vec<String>)->Param_NONBONDED{
        let mut ret:Vec<Param_NONBONDED_ATOM> = vec![];
        
        let mut cutnb:Option<f64> = None;
        let mut ctofnb:Option<f64> = None;
        let mut ctonnb:Option<f64> = None;
        let mut eps:Option<f64> = None;
        let mut e14fac:Option<f64> = None;
        let mut wmin:Option<f64> = None;

        for (_lii,line) in block.iter().enumerate(){
            if start_with(line,"!"){
                continue;
            }
            if start_with(line,"NONBONDED"){
                continue;
            }
            if start_with(line,"cutnb") || start_with(line,"     CUTNB"){
                //cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 
                let mut ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
                if ptt[0] == ""{
                    ptt.remove(0);
                }
                for pii in (0..ptt.len()).step_by(2){
                    if ptt[pii].eq_ignore_ascii_case("cutnb"){
                        cutnb = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                    if ptt[pii].eq_ignore_ascii_case("ctofnb"){
                        ctofnb = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                    if ptt[pii].eq_ignore_ascii_case("ctonnb"){
                        ctonnb = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                    if ptt[pii].eq_ignore_ascii_case("eps"){
                        eps = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                    if ptt[pii].eq_ignore_ascii_case("e14fac"){
                        e14fac = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                    if ptt[pii].eq_ignore_ascii_case("wmin"){
                        wmin = Some(ptt[pii+1].parse::<f64>().unwrap_or_else(|_|panic!("Couldn't parse {} ",ptt[pii+1])));
                    }
                }        
                continue;
            }
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() < 3{
                continue;
            }
            let aname:String = ptt[0].clone();
            let pola = ptt[1].parse::<f64>().expect(format!("error in parse {}",line).as_str());
            let welld = ptt[2].parse::<f64>().expect(format!("error in parse {}",line).as_str());
            let rmin_05 = ptt[3].parse::<f64>().expect(format!("error in parse {}",line).as_str());
            if welld < 0.0{
 
                let mut eps_1_4 = 0.0;
                let mut rmin_05_1_4 = 0.0;
                let mut flag1_4:bool = false;
                if ptt.len() > 4{
                    eps_1_4 = ptt[5].parse::<f64>().expect(format!("error in parse {}",line).as_str());
                    rmin_05_1_4 = ptt[6].parse::<f64>().expect(format!("error in parse {}",line).as_str());
                    flag1_4 = true;
                }
                let nb:Param_NONBONDED_ATOM = Param_NONBONDED_ATOM{
                    atom_type:aname,
                    polarizability:0.0,
                    number_of_effective_electrons:0.0,
                    eps:welld,
                    rmin_05:rmin_05,
                    tanford_kirkwood:false,
                    eps_1_4:eps_1_4,
                    rmin_05_1_4:rmin_05_1_4,
                    flag_1_4:flag1_4
                };
                ret.push(nb);
            }else{
                let nb:Param_NONBONDED_ATOM = Param_NONBONDED_ATOM{
                    atom_type:aname,
                    polarizability:pola,
                    number_of_effective_electrons:welld,
                    eps:0.0,
                    rmin_05:0.0,
                    tanford_kirkwood:true,
                    eps_1_4:0.0,
                    rmin_05_1_4:0.0,
                    flag_1_4:false
                };
                ret.push(nb);
            }
            if ptt.len() == 0{
                continue;
            }
        }

        return Param_NONBONDED{
            option_param:Param_NONBONDED_OPTION{
                cutnb:cutnb.unwrap_or_else(||panic!("cutnb was not defined")),
                ctofnb:ctofnb.unwrap_or_else(||panic!("ctofnb was not defined")),
                ctonnb:ctonnb.unwrap_or_else(||panic!("ctonnb was not defined")),
                eps:eps.unwrap_or_else(||panic!("eps was not defined")),
                e14fac:e14fac.unwrap_or_else(||panic!("e14fac was not defined")),
                wmin:wmin.unwrap_or_else(||panic!("wmin was not defined")),
            },
            atom_param:ret
        };
    }

    pub fn parse_cmap_block(block:&Vec<String>)->Vec<Param_CMAP>{
        let mut ret:Vec<Param_CMAP> = vec![];
        let atomline_regex:Regex = Regex::new(r"([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)").unwrap();
        let mut resname:String = "UNK".to_string();
        let mut prolineflag:bool = false;
        let mut buff:Vec<f64> = vec![];
        let mut atoms0:(String,String,String,String) = ("".to_string(),"".to_string(),"".to_string(),"".to_string());
        let mut atoms1:(String,String,String,String) = ("".to_string(),"".to_string(),"".to_string(),"".to_string());
        let mut varnum:usize = 0;
        let mut varnum_2:usize = 0;
        for line in block.iter(){
            
            if start_with(line,"CMAP"){
                continue;
            }
            if start_with(line,"!"){
                if start_with(line,"!2 adjacent prolines"){
                    resname = "PRO".to_string();
                    prolineflag = true;
                }else if start_with(line,"! alanine map"){
                    resname = "ALA".to_string();
                    prolineflag = false;
                }else if start_with(line,"! alanine before proline map"){
                    resname = "ALA".to_string();
                    prolineflag = true;
                }else if start_with(line,"! proline"){
                    assert_eq!(line,"! proline");
                    resname = "PRO".to_string();
                    prolineflag = false;
                }else if start_with(line,"!2 adjacent prolines"){
                    resname = "PRO".to_string();
                    prolineflag = true;
                }else if start_with(line,"! glycine map"){
                    resname = "GLY".to_string();
                    prolineflag = false;
                }else if start_with(line,"! glycine before proline map"){
                    resname = "GLY".to_string();
                    prolineflag = true;
                }
                continue;
            }
            if let Some(x) = atomline_regex.captures(line.as_str()){
                if buff.len() > 0{
                    panic!("CMAP format error!");
                }
                let mut vatt:Vec<String> = vec![];
                for ii in 0..8{
                    vatt.push(x.get(ii+1).unwrap().as_str().to_string());
                }
                varnum = x.get(9).unwrap().as_str().parse::<usize>().expect(format!("parse error {}",line).as_str());
                atoms0 = (vatt[0].clone(),vatt[1].clone(),vatt[2].clone(),vatt[3].clone());
                atoms1 = (vatt[4].clone(),vatt[5].clone(),vatt[6].clone(),vatt[7].clone());
                varnum_2 = varnum*varnum;
                continue;
            }

            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            for pp in ptt.iter(){
                if pp.len() > 0{
                    buff.push(pp.parse::<f64>().expect(format!("f64 parse error {}",line).as_str()));
                } 
            }
            //フォーマットがめちゃ処理しづらい・・・
            if buff.len() == varnum_2{
                if resname.len() == 0{
                    panic!("format error?? {:?} ",block);
                }
                let step:f64 = 360.0/varnum as f64;
                let mut phi_psi_val:Vec<(f64,f64,f64)> = vec![];
                for (bii,bb) in buff.iter().enumerate(){
                    let phh:f64 = -180.0+step*((bii/varnum) as f64);
                    let pss:f64 = -180.0+step*((bii%varnum) as f64);
                    phi_psi_val.push((phh,pss,*bb));
                }
                let rr = Param_CMAP{
                    residue_name:resname,
                    before_proline:prolineflag,
                    atoms0:atoms0.clone(),
                    atoms1:atoms1.clone(),
                    phi_psi_val:phi_psi_val
                };
                ret.push(rr);
                
                buff = vec![];
                resname = "".to_string();
                atoms0 = ("".to_string(),"".to_string(),"".to_string(),"".to_string());
                atoms1 = ("".to_string(),"".to_string(),"".to_string(),"".to_string());
            }
        }
        if buff.len() > 0{
            panic!("???CMAP format error? {:?}",block);
        }
        return ret;
    }

    pub fn parse_improper_block(block:&Vec<String>)->Vec<Param_IMPROPER>{
        let mut ret:Vec<Param_IMPROPER> = vec![];
        for line in block.iter(){
            if start_with(line,"!"){
                continue;
            }
            let ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() < 4{//なんか変なデータがあったらエラーになるように
                continue;
            }
            let impr:Param_IMPROPER = Param_IMPROPER{
                atoms:(ptt[0].clone(),ptt[1].clone(),ptt[2].clone(),ptt[3].clone()),
                kpsi:ptt[4].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                n:ptt[5].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
                psi0:ptt[5].as_str().parse::<f64>().expect(format!("f64 parse error {} ",line).as_str()),
            };
            ret.push(impr);
            if ptt.len() == 0{
                continue;
            }
        }
        return ret;
    }

    
    
    pub fn split_file(filename:&str)->Vec<Vec<String>>{
        //par_all22_propt_prm.inp 以外未確認
        let mut ret:Vec<Vec<String>> = vec![];
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut buff:Vec<String> = vec![];
        let labels:Vec<&str> = vec![
            "ATOMS","BONDS","ANGLES","DIHEDRALS","IMPROPER","CMAP","NBFIX"
        ];
        for (_lcount,line_) in reader.lines().enumerate() {
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT.replace_all(&line_, "")).to_string();

            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }
            for ll in labels.iter(){
                if line == *ll{
                    ret.push(buff);
                    buff = vec![];
                    break;
                }
            }

            if start_with(&line,"NONBONDED "){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"HBOND "){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"END"){
                ret.push(buff);
                buff = vec![];
            }
            buff.push(line);
        }
        return ret;

    }
    
    pub fn split_file_param19_inp(filename:&str)->Vec<Vec<String>>{
        //param19.inp 以外未確認
        let mut ret:Vec<Vec<String>> = vec![];
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut buff:Vec<String> = vec![];
        let labels:Vec<&str> = vec![
            "BOND","THETAS","PHI","IMPHI","NBFIX"
        ];
        for (_lcount,line_) in reader.lines().enumerate() {
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT_STRICT.replace_all(&line_, "")).to_string();

            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }
            for ll in labels.iter(){
                if line == *ll{
                    ret.push(buff);
                    buff = vec![];
                    break;
                }
            }

            if start_with(&line,"NONBONDED  NBXMOD 5"){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"HBOND AEXP 4 REXP 6"){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"END"){
                ret.push(buff);
                buff = vec![];
            }
            buff.push(line);
        }
        return ret;
    }

    pub fn load_top_all22_inp(filename:&str)->CHARMMParam{
        let blocks = CHARMMParam::split_file_toph19_inp(filename);
        let mut resi:Vec<Param_RESI> = vec![];
        let mut mass:Vec<Param_MASS> = vec![];
        for bb in blocks.iter(){
            if start_with(&bb[0],"RESI ") || start_with(&bb[0],"PRES "){
                let ress:Param_RESI = CHARMMParam::parse_resi_block(bb);
                resi.push(ress);
            }else if start_with(&bb[0],"MASS "){
                mass.append(&mut CHARMMParam::parse_mass_block(bb));
            }
        }
        let bonds:Vec<Param_BOND> = vec![];
        let angles:Vec<Param_ANGLE> = vec![];
        let dihedrals:Vec<Param_DIHEDRAL> = vec![];
        let impropers:Vec<Param_IMPROPER> = vec![];
        let cmaps:Vec<Param_CMAP> = vec![];
        let nonbonded:Param_NONBONDED = Param_NONBONDED{
            option_param:Param_NONBONDED_OPTION{ 
                cutnb:0.0,
                ctofnb:0.0,
                ctonnb:0.0,
                eps:0.0,
                e14fac:0.0,
                wmin:0.0
            },
            atom_param:vec![]
        };
        return CHARMMParam{
            resi,
            mass,
            bonds,
            angles,
            dihedrals,
            impropers,
            cmaps,
            nonbonded
        };

    }

    pub fn load_toph19_inp(filename:&str)->CHARMMParam{
        let blocks = CHARMMParam::split_file_toph19_inp(filename);
        let mut resi:Vec<Param_RESI> = vec![];
        let mut mass:Vec<Param_MASS> = vec![];
        for bb in blocks.iter(){
            if start_with(&bb[0],"RESI ") || start_with(&bb[0],"PRES "){
                let ress:Param_RESI = CHARMMParam::parse_resi_block_19(bb);
                resi.push(ress);
            }else if start_with(&bb[0],"MASS "){
                mass.append(&mut CHARMMParam::parse_mass_block(bb));
            }
        }
        let bonds:Vec<Param_BOND> = vec![];
        let angles:Vec<Param_ANGLE> = vec![];
        let dihedrals:Vec<Param_DIHEDRAL> = vec![];
        let impropers:Vec<Param_IMPROPER> = vec![];
        let cmaps:Vec<Param_CMAP> = vec![];
        let nonbonded:Param_NONBONDED = Param_NONBONDED{
            option_param:Param_NONBONDED_OPTION{ 
                cutnb:0.0,
                ctofnb:0.0,
                ctonnb:0.0,
                eps:0.0,
                e14fac:0.0,
                wmin:0.0
            },
            atom_param:vec![]
        };
        return CHARMMParam{
            resi,
            mass,
            bonds,
            angles,
            dihedrals,
            impropers,
            cmaps,
            nonbonded
        };

    }
    pub fn split_file_toph19_inp(filename:&str)->Vec<Vec<String>>{
        //toph19.inp 以外未確認
        let mut ret:Vec<Vec<String>> = vec![];
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut buff:Vec<String> = vec![];
        let mut massline:Vec<String> = vec![];
        for (_lcount,line_) in reader.lines().enumerate() {
            let line_ = line_.unwrap();
            let line_ =  (*REGEX_TAILBLANK.replace_all(&line_, "")).to_string();
            let line =  (*REGEX_COMMENT.replace_all(&line_, "")).to_string();

            if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                continue;
            }
            if start_with(&line,"MASS "){
                massline.push(line);
                continue;
            }
            if start_with(&line,"RESI "){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"PRES "){
                ret.push(buff);
                buff = vec![];
            }else if start_with(&line,"END"){
                ret.push(buff);
                buff = vec![];
            }
            buff.push(line);
        }
        ret.push(massline);
        return ret;
    }


    
	pub fn parse_resi_block(block:&Vec<String>)->Param_RESI{
        let mut ret:Param_RESI = Param_RESI::new();
        for line in block.iter(){
            let mut ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() == 0{
                continue;
            }
            let head:String = ptt.remove(0);

            if head == "RESI" || head == "PRES"{
                if head == "PRES"{
                    ret.pres = true;
                }
                let fto:f64 = ptt[1].parse::<f64>().expect(format!("Error in parsing resi line. {} ",line).as_str());
                ret.charge = fto;
                ret.residue_name = ptt[0].clone();
            }else if head == "ATOM"{
                let att = AtomResiParam{
                    atom_name:ptt[0].to_string(),
                    atom_type:ptt[1].to_string(),
                    partial_charge:ptt[2].parse::<f64>().expect("Error in parsing atom line."),
                };
                ret.atom_params.push(att);
            }else if head == "DELETE"{
                ret.delete.push(ptt);
            }else if head == "CMAP"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.cmap.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "IMPR"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.impr.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "IMPH"{//IMPR と同じっぽい
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.impr.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "DIHE"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.dihe.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "BOND"{
                for ii in (0..ptt.len()).step_by(2){
                    ret.bond.push((PosAtom::sep_code(&ptt[ii]),PosAtom::sep_code(&ptt[ii+1]),1));
                }
            }else if head == "DOUBLE"{
                for ii in (0..ptt.len()).step_by(2){
                    ret.bond.push((PosAtom::sep_code(&ptt[ii]),PosAtom::sep_code(&ptt[ii+1]),2));
                }
                
            }else if head == "TRIPLE"{//ない?
                for ii in (0..ptt.len()).step_by(2){
                    let a0 = ptt[ii].to_string();
                    let a1 = ptt[ii+1].to_string();
                    for bb in ret.bond.iter_mut(){
                        if (bb.0.original == a0 && bb.1.original == a1) || (bb.1.original == a0 && bb.0.original == a1){
                            bb.2 = 3
                        }
                    }
                }
            }else if head == "DONOR"{
                let a0 = ptt[0].to_string();
                let a1 = ptt[1].to_string();
                ret.donor.push((a0,a1));
                assert!(ptt.len() == 2);
            }else if head == "ACCEPTOR"{
                if ptt.len() <= 1{
                    let a0 = ptt[0].to_string();
                    ret.acceptor.push((a0.clone(),a0));
                }else{
                    let a0 = ptt[0].to_string();
                    let a1 = ptt[1].to_string();
                    ret.acceptor.push((a0,a1));
                    assert!(ptt.len() == 2);
                }

            }else if head == "IC"{
                let a0 = ptt[0].to_string();
                let a1 = ptt[1].to_string();
                let mut a2 = ptt[2].to_string();
                let a3 = ptt[3].to_string();
                let mut impro:bool = false;
                if start_with(&a2,"*"){
                    impro = true;
                    a2 = a2[1..].to_string();
                }
                let icparam = ICResiParam{
                atoms:(PosAtom::sep_code(&a0)
                ,PosAtom::sep_code(&a1)
                ,PosAtom::sep_code(&a2)
                ,PosAtom::sep_code(&a3)),
                length01_or_02:ptt[4].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                angle012_or_021:ptt[5].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                improper:impro,    
                dehedral0123:ptt[6].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                angle123:ptt[7].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                length23:ptt[8].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str())
                };
                ret.ic_params.push(icparam);
            }
        }
		return ret;
    }
    pub fn parse_resi_block_19(block:&Vec<String>)->Param_RESI{
        let mut ret:Param_RESI = Param_RESI::new();
        for line in block.iter(){
            let mut ptt:Vec<String> = line.split_whitespace().map(|m| m.to_string()).collect();
            if ptt.len() == 0{
                continue;
            }
            let head:String = ptt.remove(0);

            if head == "RESI" || head == "PRES"{
                if head == "PRES"{
                    ret.pres = true;
                }
                let fto:f64 = ptt[1].parse::<f64>().expect(format!("Error in parsing resi line. {} ",line).as_str());
                ret.charge = fto;
                ret.residue_name = ptt[0].clone();
            }else if head == "ATOM"{
                let att = AtomResiParam{
                    atom_name:ptt[0].to_string(),
                    atom_type:ptt[1].to_string(),
                    partial_charge:ptt[2].parse::<f64>().expect("Error in parsing atom line."),
                };
                ret.atom_params.push(att);
            }else if head == "DELETE"{
                if ptt[0] == "ATOM"{
                    ptt.remove(0);
                }
                ret.delete.push(ptt);
            }else if head == "CMAP"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.cmap.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "IMPR"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.impr.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "IMPH"{//IMPR と同じっぽい
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.impr.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "DIHE"{
                let pp:Vec<PosAtom> = ptt.iter().map(|m| PosAtom::sep_code(m)).collect();
                let plen:usize = pp.len();
                for pii in (0..plen).step_by(4){
                    ret.dihe.push(vec![pp[pii+0].clone(),pp[pii+1].clone(),pp[pii+2].clone(),pp[pii+3].clone()]);
                }
            }else if head == "BOND"{
                for ii in (0..ptt.len()).step_by(2){
                    ret.bond.push((PosAtom::sep_code(&ptt[ii]),PosAtom::sep_code(&ptt[ii+1]),1));
                }
            }else if head == "DOUBLE"{
                for ii in (0..ptt.len()).step_by(2){
                    ret.bond.push((PosAtom::sep_code(&ptt[ii]),PosAtom::sep_code(&ptt[ii+1]),2));
                }
                
            }else if head == "TRIPLE"{//ない?
                for ii in (0..ptt.len()).step_by(2){
                    let a0 = ptt[ii].to_string();
                    let a1 = ptt[ii+1].to_string();
                    for bb in ret.bond.iter_mut(){
                        if (bb.0.original == a0 && bb.1.original == a1) || (bb.1.original == a0 && bb.0.original == a1){
                            bb.2 = 3
                        }
                    }
                }
            }else if head == "DONO"{
                let a0 = ptt[0].to_string();
                let a1 = ptt[1].to_string();
                ret.donor.push((a0,a1));
                assert!(ptt.len() == 2);
            }else if head == "ACCE"{
                if ptt.len() <= 1{
                    let a0 = ptt[0].to_string();
                    ret.acceptor.push((a0.clone(),a0));
                }else{
                    let a0 = ptt[0].to_string();
                    let a1 = ptt[1].to_string();
                    ret.acceptor.push((a0,a1));
                    assert!(ptt.len() == 2);
                }

            }else if head == "IC"{
                let a0 = ptt[0].to_string();
                let a1 = ptt[1].to_string();
                let mut a2 = ptt[2].to_string();
                let a3 = ptt[3].to_string();
                let mut impro:bool = false;
                if start_with(&a2,"*"){
                    impro = true;
                    a2 = a2[1..].to_string();
                }
                let icparam = ICResiParam{
                atoms:(PosAtom::sep_code(&a0)
                ,PosAtom::sep_code(&a1)
                ,PosAtom::sep_code(&a2)
                ,PosAtom::sep_code(&a3)),
                length01_or_02:ptt[4].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                angle012_or_021:ptt[5].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                improper:impro,    
                dehedral0123:ptt[6].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                angle123:ptt[7].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str()),
                length23:ptt[8].parse::<f64>().expect(format!("Erro in IC section. {} \n",line).as_str())
                };
                ret.ic_params.push(icparam);
            }else{
                if head == "GROU" || head == "GROUP" || head == "PATC"{//今の所使わない
                }else{
                    if head == "!"{
                    }else{
                        //ANGLE, DONOR, ACCEPTOR
                        eprintln!("{} in {} is ignored.",head,ret.residue_name);
                    }
                }
            }
        }

        for aa in ret.acceptor.iter_mut(){
            //Acceptor base が指定されてない時は BOND で一番最初に出てきたものを採用する
            if aa.0 == aa.1{
                for bb in ret.bond.iter(){
                    if aa.0 == bb.0.atom_name{
                        aa.1 = bb.1.atom_name.clone();
                        break;
                    }else if aa.0 == bb.1.atom_name{
                        aa.1 = bb.0.atom_name.clone();
                        break;
                    }
                }
                if aa.0 == aa.1{
                    panic!("Can not assign acceptor base atom. {} ",ret.residue_name);
                }
            }
        }

		return ret;
    }

    
	pub fn parse_mass_block(block:&Vec<String>)->Vec<Param_MASS>{
        let mut ret:Vec<Param_MASS> = vec![];
        let regg:Regex = Regex::new(r"MASS[\s]+([^\s]+)[\s]+([^\s]+)[\s]+([^\s]+)").unwrap();
        for line in block.iter(){
            if let Some(x) = regg.captures(line.as_str()){
                let ii:i64 = x.get(1).unwrap().as_str().parse::<i64>().expect("Error in MASS index.");
                let namm:String = x.get(2).unwrap().as_str().to_string();
                let vall:f64 = x.get(3).unwrap().as_str().parse::<f64>().expect("Error in MASS value.");
                ret.push(Param_MASS{
                    atom_name:namm,mass:vall,index:ii
                });
            }else{
                if let Some(_) = REGEX_NOLINE.captures(line.as_str()){
                    continue;
                }
                if start_with(line,"ATOMS"){

                }else{
                    eprintln!("{} was not parsed.",line);
                }
            }
        }
		return ret;
    }
    /*
	pub fn load_top_all22_prot_rtf(filename:&str){
        

    }
    */
    
    /*
    pub fn get_prob(&self)->f64{
		return self.prob;
	}
	*/
	
}





#[test]
fn loadparamtest(){
    //"D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\toph19.inp"
    //"D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\top_all22_prot.rtf"
    //CHARMMParam::load_toph19_inp("D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\toph19.inp");
    let topp = CHARMMParam::load_top_all22_inp((debug_env::CHARMM_DIR.to_string()+"\\top_all22_prot.rtf").as_str());
    let mut parr = CHARMMParam::load_param((debug_env::CHARMM_DIR.to_string()+"\\par_all22_prot.prm").as_str());
    parr.resi = topp.resi;
}


/*
登録されている Residue を確認する用
#[test]
fn loadparam19test(){
    //"D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\toph19.inp"
    //"D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\top_all22_prot.rtf"
    //CHARMMParam::load_toph19_inp("D:\\dummy\\vscode_projects\\rust\\rust_pdbloader\\resources\\toppar_c36_jul18\\toppar\\toph19.inp");
    let mut parr = CHARMMParam::load_chamm19((debug_env::CHARMM_DIR.to_string()+"\\toph19.inp").as_str(),(debug_env::CHARMM_DIR.to_string()+"\\param19.inp").as_str());
    println!("Loaded charmm19===");
    for rr in parr.resi.iter(){
        for aa in rr.atom_params.iter(){
            println!("{}\t{}\t{}",rr.residue_name,aa.atom_name,aa.atom_type);
        }
    }
}

*/
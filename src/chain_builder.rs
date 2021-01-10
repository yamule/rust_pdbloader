#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
use std::collections::HashMap;
#[allow(dead_code,unused_imports)]
use regex::Regex;
use super::pdbdata::*;
#[allow(dead_code,unused_imports)]
use super::pdbdata;
use super::geometry::*;
#[allow(dead_code,unused_imports)]
use super::atom_id_map;
#[allow(dead_code,unused_imports)]
use super::residue_id_map;

#[allow(dead_code,unused_imports)]
use super::debug_env;

#[allow(dead_code,unused_imports)]
use super::process_3d::*;
#[allow(dead_code,unused_imports)]
use std::fs::File;
#[allow(dead_code,unused_imports)]
use std::f64::consts::PI;
use super::backbone_sample;
use super::side_chain_sample;
use super::sequence_alignment;

#[allow(dead_code,unused_imports)]
use super::misc_util::*;

use std::sync::Mutex;


#[allow(dead_code,unused_imports,non_camel_case_types)]
pub enum ATOM_STATE{
    Full,
    Heavy,
    Backbone,
    CA
}


lazy_static! {
    static ref AA_3_TO_1:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
    static ref AA_1_TO_3:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
}

pub fn prepare_static(){
    if AA_3_TO_1.lock().unwrap().len() != 0{
        return;
    }
    if AA_1_TO_3.lock().unwrap().len() != 0{
        return;
    }
    AA_3_TO_1.lock().unwrap().insert("ALA".to_string(),"A".to_string());
    AA_3_TO_1.lock().unwrap().insert("ARG".to_string(),"R".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASN".to_string(),"N".to_string());
    AA_3_TO_1.lock().unwrap().insert("ASP".to_string(),"D".to_string());
    AA_3_TO_1.lock().unwrap().insert("CYS".to_string(),"C".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLN".to_string(),"Q".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLU".to_string(),"E".to_string());
    AA_3_TO_1.lock().unwrap().insert("GLY".to_string(),"G".to_string());
    
    AA_3_TO_1.lock().unwrap().insert("HIS".to_string(),"H".to_string());
    AA_3_TO_1.lock().unwrap().insert("HSD".to_string(),"H".to_string());
    AA_3_TO_1.lock().unwrap().insert("HSC".to_string(),"H".to_string());

    AA_3_TO_1.lock().unwrap().insert("ILE".to_string(),"I".to_string());
    AA_3_TO_1.lock().unwrap().insert("LEU".to_string(),"L".to_string());
    AA_3_TO_1.lock().unwrap().insert("LYS".to_string(),"K".to_string());
    AA_3_TO_1.lock().unwrap().insert("MET".to_string(),"M".to_string());
    AA_3_TO_1.lock().unwrap().insert("PHE".to_string(),"F".to_string());
    AA_3_TO_1.lock().unwrap().insert("PRO".to_string(),"P".to_string());
    AA_3_TO_1.lock().unwrap().insert("SER".to_string(),"S".to_string());
    AA_3_TO_1.lock().unwrap().insert("THR".to_string(),"T".to_string());
    AA_3_TO_1.lock().unwrap().insert("TRP".to_string(),"W".to_string());
    AA_3_TO_1.lock().unwrap().insert("TYR".to_string(),"Y".to_string());
    AA_3_TO_1.lock().unwrap().insert("VAL".to_string(),"V".to_string());
    AA_3_TO_1.lock().unwrap().insert("UNK".to_string(),"X".to_string());
    for a in AA_3_TO_1.lock().unwrap().iter(){
        if a.0 == "HSD" || a.0 == "HSC"{
            continue;
        }
        AA_1_TO_3.lock().unwrap().insert(a.1.to_string(),a.0.to_string());
    }
}


pub fn place_backbone_ca(target:&mut PDBComp,backbone:&backbone_sample::BackboneSample){
    let ca_samplepos = &backbone.atoms[backbone_sample::CA as usize].get_xyz();
    let ca_targetpos = target.get_CA().expect("???").get_xyz();
    let mov:(f64,f64,f64) = (
        ca_targetpos.0-ca_samplepos.0
        ,ca_targetpos.1-ca_samplepos.1
        ,ca_targetpos.2-ca_samplepos.2);
    for aa in target.iter_mut_atoms(){
        if backbone_sample::get_atom_index(&aa.atom_code) > -1{
            let pos = backbone.atoms[backbone_sample::get_atom_index(&aa.atom_code) as usize].get_xyz();
            aa.set_xyz(pos.0+mov.0,pos.1+mov.1,pos.2+mov.2);
        }
    }
}
pub fn convert_aa_3_to_1(seq:&Vec<String>)->Vec<String>{
    prepare_static();
    let mut ret:Vec<String> = vec![];
    for ss in seq.iter(){
        if ss.len() == 0{
            continue;
        }
        if ss.len() != 3{
            panic!("This function changes 3 amino acid code to 1 amino acid letter. You provided {}",ss);
        }
        ret.push(AA_3_TO_1.lock().unwrap().get(ss).unwrap_or(&("X".to_string())).clone());
    }
    return ret;
}
pub fn convert_aa_1_to_3(seq:&Vec<String>)->Vec<String>{
    prepare_static();
    let mut ret:Vec<String> = vec![];
    for ss in seq.iter(){
        if ss.len() == 0{
            continue;
        }
        if ss.len() != 1{
            panic!("This function changes 1 amino acid letter to 3 amino acid code. You provided {}",ss);
        }
        ret.push(AA_1_TO_3.lock().unwrap().get(ss).unwrap_or(&("UNK".to_string())).clone());
    }
    return ret;
}


//SideChainSample の SideChain を、Target の N-CA-C に合わせた状態の座標を、Target の SideChain にコピーする
pub fn fit_side_chain(target:&mut PDBComp,sidechain:&side_chain_sample::SideChainSample){
    let mut target_ca:Point3D = Point3D::new(0.0,0.0,0.0);
    target_ca.set_vector(target.get_CA().unwrap());
    
    let mut target_c:Point3D = Point3D::new(0.0,0.0,0.0);
    target_c.set_vector(target.get_C().unwrap());
    
    let mut target_n:Point3D = Point3D::new(0.0,0.0,0.0);
    target_n.set_vector(target.get_N().unwrap());
    
    let refpoints:HashMap<String,Point3D> = sidechain.get_point_copy();

    let ref_ca = refpoints.get("CA").unwrap();
    let ref_n = refpoints.get("N").unwrap();
    let ref_c = refpoints.get("C").unwrap();

    //ちゃんとマップされているか見るべきだろうか・・・
    //水素がある場合もありそうだが
    let mut mover:Vec<&mut dyn Vector3D> = vec![];
    for aa in target.iter_mut_atoms(){
        if refpoints.contains_key(&aa.atom_code){
            aa.set_vector(refpoints.get(&aa.atom_code).unwrap());
            mover.push(aa);
        }
    }
    fit_to_vector(ref_n,ref_ca,ref_c,
     &target_n,&target_ca,&target_c
    ,&mut mover);

}

//Target の N-CA-C を Template の N-CA-C 状態にマップする
pub fn fit_n_ca_c(target:&mut PDBComp,template:&PDBComp){
    let mut target_ca:Point3D = Point3D::new(0.0,0.0,0.0);
    target_ca.set_vector(target.get_CA().unwrap());
    
    let mut target_c:Point3D = Point3D::new(0.0,0.0,0.0);
    target_c.set_vector(target.get_C().unwrap());
    
    let mut target_n:Point3D = Point3D::new(0.0,0.0,0.0);
    target_n.set_vector(target.get_N().unwrap());
    
    let mut template_ca:Point3D = Point3D::new(0.0,0.0,0.0);
    template_ca.set_vector(template.get_CA().unwrap());
    
    let mut template_c:Point3D = Point3D::new(0.0,0.0,0.0);
    template_c.set_vector(template.get_C().unwrap());
    
    let mut template_n:Point3D = Point3D::new(0.0,0.0,0.0);
    template_n.set_vector(template.get_N().unwrap());
    
    
    let mut mover:Vec<&mut dyn Vector3D> = vec![];
    for aa in target.iter_mut_atoms(){
        mover.push(aa);
    }
    fit_to_vector(&target_n,&target_ca,&target_c,
     &template_n,&template_ca,&template_c
    ,&mut mover);
}


pub fn build_dirty_chain(aa:&Vec<String>,backbone:&backbone_sample::BackboneSet,side_chain:&side_chain_sample::SideChainSet)->Vec<PDBComp>{
    let mut ret:Vec<PDBComp> = vec![];
    for (sii,ssname) in aa.iter().enumerate(){
        let mut aacode:String = ssname.to_string();
        if !backbone.has_sample(&aacode){
            aacode = "ALA".to_string();
            eprintln!("{} was changed to ALA.",ssname);
        }
        let mut res = PDBComp::new();
        res.set_name(&aacode);
        res.set_seq_id(sii as i64);
        
        //AA に含まれる原子をとりあえず作る
        let mut atoms:HashMap<String,PDBAtom> = HashMap::new();
        let bbb = backbone.get_sample_at(&aacode,0);
        for ba in bbb.atoms.iter(){
            atoms.insert(ba.atom_code.clone(),ba.clone());
        }

        let sss = side_chain.get_sample_at(&aacode,0);
        for (_,sa) in sss.atoms.iter(){
            atoms.insert(sa.atom_code.clone(),sa.clone());
        }
        
        for (_,aa) in atoms.into_iter(){
            res.add_atom(aa);
        }
        res.get_CA_mut().as_mut().unwrap().set_xyz(5.0*sii as f64,0.0,0.0);
        place_backbone_ca(&mut res,bbb);
        fit_side_chain(&mut res,sss);
        ret.push(res);
    }
    return ret;
}

//alignment とテンプレートを与えると作成された PDBResidue と、アラインされたかどうかの Boolean 値が入った配列を返す。
//smithwaterman で並べるので template_string と template に不一致があってよい
pub fn build_from_alignment(query_string:&Vec<String>,template_string:&Vec<String>
    ,template:&Vec<&PDBComp>
    ,backbone:&backbone_sample::BackboneSet,side_chain:&side_chain_sample::SideChainSet)->Vec<(PDBComp,bool)>{
    assert_eq!(query_string.len(),template_string.len());
    prepare_static();

    let mut templatestring_to_query_map:HashMap<usize,usize> = HashMap::new();
    let mut qpos:usize = 0;
    let mut tpos:usize = 0;
    for (qq,tt) in query_string.iter().zip(template_string.iter()){
        if qq == "."{
            panic!(". is not allowed.");
        }
        if tt == "."{
            panic!(". is not allowed.");
        }
        if qq != "-" && tt != "-"{
            templatestring_to_query_map.insert(tpos,qpos);
        }
        if qq != "-"{
            qpos += 1;        
        }
        if tt != "-"{
            tpos += 1;
        }
    }

    let mut query_aa:Vec<String> = vec![];
    for qq in query_string.iter(){
        if qq == "-"{
            continue;
        }
        if AA_1_TO_3.lock().unwrap().contains_key(qq){
            query_aa.push(AA_1_TO_3.lock().unwrap().get(qq).unwrap().to_string());
        }else{
            eprintln!("{} was changed to A.",qq);
            query_aa.push("ALA".to_string());
        }
    }
    eprintln!("Building dirty chain...");
    let mut modelled_res:Vec<PDBComp> = build_dirty_chain(&query_aa, backbone, side_chain);

    let mut ss1 = sequence_alignment::SeqData::new();
    let mut ss2 = sequence_alignment::SeqData::new();
    let mut sw = sequence_alignment::SequenceAlignment::new(Box::new(sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,sequence_alignment::ALIGN_LOCAL);

    for ss in template_string.iter(){
        if ss == "-"{
            continue;
        }
        ss1.seq.push(ss.clone());
    }

    for ss in template.iter(){//PDB に登録されている残基
        ss2.seq.push(AA_3_TO_1.lock().unwrap().get(ss.get_name()).unwrap_or(&("X".to_string())).clone());
    }
    //eprintln!("{:?}",ss1.seq);
    //eprintln!("{:?}",ss2.seq);
    let res = sw.align(&ss1,&ss2,true);
    let alen:usize = res.0.len();
    let mut xpos:usize = 0;
    let mut ypos:usize = 0;
    let mut mapper:Vec<(usize,usize)> = vec![];

    let mut mappedflag:Vec<bool> = vec![false;modelled_res.len()];
    
    for ii in 0..alen{
        if res.0[ii] != "-" && res.1[ii] != "-"{
            if templatestring_to_query_map.contains_key(&xpos){
                let qqpos:usize = *templatestring_to_query_map.get(&xpos).unwrap();
                if let Some(_) = template[ypos].get_CA(){
                    if let Some(_) = template[ypos].get_N(){
                        if let Some(_) = template[ypos].get_C(){
                            mappedflag[qqpos] = true;
                            if template[ypos].get_name() != AA_1_TO_3.lock().unwrap().get(&res.0[ii]).unwrap(){
                                eprintln!("Warning: Unmatched sequence {} {}",AA_1_TO_3.lock().unwrap().get(&res.0[ii]).unwrap(),template[ypos].get_name());
                            }
                            mapper.push((qqpos,ypos));
                        }
                    }
                }
            }
        }
        if res.0[ii] != "-"{
            xpos += 1;        
        }
        if res.1[ii] != "-"{
            ypos += 1;
        }
    }
    
    //eprintln!("{:?}",mapper);
    eprintln!("Mapping on the template residues...");
    place_on_template(&mapper,&mut modelled_res,template);
    let ret:Vec<(PDBComp,bool)> = modelled_res.into_iter().zip(mappedflag).collect();

    return ret;
}
pub fn place_on_template(qery_template_map:&Vec<(usize,usize)>,query:&mut Vec<PDBComp>,template:&Vec<&PDBComp>){
    for rr in qery_template_map.iter(){
        fit_n_ca_c(&mut query[rr.0],&template[rr.1]);
    }
}


#[test]
fn chain_build_test(){
    let allaa:Vec<String> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].iter().map(|m|m.to_string()).collect();

    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);

    let mut ress:Vec<PDBComp> = build_dirty_chain(&allaa,&bset,&sset);
    let mut acount:i64 = 1;
    let mut pdbstr:Vec<String> = vec![];
    for (ii,rr) in ress.iter_mut().enumerate(){
        rr.set_seq_id((ii+1) as i64);
        for aa in rr.iter_mut_atoms(){
            aa.set_serial_number(acount);
            acount += 1;
        }
        for aa in rr.iter_atoms(){
            pdbstr.push(aa.get_pdb_atom_line_string("A",rr.get_name(),rr.get_seq_id(),""));
        }
    }
    write_to_file("test/dirty.pdb",pdbstr);
}
#[test]
fn tbm_test(){
    prepare_static();
    let allaa:Vec<String> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].iter().map(|m|m.to_string()).collect();
    let pdbb:pdbdata::PDBEntry = pdbdata::load_pdb("D:/dummy/vscode_projects/rust/rust_pdbloader/example_files/1a4w_part.pdb");

    let residues:Vec<&PDBComp> = pdbb.get_model_at(0).get_entity_at(0).get_asym_at(0).iter_comps().collect();
    let mut dummystring:Vec<String> = vec![];
    for rr in residues.iter(){
        dummystring.push(AA_3_TO_1.lock().unwrap().get(rr.get_name()).unwrap().clone());
    }
    let dummyquery:Vec<String> = vec!["K".to_string();dummystring.len()];

    eprintln!("Loading backbone samples...");
    let bset = backbone_sample::BackboneSet::new(debug_env::ROTAMER_DIR);
    for aa in allaa.iter(){
        eprintln!("{} {}",aa,bset.get_num(aa));
    }

    eprintln!("Loading sidechain samples...");
    let sset = side_chain_sample::SideChainSet::new(debug_env::ROTAMER_DIR);
    for aa in allaa.iter(){
        eprintln!("{} {}",aa,sset.get_num(aa));
    }


    let mut ress:Vec<(PDBComp,bool)> = build_from_alignment(&dummyquery,&dummystring,&residues,&bset,&sset);
    let mut acount:i64 = 1;
    let mut pdbstr:Vec<String> = vec![];
    for (ii,rr) in ress.iter_mut().enumerate(){
        rr.0.set_seq_id((ii+1) as i64);
        for aa in rr.0.iter_mut_atoms(){
            aa.set_serial_number(acount);
            acount += 1;
        }
        for aa in rr.0.iter_atoms(){
            pdbstr.push(aa.get_pdb_atom_line_string("A",rr.0.get_name(),rr.0.get_seq_id(),""));
        }
    }
    write_to_file("test/tbm_check.pdb",pdbstr);
}





#[test]
fn aa_change_test(){
    prepare_static();
    let allaa3:Vec<String> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","UNK"
    ].iter().map(|m|m.to_string()).collect();
    let allaa1:Vec<String> = vec![
            "A","R","N","D","C","Q","E","G","H","I"
            ,"L","K","M","F","P","S","T","W","Y","V","X"
    ].iter().map(|m|m.to_string()).collect();
    let c3 = convert_aa_3_to_1(&allaa3);
    let c1 = convert_aa_1_to_3(&allaa1);
    assert_eq!(allaa1,c3);
    assert_eq!(allaa3,c1);
}
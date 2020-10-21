
use std::sync::Mutex;
use std::collections::HashMap;

lazy_static! {
    static ref RRESIDUE_ID:Mutex<HashMap<String,i64>> = Mutex::new(HashMap::new());
    static ref RESIDUE_NAMES:Vec<String> = vec![
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE"
            ,"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    ].iter().map(|m|m.to_string()).collect();
}

pub fn prepare(){
    for (ii,rr) in RESIDUE_NAMES.iter().enumerate(){
        RRESIDUE_ID.lock().unwrap().insert(rr.to_string(), ii as i64);
    }
}

pub fn get_index(s:&str) -> i64 {
    if let Some(x) = RRESIDUE_ID.lock().unwrap().get(s){
        return *x;
    }else{
        return -1;
    }
}

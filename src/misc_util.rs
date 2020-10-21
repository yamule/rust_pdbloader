
use std::collections::HashMap;
use regex::Regex;
use std::io::{BufWriter,Write};
use std::fs;

lazy_static! {
    pub static ref PAT_LABEL:Regex = Regex::new(r"^###[\s]*([^\s]+)[\s]+([^\s].*)").unwrap();
    pub static ref PAT_ELEM:Regex = Regex::new(r"^([^:]+):(.+)").unwrap();

    pub static ref REGEX_WS:Regex = Regex::new(r"[\s]").unwrap();//
    pub static ref REGEX_WS_P:Regex = Regex::new(r"[\s]+").unwrap();

    pub static ref REGEX_NOLINE:Regex = Regex::new(r"^[\s]*$").unwrap();//空白行
    pub static ref REGEX_TAILBLANK:Regex = Regex::new(r"[\s]*$").unwrap();//最後に改行とか空白とかある場合の削除用
    pub static ref REGEX_POSCODE:Regex = Regex::new(r"^([+\-*])(.+)").unwrap();//前後残基
}


pub fn line_to_hash(ss:&str)-> HashMap<String,String>{
    let ppt: Vec<&str> = ss.split("\t").collect();
    let mut ret:HashMap<String,String> = HashMap::new();
    for (_,p) in ppt.iter().enumerate(){
        let mat = PAT_ELEM.captures(p);
        if let Some(x) = mat{
            ret.insert(x.get(1).unwrap().as_str().to_string(),x.get(2).unwrap().as_str().to_string());
        }else{
            
            eprintln!("{} was not parsed.",p);
            
        }
    }
    return ret;
}


pub fn start_with(target:&str,fragment:&str)-> bool{
    if let Some(x) = target.find(fragment){
        if x == 0{
            return true;
        }
    }
    return false;
}


pub fn write_to_file(filename:&str,contents:Vec<String>){
    let mut f = BufWriter::new(fs::File::create(filename).unwrap());
     for ll in contents{
        f.write_all(ll.as_bytes()).unwrap();
        f.write_all("\n".as_bytes()).unwrap();
    }
}



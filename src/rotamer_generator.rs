use std::fs;
use regex::Regex;


fn generate_intermediate_files(dirname:&str){
    let paths = fs::read_dir(dirname).unwrap();
    let exx =  Regex::new(r"").unwrap();
    for path in paths {
        
    }
}

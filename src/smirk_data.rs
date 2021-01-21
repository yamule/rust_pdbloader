

use std::sync::Mutex;
use std::collections::HashMap;


lazy_static! {
    static ref ELEMENT_NUM_TO_NAME:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
    static ref ELEMENT_NAME_TO_NUM:Mutex<HashMap<String,String>> = Mutex::new(HashMap::new());
    static ref PREPARED:Mutex<bool> = Mutex::new(false);
}

pub fn prepare(){
    //https://ja.wikipedia.org/wiki/%E5%91%A8%E6%9C%9F%E8%A1%A8
    let svec:Vec<&str> = vec!["#1","H","#2","He","#3","Li","#4","Be","#5","B","#6","C","#7","N","#8","O","#9","F","#10","Ne","#11","Na","#12","Mg","#13","Al","#14","Si","#15","P","#16","S","#17","Cl","#18","Ar","#19","K","#20","Ca","#21","Sc","#22","Ti","#23","V","#24","Cr"
    ,"#25","Mn","#26","Fe","#27","Co","#28","Ni","#29","Cu","#30","Zn","#31","Ga","#32","Ge","#33","As","#34","Se","#35","Br","#36","Kr","#55","Cs","#56","Ba","#72","Hf","#73","Ta","#74","W","#75","Re","#76","Os","#77","Ir","#78","Pt","#79","Au"
    ,"#80","Hg","#81","Tl","#82","Pb","#83","Bi","#84","Po","#85","At","#86","Rn","#87","Fr","#88","Ra","#104","Rf","#105","Db","#106","Sg","#107","Bh","#108","Hs","#109","Mt","#110","Ds","#111","Rg","#112","Cn","#113","Nh","#114","Fl","#115","Mc"
    ,"#116","Lv","#117","Ts","#118","Og","#37","Rb","#38","Sr","#39","Y","#40","Zr","#41","Nb","#42","Mo","#43","Tc","#44","Ru","#45","Rh","#46","Pd","#47","Ag","#48","Cd","#49","In","#50","Sn","#51","Sb","#52","Te","#53","I","#54","Xe","#154","Upq"
    ,"#155","Upp","#156","Uph","#157","Ups","#158","Upo","#159","Upe","#160","Uhn","#161","Uhu","#162","Uhb","#163","Uht","#164","Uhq","#165","Uhp","#166","Uhh","#167","Uhs","#168","Uho","#169","Uhe","#170","Usn","#171","Ubu","#172","Ubb","#173","Ubt"
    ,"#119","Uue","#120","Ubn","#57","La","#58","Ce","#59","Pr","#60","Nd","#61","Pm","#62","Sm","#63","Eu","#64","Gd","#65","Tb","#66","Dy","#67","Ho","#68","Er","#69","Tm","#70","Yb","#71","Lu","#89","Ac","#90","Th","#91","Pa","#92","U","#93","Np"
    ,"#94","Pu","#95","Am","#96","Cm","#97","Bk","#98","Cf","#99","Es","#100","Fm","#101","Md","#102","No","#103","Lr","#139","Ute","#140","Uqn","#141","Uqu","#142","Uqb","#143","Uqt","#144","Uqq","#145","Uqp","#146","Uqh","#147","Uqs","#148","Uqo"
    ,"#149","Uqe","#150","Upn","#151","Upu","#152","Upb","#153","Upt","#121","Ubu","#122","Ubb","#123","Ubt","#124","Ubq","#125","Ubp","#126","Ubh","#127","Ubs","#128","Ubo","#129","Ube","#130","Utn","#131","Utu","#132","Utb","#133","Utt"
    ,"#134","Utq","#135","Utp","#136","Uth","#137","Uts","#138","Uto"];
    
    let slen:usize = svec.len();
    for ii in (0..slen).step_by(2) {
        ELEMENT_NUM_TO_NAME.lock().unwrap().insert(svec[ii].to_string(), svec[ii+1].to_string());
        ELEMENT_NAME_TO_NUM.lock().unwrap().insert(svec[ii+1].to_string(), svec[ii].to_string());
    }
}

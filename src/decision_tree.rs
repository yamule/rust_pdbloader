extern crate regex;
extern crate rand;

#[allow(dead_code)]
//use std::io::prelude::*;
use std::fs::File;
use std::collections::HashSet;
use std::collections::HashMap;
use std::f64;

use std::fs;
use std::io::{BufWriter,Write,BufReader,BufRead};
use self::regex::Regex;//module ファイル内で extern crate した場合 self が必要。lib.rs で extern crate すると必要ない。
use self::rand::prelude::*;
use super::predictor::Predictor;




pub const LABEL_PREDICTOR:&str = "predictor";
const LABEL_TARGETTYPE:&str = "target_type";
const LABEL_FORMATVERSION:&str = "format_version";
const LABEL_OPTIONS:&str = "options";
const LABEL_OPTIONS_NUMCLASSES:&str = "num_classes";

const LABEL_VARIABLE:&str = "variable";
const LABEL_VARIABLE_NAME:&str = "name";
const LABEL_VARIABLE_INDEX:&str = "index";

const LABEL_NODE:&str = "node";
const LABEL_NODE_ID:&str = "id";
const LABEL_NODE_VARID:&str = "vi";
const LABEL_NODE_THRESHOLD:&str = "th";
const LABEL_NODE_CHILDLOWER:&str = "cl";
const LABEL_NODE_CHILDHIGHER:&str = "ch";
const LABEL_NODE_LEAFANSWERSTRING:&str = "ans";

const LABEL_ANSCODE_USIZE:&str = "usize";
const LABEL_ANSCODE_CLASSCOUNT:&str = "classcount";
const LABEL_ANSCODE_F64VEC:&str = "f64vec";
const LABEL_ANSCODE_F64:&str = "f64";


const COUNTCLASS_DELIMITER:&str = ",";




#[allow(non_camel_case_types)]
#[derive(Debug)]
pub enum MaxBiasFunc{
SPFUNC_ACCMIN,
SPFUNC_ACCSUM,
SPFUNC_MCC,
SPFUNC_F1_MIN,
SPFUNC_F1_SUM,
SPFUNC_PPV_MIN,
SPFUNC_PPV_SUM,
SPFUNC_FDR_MAX,
SPFUNC_FDR_SUM,
SPFUNC_TPR_MIN,
SPFUNC_TPR_SUM,
SPFUNC_PLR_MIN,
SPFUNC_PLR_SUM,
SPFUNC_NLR_MIN,
SPFUNC_NLR_SUM,
SPFUNC_DOR,
SPFUNC_INFORMEDNESS,
SPFUNC_COHEN_KAPPA,
SPFUNC_FLEISS_KAPPA,
SPFUNC_CHISQ
}

pub fn write_to_file(filename:&str,contents:Vec<String>){
    let mut f = BufWriter::new(fs::File::create(filename).unwrap());
     for ll in contents{
        f.write_all(ll.as_bytes()).unwrap();
        f.write_all("\n".as_bytes()).unwrap();
    }
}

//XX=YY<separatpr>ZZ=AA<separator>... という文字列を渡すと XX, ZZ をキーに、YY, ZZ を値にしたハッシュマップにする
pub fn line_to_hashmap(line:String,separator:& Regex)->HashMap<String,String>{
    let mut ret:HashMap<String,String> = HashMap::new();
    let v: Vec<&str> = separator.split(line.as_str()).collect();
    let key_value = Regex::new(r"^([^=]+)=(.+)$").unwrap();
    for vv in v{
        let caps = key_value.captures(vv);
        match caps{
            Some(cc)=>{
                ret.insert(
                cc.get(1).map_or("".to_string(), |m| m.as_str().to_string())
                ,cc.get(2).map_or("".to_string(), |m| m.as_str().to_string()));
            },
            None => {
                if vv.len() > 0 {
                    println!("{}\nin\n{}\nwas not parsed.",vv,line);
                }
            }
        }
    }
    return ret;
}


//ちゃんとキーがあるか調べる。
//unwrap がきちんとエラーメッセージを出さないことがあるため。
pub fn check_keys(hm:&HashMap<String,String>,k:Vec<&str>){
	let mut panicflag = false;
	for kk in k{
		if !hm.contains_key(kk){
			panicflag = true;
			eprintln!("HashMap does not have {}.",kk);
		}
	}
	if panicflag{
		panic!("Panicked in hash map check.");
	}
}


#[derive(Debug)]
pub enum SplitFunctionType{
	Gini,
	MaxBiasMCC,
	ExperimentalMaxBias,
	Entropy//現在未対応
}



pub struct DecisionTreeOptions{
	pub max_depth:Option<usize>,
	pub max_leaf_nodes:Option<usize>,
	pub min_samples_split:Option<usize>,//min_samplefraction_split が定義誰ている場合そちら優先
	pub min_samplefraction_split:Option<f64>,
	pub feature_fraction :Option<f64>,
	pub split_function_type:SplitFunctionType,
	pub random_seed:Option<u64>
}


#[derive(Debug)]
pub struct SimpleDecisionTree<'a>{
	num_classes:usize,
	pub nodes:Vec<Box<DTNode<'a>>>,
	pub var_names:Option<Vec<String>>
}

impl<'a> Predictor for SimpleDecisionTree<'a>{
	fn predict_class(&self,v:&Vec<f64>)->Vec<f64>{
		return self.predict_f64vec(v).iter().map(|m|{m.clone()}).collect();
	}
	fn predict_value(&self,_:&Vec<f64>)->f64{
		panic!("not implemented");
	}	
}

impl<'a> SimpleDecisionTree<'a>{
	pub fn push_node(&mut self,mut d:DTNode<'a>)->usize{
 		d.id = self.nodes.len();
		let rid = d.id;
		self.nodes.push(Box::new(d));
		return rid;
	}
	pub fn set_var_names(&mut self,v:&Vec<&str>){
		let mut va:Vec<String> = vec![];
		for ss in v.iter(){
			va.push(ss.to_string());
		}
		self.var_names = Some(va);
	}


	pub fn make_copy(&self,copysamples:bool)->SimpleDecisionTree{
		return SimpleDecisionTree{
			num_classes:self.num_classes.clone(),
			nodes:self.nodes.iter().map(|m|{Box::new(m.make_copy(copysamples))}).collect(),
			var_names:match self.var_names.as_ref(){
				Some(x)=>{Some(x.clone())},
				_=>{None}
			}
		};
	}
	

	//以下のように nodes からの取り出しをメソッドから行おうとするとボローチェックを通らなかった
	pub fn get_node_at(&mut self,i:usize)->&'a mut DTNode{
		return &mut self.nodes[i];
	}
	pub fn get_node_ref_at(&self,i:usize)->&'a DTNode{
		return &self.nodes[i];
	}

	pub fn build_classifier(samples:&'a mut Vec<Box<&'a SampleData_usize>>
	,opp:DecisionTreeOptions
	,vnames:Option<&Vec<&str>>)->SimpleDecisionTree<'a>{
		//枝の深さ
		let max_depth:usize = opp.max_depth.unwrap_or(samples.len());
		
		let max_leaf_nodes:usize = opp.max_leaf_nodes.unwrap_or(samples.len());
		
		let mut min_samples_split:usize = opp.min_samples_split.unwrap_or(1);
		
		if let Some(x) = opp.min_samplefraction_split{
			//min_samplefraction_split が定義されている場合そちらを使用
			min_samples_split = (x*(samples.len() as f64)).max(1.0) as usize ; 
		}


		//ノード分割の際検討される変数の割合
		let feature_fraction :f64  = opp.feature_fraction.unwrap_or(1.0);
		
		let split_function_type = opp.split_function_type;

		let mut rgen:StdRng;//乱数生成器

		match opp.random_seed{
			Some(x)=>{
				rgen = SeedableRng::seed_from_u64(x);
			},
			None => {
				rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
			}
		}

		if let SplitFunctionType::Gini = split_function_type{
		}else if let SplitFunctionType::MaxBiasMCC = split_function_type{
		}else if let SplitFunctionType::Entropy = split_function_type{
		}else{
			panic!("{:?} is not supported",split_function_type);
		}
		let mut classmax:usize = 0;//class をそのまま index に使っている
		for vv in samples.iter(){
			if vv.target_usize > classmax{
				classmax = vv.target_usize;
			}
		}
		let mut ret:SimpleDecisionTree = SimpleDecisionTree{
			num_classes:classmax+1,
			nodes:vec![],
			var_names:None
		};
		if let Some(x)=vnames{
			let nn:Vec<String> = x.iter().map(|m|{m.to_string()}).collect();
			ret.var_names = Some(nn);
		}
		let varnum:usize = samples[0].values.len();
		let rsamples:Vec<Box<&SampleData_usize>>
		 = samples.iter().map(|m|{m.clone()}).collect();
		
		let rootnode:DTNode = DTNode{
			splitting_var:None,
			xsamples:rsamples,
			child_lower:None,
			child_higher_or_eq:None,
			id:0,
			depth:0,
			dummy:false,
			answer_string:None,
			answer_usize:None,
			answer_f64:None,
			answer_f64vec:None,
		};
		let rootnodeid = ret.push_node(rootnode);
		let mut node_qued:Vec<usize> = vec![];
		node_qued.push(rootnodeid);
		
		while node_qued.len() > 0{
			//let cnodeid = node_qued.remove(0);
            //最多メンバーを持つクラス以外のメンバーが多いノードから分割する
            let mut targetnodeindex = 0;
            let mut maxpenal:f64 = 0.0;
            let mut zcou:Vec<f64> = vec![0.0;ret.num_classes];
            for (nii,nn) in node_qued.iter().enumerate(){
                let mut zsum = 0.0;
                let mut zmax:f64 = 0.0;
                let mut _zmaxindex:usize = 0;
                for zz in ret.nodes[*nn].xsamples.iter(){
                    zcou[(*zz).target_usize] += 1.0;
                    zsum += 1.0;
                    if zcou[(*zz).target_usize] > zmax{
                        zmax = zcou[(*zz).target_usize];
                        _zmaxindex = (*zz).target_usize;
                    }
                }
                if maxpenal < zsum - zmax{
                    maxpenal = zsum - zmax;
                    targetnodeindex = nii;
                }
            }
            let cnodeid = node_qued.remove(targetnodeindex);


			let original_loss = match split_function_type{
				SplitFunctionType::Gini =>{
					SampleData_usize::calc_gini(classmax,&ret.nodes[cnodeid].xsamples)
				},
				SplitFunctionType::Entropy =>{
					SampleData_usize::calc_entropy(classmax,&ret.nodes[cnodeid].xsamples)
				},
				SplitFunctionType::MaxBiasMCC =>{
					let mut samplecount_class:Vec<usize> = vec![0;classmax+1];
					for ss in ret.nodes[cnodeid].xsamples.iter(){
						samplecount_class[ss.target_usize]+=1;
					}
					let mut classcount = 0;
					for ii in 0..samplecount_class.len(){
						if samplecount_class[ii] > 0{
							classcount += 1;
						}
					}
					if classcount == 1{
						0.0
					}else{
						f64::INFINITY
					}
				},
				_=>{
					panic!("{:?} is not Supported.",split_function_type);
				}
			};
			if original_loss == 0.0{
				continue;
			}
			if ret.nodes[cnodeid].depth >= max_depth{
				continue;
			}
			if ret.nodes.len() >= max_leaf_nodes{
				break;
			}
			let mut dvec:Vec<usize> = (0..varnum).collect();
			let mut var_check:Vec<usize> = vec![];
			dvec.shuffle(&mut rgen);//連番をランダムにシャッフルする
			while dvec.len() > 0{//fraction 分のシャッフルされた連番を前からとる
				let cnum :usize = dvec.remove(0);
				if var_check.len() as f64/ varnum as f64 > feature_fraction{
					break;
				
				}
				var_check.push(cnum); 
			}
			let mut minres:Option<SplitResult> = None;
			for ii in var_check.iter(){
				
				let sres:Result<SplitResult,String> = match split_function_type{
				SplitFunctionType::Gini =>{
					SampleData_usize::calc_threshold_gini(classmax,&mut ret.nodes[cnodeid].xsamples,*ii)
				},
				SplitFunctionType::Entropy =>{
					SampleData_usize::calc_threshold_entropy(classmax,&mut ret.nodes[cnodeid].xsamples,*ii)
				},
				SplitFunctionType::MaxBiasMCC =>{
					SampleData_usize::calc_threshold_maxbias(classmax,&mut ret.nodes[cnodeid].xsamples,*ii,&MaxBiasFunc::SPFUNC_MCC)
				},
				_=>{panic!("undefined split_function_type {:?}",split_function_type);}
				};

				if let Err(_) = sres{
					continue;
				}

				let sres = sres.unwrap();
				if sres.num_lower < min_samples_split ||  sres.num_higher < min_samples_split {
					continue;
				}
				if sres.loss_weighted >= original_loss{//減少が少ない場却下してもいいかも
					continue;
				}
				if let None = minres{
					minres = Some(sres);
				}else{
					if minres.as_ref().unwrap().loss_weighted > sres.loss_weighted{
						minres = Some(sres);
					}
				}
			}
			
			match minres{
				Some(res)=>{
					ret.nodes[cnodeid].set_split(SplittingVar{
						var_index:res.varindex,
						threshold:res.threshold
					});
				
					//計算結果確認用println!("{} {} ",res.loss_higher,res.loss_lower);
					let mut lowers:Vec<Box<& SampleData_usize>> = Vec::new();
					let mut highers:Vec<Box<& SampleData_usize>> = Vec::new();
					
					//let dvec = ret.get_node_ref_at(cnodeid).xsamples.clone();
					let dvec = ret.nodes[cnodeid].xsamples.clone();
					for sxx in dvec.into_iter(){
						if res.threshold > sxx.get_value_at(res.varindex){
							lowers.push(Box::new(*sxx));
							//lowers.push(sxx.clone());
							//for kkk in samples.iter(){//ret に入れているとこうしないとボローチェック通らない。。。
							//samples 内全部回すと当然遅いので不許可
							//	if (**kkk)  as *const _ == (*sxx)  as *const _ {
							//		lowers.push(kkk.clone());
							//	}
							//}
						}else{
							highers.push(Box::new(*sxx));
						}
					}
				
					//threshold が minthreshold_start と同じ場合は通るが、変更するバージョンも入れたのでやめた
					//assert_eq!(lowers.len(),res.num_lower);
					//assert_eq!(highers.len(),res.num_higher);
					
					let lnode:DTNode = DTNode{
						splitting_var:None,
						xsamples:lowers,
						child_lower:None,
						child_higher_or_eq:None,
						id:0,
						dummy:false,
						depth:ret.nodes[cnodeid].depth+1,
						answer_string:None,
						answer_usize:None,
						answer_f64:None,
						answer_f64vec:None,

					};

					let hnode:DTNode = DTNode{
						splitting_var:None,
						xsamples:highers,
						child_lower:None,
						child_higher_or_eq:None,
						id:0,
						depth:ret.nodes[cnodeid].depth+1,
						dummy:false,
						answer_string:None,
						answer_usize:None,
						answer_f64:None,
						answer_f64vec:None,
					};
					let lid = ret.push_node(lnode);
					let hid = ret.push_node(hnode);
					ret.nodes[cnodeid].child_lower = Some(lid);
					ret.nodes[cnodeid].child_higher_or_eq = Some(hid);
					node_qued.push(lid);
					node_qued.push(hid);
				
				},
				_=>{
				}
			}
		}


		for nn in ret.nodes.iter_mut(){
			nn.set_answer_string((LABEL_ANSCODE_CLASSCOUNT.to_string()+"="+&nn.get_answer_count_string()).as_str());
		}

		return ret;
	}
	pub fn build_experimental_maxbias_classifier(samples:&'a mut Vec<Box<&'a SampleData_usize>>
	,opp:DecisionTreeOptions
	,vnames:Option<&Vec<&str>>
	,spfunc:&MaxBiasFunc)->SimpleDecisionTree<'a>{
		//枝の深さ
		let max_depth:usize = opp.max_depth.unwrap_or(samples.len());
		
		//後入れ先出し処理なのでまんべんなく分割されるはずだが不明。
		let max_leaf_nodes:usize = opp.max_leaf_nodes.unwrap_or(samples.len());
		
		let mut min_samples_split:usize = opp.min_samples_split.unwrap_or(1);
		
		if let Some(x) = opp.min_samplefraction_split{
			//min_samplefraction_split が定義されている場合そちらを使用
			min_samples_split = (x*(samples.len() as f64)).max(1.0) as usize ; 
		}


		//ノード分割の際検討される変数の割合
		let feature_fraction :f64  = opp.feature_fraction.unwrap_or(1.0);
		
		//Max Bias 用 以外対応してない
		let split_function_type = SplitFunctionType::ExperimentalMaxBias;
		if let SplitFunctionType::ExperimentalMaxBias = opp.split_function_type{

		}else{
			//Experimental 用なので
			panic!("using wrong function.");
		}
		let mut rgen:StdRng;//乱数生成器

		match opp.random_seed{
			Some(x)=>{
				rgen = SeedableRng::seed_from_u64(x);
			},
			None => {
				rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
			}
		}

		let mut classmax:usize = 0;//class をそのまま index に使っている
		for vv in samples.iter(){
			if vv.target_usize > classmax{
				classmax = vv.target_usize;
			}
		}
		let mut ret:SimpleDecisionTree = SimpleDecisionTree{
			num_classes:classmax+1,
			nodes:vec![],
			var_names:None
		};
		if let Some(x)=vnames{
			let nn:Vec<String> = x.iter().map(|m|{m.to_string()}).collect();
			ret.var_names = Some(nn);
		}
		let varnum:usize = samples[0].values.len();
		let rsamples:Vec<Box<&SampleData_usize>>
		 = samples.iter().map(|m|{m.clone()}).collect();
		
		let rootnode:DTNode = DTNode{
			splitting_var:None,
			xsamples:rsamples,
			child_lower:None,
			child_higher_or_eq:None,
			id:0,
			depth:0,
			dummy:false,
			answer_string:None,
			answer_usize:None,
			answer_f64:None,
			answer_f64vec:None,
		};
		let rootnodeid = ret.push_node(rootnode);
		let mut node_qued:Vec<usize> = vec![];
		node_qued.push(rootnodeid);
		
		while node_qued.len() > 0{
			let cnodeid = node_qued.remove(0);
			let original_loss = match split_function_type{
				
				SplitFunctionType::ExperimentalMaxBias =>{
					let mut samplecount_class:Vec<usize> = vec![0;classmax+1];
					for ss in ret.nodes[cnodeid].xsamples.iter(){
						samplecount_class[ss.target_usize]+=1;
					}
					let mut classcount = 0;
					for ii in 0..samplecount_class.len(){
						if samplecount_class[ii] > 0{
							classcount += 1;
						}
					}

					if classcount == 1{
						0.0
					}else{
						f64::INFINITY
					}

				},
				_=>{
					panic!("{:?} is not Supported.",split_function_type);
				}
			};
			if original_loss == 0.0{
				continue;
			}
			if ret.nodes[cnodeid].depth >= max_depth{
				continue;
			}
			if ret.nodes.len() >= max_leaf_nodes{
				break;
			}
			let mut dvec:Vec<usize> = (0..varnum).collect();
			let mut var_check:Vec<usize> = vec![];
			dvec.shuffle(&mut rgen);//連番をランダムにシャッフルする
			while dvec.len() > 0{//fraction 分のシャッフルされた連番を前からとる
				let cnum :usize = dvec.remove(0);
				if var_check.len() as f64/ varnum as f64 > feature_fraction{
					break;
				
				}
				var_check.push(cnum); 
			}
			let mut minres:Option<SplitResult> = None;
			for ii in var_check.iter(){
				
				let sres:Result<SplitResult,String> = match split_function_type{
				SplitFunctionType::ExperimentalMaxBias =>{
					SampleData_usize::calc_threshold_maxbias(classmax,&mut ret.nodes[cnodeid].xsamples,*ii,&spfunc)
				},
				_=>{panic!("undefined split_function_type {:?}",split_function_type);}
				};

				if let Err(_) = sres{
					continue;
				}

				let sres = sres.unwrap();
				if sres.num_lower < min_samples_split ||  sres.num_higher < min_samples_split {
					continue;
				}
				if sres.loss_weighted >= original_loss{//減少が少ない場却下してもいいかも
					continue;
				}
				if let None = minres{
					minres = Some(sres);
				}else{
					if minres.as_ref().unwrap().loss_weighted > sres.loss_weighted{
						minres = Some(sres);
					}
				}
			}
			
			match minres{
				Some(res)=>{
					ret.nodes[cnodeid].set_split(SplittingVar{
						var_index:res.varindex,
						threshold:res.threshold
					});

					let mut lowers:Vec<Box<& SampleData_usize>> = Vec::new();
					let mut highers:Vec<Box<& SampleData_usize>> = Vec::new();
					
					//let dvec = ret.get_node_ref_at(cnodeid).xsamples.clone();
					let dvec = ret.nodes[cnodeid].xsamples.clone();
					for sxx in dvec.into_iter(){
						if res.threshold > sxx.get_value_at(res.varindex){
							lowers.push(Box::new(*sxx));
						}else{
							highers.push(Box::new(*sxx));
						}
					}
					let lnode:DTNode = DTNode{
						splitting_var:None,
						xsamples:lowers,
						child_lower:None,
						child_higher_or_eq:None,
						id:0,
						dummy:false,
						depth:ret.nodes[cnodeid].depth+1,
						answer_string:None,
						answer_usize:None,
						answer_f64:None,
						answer_f64vec:None,

					};

					let hnode:DTNode = DTNode{
						splitting_var:None,
						xsamples:highers,
						child_lower:None,
						child_higher_or_eq:None,
						id:0,
						depth:ret.nodes[cnodeid].depth+1,
						dummy:false,
						answer_string:None,
						answer_usize:None,
						answer_f64:None,
						answer_f64vec:None,
					};
					let lid = ret.push_node(lnode);
					let hid = ret.push_node(hnode);
					ret.nodes[cnodeid].child_lower = Some(lid);
					ret.nodes[cnodeid].child_higher_or_eq = Some(hid);
					node_qued.push(lid);
					node_qued.push(hid);
				
				},
				_=>{
				}
			}
		}

		
		for nn in ret.nodes.iter_mut(){
			nn.set_answer_string((LABEL_ANSCODE_CLASSCOUNT.to_string()+"="+&nn.get_answer_count_string()).as_str());
		}
		return ret;
	}

	pub fn remap_samples(&mut self,samples:&'a Vec<Box<&'a SampleData_usize>>,append:bool){
		if ! append{
			for nn in self.nodes.iter_mut(){
				nn.xsamples = vec![];
			}
		}

		//predict と同じだが、オーバーヘッドありそうなので
		for ss in samples.iter(){
			let mut currentid:usize = 0;
			let stopcount:usize = self.nodes.len();
			let mut lcou:usize = 0;
			loop{
				self.nodes[currentid].xsamples.push(Box::new(ss));
				if self.nodes[currentid].is_leaf(){
					break;
				}
				let splitv = &self.nodes[currentid].splitting_var.as_ref().unwrap();
				if ss.values[splitv.var_index] < splitv.threshold{
					currentid = self.nodes[currentid].child_lower.unwrap();
				}else{
					currentid = self.nodes[currentid].child_higher_or_eq.unwrap();
				}
				if lcou >= stopcount{
					panic!("This tree may have circular reference.");
				} 
				lcou += 1;
			}
		}
	}	
    
	pub fn from_file(filename:&str)->Vec<SimpleDecisionTree<'a>>{
		
        let file = File::open(filename).unwrap();
		let reader = BufReader::new(file);
		
        let mut ret:Vec<SimpleDecisionTree> =vec![];
		let mut buff:Vec<String> =vec![];
		for line in reader.lines(){
			
            let line = line.unwrap();
			match line.find("<<"){
				Some(x)=>{
					if x == 0{
						if buff.len() > 0 {
							buff.push(line);
							ret.push(SimpleDecisionTree::from_string(buff));
							buff = vec![];
						}
					}else{
						buff.push(line);
					}
				},
				_=>{
					buff.push(line);
				}
			}

		}
		return ret;	
	}
	pub fn predict_f64vec(&self,vars:&Vec<f64>)->&Vec<f64>{
		let rnode = self.predict_leafnode(vars);
		return &rnode.get_answer_f64vec();
	}
	pub fn predict_leafnode(&self,vars:&Vec<f64>)->&DTNode<'a>{
		let mut currentid:usize = 0;
		let stopcount:usize = self.nodes.len();
		let mut lcou:usize = 0;
		loop{
			if self.nodes[currentid].is_leaf(){
				return &self.nodes[currentid];
			}
			let splitv = &self.nodes[currentid].splitting_var.as_ref().unwrap();
			if vars[splitv.var_index] < splitv.threshold{
				currentid = self.nodes[currentid].child_lower.unwrap();
			}else{
				currentid = self.nodes[currentid].child_higher_or_eq.unwrap();
			}
			if lcou >= stopcount{
				panic!("This tree may have circular reference.");
			} 
			lcou += 1;
		}
	}


	pub fn normalize_usizevec(vv:Vec<usize>)->Vec<f64>{
		let mut ssum:f64 = vv.iter().fold(0, |sum, a| sum + a) as f64;
		if ssum == 0.0{
			ssum = 1.0;
		}
		return vv.iter().map(|m|{(*m as f64)/ssum}).collect();
	}
	pub fn oob_pruning(source:&'a SimpleDecisionTree,oobsample:&'a Vec<Box<&'a SampleData_usize>>)->SimpleDecisionTree<'a>{
		let mut ret = source.make_copy(false);
		ret.remap_samples(oobsample,false);
		let mut parentnode:Vec<usize> = (0..ret.nodes.len()).map(|m|{m}).collect();
		for (ii,nn) in ret.nodes.iter().enumerate(){	
			if ii != nn.id {
				panic!("index in vec:{} id:{} must be the same.",ii,nn.id);
			}
			if let Some(x) = nn.child_lower{
				parentnode[x] = ii;
			}
			if let Some(x) = nn.child_higher_or_eq{
				parentnode[x] = ii;
			}
		}
		let mut leaf_parent:Vec<usize> = vec![];
		let mut processed:HashSet<usize> = HashSet::new();
		for nn in ret.nodes.iter(){

			if nn.is_leaf(){
				if !processed.contains(&parentnode[nn.id]){
					leaf_parent.push(parentnode[nn.id]);
					processed.insert(parentnode[nn.id]);
				}
			}
		}
//ここから
//Leafだけでなく全部見たほうがいいか

		while leaf_parent.len() > 0{
			let tt = leaf_parent.remove(0);
			if ret.nodes[tt].is_leaf(){
				continue;
			}
			let sparent:Vec<f64> = SimpleDecisionTree::normalize_usizevec(source.nodes[tt].get_answer_count(source.num_classes));
			let dparent:Vec<f64> = SimpleDecisionTree::normalize_usizevec(ret.nodes[tt].get_answer_count(source.num_classes));
			
			/*
			let lid:usize = *ret.nodes[tt].child_lower.as_ref().unwrap();
			let hid:usize = *ret.nodes[tt].child_higher_or_eq.as_ref().unwrap();

			let s_lower:Vec<f64> = SimpleDecisionTree::normalize_usizevec(source.nodes[lid].get_answer_count(source.num_classes));
			let d_lower:Vec<f64> = SimpleDecisionTree::normalize_usizevec(ret.nodes[lid].get_answer_count(source.num_classes));
			
			let s_higher:Vec<f64> = SimpleDecisionTree::normalize_usizevec(source.nodes[hid].get_answer_count(source.num_classes));
			let d_higher:Vec<f64> = SimpleDecisionTree::normalize_usizevec(ret.nodes[hid].get_answer_count(source.num_classes));
			*/
			

			if sparent.iter().fold(0.0,|s,a|{s+a}) == 0.0{
				//こうなるのはおかしい
				eprintln!("Didn't you perform Sample mapping?");
			}

			if dparent.iter().fold(0.0,|s,a|{s+a}) == 0.0{
				if !processed.contains(&parentnode[tt]){
					leaf_parent.push(parentnode[tt]);
					processed.insert(parentnode[tt]);
				}
				continue;
				//dsum = 1.0;
			}


			

			
		}
		return ret;
	}



	pub fn predict(&self,vars:&Vec<f64>)->&str{
		let rnode = self.predict_leafnode(vars);
		return rnode.get_answer_string();
	}
    pub fn from_string(vlines:Vec<String>)
    -> SimpleDecisionTree<'a>{
        let labreg = Regex::new(r"^([^:]+):(.+)$").unwrap();
        let whitespacelike = Regex::new(r"[\t \n\r]").unwrap();
        let tab = Regex::new(r"\t").unwrap();
		
		let mut maxvarindex:i64 = -1;
		let mut maxnodeindex:i64 = -1;
		let mut varnamemap:HashMap<usize,String> = HashMap::new();
		let mut nodemap:HashMap<usize,DTNode> = HashMap::new();
		let mut numclass = 0_usize;
		let mut ret = SimpleDecisionTree{num_classes:0,
			nodes:vec![],
			var_names:None
		};
        for ll in vlines{
            match ll.find("<<"){//最初と最後にゴミがついていることがあるので
                Some(s) => {
                    if s == 0 {
                        continue;
                    }
                },
                _ => {}
            };
			
            match ll.find("#"){
                Some(s) => {
					if s == 0{
						continue;
					}
				},
				_=>{}
			}
			let duplabel:bool = false;
            match ll.find(">>"){
                Some(s) => {
                    if s == 0 {
						if duplabel{
							panic!("Multiple tree params found.");
						}
						let hm = line_to_hashmap(ll,&whitespacelike);
						let plabel = ">>".to_string()+LABEL_PREDICTOR;
						check_keys(&hm,vec![&plabel
						,LABEL_TARGETTYPE
						,LABEL_FORMATVERSION]);
						if hm.get(&plabel).unwrap() != "simple_decision_tree"
						|| hm.get(LABEL_TARGETTYPE).unwrap() != "usize_classification"
						|| hm.get(LABEL_FORMATVERSION).unwrap() != "0.2" {
							panic!("error in code {:?}",hm);
						}
                        continue;
                    }
                },
                _ => {}
            };

            let caps = labreg.captures(ll.as_str());
            match caps{
                Some(cc)=>{
                    let lab = cc.get(1).map_or("", |m| m.as_str());
                    let contents = cc.get(2).map_or("", |m| m.as_str());
                    match lab{
                        LABEL_VARIABLE=>{
							//カラム名のチェックに利用する。なくても可
                            let hm = line_to_hashmap(contents.to_string(),&tab);
                            let vvindex = hm.get(LABEL_VARIABLE_INDEX).unwrap().parse::<usize>().unwrap();
                            let name = hm.get(LABEL_VARIABLE_NAME).unwrap();
							if vvindex as i64 >  maxvarindex{
								maxvarindex = vvindex as i64;
							}
							if varnamemap.contains_key(&vvindex){
								panic!("Index duplication {} .",vvindex);
							}
							varnamemap.insert(vvindex,name.to_string());
                        },
                        LABEL_NODE=>{
                            let hm = line_to_hashmap(contents.to_string(),&tab);
							let nindex:usize = hm.get(LABEL_NODE_ID).unwrap().parse::<usize>().unwrap();
							if nindex as i64 >  maxnodeindex{
								maxnodeindex = nindex as i64;
							}
							if nodemap.contains_key(&nindex){
								panic!("duplicated node id {} .",nindex);
							}
							let mut nnode = DTNode::new(nindex);
							if hm.contains_key(LABEL_NODE_LEAFANSWERSTRING){
								nnode.set_answer_string(hm.get(LABEL_NODE_LEAFANSWERSTRING).unwrap());
							}else{
								let child_l:usize = hm.get(LABEL_NODE_CHILDLOWER).unwrap().parse::<usize>().unwrap();
								let child_h:usize = hm.get(LABEL_NODE_CHILDHIGHER).unwrap().parse::<usize>().unwrap();
								let vii:usize = hm.get(LABEL_NODE_VARID).unwrap().parse::<usize>().unwrap();
								let v:f64 = hm.get(LABEL_NODE_THRESHOLD).unwrap().parse::<f64>().unwrap();
								nnode.child_lower = Some(child_l);
								nnode.child_higher_or_eq = Some(child_h);
								nnode.set_split(SplittingVar{var_index:vii,threshold:v});
							}
							nodemap.insert(nindex,nnode);
                        },
                        LABEL_OPTIONS=>{
                            let hm = line_to_hashmap(contents.to_string(),&tab);
							if hm.contains_key(LABEL_OPTIONS_NUMCLASSES){
								numclass = hm.get(LABEL_OPTIONS_NUMCLASSES).unwrap().parse::<usize>().unwrap();
							}
                        },
                        _ =>{
                            println!("{} is not defined.",lab);
                        }
                    }
                },
                _ =>{
                  println!("{} was not parsed.",ll);
                }
            }
        }
		if maxvarindex > -1{
			//カラム名指定がある場合
			let mut var_names:Vec<String> = vec![];
			for ii in 0..((maxvarindex+1) as usize){
				if varnamemap.contains_key(&ii){
					var_names.push(varnamemap.get(&ii).unwrap().to_string());;
				}else{
					var_names.push("".to_string());
				}
			}
			ret.var_names = Some(var_names);
		}
		
		for ii in 0..((maxnodeindex+1) as usize){
			if nodemap.contains_key(&ii){
				ret.nodes.push(Box::new(nodemap.remove(&ii).unwrap()));
			}else{
				ret.nodes.push(Box::new(DTNode::new(ii)));
			}
		}
		ret.num_classes = numclass;
		let rreg =  Regex::new(COUNTCLASS_DELIMITER).unwrap();

		for nn in ret.nodes.iter_mut(){
			if nn.is_leaf(){
				let mut ansvec:Vec<f64> = vec![0.0;ret.num_classes];
				let ss = nn.get_answer_string();

				let pcode = Regex::new(r"^([^=]+)=(.+)$").unwrap();
				let caps = pcode.captures(ss);
				//println!("{:?}",caps);
				match caps{
					Some(c)=>{
						
 	  	                let scode:&str = c.get(1).map_or("", |m| m.as_str());
 	  	                let sval:&str = c.get(2).map_or("", |m| m.as_str());
						match scode{
							LABEL_ANSCODE_CLASSCOUNT=>{
								let hm:HashMap<String,String> = line_to_hashmap(sval.to_string(),&rreg);
								let mut ssum:f64 = 0.0;
								for (kk,vv) in &hm{
									let k = kk.parse::<usize>().unwrap();
									let v = vv.parse::<usize>().unwrap(); 
									ssum += v as f64;
									ansvec[k] = v as f64;
								}
								if ssum > 0.0 {
									for ii in 0..ansvec.len(){
										ansvec[ii] /= ssum;
									}
								}
								nn.set_answer_f64vec(&ansvec);
							},
							LABEL_ANSCODE_F64VEC=>{
								let hm:HashMap<String,String> = line_to_hashmap(sval.to_string(),&rreg);
								for (kk,vv) in &hm{
									let k = kk.parse::<usize>().unwrap();
									let v = vv.parse::<usize>().unwrap(); 
									ansvec[k] = v as f64;
								}
								nn.set_answer_f64vec(&ansvec);
							},
							LABEL_ANSCODE_USIZE=>{
								nn.set_answer_usize(sval.parse::<usize>().unwrap());
							},
							LABEL_ANSCODE_F64=>{
								nn.set_answer_f64(sval.parse::<f64>().unwrap());
							},
							_=>{eprintln!("{} was considered as string prediction.",ss);}
						}
					},
					_=>{
						//eprintln!("??");
					}

				}
				//let v: Vec<&str> = ss.split(
			}
		}

        return ret;
    }
	pub fn to_string(&self)->String{
		let mut ret:Vec<String> = vec![];
		ret.push(format!(">>{}={}\t{}={}\t{}={}"
		,LABEL_PREDICTOR
		,"simple_decision_tree"
		,LABEL_TARGETTYPE
		,"usize_classification"
		,LABEL_FORMATVERSION
		,"0.2"
		));
		ret.push(format!("{}:{}={}"
		,LABEL_OPTIONS
		,LABEL_OPTIONS_NUMCLASSES
		,self.num_classes));
		if let Some(x) = &self.var_names{
			ret.push("#var_names".to_string());
			for (ii,xx) in x.iter().enumerate(){
				ret.push(format!("{}:{}={}\t{}={}"
				,LABEL_VARIABLE
				,LABEL_VARIABLE_INDEX
				,ii
				,LABEL_VARIABLE_NAME
				,xx));
			}
		}
		ret.push("#nodes".to_string());
		for nn in self.nodes.iter(){
			if !nn.is_leaf(){

				ret.push(format!("{}:\t{}={}\t{}={}\t{}={}\t{}={}\t{}={}"
				,LABEL_NODE
				,LABEL_NODE_ID
				,nn.id
				,LABEL_NODE_CHILDLOWER
				,nn.child_lower.unwrap()
				,LABEL_NODE_CHILDHIGHER
				,nn.child_higher_or_eq.unwrap()
				,LABEL_NODE_VARID
				,nn.splitting_var.as_ref().unwrap().var_index
				,LABEL_NODE_THRESHOLD
				,nn.splitting_var.as_ref().unwrap().threshold));
			}else{
				
				ret.push(format!("{}:\t{}={}\t{}={}"
				,LABEL_NODE
				,LABEL_NODE_ID
				,nn.id
				,LABEL_NODE_LEAFANSWERSTRING
				,nn.get_answer_string()));
				
			}
		}
		ret.push("<<".to_string());
		let mut sret = "".to_string();
		for vv in ret{
			sret += vv.as_str();
			sret += "\n";
		}
		return sret;	
	}
}




#[allow(dead_code,non_camel_case_types)]
#[derive(Debug)]
pub struct SampleData_usize{
	pub values:Vec<f64>,
	pub name:String,
	target_usize:usize,
}

#[allow(dead_code)]
impl SampleData_usize{
	//usize と f64 を作る予定。サンプルの種類により分割方法を変えるので、分割方法はサンプルに impl する
    pub fn new(name:String,target_usize:usize,values:Vec<f64>)->SampleData_usize{
        return SampleData_usize{name,target_usize,values};
    }
    pub fn sort_by(varindex:usize,vec:&mut Vec<Box<&SampleData_usize>>){
        vec.sort_by(|a, b| match a.values.get(varindex).unwrap()
        .partial_cmp(b.values.get(varindex).unwrap()){
            Some(x)=>{
                x
            },
            _=>{
                panic!("Error in sort.")
            }
        });
    }

	pub fn get_target(&self)->usize{
		return self.target_usize
	}
	pub fn get_value_at(&self,i:usize)->f64{
		return self.values[i];
	}

	pub fn calc_gini(max_index:usize,samples:&Vec<Box<&SampleData_usize>>)->f64{
		
		let cdev:usize = max_index+1;
		let mut samplecount_class:Vec<usize> = vec![0;cdev];//クラスに属するサンプル数
		
		for ss in samples.iter(){
			samplecount_class[ss.target_usize]+=1;
		}
		
		let mut orig_gini = 0.0;
		for ii in 0..cdev{
			//https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity
			//let ratio:f64 = (samplecount_class_lower[ii] as f64)/(jj as f64 + 1.0)/(samplecount_class[ii] as f64/d.size() as f64);
			let ratio:f64 = (samplecount_class[ii] as f64)/(samples.len() as f64);
			orig_gini += ratio*ratio;
		}
		orig_gini = 1.0-orig_gini;
		return orig_gini;
	}
	
	pub fn calc_entropy(max_index:usize,samples:&Vec<Box<&SampleData_usize>>)->f64{
		
		let cdev:usize = max_index+1;
		let mut samplecount_class:Vec<usize> = vec![0;cdev];//クラスに属するサンプル数
		
		for ss in samples.iter(){
			samplecount_class[ss.target_usize]+=1;
		}
		
		let mut orig_entropy = 0.0;
		for ii in 0..cdev{
			let ratio:f64 = (samplecount_class[ii] as f64)/(samples.len() as f64);
			if ratio <= 0.0{
				continue;
			}
			orig_entropy -= ratio*(ratio.log2());
		}
		return orig_entropy;
	}

	/**
	 * gini 係数を計算し、最も小さくなる際の Threshold とその際の gini 係数を返す。
	 * target_usize をそのまま配列の index に使用しているので、最大の index を与える必要がある。
	 */
	pub fn calc_threshold_gini(max_index:usize,samples:&mut Vec<Box<&SampleData_usize>>,varindex:usize)->Result<SplitResult,String>{
		
		let cdev:usize = max_index+1;

		let mut samplecount_class:Vec<usize> = vec![0;cdev];//クラスに属するサンプル数
		let mut samplecount_class_lower:Vec<usize> = vec![0;cdev];//Threshold 以下にあるサンプル数
		
		SampleData_usize::sort_by(varindex,samples);
		let mut nolossflag = false;
		let mut minthreshold_start:f64 = 0.0;
		let mut minthreshold_end:f64 = 0.0;
		let mut mingini = f64::INFINITY;


		for ss in samples.iter(){
			samplecount_class[ss.target_usize]+=1;
		}
		
		//使ってない
		let mut _orig_gini = 0.0;
		for ii in 0..cdev{
			//https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity
			//let ratio:f64 = (samplecount_class_lower[ii] as f64)/(jj as f64 + 1.0)/(samplecount_class[ii] as f64/d.size() as f64);
			let ratio:f64 = (samplecount_class[ii] as f64)/(samples.len() as f64);
			_orig_gini += ratio*ratio;
		}
		_orig_gini = 1.0-_orig_gini;


		let ssiz:f64 = samples.len() as f64;
		let mut losses:(f64,f64) = (0.0,0.0);
		let mut nums:(usize,usize) = (0,0);
		for jj in 0..(samples.len()-1){
			samplecount_class_lower[samples[jj].get_target()]+=1;
			if samples[jj].get_value_at(varindex) != samples[jj+1].get_value_at(varindex) {
				let mut e1:f64 = 0.0;
				let mut e2:f64 = 0.0;
				for ii in 0..cdev{
					//https://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity
					//let ratio:f64 = (samplecount_class_lower[ii] as f64)/(jj as f64 + 1.0)/(samplecount_class[ii] as f64/d.size() as f64);
					let ratio:f64 = (samplecount_class_lower[ii] as f64)/(jj as f64 + 1.0);
					e1 += ratio*ratio;
				
					//double ratio = (samplecount_class[ii]-samplecount_class_lower[ii])/(double)(d.size()-jj-1)*samplecount_class[ii]/d.size();
					//let ratio:f64 = ((samplecount_class[ii]-samplecount_class_lower[ii]) as f64)/((samples.len()-jj-1) as f64)*(samplecount_class[ii] as f64/samples.len() as f64 );
					let ratio:f64 = ((samplecount_class[ii]-samplecount_class_lower[ii]) as f64)/(ssiz+(-(jj as f64)-1.0));
					e2 += ratio*ratio;
					
				}
				let ze =  (1.0-e1)*((jj+1) as f64/ssiz)+(1.0-e2)*((ssiz+(-(jj as f64)-1.0))/ssiz);
				if ze == mingini && nolossflag{
					minthreshold_end = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
				}else if ze < mingini {
					mingini = ze;
					losses.0 = 1.0-e1;
					losses.1 = 1.0-e2;
					nums.0 = jj+1;
					nums.1 = samples.len()-jj-1;
					minthreshold_start = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
					minthreshold_end = minthreshold_start;
					nolossflag = true;
				}else{
					nolossflag = false;
				}
			}
		}
		return Ok(SplitResult{
			threshold:(minthreshold_start),
			threshold_range:(minthreshold_start,minthreshold_end),
			num_lower:nums.0,
			num_higher:nums.1,
			loss_lower:losses.0,
			loss_higher:losses.1,
			loss_weighted:mingini,
			varindex:varindex});
	}
	
	/**
	 * entropy を計算し、最も小さくなる際の Threshold とその際の gini 係数を返す。
	 * target_usize をそのまま配列の index に使用しているので、最大の index を与える必要がある。
	 */
	pub fn calc_threshold_entropy(max_index:usize,samples:&mut Vec<Box<&SampleData_usize>>,varindex:usize)->Result<SplitResult,String>{
		
		let cdev:usize = max_index+1;

		let mut samplecount_class:Vec<usize> = vec![0;cdev];//クラスに属するサンプル数
		let mut samplecount_class_lower:Vec<usize> = vec![0;cdev];//Threshold 以下にあるサンプル数
		
		SampleData_usize::sort_by(varindex,samples);
		let mut nolossflag = false;
		let mut minthreshold_start:f64 = 0.0;
		let mut minthreshold_end:f64 = 0.0;
		let mut minentropy = f64::INFINITY;


		for ss in samples.iter(){
			samplecount_class[ss.target_usize]+=1;
		}
		
		let ssiz:f64 = samples.len() as f64;
		let mut losses:(f64,f64) = (0.0,0.0);
		let mut nums:(usize,usize) = (0,0);
		for jj in 0..(samples.len()-1){
			samplecount_class_lower[samples[jj].get_target()]+=1;
			if samples[jj].get_value_at(varindex) != samples[jj+1].get_value_at(varindex) {
				let mut e1:f64 = 0.0;
				let mut e2:f64 = 0.0;
				for ii in 0..cdev{
					let ratio:f64 = (samplecount_class_lower[ii] as f64)/(jj as f64 + 1.0);
					if ratio > 0.0{
						e1 -= ratio*(ratio.log2());
					}
					let ratio:f64 = ((samplecount_class[ii]-samplecount_class_lower[ii]) as f64)/(ssiz+(-(jj as f64)-1.0));
					if ratio > 0.0{
						e2 -= ratio*(ratio.log2());
					}
				}
				let ze =  e1*((jj+1) as f64/ssiz)+e2*((ssiz+(-(jj as f64)-1.0))/ssiz);
				if ze == minentropy && nolossflag{
					minthreshold_end = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
				}else if ze < minentropy {
					minentropy = ze;
					losses.0 = e1;
					losses.1 = e2;
					nums.0 = jj+1;
					nums.1 = samples.len()-jj-1;
					minthreshold_start = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
					minthreshold_end = minthreshold_start;
					nolossflag = true;
				}else{
					nolossflag = false;
				}
			}
		}
		return Ok(SplitResult{
			threshold:(minthreshold_start),
			threshold_range:(minthreshold_start,minthreshold_end),
			num_lower:nums.0,
			num_higher:nums.1,
			loss_lower:losses.0,
			loss_higher:losses.1,
			loss_weighted:minentropy,
			varindex:varindex});
	}

	/**
	 * lower に行く割合と Higher に行く割合の差の合計が最も大きい分割点を採用する
	 */
	pub fn calc_threshold_maxbias(max_index:usize,samples:&mut Vec<Box<&SampleData_usize>>,varindex:usize,spfunc:&MaxBiasFunc)->Result<SplitResult,String>{
		let cdev:usize = max_index+1;

		let mut samplecount_class:Vec<usize> = vec![0;cdev];//クラスに属するサンプル数
		let mut samplecount_class_lower:Vec<usize> = vec![0;cdev];//Threshold 以下にあるサンプル数
		
		SampleData_usize::sort_by(varindex,samples);
		let mut nolossflag = false;
		let mut minthreshold_start:f64 = 0.0;
		let mut minthreshold_end:f64 = 0.0;
		let mut minloss = f64::INFINITY;
		//let mut minloss_average = f64::INFINITY;

		for ss in samples.iter(){
			samplecount_class[ss.target_usize]+=1;
		}
		let mut classcount = 0;
		for ii in 0..samplecount_class.len(){
			if samplecount_class[ii] > 0{
				classcount += 1;
			}
		}
		if classcount == 1{
			return Err("Can not divide this node.".to_owned());
		}

		let mut losses:(f64,f64) = (0.0,0.0);
		let mut nums:(usize,usize) = (0,0);
		for jj in 0..(samples.len()-1){
			samplecount_class_lower[samples[jj].get_target()]+=1;
			if samples[jj].get_value_at(varindex) != samples[jj+1].get_value_at(varindex) {
				
				let mut samplecount_tp_high = 0_usize;//high の方に多かった時に high に入った数
				let mut samplecount_fn_high = 0_usize;//high の方に多かった時に low に入った数
				let mut samplecount_tp_low = 0_usize;//上と逆
				let mut samplecount_fn_low = 0_usize;
				

				let mut _samplecount_tp_high_normalized = 0_f64;//high の方に多かった時に high に入った数
				let mut _samplecount_fn_high_normalized = 0_f64;//high の方に多かった時に low に入った数
				let mut _samplecount_tp_low_normalized = 0_f64;//上と逆
				let mut _samplecount_fn_low_normalized = 0_f64;
				
				
				let mut classcount_low = 0_usize;
				let mut classcount_high = 0_usize;
				
				for ii in 0..cdev{
					if samplecount_class[ii] == 0{
						continue;
					}


					if samplecount_class_lower[ii] > samplecount_class[ii]-samplecount_class_lower[ii]{
					
						classcount_low += 1;
						samplecount_tp_low += samplecount_class_lower[ii];
						samplecount_fn_low += samplecount_class[ii]-samplecount_class_lower[ii];

						_samplecount_tp_low_normalized += (samplecount_class_lower[ii] as f64)/(samplecount_class[ii] as f64) ;
						_samplecount_fn_low_normalized += (samplecount_class[ii] as f64 -samplecount_class_lower[ii] as f64)/(samplecount_class[ii] as f64);



					}else{
						classcount_high += 1;
						samplecount_fn_high += samplecount_class_lower[ii];
						samplecount_tp_high += samplecount_class[ii]-samplecount_class_lower[ii];
						
						_samplecount_fn_high_normalized += (samplecount_class_lower[ii] as f64)/(samplecount_class[ii] as f64) ;
						_samplecount_tp_high_normalized += (samplecount_class[ii] as f64 -samplecount_class_lower[ii] as f64)/(samplecount_class[ii] as f64);

					}
				}



				if classcount_low*classcount_high  <= 0{
					continue;
				}
				
				let ze:f64;
				
				let samplecount_tp_high = samplecount_tp_high as f64;
				let samplecount_tp_low = samplecount_tp_low as f64;
				
				let samplecount_fn_high = samplecount_fn_high as f64;
				let samplecount_fn_low = samplecount_fn_low as f64;
				
				/*normalize すると成績悪い
				let samplecount_tp_high = samplecount_tp_high_normalized as f64;
				let samplecount_tp_low = samplecount_tp_low_normalized as f64;
				
				let samplecount_fn_high = samplecount_fn_high_normalized as f64;
				let samplecount_fn_low = samplecount_fn_low_normalized as f64;
				*/

				let classcount_high = classcount_high as f64;
				let classcount_low = classcount_low as f64;



				//let sptype = MaxBiasFunc::SPFUNC_INFORMEDNESS;

				match spfunc{
					MaxBiasFunc::SPFUNC_ACCMIN=>{
						//0.76 くらい
						//ACC min
						ze = -1.0*((samplecount_tp_low as f64)/((samplecount_tp_low+samplecount_fn_low) as f64))
						.min(samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high)); 
					}, 
					MaxBiasFunc::SPFUNC_ACCSUM=>{
						//ACC sum
						ze = -1.0*((samplecount_tp_low as f64)/((samplecount_tp_low+samplecount_fn_low) as f64)
						+ samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high)); 
					}, 
					MaxBiasFunc::SPFUNC_MCC=>{
						//0.82 くらい
						//MCC
						if (samplecount_tp_high+samplecount_fn_high)*
							(samplecount_tp_high+samplecount_fn_low)*
							(samplecount_tp_low+samplecount_fn_high)*
							(samplecount_tp_low+samplecount_fn_low)
							<= 0.0{
							panic!("???");
						}else{
							let mcc:f64 =
							(samplecount_tp_high*samplecount_tp_low
							-  samplecount_fn_high*samplecount_fn_low)
							/((samplecount_tp_high+samplecount_fn_high)*
							(samplecount_tp_high+samplecount_fn_low)*
							(samplecount_tp_low+samplecount_fn_high)*
							(samplecount_tp_low+samplecount_fn_low)
							).sqrt();
							ze = -1.0*mcc;
						}
					}, 
					MaxBiasFunc::SPFUNC_F1_MIN=>{
						//0.81 くらい
						if ((samplecount_tp_high as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64)
						*((samplecount_tp_low as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64) == 0.0{
							
							panic!("???");
						}else{
							let f1_high = (samplecount_tp_high as f64*2.0)/((samplecount_tp_high as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64);
							let f1_low = (samplecount_tp_low as f64*2.0)/((samplecount_tp_low as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64);
							ze = -1.0*f1_high.min(f1_low);				
						}
					}, 
					MaxBiasFunc::SPFUNC_F1_SUM=>{
						//0.81 くらい
						if ((samplecount_tp_high as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64)
						*((samplecount_tp_low as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64) == 0.0{
							ze = 100.0;
						}else{
							let f1_high = (samplecount_tp_high as f64*2.0)/((samplecount_tp_high as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64);
							let f1_low = (samplecount_tp_low as f64*2.0)/((samplecount_tp_low as f64*2.0)+ samplecount_fn_high as f64+samplecount_fn_low as f64);
							ze = -1.0*(f1_high+f1_low);				
						}
					}, 
					MaxBiasFunc::SPFUNC_PPV_MIN=>{

						//0.7996
						if (samplecount_tp_high+samplecount_fn_low)*(samplecount_tp_low+samplecount_fn_high) <= 0.0{
							panic!("???");
						}else{
							let ppv_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_low);
							let ppv_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_high);
							ze = -1.0*(ppv_high.min(ppv_low));
						}
					}, 
					MaxBiasFunc::SPFUNC_PPV_SUM=>{

						//0.80 くらい
						if (samplecount_tp_high+samplecount_fn_low)*(samplecount_tp_low+samplecount_fn_high) <= 0.0 {
							ze = 100.0;
						}else{
							let ppv_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_low);
							let ppv_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_high);
							ze = -1.0*(ppv_high + ppv_low);
						}
					}, 
					MaxBiasFunc::SPFUNC_FDR_MAX=>{

						//0.7996 1.0-ppv なので上と同じ
						if (samplecount_tp_high+samplecount_fn_low)*(samplecount_tp_low+samplecount_fn_high) <= 0.0 {
							ze = 100.0;
						}else{
							let fdr_high = samplecount_fn_low/(samplecount_tp_high+samplecount_fn_low);
							let fdr_low = samplecount_fn_high/(samplecount_tp_low+samplecount_fn_high);
							ze = fdr_high.max(fdr_low);
						}
					}, 
					MaxBiasFunc::SPFUNC_FDR_SUM=>{
						//0.8026
						if (samplecount_tp_high+samplecount_fn_low)*(samplecount_tp_low+samplecount_fn_high) <= 0.0 {
							ze = 100.0;
						}else{
							let fdr_high = samplecount_fn_low/(samplecount_tp_high+samplecount_fn_low);
							let fdr_low = samplecount_fn_high/(samplecount_tp_low+samplecount_fn_high);
							ze = fdr_high + fdr_low;
						}
					}, 
					MaxBiasFunc::SPFUNC_TPR_MIN=>{
						//0.7514
						if (samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0 {
							ze = 100.0;
						}else{	
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							ze = -1.0*tpr_high.min(tpr_low);
						}
					}, 
					MaxBiasFunc::SPFUNC_TPR_SUM=>{
						//0.7632
						if (samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0{
							panic!("???");
						}else{	
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							ze = -1.0*(tpr_high + tpr_low);
						}
					},
					MaxBiasFunc::SPFUNC_PLR_MIN=>{
						//0.7732
						if samplecount_fn_low*samplecount_fn_high*(samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0
						 {
							ze = 100.0;
						}else{
							let mut lrp_high:f64 = samplecount_tp_high*(samplecount_tp_low+samplecount_fn_low)
							/ (samplecount_fn_low*(samplecount_tp_high+samplecount_fn_high));
							let mut lrp_low:f64 = samplecount_tp_low*(samplecount_tp_high+samplecount_fn_high)
							/(samplecount_fn_high*(samplecount_tp_low+samplecount_fn_low));
							
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
						
							lrp_low = lrp_low.min(10000.0+tpr_high+tpr_low);
							lrp_high = lrp_high.min(10000.0+tpr_high+tpr_low);
							ze = -1.0*lrp_high.min(lrp_low);
						}
					}, 
					MaxBiasFunc::SPFUNC_PLR_SUM=>{
						//0.7601
						if samplecount_fn_low*samplecount_fn_high*(samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0 {
							ze = 100.0;
						}else{
							let mut lrp_high:f64 = samplecount_tp_high*(samplecount_tp_low+samplecount_fn_low)
							/ (samplecount_fn_low*(samplecount_tp_high+samplecount_fn_high));
							let mut lrp_low:f64 = samplecount_tp_low*(samplecount_tp_high+samplecount_fn_high)
							/(samplecount_fn_high*(samplecount_tp_low+samplecount_fn_low));
							
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
					
							lrp_low = lrp_low.min(10000.0+tpr_high+tpr_low);
							lrp_high = lrp_high.min(10000.0+tpr_high+tpr_low);
							
							ze = -1.0*(lrp_high + lrp_low);
						}	
					}, 
					MaxBiasFunc::SPFUNC_NLR_MIN=>{
						
						//0.7668
					
						if samplecount_tp_low*samplecount_tp_high*(samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0{
							ze = 100.0;
						}else{
							let mut lrm_high:f64 = samplecount_fn_high*(samplecount_tp_low+samplecount_fn_low)
							/ (samplecount_tp_low*(samplecount_tp_high+samplecount_fn_high));
							let mut lrm_low:f64 = samplecount_fn_low*(samplecount_tp_high+samplecount_fn_high)
							/(samplecount_tp_high*(samplecount_tp_low+samplecount_fn_low));
							
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
						
							lrm_low = lrm_low.min(10000.0-tpr_high-tpr_low);
							lrm_high = lrm_high.min(10000.0-tpr_high-tpr_low);
							ze = lrm_high.max(lrm_low);
						}
					}, 
					MaxBiasFunc::SPFUNC_NLR_SUM=>{
						//0.7668
						if samplecount_tp_low*samplecount_tp_high*(samplecount_tp_high+samplecount_fn_high)*(samplecount_tp_low+samplecount_fn_low) <= 0.0 {
							ze = 100.0;
						}else{
							let mut lrm_high:f64 = samplecount_fn_high*(samplecount_tp_low+samplecount_fn_low)
							/ (samplecount_tp_low*(samplecount_tp_high+samplecount_fn_high));
							let mut lrm_low:f64 = samplecount_fn_low*(samplecount_tp_high+samplecount_fn_high)
							/(samplecount_tp_high*(samplecount_tp_low+samplecount_fn_low));
							
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
						
							lrm_low = lrm_low.min(10000.0-tpr_high-tpr_low);
							lrm_high = lrm_high.min(10000.0-tpr_high-tpr_low);

							ze = lrm_high+lrm_low;
						}
					},
					MaxBiasFunc::SPFUNC_DOR=>{
						//0.7532
						if samplecount_fn_low*samplecount_fn_high <= 0.0 {
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
							ze = -1.0*((samples.len() as f64).powf(2.0)+tpr_high+tpr_low);
						}else{
							let mut dor = (samplecount_tp_low*samplecount_tp_high)/(samplecount_fn_low*samplecount_fn_high );
							
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							
							dor = dor.min((samples.len() as f64).powf(2.0)+tpr_high+tpr_low);
							ze = dor*-1.0;
						}
					}, 
					MaxBiasFunc::SPFUNC_INFORMEDNESS=>{
						//https://en.wikipedia.org/wiki/Informedness
						if classcount_high*classcount_low <= 0.0{
							ze = 1000.0;
						}else{
							let tpr_high = samplecount_tp_high/(samplecount_tp_high+samplecount_fn_high);
							let tpr_low = samplecount_tp_low/(samplecount_tp_low+samplecount_fn_low);
							ze = (tpr_low+tpr_high-1.0)*-1.0;
						}
						
					}, 
					MaxBiasFunc::SPFUNC_COHEN_KAPPA=>{
					
						//0.8096
						//https://en.wikipedia.org/wiki/Cohen%27s_kappa
						if classcount_high*classcount_low <= 0.0{
							ze = 1000.0;
						}else{
							let po = (samplecount_tp_high+samplecount_tp_low)/(samples.len() as f64);
							let pe1 = (samplecount_tp_high+samplecount_fn_low)/(samples.len() as f64)
							*(samplecount_tp_high+samplecount_fn_high)/(samples.len() as f64);
							let pe2 = (samplecount_tp_low+samplecount_fn_high)/(samples.len() as f64)
							*(samplecount_tp_low+samplecount_fn_low)/(samples.len() as f64);
							ze = -1.0*(po-pe1-pe2)/(1.0-pe1-pe2);
						}
					}, 
					MaxBiasFunc::SPFUNC_FLEISS_KAPPA=>{
					
						//0.8118
						//https://en.wikipedia.org/wiki/Fleiss%27_kappa
						if classcount_high*classcount_low <= 0.0{
							ze = 1000.0;
						}else{
							let po = (samplecount_tp_high+ samplecount_tp_low)/(samples.len() as f64 );
							let pe1 = (samplecount_tp_high*2.0+samplecount_fn_low)
							/(samples.len() as f64 *2.0);					 
							let pe2 = (samplecount_tp_low*2.0+samplecount_fn_high)
							/(samples.len() as f64 *2.0);					 
							ze = -1.0*(po - pe1*pe1 - pe2*pe2)/(1.0 - pe1*pe1 - pe2*pe2);
						}
					}, 
					MaxBiasFunc::SPFUNC_CHISQ=>{
					
						//0.8118
						//https://en.wikipedia.org/wiki/Fleiss%27_kappa
						if classcount_high*classcount_low <= 0.0{
							ze = 1000.0;
						}else{
							let spp = samples.len() as f64;
							let k1 = (samplecount_tp_high+samplecount_fn_high)/spp;
							let k2 = (samplecount_tp_low+samplecount_fn_low)/spp;
							let k3 = (samplecount_tp_high+samplecount_fn_low)/spp;
							let k4 = (samplecount_tp_low+samplecount_fn_high)/spp;
							
							let l1 = (k1*k3*spp-samplecount_tp_high).powf(2.0);
							let l2 = (k1*k4*spp-samplecount_fn_high).powf(2.0);
							let l3 = (k2*k3*spp-samplecount_fn_low).powf(2.0);
							let l4 = (k2*k4*spp-samplecount_tp_low).powf(2.0);

							ze = -1.0*(l1+l2+l3+l4);
						}
					}
				}

				if ze == minloss && nolossflag{
					minthreshold_end = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
				}else if ze < minloss {
					minloss = ze;
					losses.0 = ze;
					losses.1 = ze;
					nums.0 = jj+1;
					nums.1 = samples.len()-jj-1;
					minthreshold_start = samples[jj].get_value_at(varindex)/2.0+samples[jj+1].get_value_at(varindex)/2.0;
					minthreshold_end = minthreshold_start;
					nolossflag = true;
				}else{
					nolossflag = false;
				}

			}
		}
		return Ok(SplitResult{
			//threshold:(minthreshold_start*0.5+minthreshold_end*0.5),
			threshold:minthreshold_start,
			threshold_range:(minthreshold_start,minthreshold_end),
			num_lower:nums.0,
			num_higher:nums.1,
			loss_lower:losses.0,
			loss_higher:losses.1,
			loss_weighted:minloss,
			varindex:varindex});
	}

}

#[allow(dead_code)]
#[derive(Debug)]
pub struct SplitResult{
	threshold:f64,
	threshold_range:(f64,f64),
	num_lower:usize,
	num_higher:usize,
	loss_lower:f64,
	loss_higher:f64,
	loss_weighted:f64,
	varindex:usize,
}
#[derive(Debug)]
#[allow(dead_code)]
pub struct DTNode<'a>{
    splitting_var:Option<SplittingVar>,
    xsamples:Vec<Box<&'a SampleData_usize>>,
	child_lower:Option<usize>,
	child_higher_or_eq:Option<usize>,
	id:usize,
	depth:usize,
	dummy:bool,
	answer_string:Option<String>,
	answer_usize:Option<usize>,
	answer_f64:Option<f64>,
	answer_f64vec:Option<Vec<f64>>,
}

impl<'a> DTNode<'a>{
	pub fn make_copy(&self,copysamples:bool)->DTNode<'a>{
		let mut ret = DTNode{
		splitting_var:match self.splitting_var.as_ref(){
			Some(x)=>{Some(x.make_copy())},
			_=>{None}
		},
		xsamples:vec![],
		child_lower:match self.child_lower{
		Some(x)=>{Some(x.clone())},
		_=>{None}},
		child_higher_or_eq:match self.child_higher_or_eq{
		Some(x)=>{Some(x.clone())},
		_=>{None}},
		id:self.id.clone(),
		depth:self.depth.clone(),
		dummy:self.dummy.clone(),
		answer_string:match self.answer_string.as_ref(){
		Some(x)=>{Some(x.clone())},
		_=>{None}},
		answer_usize:match self.answer_usize.as_ref(){
		Some(x)=>{Some(x.clone())},
		_=>{None}},
		answer_f64:match self.answer_f64.as_ref(){
		Some(x)=>{Some(x.clone())},
		_=>{None}},
		answer_f64vec:match self.answer_f64vec.as_ref(){
		Some(x)=>{Some(x.clone())},
		_=>{None}
		}
		};

		if copysamples {
			ret.xsamples = self.xsamples.iter().map(|m| m.clone()).collect();
		}
		return ret;
	}
	pub fn new(ii:usize)->DTNode<'a>{
		return DTNode{
			splitting_var:None,
			xsamples:vec![],
			child_lower:None,
			child_higher_or_eq:None,
			id:ii,
			depth:0,
			dummy:true,
			answer_string:None,
			answer_usize:None,
			answer_f64:None,
			answer_f64vec:None
		};
	}


	pub fn is_dummy(&self)->bool{
		return self.dummy;
	}
	pub fn set_dummy_flag(&mut self,d:bool){
		self.dummy = d;
	}
	pub fn set_split(&mut self,s:SplittingVar){
		self.splitting_var = Some(s);
	}
	pub fn set_answer_string(&mut self,s:&str){
		self.answer_string = Some(s.to_string());
	}
	pub fn set_answer_f64(&mut self,f:f64){
		self.answer_f64 = Some(f);
	}
	
	pub fn set_answer_usize(&mut self,u:usize){
		self.answer_usize = Some(u);
	}
	
	pub fn set_answer_f64vec(&mut self,f:&Vec<f64>){
		self.answer_f64vec = Some(f.iter().map(|m|{m.clone()}).collect());
	}

	pub fn get_answer_string(&self)->&str{
		if let Some(x) = self.answer_string.as_ref(){
			return x.as_str();
		}else{
			panic!("This node does not have string data.");
		}
	}
	
	pub fn get_answer_f64(&self)->&f64{
		if let Some(x) = self.answer_f64.as_ref(){
			return &x;
		}else{
			panic!("This node does not have f64 data.");
		}
	}
	
	
	pub fn get_answer_usize(&self)->&usize{
		if let Some(x) = self.answer_usize.as_ref(){
			return &x;
		}else{
			panic!("This node does not have usize data.");
		}
	}
	
	pub fn get_answer_f64vec(&self)->&Vec<f64>{
		if let Some(x) = self.answer_f64vec.as_ref(){
			return &x;
		}else{
			panic!("This node does not have f64vec data.");
		}
	}

	/**
	 * Vec に含まれる Sample が持っている target の数を数えて (target, count) のタプルの Vec として返す。 
	 */
	pub fn count_classes(samplevec:&Vec<Box<&SampleData_usize>>)->Vec<(usize,usize)>{
		let mut hm:HashMap<usize,usize> = HashMap::new();
		for xx in samplevec.iter(){
			*hm.entry(xx.target_usize).or_insert(0) += 1;
		}
		let mut vvec: Vec<_> = hm.into_iter().collect();
		vvec.sort_by(|a, b| a.1.cmp(&b.1).reverse());
		return vvec;
	}


	/**
	 * クラスに属するサンプルの数を返す。
	 */
	pub fn get_answer_count(&self,vecsize:usize)->Vec<usize>{
		let mut ret:Vec<usize> = vec![0;vecsize];
		
		let vvec:Vec<(usize,usize)> = DTNode::count_classes(&self.xsamples);
		for vv in vvec.iter(){
			ret[vv.0]= vv.1;
		}
		return ret;
	}

	/**
	 * ID と 数のペアを文字列にして返す。
	 */
	
	pub fn get_answer_count_string(&self)->String{
		if self.xsamples.len() == 0{
			return "".to_string();
		}
		let mut ret:String = "".to_string();
		let vvec:Vec<(usize,usize)> = DTNode::count_classes(&self.xsamples);
		for (ii,vv) in vvec.iter().enumerate(){
			ret += format!("{}={}",vv.0.to_string(),vv.1.to_string()).as_str();
			if ii != vvec.len() -1{
				ret += COUNTCLASS_DELIMITER;
			}
		}
		return ret;
	}
	
	pub fn get_answer_majority(&self)->String{
		let vvec:Vec<(usize,usize)> = DTNode::count_classes(&self.xsamples);
		return vvec[0].0.to_string(); 
	}


	pub fn is_leaf(&self)->bool{
		if let None = self.child_lower{
			if let None = self.child_higher_or_eq{
				return true; 
			}
		}
		return false;
	}
}


#[derive(Debug)]
#[allow(dead_code)]
pub struct SplittingVar{
    pub var_index:usize,
    pub threshold:f64
}

#[allow(dead_code)]
impl SplittingVar{

	pub fn make_copy(&self)->SplittingVar{
		return SplittingVar{var_index:self.var_index.clone(),threshold:self.threshold.clone()};
	}
    pub fn new(var_index:usize,threshold:f64)->SplittingVar{
        return SplittingVar{var_index,threshold};
    }
    pub fn set_var_index(&mut self,i:usize){
        self.var_index = i;
    }
    pub fn set_threshold(&mut self,t:f64){
        self.threshold = t;
    }
}

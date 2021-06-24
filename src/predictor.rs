use std::collections::HashMap;

#[allow(dead_code)]
pub struct PredictionMachine{
    pub usize_to_classname:HashMap<usize,String>,
    pub classname_to_usize:HashMap<String,usize>,
    pub predictors:Vec<Box<Predictor>>
}


pub trait Predictor {
    fn predict_value(&self,v:&Vec<f64>) -> f64;//regression
    fn predict_class(&self,v:&Vec<f64>) -> Vec<f64>;//classification
}

#[derive(Debug)]
pub enum AnswerType{
	USIZE,
	F64,
    F64VEC,
    STRING,
}

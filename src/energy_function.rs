
use super::ff_env;

//他言語で作った場合に丁度の値の場合誤差になるので
//PDB は精度小数点以下 3 くらいなのでこれくらいで十分と思う
pub const ROUNDING_EPSILON:f64 = 0.00001;

pub trait EnergyFunction{
    fn calc_energy(&self,mdenv:&ff_env::FFEnv,atom_level_energy:&mut Vec<f64>,weight:f64)->f64;
}

pub trait Binning{
    fn get_div_unit(&self)->f64;
    fn get_num_bins(&self)->usize;
    fn get_start_point(&self)->f64;
}

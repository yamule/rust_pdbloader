
#[allow(unused_imports)]
use super::geometry;
use super::geometry::Vector3D;
use super::geometry::Point3D;
use super::process_3d;
use super::energy_function;

use std::f64::consts::PI;
const DEGREE_TO_RADIAN2:f64 = PI/180.0*PI/180.0;
pub const EPSILON:f64 = 1.0e-20;


impl Vector3D for MDAtom{

    fn get_xyz(&self)-> (f64,f64,f64){
        return (self.get_x(),self.get_y(),self.get_z());
    }

    fn set_vector(&mut self,v:&dyn Vector3D){
        self.set_x(v.get_x());
        self.set_y(v.get_y());
        self.set_z(v.get_z());
    }

    fn set_x(&mut self,xx:f64){
        self.x = xx;
    }
    fn set_y(&mut self,yy:f64){
        self.y = yy;
    }
    fn set_z(&mut self,zz:f64){
        self.z = zz;
    }
    fn get_x(&self)->f64{
        return self.x;
    }
    fn get_y(&self)->f64{
        return self.y;
    }
    fn get_z(&self)->f64{
        return self.z;
    }
    
    fn set_xyz(&mut self,xx:f64,yy:f64,zz:f64){
        self.set_x(xx);
        self.set_y(yy);
        self.set_z(zz);
    }
    fn distance(&self,point:&dyn Vector3D)->f64{
        let dd = (self.get_x() - point.get_x()).powf(2.0)+
        (self.get_y() - point.get_y()).powf(2.0)+
        (self.get_z() - point.get_z()).powf(2.0);
        if dd == 0.0{
            return 0.0;
        }
        return dd.sqrt();
    }
}


#[derive(Debug,Clone)]
pub struct MDAtom{
    pub atom_type:String,
    pub atom_name:String,
    pub chain_name:String,
    pub residue_name:String,
    pub residue_number:i64,
    pub residue_ins_code:String,
    
    pub residue_index_in_chain:i64,
    pub atom_index:usize,
    
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub mass:f64,
    pub charge:f64,
    pub created:bool,
    pub unplaced:bool,

    pub nb_epsilon:f64,
    pub nb_r1_2:f64,//Rmin1/2

    pub nb_14flag:bool,//無いの多い
    pub nb_14_epsilon:f64,
    pub nb_14_r1_2:f64,
    //pub hybrid_orbit:HybridOrbit,

    pub nterminal:bool,
    pub cterminal:bool,
}



pub struct MDEnv{
    pub atoms:Vec<MDAtom>,
    pub dist:Vec<Vec<f64>>,
    pub num_edges:Vec<Vec<u64>>,//ある原子とある原子の最短距離（存在するエッジの数）。一定数以上についてはカウントしない。今はその一定数は 5。
}
impl MDEnv{
}



pub fn calc_dihedral_angle_radian(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D)->f64{
    let norm0:(f64,f64,f64) = process_3d::calc_norm(
        v0.get_x() - v1.get_x(),
        v0.get_y() - v1.get_y(),
        v0.get_z() - v1.get_z(),
        v2.get_x() - v1.get_x(),
        v2.get_y() - v1.get_y(),
        v2.get_z() - v1.get_z()
    );
    let norm1:(f64,f64,f64) = process_3d::calc_norm(
        v1.get_x() - v2.get_x(),
        v1.get_y() - v2.get_y(),
        v1.get_z() - v2.get_z(),
        v3.get_x() - v2.get_x(),
        v3.get_y() - v2.get_y(),
        v3.get_z() - v2.get_z()
    );
    let norm2:(f64,f64,f64) = process_3d::calc_norm(
        v2.get_x() - v1.get_x(),
        v2.get_y() - v1.get_y(),
        v2.get_z() - v1.get_z(),
        norm0.0,
        norm0.1,
        norm0.2,
    );
    let mut direc:f64 = 1.0;
    if process_3d::distance(&(norm1.0,norm1.1,norm1.2),&(norm2.0,norm2.1,norm2.2))
    >  (2.0_f64).sqrt(){
        direc = -1.0;
    }
    let drad = calc_angle_radian(
        &Point3D{x:norm0.0,y:norm0.1,z:norm0.2},
        &Point3D{x:0.0,y:0.0,z:0.0},
        &Point3D{x:norm1.0,y:norm1.1,z:norm1.2}
    )*direc;

    return drad;
}

pub fn calc_dihedral_angle(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D,v3:&dyn Vector3D)->f64{
    return calc_dihedral_angle_radian(v0, v1, v2, v3)*180.0/PI;
}


//v0-v1-v2 の angle を radian で返す。
//あらかじめ standerdize しておく必要はない。
pub fn calc_angle_radian(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D)->f64{
    let mut vx0 = v0.get_x() - v1.get_x();
    let mut vy0 = v0.get_y() - v1.get_y();
    let mut vz0 = v0.get_z() - v1.get_z();

    let mut vx2 = v2.get_x() - v1.get_x();
    let mut vy2 = v2.get_y() - v1.get_y();
    let mut vz2 = v2.get_z() - v1.get_z();

    let mut dis0 = vx0*vx0+vy0*vy0+vz0*vz0;
    let mut dis2 = vx2*vx2+vy2*vy2+vz2*vz2;
    if dis0 == 0.0 || dis2 == 0.0{
        return 0.0;
    }
    dis0 = dis0.sqrt();
    dis2 = dis2.sqrt();
    vx0 /= dis0;
    vy0 /= dis0;
    vz0 /= dis0;
    
    vx2 /= dis2;
    vy2 /= dis2;
    vz2 /= dis2;

    let dz = process_3d::distance(&(vx0,vy0,vz0),&(vx2,vy2,vz2));
    
    return ((2.0-dz*dz)/2.0).max(-1.0).min(1.0).acos();
}

//v0-v1-v2 の angle を 360 degree で返す。
pub fn calc_angle(v0:&dyn Vector3D,v1:&dyn Vector3D,v2:&dyn Vector3D)->f64{
    return calc_angle_radian(v0,v1,v2)*180.0/PI;
}



#[derive(Debug,Clone)]
pub struct DihedralVars{
    pub atoms:(usize,usize,usize,usize),
    
    kchi:f64,
    n:f64,//usize だが変換が面倒なので
    delta:f64,
    debug_string:String
}
impl energy_function::EnergyFunction for DihedralVars{
    fn calc_energy(&self,mdenv:&MDEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{

        let phi:f64 = calc_dihedral_angle(&mdenv.atoms[self.atoms.0], &mdenv.atoms[self.atoms.1],&mdenv.atoms[self.atoms.2],&mdenv.atoms[self.atoms.3]);
        let dsc:f64 = self.kchi*(1.0+((self.n*phi-self.delta)/180.0*PI).cos())*weight;
        
        atom_level_energy[self.atoms.0] += dsc/4.0;
        atom_level_energy[self.atoms.1] += dsc/4.0;
        atom_level_energy[self.atoms.2] += dsc/4.0;
        atom_level_energy[self.atoms.3] += dsc/4.0;
        
        return dsc;
    }
}

#[derive(Debug,Clone)]
pub struct IMPRVars{
    pub atoms:(usize,usize,usize,usize),
    kpsi:f64,
    psi0:f64
}

impl energy_function::EnergyFunction for IMPRVars{
    fn calc_energy(&self,mdenv:&MDEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        let psi = calc_dihedral_angle(
            &mdenv.atoms[self.atoms.0]
            ,&mdenv.atoms[self.atoms.1]
            ,&mdenv.atoms[self.atoms.2]
            ,&mdenv.atoms[self.atoms.3])*-1.0;
        let dsc = self.kpsi*(psi-self.psi0)*(psi-self.psi0) *DEGREE_TO_RADIAN2*weight;
        atom_level_energy[self.atoms.0] += dsc/4.0;
        atom_level_energy[self.atoms.1] += dsc/4.0;
        atom_level_energy[self.atoms.2] += dsc/4.0;
        atom_level_energy[self.atoms.3] += dsc/4.0;
        return dsc;
    }
}


#[derive(Debug,Clone)]
#[allow(dead_code)]
pub struct CMAPVars{
    atoms:(usize,usize,usize,usize),
    matrixindex:usize,//Matrix 全部コピーするのは冗長なので
}

impl energy_function::EnergyFunction for CMAPVars{
    
#[allow(unused_variables)]
    fn calc_energy(&self,mdenv:&MDEnv, atom_level_energy: &mut Vec<f64>,weight:f64)->f64{
        //CMAP は別の方法使った方が良いと思う
        panic!("not implemented yet");
    }
}



use core::f64;
use std::fs;
use std::cmp::Ordering;
use std::io::{BufWriter,Write};
use regex::Regex;
use std::collections::{HashMap,HashSet};
use crate::{geometry::{Point3D, Vector3D}, mmcif_process, pdbdata::{PDBAsym, PDBAtom,PDBComp, PDBEntry}, process_3d, structural_alignment::align};
use super::structural_alignment;
use super::matrix_process;
use chrono::{Local};


pub struct PseudoAtom{
    pub x:f64,
    pub y:f64,
    pub z:f64,
    pub radius:f64
}

impl Vector3D for PseudoAtom{

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

pub struct SideChainCylinder<'a>{
    pub n:&'a PseudoAtom,
    pub ca:&'a PseudoAtom,
    pub c:&'a PseudoAtom,
    pub endpoint:Point3D,
    pub _x:Point3D,//n ca c の norm
    pub _y:Point3D,//
    pub _z:Point3D,//
    pub length:f64,//CA から最終点の距離
    pub num_sep:usize,//length を何個に分割するか//90 度ごと 4 つに分割するので、領域はこれ x4
}
impl<'a> SideChainCylinder<'a>{
    pub fn update(&mut self){
        let mut nn = self.n.get_xyz();
        let caa = self.ca.get_xyz();
        let mut cc = self.c.get_xyz();
        nn = process_3d::standardize(nn.0-caa.0,nn.1-caa.1,nn.2-caa.2);
        cc = process_3d::standardize(cc.0-caa.0,cc.1-caa.1,cc.2-caa.2);


        
        let yy:(f64,f64,f64) = process_3d::calc_norm_t(
            &nn,
            &(0.0,0.0,0.0),
            &cc
        );
        let zz = process_3d::standardize(
            -1.0*(nn.0/2.0+cc.0/2.0)
            ,-1.0*(nn.1/2.0+cc.1/2.0)
            ,-1.0*(nn.2/2.0+cc.2/2.0)
            );

        let xx:(f64,f64,f64) = process_3d::calc_norm_t(
            &zz,
            &(0.0,0.0,0.0),
            &yy
        );
        

         


    }
    pub fn is_in(&self,atom:&dyn Vector3D){
        
    }
}






use core::f64;
use std::fs;
use std::cmp::Ordering;
use std::io::{BufWriter,Write};
use regex::Regex;
use std::collections::{HashMap,HashSet};
use crate::{geometry::{Point3D, Vector3D}, mmcif_process, pdbdata::{PDBAsym, PDBAtom,PDBComp, PDBEntry}, process_3d, structural_alignment::align};
use super::structural_alignment;
use super::matrix_process;
use super::geometry;
use chrono::{Local};
use std::f64::consts::PI;

const to_radian:f64 = PI/180.0;
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
impl PseudoAtom{
    pub fn new(xx:f64,yy:f64,zz:f64)->PseudoAtom{
        return PseudoAtom{
            x:xx,y:yy,z:zz,radius:1.0
        };
    }
    pub fn get_radius(&self)->f64{
        return self.radius;
    }
    pub fn set_radius(&mut self,r:f64){
        self.radius = r;
    }
}
pub struct SideChainCylinder<'a>{
    pub n:&'a PseudoAtom,
    pub ca:&'a PseudoAtom,
    pub c:&'a PseudoAtom,
    pub endpoint:Point3D,
    pub rx:Point3D,//n ca c の norm と rz 
    pub ry:Point3D,//n ca c の norm
    pub rz:Point3D,//n と c を 単位ベクトル化して中点から ca を通る線。
    pub length:f64,//CA から最終点の距離
    pub num_sep:usize,//length を何個に分割するか//90 度ごと 4 つに分割するので、領域はこれ x4
}
impl<'a> SideChainCylinder<'a>{
    pub fn new(n:&'a PseudoAtom,ca:&'a PseudoAtom,c:&'a PseudoAtom)->SideChainCylinder<'a>{
        let mut ret = SideChainCylinder{
        n:n,
        ca:ca,
        c:c,
        endpoint:Point3D::new(0.0,0.0,0.0),
        rx:Point3D::new(0.0,0.0,0.0),//n ca c の norm と rz 
        ry:Point3D::new(0.0,0.0,0.0),//n ca c の norm
        rz:Point3D::new(0.0,0.0,0.0),//n と c を 単位ベクトル化して中点から ca を通る線。
        length:8.0,//CA から最終点の距離
        num_sep:4
        };
        ret.update();
        return ret;
    }
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
        self.rz.set_xyz(zz.0*self.length,zz.1*self.length,zz.2*self.length);
        //C-CA-CB の ANGLE は charmm19  (param19)で 
        //C    CH1E CH3E    70.0     106.5
        //N-CA-CB の ANGLE は 
        //CH3E CH1E NH1     65.0     108.5
        //N-CA-C は
        //C    CH1E NH1     45.0     111.6
        //ALA の IMPROPER (toph19)は 
        //IC   N    C    *CA  CB     0.0000    0.00  120.00    0.00   0.0000
        //Excel で値を入れて計算してみると少しずれる。。。
        let xx:(f64,f64,f64) = process_3d::calc_norm_t(
            &zz,
            &(0.0,0.0,0.0),
            &yy
        );

        let mut cbb = Point3D::from_tuple(&cc);
        process_3d::rotate_3d(&mut vec![&mut cbb],&Point3D::from_tuple(&yy),-106.5*to_radian);
        process_3d::rotate_3d(&mut vec![&mut cbb],&Point3D::from_tuple(&cc),60.0*to_radian);
        if true{
            let mut ggeo = geometry::Geometry::new();
            for a in vec![&nn,&caa,&cc]{
                let mut mm = geometry::Geometry::generate_sphere(a,0.25,8,8);
                for f in mm.1.iter_mut(){
                    f.set_color(&vec![255,0,255]);
                }
                ggeo.add_objects(mm);
            }

            let mut mm = geometry::Geometry::generate_sphere(&cbb.get_xyz(),0.25,8,8);
            println!("{:?}",&cbb);
            for f in mm.1.iter_mut(){
                f.set_color(&vec![255,0,0]);
            }
            ggeo.add_objects(mm);


            for a in vec![&(0.822,-1.200,-0.000),&(0.000,0.000,0.000),&(-0.949,-0.005,1.177),&(0.860,1.255,0.000)]{
                let mut mm = geometry::Geometry::generate_sphere(a,0.25,8,8);
                for f in mm.1.iter_mut(){
                    f.set_color(&vec![255,255,255]);
                }
                ggeo.add_objects(mm);
            }


            ggeo.calc_all_norms();
            ggeo.add_colortile_material();
            ggeo.calc_all_norms();
            let mut gv:Vec<geometry::Geometry> = vec![ggeo];
        
            geometry::Geometry::save("test/anglecheck.obj",&mut gv);
        }


    }
    pub fn is_in(&self,atom:&dyn Vector3D){
        
    }
}




#[test]
fn coarse_grained_test(){
    let n:PseudoAtom = PseudoAtom::new(0.822,-1.200,-0.000);
    let ca:PseudoAtom = PseudoAtom::new(0.000,0.000,0.000);
    let cb:PseudoAtom = PseudoAtom::new(-0.949,-0.005,1.177);
    let c:PseudoAtom = PseudoAtom::new(0.860,1.255,0.000);
    let cyl = SideChainCylinder::new(&n,&ca,&c);   
}
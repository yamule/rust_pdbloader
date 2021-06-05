use core::f64;
use std::fs;
use std::cmp::Ordering;
use std::io::{BufWriter,Write};
use regex::Regex;
use std::collections::{HashMap,HashSet};
use crate::process_3d::standardize;
use crate::{geometry::{Point3D, Vector3D}, mmcif_process, pdbdata::{PDBAsym, PDBAtom,PDBComp, PDBEntry}, process_3d, structural_alignment::align};
use super::structural_alignment;
use super::matrix_process;
use super::geometry;
use chrono::{Local};
use std::f64::consts::PI;
use rand::prelude::*;

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
    pub radius:f64,
    pub rx:Point3D,//ry と rz の法線であり、nc center->c のベクトル
    pub ry:Point3D,//n ca c の norm を CB 分回転したベクトル 
    pub rz:Point3D,//CA →見積もった CB へのベクトル
    pub length:f64,//CA から最終点の距離
    pub rz_length:Point3D,
    pub num_sep:usize,//length を何個に分割するか//90 度ごと 4*2 つに分割するので、領域はこれ x8
}
impl<'a> SideChainCylinder<'a>{
    pub fn new(n:&'a PseudoAtom,ca:&'a PseudoAtom,c:&'a PseudoAtom,llen:f64,ssep:usize,radi:f64)->SideChainCylinder<'a>{
        let mut ret = SideChainCylinder{
        n:n,
        ca:ca,
        c:c,
        radius:radi,
        rx:Point3D::new(0.0,0.0,0.0),//n ca c の norm と rz の norm
        ry:Point3D::new(0.0,0.0,0.0),//n ca c の norm
        rz:Point3D::new(0.0,0.0,0.0),//ca->cb の単位ベクトル
        rz_length:Point3D::new(0.0,0.0,0.0),//ca->cb を length 分伸ばした線分。他原子の評価の際に毎回計算する必要があるのであらかじめ計算しておく
        length:llen,//CA から最終点の距離
        num_sep:ssep
        };
        ret.update();
        return ret;
    }

    pub fn generate_obj(&self)->Vec<(Vec<Point3D>,Vec<geometry::Face>)>{
        let st:(f64,f64,f64) = self.ca.get_xyz();
        let en:(f64,f64,f64) = (
            st.0+self.rz.get_x()*self.length
            ,st.1+self.rz.get_y()*self.length
            ,st.2+self.rz.get_z()*self.length
        );
        let mut ret:Vec<(Vec<Point3D>,Vec<geometry::Face>)> = vec![];
        ret.push(geometry::Geometry::generate_cylinder(&st,&en,self.radius,8,false));
        ret.push(geometry::Geometry::generate_sphere(&self.ca.get_xyz(),0.5,8,8));
        ret.push(geometry::Geometry::generate_sphere(&self.n.get_xyz(),0.5,8,8));
        ret.push(geometry::Geometry::generate_sphere(&self.c.get_xyz(),0.5,8,8));
        return ret;        
    }

    pub fn update(&mut self){
        let mut nn = self.n.get_xyz();
        let caa = self.ca.get_xyz();
        let mut cc = self.c.get_xyz();
        nn = process_3d::standardize(nn.0-caa.0,nn.1-caa.1,nn.2-caa.2);
        cc = process_3d::standardize(cc.0-caa.0,cc.1-caa.1,cc.2-caa.2);


        
        //CHARMM の角度と二面角と結合長から計算したら微妙に違ったので
        //もう適当でいいやという気分になった
        //C-CA-N から CB 位置を計算する
        //法線方向の距離
        let pyy = (18.5*to_radian).cos()*((60.0*to_radian).sin());
        //standardize された nc 中心点逆方向ベクトルの距離
        let pzz =  (1.0-pyy*pyy).sqrt();


        let yy:(f64,f64,f64) = process_3d::calc_norm_t(
            &nn,
            &(0.0,0.0,0.0),
            &cc
        );
        let nccenter = (nn.0/2.0+cc.0/2.0,nn.1/2.0+cc.1/2.0,nn.2/2.0+cc.2/2.0);
        self.rx.set_xyz(cc.0-nccenter.0,cc.1-nccenter.1,cc.2-nccenter.2);
        self.rx.standardize();
        

        let mut pcx:(f64,f64,f64) = standardize(nn.0/2.0+cc.0/2.0,nn.1/2.0+cc.1/2.0,nn.2/2.0+cc.2/2.0);
        pcx.0 *= -1.0;
        pcx.1 *= -1.0;
        pcx.2 *= -1.0;

        self.rz.set_xyz(
            yy.0*-pyy+pcx.0*pzz
            ,yy.1*-pyy+pcx.1*pzz
            ,yy.2*-pyy+pcx.2*pzz
        );
        self.rz.standardize();

        let ynorm = process_3d::calc_norm(
            self.rz.get_x(),
            self.rz.get_y(),
            self.rz.get_z(),
            self.rx.get_x(),
            self.rx.get_y(),
            self.rx.get_z()
        );
        self.ry.set_xyz(ynorm.0,ynorm.1,ynorm.2);
        

        self.rz_length.set_xyz(
            self.rz.get_x()*self.length,
            self.rz.get_y()*self.length,
            self.rz.get_z()*self.length
        );
        if false{//デバッグ用コード
            
            
            let candidate1:(f64,f64,f64)
            = (
                (yy.0*-pyy+pcx.0*pzz)*1.52
                ,
                (yy.1*-pyy+pcx.1*pzz)*1.52
                ,
                (yy.2*-pyy+pcx.2*pzz)*1.52
            );
        
            //CHARMM の Force field param から CB 位置を計算する
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
            //bond は
            //CH1E CH3E   225.0       1.52
            //N-CA-C は
            //C    CH1E NH1     45.0     111.6
            //ALA の IMPROPER (toph19)は 
            //IC   N    C    *CA  CB     0.0000    0.00  120.00    0.00   0.0000
            //値を入れて計算してみると少しずれる。。。
            let xx:(f64,f64,f64) = process_3d::calc_norm_t(
                &zz,
                &(0.0,0.0,0.0),
                &yy
            );

            let mut cbb = Point3D::from_tuple(&cc);
            cbb.standardize();
            cbb.set_x(cbb.get_x()* 1.52);
            cbb.set_y(cbb.get_y()* 1.52);
            cbb.set_z(cbb.get_z()* 1.52);

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
                println!("cand1: {} {:?}",process_3d::distance(
                    &(-0.949,-0.005,1.177)
                    ,&candidate1
                ),candidate1);
                println!("cand2: {} {:?}",process_3d::distance(
                    &(-0.949,-0.005,1.177)
                    ,&cbb.get_xyz()
                ),&cbb.get_xyz());

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
    }

    //分割中の近い位置を 0 として何番目の分割内にあるか
    //,
    //方向
    //7 3 0 4
    //6 2 1 5
    //で返す。
    pub fn get_position_of(&self,atom:&dyn Vector3D)->Option<(u8,u8)>{
        if atom.distance(self.ca).powf(2.0) > self.length*self.length + self.radius*self.radius{
            return None;
        }
        let mut ppos = atom.get_xyz();
        ppos.0 -= self.ca.get_x();
        ppos.1 -= self.ca.get_y();
        ppos.2 -= self.ca.get_z();
        
        let a:f64 = process_3d::distance(&ppos,&self.rz_length.get_xyz());
        let b:f64 = process_3d::distance(&ppos,&(0.0,0.0,0.0));
        if b*b >= a*a+self.length*self.length || a*a >= b*b+self.length*self.length{            
            return None;
        }

        let d = (b*b-a*a+self.length*self.length)/(2.0*self.length);
        assert!(d > 0.0);
        let ratio = d/self.length;
        let rpos:(f64,f64,f64) = (
            self.rz_length.get_x()*ratio,
            self.rz_length.get_y()*ratio,
            self.rz_length.get_z()*ratio
        );

        //ppos を原点中心の相対的な位置に移動
        ppos.0 -= rpos.0;
        ppos.1 -= rpos.1;
        ppos.2 -= rpos.2;

        let mut ddist = ppos.0*ppos.0+ppos.1*ppos.1+ppos.2*ppos.2;//中心点からの距離
        
        if ddist > 0.0{
            ddist = ddist.sqrt();
            ppos.0 /= ddist;
            ppos.1 /= ddist;
            ppos.2 /= ddist;
        }//ppos は単位ベクトルに直した
        if ddist > self.radius {
            return None;
        }
        let bdis = process_3d::distance(
            &ppos
            ,&(self.rx.get_xyz())
        );
        let mut rad_1 = (1.0-bdis*bdis/2.0).acos();
        
        let cdis = process_3d::distance(&ppos,&self.ry.get_xyz());

        if cdis > (2.0_f64).sqrt(){
            rad_1 *= -1.0;
        }
        let direc_code:u8 = 
        if rad_1 < PI/2.0 && rad_1 >= 0.0{
            0
        }else if rad_1 >= -1.0*PI/2.0 && rad_1 < 0.0{
            1
        }else if rad_1 < -1.0*PI/2.0 && rad_1 >= -1.0*PI{
            2
        }else{
            3
        };
        let direc_code:u8 = if ddist < self.radius/2.0{
            direc_code
        }else{
            direc_code +3
        };
        let rstep = 1.0/(self.num_sep as f64);
        let mut pcode:u8 = 0;
        for ii in 1..=self.num_sep{
            if rstep*(ii as f64) > ratio{
                pcode = ii as u8;
                break;
            }
        }
        return Some((pcode,direc_code));
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
    let cyl = SideChainCylinder::new(&n,&ca,&c,8.0,3,6.0);
    let obj = cyl.generate_obj();
    let mut geom:geometry::Geometry = geometry::Geometry::new();
    for bb in obj.into_iter(){
        geom.add_objects(bb);
    }
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(10);
    
    for _ in 0..500{
        let spos:(f64,f64,f64) = (rgen.gen_range(-10.0,10.0),rgen.gen_range(-10.0,10.0),rgen.gen_range(-10.0,10.0));
        let mut spp = 
        geometry::Geometry::generate_sphere(&spos,0.1,8,8);
        let res = cyl.get_position_of(&Point3D::new(spos.0,spos.1,spos.2));
        if let Some(x) = res{
            //let r = 80*x.0;
            let g = 36*x.1;
            g が間違っていそう
            let r = 0;
            geometry::Face::color_faces(&mut spp.1,&vec![r,g,0]);
        }
        geom.add_objects(spp);
    }


    geom.add_objects(geometry::Geometry::generate_sphere(&cb.get_xyz(),0.5,8,8));

    geom.calc_all_norms();
    geom.add_colortile_material();
    geom.calc_all_norms();

    geometry::Geometry::save("test/cylindercheck.obj",&mut vec![geom]);    
}

#[test]
fn coarse_grained_test2(){
    let x = (18.5*to_radian).cos()*((60.0*to_radian).cos());
    let y = (18.5*to_radian).sin();
    let z = (18.5*to_radian).cos()*((60.0*to_radian).sin());
    let p =  (1.0-z*z).sqrt();
    println!("{} {} {}"
    ,x*(PI+(180.0-111.6)*to_radian/-2.0).cos()-y*(PI+(180.0-111.6)*to_radian/-2.0).sin()
    ,y*(PI+(180.0-111.6)*to_radian/-2.0).cos()+x*(PI+(180.0-111.6)*to_radian/-2.0).sin()
    ,z);
    println!("{} {} {}"
    ,x
    ,y
    ,z);
}
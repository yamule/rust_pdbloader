
use super::geometry::*;

#[allow(unused_imports)]
use std::f64::consts::PI;


#[allow(unused_imports)]
use rand::SeedableRng;
#[allow(unused_imports)]
use rand::rngs::StdRng;
#[allow(unused_imports)]
use rand::Rng;


/**
 * 2ベクトルの外積を返す
    * @param v1x
    * @param v1y
    * @param v1z
    * @param v2x
    * @param v2y
    * @param v2z
    * @return 
    */
pub fn calc_norm(v1x:f64,v1y:f64,v1z:f64,v2x:f64,v2y:f64,v2z:f64)->(f64,f64,f64){
    return standardize(
    (v1y*v2z-v1z*v2y)*-1.0,
    (v1z*v2x-v1x*v2z)*-1.0,
    (v1x*v2y-v1y*v2x)*-1.0);
}
//_t は tuple の意味
pub fn calc_norm_t(
    p0:&(f64,f64,f64),
    p1:&(f64,f64,f64),
    p2:&(f64,f64,f64)
)->(f64,f64,f64){
    return calc_norm(
        p0.0-p1.0,
        p0.1-p1.1,
        p0.2-p1.2,
        p2.0-p1.0,
        p2.1-p1.1,
        p2.2-p1.2
    );
}


pub fn standardize(xx:f64,yy:f64,zz:f64) -> (f64,f64,f64){
    let mut llen:f64 = xx*xx+yy*yy+zz*zz;
    if llen > 0.0{
        llen = llen.sqrt();
        return (xx/llen,yy/llen,zz/llen);
    }
    return (0.0,0.0,0.0);
}

pub fn standardize_e(ll:&mut dyn Vector3D){
    let r = standardize(ll.get_x(),ll.get_y(),ll.get_z());
    ll.set_xyz(r.0,r.1,r.2);
}

pub fn distance(xyz_a:&(f64,f64,f64),xyz_b:&(f64,f64,f64))->f64{
    let x = xyz_a.0-xyz_b.0;
    let y = xyz_a.1-xyz_b.1;
    let z = xyz_a.2-xyz_b.2;
    let r:f64 = x*x+y*y+z*z;
    if r > 0.0{
        return r.sqrt();
    }
    return 0.0;
}

pub fn calc_norm_v(v0:&Point3D,v1:&Point3D,v2:&Point3D)->(f64,f64,f64){
    return calc_norm(
        v0.get_x()-v1.get_x(),
        v0.get_y()-v1.get_y(),
        v0.get_z()-v1.get_z(),
        v2.get_x()-v1.get_x(),
        v2.get_y()-v1.get_y(),
        v2.get_z()-v1.get_z());
}



pub fn subtracted_t(a:&(f64,f64,f64),b:&(f64,f64,f64))->(f64,f64,f64){
    return (a.0-b.0,a.1-b.1,a.2-b.2);
}

pub fn added_t(a:&(f64,f64,f64),b:&(f64,f64,f64))->(f64,f64,f64){
    return (a.0+b.0,a.1+b.1,a.2+b.2);
}



/*
* get_xzero_rot_angle の結果を rotateTo に使えるように目的のベクトルが X 座標以外 0 となるような回転行列を返す
*/
pub fn rotate_xzero(xx:f64,yy:f64,zz:f64,angle:(f64,f64,f64,f64))->(f64,f64,f64){
    let rx = xx*angle.0+zz*angle.1;
    let rz = -xx*angle.1+zz*angle.0;

    let px = rx*angle.2+yy*angle.3;
    let py = -rx*angle.3+yy*angle.2;
    return (px,py,rz);
}


/*
* get_xzero_rot_angle で取られた結果をある三次元上の点に適用した結果を返す
*/
pub fn rotate_with(xx:f64,yy:f64,zz:f64,angle:(f64,f64,f64,f64))->(f64,f64,f64){
    let px = xx*angle.2-yy*angle.3;
    let py = xx*angle.3+yy*angle.2;

    let rx = px*angle.0-zz*angle.1;
    let rz = px*angle.1+zz*angle.0;
    return (rx,py,rz);
}


/**
 * 
    * t1<-t2->t3 というベクトルが作る三角形を s1<-s2->s3 というベクトルが作る三角形に合うように vset を回転する
    * s1 が X 方向になるように計算している。
    * s2 と t2 が同じ位置になっていると思う。
    * t1 や s1 等が vset に入っていると結果がおかしくなるが、その場合ボローチェッカーを通らないはず。
    * @param d1
    * @param d2
    * @param d3
    * @param s1
    * @param s2
    * @param s3
    * @param vset 
    */

pub fn fit_to_vector(t1:&dyn Vector3D,t2:&dyn Vector3D,t3:&dyn Vector3D
    ,s1:&dyn Vector3D,s2:&dyn Vector3D,s3:&dyn Vector3D,
    vset_:&mut Vec<&mut dyn Vector3D>){


    let angle:(f64,f64,f64,f64) = get_xzero_rot_angle(
    s1.get_x()-s2.get_x()
    ,s1.get_y()-s2.get_y()
    ,s1.get_z()-s2.get_z());
    let sst:(f64,f64,f64) = rotate_xzero(
        s3.get_x()-s2.get_x()
        ,s3.get_y()-s2.get_y()
        ,s3.get_z()-s2.get_z()
        ,angle);
    
    let mut vc = Point3D::new(0.0,sst.1,sst.2);
    vc.standardize();

    let cos1 = vc.get_z();
    let sin1 = vc.get_y();

    
    let angle2 = get_xzero_rot_angle(
    t1.get_x()-t2.get_x()
    ,t1.get_y()-t2.get_y()
    ,t1.get_z()-t2.get_z());
    
    let dd = rotate_xzero(
    t3.get_x()-t2.get_x()
    ,t3.get_y()-t2.get_y()
    ,t3.get_z()-t2.get_z(),angle2);
    
    let mut vcd = Point3D::new(0.0,dd.1,dd.2);
    vcd.standardize();
    
    let cos2 = vcd.get_z();
    let sin2 = vcd.get_y();
    
    let mut c1 = Point3D::from_vector3d(t1);
    let mut c2 = Point3D::from_vector3d(t2);
    let mut c3 = Point3D::from_vector3d(t3);
    let mut vset:Vec<&mut dyn Vector3D> = vec![];
    for v in vset_.iter_mut(){
        vset.push(*v);
    }
    vset.push(&mut c1);
    vset.push(&mut c2);
    vset.push(&mut c3);
    for v in vset.iter_mut(){
        v.set_xyz(v.get_x()-t2.get_x(),v.get_y()-t2.get_y(),v.get_z()-t2.get_z());

        let mut dtmp1 = rotate_xzero(v.get_x(),v.get_y(),v.get_z(),angle2);
        let mut vcdd = Point3D::new(dtmp1.0,dtmp1.1,dtmp1.2);
        //println!("{:?}",vcdd);
        dtmp1.1 = vcdd.get_y()*cos2-vcdd.get_z()*sin2;
        dtmp1.2 = vcdd.get_y()*sin2+vcdd.get_z()*cos2;

        //println!("{:?}",dtmp1);
        vcdd.set_y(dtmp1.1);
        vcdd.set_z(dtmp1.2);
        
        dtmp1.1 = vcdd.get_y()*cos1+vcdd.get_z()*sin1;
        dtmp1.2 = vcdd.get_z()*cos1-vcdd.get_y()*sin1;
        
        let l = rotate_with(dtmp1.0,dtmp1.1,dtmp1.2,angle);
        //v.set_xyz(l.0+s2.get_x(),l.1+s2.get_y(),l.2+s2.get_z());
        v.set_xyz(l.0,l.1,l.2);
    }

    let snorm0 = calc_norm_t(
        &s1.get_xyz(),
        &s2.get_xyz(),
        &s3.get_xyz(),
    );
    let ps1 = standardize(
    s1.get_x()-s2.get_x(),
    s1.get_y()-s2.get_y(),
    s1.get_z()-s2.get_z());
    
    let snorm1 = calc_norm_t(  
    &ps1,
        &(0.0,0.0,0.0),
        &snorm0
    );


    //全部 0 起点で standardize されている必要がある
    let get_radian = |
    p:&(f64,f64,f64)
    ,xaxis:&(f64,f64,f64)
    ,yaxis:&(f64,f64,f64)
    |->f64{
        let d:f64 = distance(xaxis, p);
        let radp = ((2.0-d*d)/2.0).acos();
        let zd:f64 = distance(yaxis, p);
        if zd > (2.0_f64).sqrt(){
            return radp*-1.0;
        }
        return radp;
    };
    c1.subtract(&c2.get_xyz());
    c3.subtract(&c2.get_xyz());
    c1.standardize();
    c3.standardize();
    let mut sc1 = Point3D::from_vector3d(s1);
    sc1.subtract(&s2.get_xyz());
    sc1.standardize();
    let mut sc3 = Point3D::from_vector3d(s3);
    sc3.subtract(&s2.get_xyz());
    sc3.standardize();
    
    let rad1 = get_radian(&c1.get_xyz(),&ps1,&snorm1);
    let rad1b = get_radian(&sc1.get_xyz(),&ps1,&snorm1);

    let rad3 = get_radian(&c3.get_xyz(),&ps1,&snorm1);
    let rad3b = get_radian(&sc3.get_xyz(),&ps1,&snorm1);
    
    let rot = (rad3b-rad3)/2.0+(rad1b-rad1)/2.0;
    
    let mut vset:Vec<&mut dyn Vector3D> = vec![];
    for v in vset_.iter_mut(){
        //v.set_xyz(v.get_x()-s2.get_x(), v.get_y()-s2.get_y(), v.get_z()-s2.get_z());
        vset.push(*v);
    }
    rotate_3d(&mut vset,&(Point3D::from_tuple(&snorm0)),rot);
    for v in vset_.iter_mut(){
        v.set_xyz(v.get_x()+s2.get_x(), v.get_y()+s2.get_y(), v.get_z()+s2.get_z());
    }
    
    
    //println!("{} {} {} {:?} {:?} {:?}",rad1,rad3,rad3b,sc3,ps1,snorm1);
    
}


//t1->t2 と s1->s2 が合う昔のコード
pub fn fit_to_vector_old(t1:&dyn Vector3D,t2:&dyn Vector3D,t3:&dyn Vector3D
    ,s1:&dyn Vector3D,s2:&dyn Vector3D,s3:&dyn Vector3D,
    vset_:&mut Vec<&mut dyn Vector3D>){


    let angle:(f64,f64,f64,f64) = get_xzero_rot_angle(
    s1.get_x()-s2.get_x()
    ,s1.get_y()-s2.get_y()
    ,s1.get_z()-s2.get_z());
    let sst:(f64,f64,f64) = rotate_xzero(
        s3.get_x()-s2.get_x()
        ,s3.get_y()-s2.get_y()
        ,s3.get_z()-s2.get_z()
        ,angle);
    
    let mut vc = Point3D::new(0.0,sst.1,sst.2);
    vc.standardize();

    let cos1 = vc.get_z();
    let sin1 = vc.get_y();

    
    let angle2 = get_xzero_rot_angle(
    t1.get_x()-t2.get_x()
    ,t1.get_y()-t2.get_y()
    ,t1.get_z()-t2.get_z());
    
    let dd = rotate_xzero(
    t3.get_x()-t2.get_x()
    ,t3.get_y()-t2.get_y()
    ,t3.get_z()-t2.get_z(),angle2);
    
    let mut vcd = Point3D::new(0.0,dd.1,dd.2);
    vcd.standardize();
    
    let cos2 = vcd.get_z();
    let sin2 = vcd.get_y();
    
    let mut c1 = Point3D::from_vector3d(t1);
    let mut c2 = Point3D::from_vector3d(t2);
    let mut c3 = Point3D::from_vector3d(t3);
    let mut vset:Vec<&mut dyn Vector3D> = vec![];
    for v in vset_.iter_mut(){
        vset.push(*v);
    }
    vset.push(&mut c1);
    vset.push(&mut c2);
    vset.push(&mut c3);
    for v in vset.iter_mut(){
        v.set_xyz(v.get_x()-t2.get_x(),v.get_y()-t2.get_y(),v.get_z()-t2.get_z());

        let mut dtmp1 = rotate_xzero(v.get_x(),v.get_y(),v.get_z(),angle2);
        let mut vcdd = Point3D::new(dtmp1.0,dtmp1.1,dtmp1.2);
        //println!("{:?}",vcdd);
        dtmp1.1 = vcdd.get_y()*cos2-vcdd.get_z()*sin2;
        dtmp1.2 = vcdd.get_y()*sin2+vcdd.get_z()*cos2;

        //println!("{:?}",dtmp1);
        vcdd.set_y(dtmp1.1);
        vcdd.set_z(dtmp1.2);
        
        dtmp1.1 = vcdd.get_y()*cos1+vcdd.get_z()*sin1;
        dtmp1.2 = vcdd.get_z()*cos1-vcdd.get_y()*sin1;
        
        let l = rotate_with(dtmp1.0,dtmp1.1,dtmp1.2,angle);
        v.set_xyz(l.0+s2.get_x(),l.1+s2.get_y(),l.2+s2.get_z());
    }
    //println!("{} {} {} {:?} {:?} {:?}",rad1,rad3,rad3b,sc3,ps1,snorm1);
    
}
    
/**
 * query1->query2 というベクトルを template1->template2 というベクトルに合うように vset を回転する
 * query1 等が vset に入っていると結果がおかしくなるが、その場合ボローチェッカーを通らないはず。
 * @param dest1
 * @param dest2
 * @param src1
 * @param src2
 * @param vset 
*/
pub fn fit_to_vector_2(query1:&dyn Vector3D,query2:&dyn Vector3D
    ,template1:&dyn Vector3D,template2:&dyn Vector3D
    ,vset:&mut Vec<&mut dyn Vector3D>){
    let angle:(f64,f64,f64,f64) = get_xzero_rot_angle(template2.get_x()-template1.get_x(),template2.get_y()-template1.get_y(),template2.get_z()-template1.get_z());
    let angle2:(f64,f64,f64,f64) = get_xzero_rot_angle(query2.get_x()-query1.get_x(),query2.get_y()-query1.get_y(),query2.get_z()-query1.get_z());

    
    for v in vset.iter_mut(){
        v.set_xyz(v.get_x()-query1.get_x(),v.get_y()-query1.get_y(),v.get_z()-query1.get_z());
        let d1:(f64,f64,f64) = rotate_xzero(v.get_x(),v.get_y(),v.get_z(),angle2);
        let dd:(f64,f64,f64) = rotate_with(d1.0,d1.1,d1.2,angle);
        v.set_xyz(dd.0,dd.1,dd.2);
        v.set_xyz(v.get_x()-template1.get_x(),v.get_y()-template1.get_y(),v.get_z()-template1.get_z());
    }
}

    
/**
 * あるベクトルの Y 軸中心の回転 sin, cos, Z 軸中心の回転 sin, cos を得る
* x 軸方向にのみ値が入っているベクトルを掛けると元の線分と同じベクトルになる
* 
*/
pub fn get_xzero_rot_angle(xx:f64,yy:f64,zz:f64)->(f64,f64,f64,f64){
    let mut lenxz:f64 = xx*xx+zz*zz;
    let mut ysin:f64 = 0.0;
    let mut ycos:f64 = 1.0;
    
    
    let mut px:f64 = xx;
    let py:f64 = yy;
    let mut _pz:f64 = zz;
    if lenxz > 0.0{
        lenxz = lenxz.sqrt();
        ycos = xx/lenxz;
        ysin = zz/lenxz;
        px = xx*ycos+zz*ysin;
        _pz = 0.0;
    }
    
        
    let mut zsin:f64 = 0.0;
    let mut zcos:f64 = 1.0;
    
    let mut lenyy:f64 = px*px+py*py;
    if lenyy > 0.0{
        lenyy = lenyy.sqrt();
        zsin = py/lenyy;
        zcos = px/lenyy;
    }
    return (ycos,ysin,zcos,zsin);
}



//axisvec を軸として radian 度 回転する
//axisvec は standardize されている必要がある
//https://ja.wikipedia.org/wiki/%E3%83%AD%E3%83%89%E3%83%AA%E3%82%B2%E3%82%B9%E3%81%AE%E5%9B%9E%E8%BB%A2%E5%85%AC%E5%BC%8F
pub fn rotate_3d(target:&mut Vec<&mut dyn Vector3D>,axisvec:&dyn Vector3D,radian:f64){
    let cc = radian.cos();
    let ss = radian.sin();
    let nx = axisvec.get_x();
    let ny = axisvec.get_y();
    let nz = axisvec.get_z();

    let vv:Vec<Vec<f64>> = vec![
         vec![cc+nx*nx*(1.0-cc)     ,nx*ny*(1.0-cc)-nz*ss   ,nx*nz*(1.0-cc)+ny*ss   ]
        ,vec![ny*nx*(1.0-cc)+nz*ss  ,cc+ny*ny*(1.0-cc)      ,ny*nz*(1.0-cc)-nx*ss   ]
        ,vec![nz*nx*(1.0-cc)-ny*ss  ,nz*ny*(1.0-cc)+nx*ss   ,cc+nz*nz*(1.0-cc)      ]
    ];
    
    for tt in target.iter_mut(){
        let x = tt.get_x();
        let y = tt.get_y();
        let z = tt.get_z();
        tt.set_xyz(
            x*vv[0][0] + y*vv[0][1] + z*vv[0][2],
            x*vv[1][0] + y*vv[1][1] + z*vv[1][2],
            x*vv[2][0] + y*vv[2][1] + z*vv[2][2]
        );
    }
}




#[test]
fn testrotate(){
    let mut tmpp = Point3D::new((PI/6.0).cos(),0.5,0.1);
    let axx = Point3D::new(0.0,0.0,1.0);
    rotate_3d(&mut vec![&mut tmpp], &axx,PI/3.0);
    
    assert_eq!((tmpp.get_x()*100.0+0.1) as i64,0 as i64);
    assert_eq!((tmpp.get_y()*100.0+0.1) as i64,100 as i64);
    assert_eq!((tmpp.get_z()*100.0+0.1) as i64,10 as i64);

    let mut tmpp = Point3D::new((PI/6.0).cos(),0.1,0.5);
    let axx = Point3D::new(0.0,1.0,0.0);
    rotate_3d(&mut vec![&mut tmpp], &axx,-1.0*PI/3.0);//右手系なので Y だけ逆？
    
    assert_eq!((tmpp.get_x()*100.0+0.1) as i64,0 as i64);
    assert_eq!((tmpp.get_y()*100.0+0.1) as i64,10 as i64);
    assert_eq!((tmpp.get_z()*100.0+0.1) as i64,100 as i64);
    
    let mut tmpp = Point3D::new(0.1,(PI/6.0).cos(),0.5);
    let axx = Point3D::new(1.0,0.0,0.0);
    rotate_3d(&mut vec![&mut tmpp], &axx,PI/3.0);
    assert_eq!((tmpp.get_x()*100.0+0.1) as i64,10 as i64);
    assert_eq!((tmpp.get_y()*100.0+0.1) as i64,0 as i64);
    assert_eq!((tmpp.get_z()*100.0+0.1) as i64,100 as i64);
    
}

#[test]
fn testrotate2(){
    let mut rgen:StdRng = SeedableRng::seed_from_u64(123);
    for _ in 0..100{
        let v1 = Point3D::new(rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        let v2 = Point3D::new(rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        let v3 = Point3D::new(rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        let v4 = Point3D::new(rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));

        let mut v1b = v1.get_copy();
        let mut v2b = v2.get_copy();
        let mut v3b = v3.get_copy();
        let mut v4b = v4.get_copy();
        let mut rott = Point3D::new(rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0),rgen.gen_range(-100.0,100.0));
        rott.standardize();
        rotate_3d(&mut vec![&mut v1b,&mut v2b,&mut v3b,&mut v4b],&rott,rgen.gen_range(0.0,10.0));


        
        assert!((v1b.distance(&v2b) - v1.distance(&v2)).abs() < 0.0001);
        assert!((v1b.distance(&v3b) - v1.distance(&v3)).abs() < 0.0001);
        assert!((v3b.distance(&v2b) - v3.distance(&v2)).abs() < 0.0001);


        let mut v1c = v1b.get_copy();
        let mut v2c = v2b.get_copy();
        let mut v3c = v3b.get_copy();
        let mut v4c = v4b.get_copy();
    
        
        assert!((v1b.distance(&v2b) - v1c.distance(&v2c)).abs() < 0.0001);
        assert!((v1b.distance(&v3b) - v1c.distance(&v3c)).abs() < 0.0001);
        assert!((v3b.distance(&v2b) - v3c.distance(&v2c)).abs() < 0.0001);

        fit_to_vector(&v1b,&v2b,&v3b,
        &v1,&v2,&v3,
        &mut vec![&mut v3c,&mut v4c,&mut v2c,&mut v1c]);
        
        //println!("{:?} {:?} {:?} {} ",v1,v1b,v1c,v1c.distance(&v1));
        //println!("{:?} {:?} {:?} {} ",v2,v2b,v2c,v2c.distance(&v2));
        //println!("{:?} {:?} {:?} {} ",v3,v3b,v3c,v3c.distance(&v3));
        //println!("{:?} {:?} {:?} {} ",v4,v4b,v4c,v4c.distance(&v4));
        
        assert!((v1b.distance(&v2b) - v1c.distance(&v2c)).abs() < 0.0001);
        assert!((v1b.distance(&v3b) - v1c.distance(&v3c)).abs() < 0.0001);
        assert!((v3b.distance(&v2b) - v3c.distance(&v2c)).abs() < 0.0001);

        assert!(v1c.distance(&v1) < 0.0001);
        assert!(v2c.distance(&v2) < 0.0001);
        assert!(v3c.distance(&v3) < 0.0001);
        assert!(v4c.distance(&v4) < 0.0001);
    }

}
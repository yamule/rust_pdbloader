#[allow(dead_code,unused_imports)]
use std::io::{BufWriter,Write,BufReader,BufRead};
#[allow(dead_code,unused_imports)]
use std::collections::HashMap;
#[allow(dead_code,unused_imports)]
use std::fs;
#[allow(dead_code,unused_imports)]
use super::process_3d::*;

extern crate regex;
use regex::Regex;

fn write_to_file(filename:&str,contents:Vec<String>){
    let mut f = BufWriter::new(fs::File::create(filename).unwrap());
     for ll in contents{
        f.write_all(ll.as_bytes()).unwrap();
        f.write_all("\n".as_bytes()).unwrap();
    }
}

#[derive(Debug)]
pub struct Point2D{
    pub x:f64,
    pub y:f64
}

impl Point2D{
    
    fn to_string(&self,flag32:bool)->String{
        if flag32{
            return format!("{} {}",self.x as f32,self.y as f32);

        }else{
            return format!("{} {}",self.x,self.y);
        }
    }
}
pub trait Vector3D{
    fn get_x(&self)->f64;
    fn get_y(&self)->f64;
    fn get_z(&self)->f64;
    fn get_xyz(&self)->(f64,f64,f64);

    fn set_x(&mut self,xx:f64);
    fn set_y(&mut self,yy:f64);
    fn set_z(&mut self,zz:f64);
    fn set_xyz(&mut self,xx:f64,yy:f64,zz:f64);
    fn set_vector(&mut self,v:&dyn Vector3D);
    fn distance(&self,p:&dyn Vector3D)->f64;
}

#[derive(Debug,Clone,Copy)]
pub struct Point3D{
    pub x:f64,
    pub y:f64,
    pub z:f64
}

impl Point3D{
    pub fn new(x:f64,y:f64,z:f64)->Point3D{
        return Point3D{x,y,z};
    }
    pub fn new_tup(xyz:&(f64,f64,f64))->Point3D{
        return Point3D::new(xyz.0,xyz.1,xyz.2);
    }

    pub fn from_vector3d(v:&dyn Vector3D)->Point3D{
        return Point3D{x:v.get_x(),y:v.get_y(),z:v.get_z()}; 
    }

    pub fn from_tuple(v:&(f64,f64,f64))->Point3D{
        return Point3D{x:v.0,y:v.1,z:v.2}; 
    }
    pub fn zero()->Point3D{
        return Point3D::new(0.0,0.0,0.0);
    }
    
    pub fn get_copy(&self)->Point3D{
        return Point3D::new(self.get_x(),self.get_y(),self.get_z());
    }

    pub fn add(&mut self,v:&(f64,f64,f64)){
        self.x += v.0;
        self.y += v.1;
        self.z += v.2;
    }

    pub fn subtract(&mut self,v:&(f64,f64,f64)){
        self.x -= v.0;
        self.y -= v.1;
        self.z -= v.2;
    }

    pub fn multiply(&mut self,m:f64){
        self.x *= m;
        self.y *= m;
        self.z *= m;
    }
    pub fn standardize(&mut self){
        standardize_e(self);
    }
    fn to_string(&self,flag32:bool)->String{
        if flag32{
            return format!("{} {} {}",self.x as f32,self.y as f32,self.z as f32);

        }else{
            return format!("{} {} {}",self.x,self.y,self.z);
        }
    }
} 
impl Vector3D for Point3D{
    fn set_xyz(&mut self,xx:f64,yy:f64,zz:f64){
        self.set_x(xx);
        self.set_y(yy);
        self.set_z(zz);
    }
    
    fn set_vector(&mut self,v:&dyn Vector3D){
        self.set_x(v.get_x());
        self.set_y(v.get_y());
        self.set_z(v.get_z());
    }

    fn get_xyz(&self)-> (f64,f64,f64){
        return (self.get_x(),self.get_y(),self.get_z());
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
    
    fn distance(&self,p:&dyn Vector3D)->f64{
        return distance(
        &(self.get_x(),self.get_y(),self.get_z())
        ,&(p.get_x(),p.get_y(),p.get_z())
        );
    }
}
pub struct Face{
    index_v:Vec<usize>,
    index_vt:Vec<usize>,
    index_vn:Vec<usize>,
    norm:Option<Point3D>,
    dead_flag:bool
}

impl Face{
    pub fn new()->Face{
        return Face{index_v:vec![],
        index_vt:vec![],
        index_vn:vec![],
        norm:None,
        dead_flag:false,
        };
    }
    pub fn calc_center_point(&self,parentgeom:&Geometry)->Point3D{
        let v0:&Point3D = &parentgeom.vertices[self.index_v[0]];
        let v1:&Point3D = &parentgeom.vertices[self.index_v[1]];
        let v2:&Point3D = &parentgeom.vertices[self.index_v[2]];
        let mut ret:Point3D = Point3D::zero();
        ret.add(&v0.get_xyz());
        ret.add(&v1.get_xyz());
        ret.add(&v2.get_xyz());
        ret.set_xyz(ret.get_x()/3.0,ret.get_y()/3.0,ret.get_z()/3.0);
        return ret;
    }
    pub fn calc_norm(&self,parentgeom:&Geometry)->Point3D{
        let v0:&Point3D = &parentgeom.vertices[self.index_v[0]];
        let v1:&Point3D = &parentgeom.vertices[self.index_v[1]];
        let v2:&Point3D = &parentgeom.vertices[self.index_v[2]];
        let tt = calc_norm_v(v0,v1,v2);
        return Point3D{x:tt.0,y:tt.1,z:tt.2};
    }
    pub fn set_norm(&mut self,v:Point3D){
        self.norm = Some(v);
    }
    pub fn get_norm(&self)->&Point3D{
        if let None = self.norm{
            panic!("The snorms must be calculated at first.");
        }
        return self.norm.as_ref().unwrap();
    }
    pub fn is_dead(&self)->bool{
        return self.dead_flag;
    }
    pub fn to_string(&self,voffset:usize,vtoffset:usize,vnoffset:usize)->String{
        let mut ret:String = "".to_string();
		for ii in 0..self.index_v.len(){
			ret += format!("{}/",self.index_v[ii]+voffset).as_str();
			if self.index_vt.len() > ii{
				ret += format!("{}",self.index_vt[ii]+vtoffset).as_str();
			}
			ret += "/";
			if self.index_vn.len() > ii {
				ret += format!("{}",self.index_vn[ii]+vnoffset).as_str();
			}
			//if(ii != vNum()-1){
			if ii != self.index_v.len()-1 {
				ret += " ";
			}
		}
        return ret;
    }
}


pub struct Material{//途中
    name:String
}

pub struct Geometry{
    pub faces:Vec<Face>,
    pub vertices:Vec<Point3D>,
    norms:Vec<Point3D>,
    texture_vertices:Vec<Point2D>,
    v_to_f:Vec<Vec<usize>>,//vertices->face へのマップ
    material:Option<Material>,
    id:String
}



impl Geometry{
    pub fn new()->Geometry{
        return Geometry{
    faces:vec![],
    vertices:vec![],
    norms:vec![],
    texture_vertices:vec![],
    v_to_f:vec![],//vertices->face へのマップ
    material:None,
    id:"dummy".to_string()
        };
    }
    pub fn print(&self,f:&mut BufWriter<fs::File>,voffset:usize,vtoffset:usize,vnoffset:usize,use_f32:bool){

        f.write_all(("o ".to_string()+self.id.as_str()+"\n").as_bytes()).unwrap();
		for v in self.vertices.iter(){
			f.write_all(("v ".to_string()+v.to_string(use_f32).as_str()+"\n").as_bytes()).unwrap();
		}

		for v in self.texture_vertices.iter(){
			f.write_all(("vt ".to_string()+v.to_string(use_f32).as_str()+"\n").as_bytes()).unwrap();
		}
		for v in self.norms.iter(){
			f.write_all(("vn ".to_string()+v.to_string(use_f32).as_str()+"\n").as_bytes()).unwrap();
		}

		
		if let Some(x) = self.material.as_ref(){
			f.write_all(("usemtl ".to_string()+x.name.as_str()+"\n").as_bytes()).unwrap();
		}else{
			f.write_all("usemtl None\n".as_bytes()).unwrap();
		}
		for fa in self.faces.iter(){
			if fa.is_dead() {
			}else{
				f.write_all(("f ".to_string()+fa.to_string(voffset,vtoffset,vnoffset).as_str()+"\n").as_bytes()).unwrap();
			}
		}
    }
#[allow(non_snake_case)]
    pub fn V(&self,i:usize)->&Point3D{
        return &self.vertices[i];
    }

    pub fn add_vertex(&mut self,v:Point3D)->usize{
        self.vertices.push(v);
        return self.vertices.len()-1;
    }

    pub fn add_face(& mut self,a:usize,b:usize,c:usize){
        let mut fa = Face::new();
        fa.index_v.push(a);
        fa.index_v.push(b);
        fa.index_v.push(c);
        self.faces.push(fa);
    }
    pub fn add_reversi_face(&mut self){
        let flen:usize = self.faces.len();
        for ii in 0..flen{
            let f:&Face = &self.faces[ii];
            let mut fv:Vec<usize> = f.index_v.clone();
            assert!(fv.len() == 3);//3 しか対応してない
            fv.reverse();
            self.add_face(fv[0],fv[1],fv[2])
        }
    }

    /**
     * ある Vertex をどの Face が保有しているかを示すマップを作成する。
     */
    pub fn calc_all_v_to_f(&mut self){
        self.v_to_f = vec![vec![];self.vertices.len()];
        for (fii,ff) in self.faces.iter().enumerate(){
            for ii in ff.index_v.iter(){
                self.v_to_f[*ii].push(fii);
            }
        }
    }
    pub fn calc_all_norms(&mut self){
        if self.v_to_f.len() != self.vertices.len(){
            self.calc_all_v_to_f();
        }
        let mut vns:Vec<Point3D> = vec![];
        for ff in self.faces.iter(){
            let vn = ff.calc_norm(self);
            vns.push(vn);
        }
        for (vni,vv) in vns.into_iter().enumerate(){
            self.faces[vni].set_norm(vv);
        }

        self.norms = vec![];
        for vii in 0..self.vertices.len(){
            let fv:&Vec<usize> = &self.v_to_f[vii];
            let mut vn:Point3D = Point3D::zero();
            for fii in fv{

                let v = self.faces[*fii].get_norm();
                vn.add(&v.get_xyz());
            }
            if fv.len() > 0{
                vn.standardize();
            }
            self.norms.push(vn);
        }
        for ff in self.faces.iter_mut(){
            ff.index_vn = vec![];
            for ii in ff.index_v.iter(){
                ff.index_vn.push(*ii);
            }
        }
    }
    pub fn add_box(&mut self,center:&dyn Vector3D,size:f64){
        let bb = generate_box(center, size,size,size);
        let mut vmap:Vec<usize> = vec![];
        for v in bb.0.iter(){
            vmap.push(self.add_vertex((*v).clone()));
        }
        for f in bb.1.iter(){
            self.add_face(vmap[f.index_v[0]]
                ,vmap[f.index_v[1]]
                ,vmap[f.index_v[2]]
                );
        }
    }
	pub fn save(filename:&str,geoms:&Vec<Geometry>){
		//StringBuffer sb = new StringBuffer();
		let mut voffset:usize = 1;
		let mut vnoffset:usize = 1;
		let mut vtoffset:usize = 1;
		let materialfile:String = filename.to_string()+".mtl";


		//saveDummyMaterial(this.materialFile);
		//仮
        let mut matelsb:Vec<String> = vec![];
        matelsb.push("newmtl material1".to_string());
        matelsb.push("Ka 1.0 1.0 1.0".to_string());
        matelsb.push("Kd 0.8 0.8 0.8".to_string());
        matelsb.push("Ks 0.0 0.0 0.0".to_string());
        matelsb.push("Ke 0.2 0.2 0.2".to_string());
        matelsb.push("Ns 0.0".to_string());
        matelsb.push("Ni 1.0".to_string());
        matelsb.push("d 1.0".to_string());
        matelsb.push("illum 2".to_string());
        write_to_file(&materialfile,matelsb);


        let mut f = BufWriter::new(fs::File::create(filename).unwrap());
        
        let re = Regex::new(r".*[/\\]").unwrap();
        let filep:String = re.replace_all(materialfile.as_str(), "").to_string();
        
        f.write_all(("mtllib ".to_string()+filep.as_str()+"\n").as_bytes()).unwrap();

        for geom in geoms{
            //texture は 0.0～1.0 にスケーリングしないといけない
            geom.print(&mut f,voffset,vtoffset,vnoffset,true);
            voffset += geom.vertices.len();
            vtoffset += geom.texture_vertices.len();
            vnoffset += geom.norms.len();
        }
	}
}

fn generate_box(center:&dyn Vector3D,wid:f64,hei:f64,dep:f64)->(Vec<Point3D>,Vec<Face>){
    let mut vret:Vec<Point3D> = vec![];
    let mut fret:Vec<Face> = vec![];
    let w2 = wid/2.0;
    let h2 = hei/2.0;
    let d2 = dep/2.0;
    let xx = center.get_x();
    let yy = center.get_y();
    let zz = center.get_z();
    vret.push(Point3D::new(-w2,-h2,-d2));
    vret.push(Point3D::new( w2,-h2,-d2));
    vret.push(Point3D::new( w2,-h2, d2));
    vret.push(Point3D::new(-w2,-h2, d2));
    
    vret.push(Point3D::new(-w2, h2,-d2));
    vret.push(Point3D::new( w2, h2,-d2));
    vret.push(Point3D::new( w2, h2, d2));
    vret.push(Point3D::new(-w2, h2, d2));
    for vv in vret.iter_mut(){
        let tv = vv.get_xyz();
        vv.set_xyz(tv.0+xx,tv.1+yy,tv.2+zz);
    }

    for ii in 0..4{
        fret.push(
            Face{
                index_v:vec![(ii+1)%4+4,(ii+1)%4,ii],
                index_vt:vec![],
                index_vn:vec![],
                norm:None,
                dead_flag:false 
            }
        );
        fret.push(
            Face{
                index_v:vec![ii,ii+4,(ii+1)%4+4],
                index_vt:vec![],
                index_vn:vec![],
                norm:None,
                dead_flag:false 
            }
        );
    }
    let top_bottom:Vec<Vec<usize>> 
    = vec![vec![0,1,2],vec![2,3,0],vec![6,5,4],vec![4,7,6]];
    for tt in top_bottom.into_iter(){
        fret.push(
            Face{
                index_v:tt,
                index_vt:vec![],
                index_vn:vec![],
                norm:None,
                dead_flag:false 
            }
        );
    }
    
    return (vret,fret);
}







#[test]
fn geomtest(){
    let mut ggeo = Geometry::new();
    ggeo.add_vertex(Point3D{x:0.0,y:0.0,z:0.0});
    ggeo.add_vertex(Point3D{x:0.0,y:1.0,z:0.0});
    ggeo.add_vertex(Point3D{x:0.0,y:0.0,z:1.0});
    ggeo.add_vertex(Point3D{x:1.0,y:1.0,z:0.0});
    ggeo.add_face(0,1,2);
    ggeo.add_face(0,1,3);
    ggeo.add_face(0,2,1);
    ggeo.add_face(0,3,1);
    ggeo.calc_all_norms();
    let gv:Vec<Geometry> = vec![ggeo];

    Geometry::save("test/testgeom.obj",&gv);
}




#[test]
fn boxtest(){
    let mut ggeo = Geometry::new();
    ggeo.add_box(&Point3D::new(0.0,0.0,0.0),4.0);

    ggeo.calc_all_norms();
    let gv:Vec<Geometry> = vec![ggeo];

    Geometry::save("test/boxtestgeom.obj",&gv);
}

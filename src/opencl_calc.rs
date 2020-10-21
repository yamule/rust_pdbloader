//extern crate ocl;
/*
use ocl;
use ocl_core;



extern crate rand;
use super::pdbdata::*;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
use self::rand::prelude::*;
use super::geometry::Vector3D;

const EPSILON:f32 = 0.000000001;
const EPSILON64:f64 = EPSILON as f64;

pub fn drmsd(atoms:&Vec<PDBAtom>,reference_dist:&Vec<Vec<f64>>)->f64{
    let alen:usize = atoms.len();
    let mut ret:f64 = 0.0;
    for ii in 0..(alen-1){
        for jj in (ii+1)..alen{
            if reference_dist[ii][jj] > 0.0{
                let d:f64 = atoms[ii].distance(&atoms[jj]);
                ret += (d-reference_dist[ii][jj]).powf(2.0);
            }
        }
    }
    if ret > 0.0{
        return (ret/(alen as f64)).sqrt();
    }
    return 0.0;
}

pub fn random_movement(atoms:&mut Vec<PDBAtom>,maxval:f64,targetnum:usize,rgen:&mut StdRng){
    let alen = atoms.len();
    let mut indices:Vec<usize> = (0..alen).collect();
    indices.shuffle(rgen);
    for ii in 0..targetnum{
        let aa:&mut PDBAtom = &mut atoms[indices[ii]];
        let x = aa.get_x()+rgen.gen_range(0.0,maxval);
        let y = aa.get_y()+rgen.gen_range(0.0,maxval);
        let z = aa.get_z()+rgen.gen_range(0.0,maxval);
        aa.set_xyz(x,y,z);
    }
}
pub fn copy_xyz(src:&Vec<PDBAtom>,dest:&mut Vec<PDBAtom>){

    for (ii,aa) in src.iter().enumerate(){
        dest[ii].set_x(aa.get_x());
        dest[ii].set_y(aa.get_y());
        dest[ii].set_z(aa.get_z());
    }
}

pub struct VecAndBuffer<T:ocl_core::OclPrm>{
    pub vector:Vec<T>,
    pub buffer:ocl::Buffer<T>,
}
impl<T:ocl_core::OclPrm> VecAndBuffer<T>{
    pub fn new(vvec:Vec<T>,pro_que:&ocl::ProQue)->VecAndBuffer<T>{
        let buffer =  match 
        ocl::Buffer::<T>::builder()
        .queue(pro_que.queue().clone())
        .flags(ocl::MemFlags::new().read_write().copy_host_ptr())
        .len(vvec.len())
        .copy_host_slice(&vvec.as_slice())
        .build(){
            Ok(x) => {x},
            Err(e) =>{
                panic!("{:?}",e);
            }
        };
        
        if let Err(x) = buffer.write(&vvec).enq(){
            panic!("{:?}",x);
        }
        return VecAndBuffer{vector:vvec,buffer:buffer};
    }

    pub fn read(&mut self){
        if let Err(x) = self.buffer.read(&mut self.vector).enq(){
            panic!("error when read from buff {:?}",x);
        }
    }
    pub fn write(&mut self){
        if let Err(x) = self.buffer.write(&self.vector).enq(){
            panic!("error when read from buff {:?}",x);
        }
    }

}

pub struct MCDistGpuBuffers{
    pub pro_que:ocl::ProQue,
    pub distresult:VecAndBuffer<f32>,
    pub ref_distresult:VecAndBuffer<f32>,
    pub ref_weights:VecAndBuffer<f32>,
    pub coordinates:VecAndBuffer<f32>,
    pub coordinates_prev:VecAndBuffer<f32>,
    pub acceptflag:VecAndBuffer<i32>,
    pub random_offset:VecAndBuffer<i32>,
    pub movement_factor:VecAndBuffer<f32>,
    pub kernels:Vec<ocl::Kernel>,
    pub rgen:StdRng,
    pub prev_score:f64



}

impl MCDistGpuBuffers{


    pub fn new(atoms:&Vec<&PDBAtom>
        ,ref_dist:&Vec<Vec<f64>>
        ,ref_wei:Option<&Vec<Vec<f64>>>
        ,random_seed:Option<u64>)->MCDistGpuBuffers{
        
        let mut rgen:StdRng;//乱数生成器
        match random_seed{
            Some(x)=>{
                rgen = SeedableRng::seed_from_u64(x);
            },
            None => {
                rgen =   SeedableRng::from_rng(rand::thread_rng()).unwrap();
            }
        }
        let src = format!(r#"
        //atom-atom 間の距離を計算する
        __kernel void calc_dist(
            __global float* coordinates
            ,__global float* prev_coordinates
            ,__global float* current_dist
            ,__global float* ref_dist
            ,__global float* ref_weight
            ,__global float* prev_score
            ,__global int* accept
            ,int num_atoms
            ,int num_iter
            ,__global float* random_buffer
            ,int random_buffer_size) {{
            
            int gall = get_global_id(0)+1;

            int pp = (int)(sqrt((float)gall*2.0f)+0.25f);
            int rr = pp*(pp+1)/2;
            int g0 = pp;
            int g1 = gall-rr-1;
            if(gall == rr){{
                g0 = pp-1;
                g1 = pp-1;
            }}
            if(gall < rr){{
                pp --;
                rr =  pp*(pp+1)/2;
                g0 = pp;
                g1 = gall-rr-1;
            }}


            float EPSILON = {};
            if(g1 >= num_atoms || g0 >= num_atoms ){{
            }}else{{

                current_dist[gall-1] = sqrt(
                (coordinates[g0*3+0]-coordinates[g1*3+0])*(coordinates[g0*3+0]-coordinates[g1*3+0])
                +(coordinates[g0*3+1]-coordinates[g1*3+1])*(coordinates[g0*3+1]-coordinates[g1*3+1])
                +(coordinates[g0*3+2]-coordinates[g1*3+2])*(coordinates[g0*3+2]-coordinates[g1*3+2])
                +EPSILON) ;
            }}
        }}

        //ref dist との差異を計算して score を計算し、accept か reject かを決定する
        __kernel void calc_diff(
            __global float* coordinates
            ,__global float* prev_coordinates
            ,__global float* current_dist
            ,__global float* ref_dist
            ,__global float* ref_weight
            ,__global float* prev_score
            ,__global int* accept
            ,__global int* random_offset
            ,int num_atoms
            ,int num_iter
            ,__global float* random_buffer
            ,int random_buffer_size) {{
                
            int gall = get_global_id(0);
            float EPSILON = {};
            if(gall < num_atoms){{
                //printf(";;;%d %d\n",g1,g0);
                float diff = 0.0;
                float weight_sum = 0.0;
                
                for(int ii = 0;ii < num_atoms;ii++){{
                    int xx = max(ii,gall);
                    int yy = min(ii,gall);
                    if(xx == yy){{
                        continue;
                    }}
                    int jj = xx*(xx+1)/2 + yy;
                    if(ref_weight[jj] <= 0.0){{
                        continue;
                    }}
                    weight_sum += ref_weight[jj];
                    diff += ref_weight[jj]*(current_dist[jj]-ref_dist[jj])*(current_dist[jj]-ref_dist[jj]);
                    
                }}

                diff /= weight_sum+EPSILON;
                diff = sqrt(diff+EPSILON);
                accept[gall] = 0;
                float score = 1.0/(diff+EPSILON);
                //printf("scores: %f %f %f\n",diff,score,prev_score[gall]);
                if(score > prev_score[gall]){{
                    prev_score[gall] = score;
                    accept[gall] = 1;
                }}else if(1.0/(1.0+exp((float)(prev_score[gall]/score))*10.0)
                 > random_buffer[(random_offset[0]+gall)%random_buffer_size]){{
                    prev_score[gall] = score;
                    accept[gall] = 1;
                }}
            }}
            /*
            結果確認用
            if(gall == num_atoms){{
                double diff = 0.0;
                float weight_sum = 0.0;
                int buffsiz = (num_atoms+1)*num_atoms/2;
                for(int jj = 0;jj < buffsiz;jj++){{
                    if(ref_weight[jj] > 0.0){{
                        weight_sum += ref_weight[jj];
                        diff += ref_weight[jj]*(current_dist[jj]-ref_dist[jj])*(current_dist[jj]-ref_dist[jj]);
                    }}
                }}
                diff /= weight_sum+EPSILON;
                diff = sqrt(diff+EPSILON);
                float score = 1.0/(diff+EPSILON);
                printf("scores: %f %f %f %f\n",diff,score,prev_score[gall],weight_sum);
                prev_score[gall] = score;
            }}
            */
            
        }}
        __kernel void random_movement(
            __global float* coordinates
            ,__global float* prev_coordinates
            ,__global float* current_dist
            ,__global float* ref_dist
            ,__global float* ref_weight
            ,__global float* prev_score
            ,__global int* accept
            ,__global int* random_offset
            ,__global float* movement_factor
            ,int num_atoms
            ,int num_iter
            ,__global float* random_buffer
            ,int random_buffer_size) {{
        
            float EPSILON = {};
            
            /*
            一部だけ動かした方が良かったりもする
            //int gid = (int)(random_buffer[(get_global_id(0)+random_offset[0]*2)%random_buffer_size]*num_atoms);
            //int ggid = get_global_id(0);
            //if(ggid  < 300){{
            */
            int gid = get_global_id(0);
            if(gid  < num_atoms){{
                //printf(";;;%d %d\n",g1,g0);

                if(accept[gid] == 1){{
                    prev_coordinates[gid*3] = coordinates[gid*3];
                    prev_coordinates[gid*3+1] = coordinates[gid*3+1];
                    prev_coordinates[gid*3+2] = coordinates[gid*3+2];
                }}else{{
                    coordinates[gid*3] = prev_coordinates[gid*3];
                    coordinates[gid*3+1] = prev_coordinates[gid*3+1];
                    coordinates[gid*3+2] = prev_coordinates[gid*3+2];
                }}
                int buffpos = (gid+random_offset[0])*3;
                float mfac = movement_factor[0];
                coordinates[gid*3] += random_buffer[buffpos%random_buffer_size]*mfac-mfac/2.0;
                buffpos++;
                coordinates[gid*3+1] += random_buffer[buffpos%random_buffer_size]*mfac-mfac/2.0;
                buffpos++;
                coordinates[gid*3+2] += random_buffer[buffpos%random_buffer_size]*mfac-mfac/2.0;
            }}
        }}
        "#,EPSILON,EPSILON,EPSILON);

        //let dimss = ocl::SpatialDims::new(Some(atoms.len())
         //   ,Some(atoms.len()),Some(1)).unwrap();
         let buff_size:usize = 400;//最大原子数
         let dimss = ocl::SpatialDims::new(Some(buff_size*(buff_size/2+1)),None,None).unwrap();
         let dimss2 = ocl::SpatialDims::new(Some(buff_size),None,None).unwrap();
        let pro_que = match ocl::ProQue::builder()
        .src(src).dims(dimss)
        .build(){
            Ok(x)=>{x},
            Err(e)=>{panic!("{:?}",e);}
        };
        
        let mut ref_distresult_:Vec<f32> = vec![];
        let mut ref_weights_:Vec<f32> = vec![];
        if let Some(x) = ref_wei{
            for xx in 0..atoms.len(){
                for yy in 0..(xx+1){
                    if xx == yy {
                        ref_weights_.push(0.0 as f32);
                    }else{
                        ref_weights_.push( x[xx][yy] as f32);
                    }
                }
            }
        }else{
            for xx in 0..atoms.len(){
                for _ in 0..(xx+1){
                    ref_weights_.push(1.0 as f32);
                }
            }
        }
        for xx in 0..atoms.len(){
            for yy in 0..(xx+1){
                ref_distresult_.push(ref_dist[xx][yy] as f32);
                if ref_dist[xx][yy] <= 0.0 {
                    ref_weights_[ref_distresult_.len() -1] = 0.0;
                }
            }
        }

        let ref_distresult = VecAndBuffer::<f32>::new(ref_distresult_,&pro_que);
        let ref_weights = VecAndBuffer::<f32>::new(ref_weights_,&pro_que);
        
        let acceptflag = VecAndBuffer::<i32>::new(vec![0_i32;buff_size],&pro_que);
        let random_offset = VecAndBuffer::<i32>::new(vec![0_i32],&pro_que);
        let movement_factor = VecAndBuffer::<f32>::new(vec![1.0_f32],&pro_que);
        let prev_score = VecAndBuffer::<f32>::new(vec![std::f32::NEG_INFINITY;buff_size],&pro_que);
        


        let mut coordinates:Vec<f32> = vec![];
        for _aa in atoms.iter(){
            //coordinates.push(_aa.get_x() as f32);
            //coordinates.push(_aa.get_y() as f32);
            //coordinates.push(_aa.get_z() as f32);
            coordinates.push(rgen.gen_range(0.0_f32,1.0_f32));
            coordinates.push(rgen.gen_range(0.0_f32,1.0_f32));
            coordinates.push(rgen.gen_range(0.0_f32,1.0_f32));
        }
        let coordinates = VecAndBuffer::<f32>::new(coordinates,&pro_que);
        let coordinates_prev = VecAndBuffer::<f32>::new(coordinates.vector.clone(),&pro_que);
        let distresult = VecAndBuffer::<f32>::new(vec![0.0f32; atoms.len()*atoms.len()],&pro_que);
        
        
        let mut random_buffer_:Vec<f32> = vec![];

        for _ in 0..100000{
            random_buffer_.push(rgen.gen_range(0.0,1.0));
        }
        let random_buffer = VecAndBuffer::<f32>::new(random_buffer_,&pro_que);

        /* __global float* coordinates
            ,__global float* prev_coordinates
            ,__global float* current_dist
            ,__global float* ref_dist
            ,__global float* ref_weight
            ,__global float prev_score
            ,__global int accept
            ,int num_atoms
            ,int num_iter
            ,__global float* random_buffer
            ,int random_buffer_size) {{
         */

        let mut kernel = match
        pro_que.kernel_builder("calc_dist").global_work_size(&dimss)
        .arg(&coordinates.buffer)
        .arg(&coordinates_prev.buffer)
        .arg(&distresult.buffer)
        .arg(&ref_distresult.buffer)
        .arg(&ref_weights.buffer)
        .arg(&prev_score.buffer)
        .arg(&acceptflag.buffer)
        .arg(atoms.len() as i32)
        .arg(1 as i32)
        .arg(&random_buffer.buffer)
        .arg(random_buffer.vector.len() as i32)
        .build(){
            Ok(x)=>{x},
            Err(e)=>{
                panic!("{:?}",e);
            }
        };
        kernel.set_default_local_work_size(dimss2);

        let mut kernel2 = match
        pro_que.kernel_builder("calc_diff").global_work_size(&dimss)
        .arg(&coordinates.buffer)
        .arg(&coordinates_prev.buffer)
        .arg(&distresult.buffer)
        .arg(&ref_distresult.buffer)
        .arg(&ref_weights.buffer)
        .arg(&prev_score.buffer)
        .arg(&acceptflag.buffer)
        .arg(&random_offset.buffer)
        .arg(atoms.len() as i32)
        .arg(1 as i32)
        .arg(&random_buffer.buffer)
        .arg(random_buffer.vector.len() as i32)
        .build(){
            Ok(x)=>{x},
            Err(e)=>{
                panic!("{:?}",e);
            }
        };
        kernel2.set_default_local_work_size(dimss2);

        
        let mut kernel3 = match
        pro_que.kernel_builder("random_movement").global_work_size(&dimss)
        .arg(&coordinates.buffer)
        .arg(&coordinates_prev.buffer)
        .arg(&distresult.buffer)
        .arg(&ref_distresult.buffer)
        .arg(&ref_weights.buffer)
        .arg(&prev_score.buffer)
        .arg(&acceptflag.buffer)
        .arg(&random_offset.buffer)
        .arg(&movement_factor.buffer)
        .arg(atoms.len() as i32)
        .arg(1 as i32)
        .arg(&random_buffer.buffer)
        .arg(random_buffer.vector.len() as i32)
        .build(){
            Ok(x)=>{x},
            Err(e)=>{
                panic!("{:?}",e);
            }
        };
        kernel3.set_default_local_work_size(dimss2);

        let kernels:Vec<ocl::Kernel> = vec![kernel,kernel2,kernel3];
        let prev_score:f64 = std::f64::NEG_INFINITY;
        let ret = MCDistGpuBuffers{
            pro_que,
            distresult,
            ref_distresult,
            ref_weights,
            coordinates,
            coordinates_prev,
            acceptflag,
            random_offset,
            movement_factor,
            kernels,
            rgen,
            prev_score,
        };
        return ret;
    }
    pub fn calc(&mut self){
        unsafe {
            for ii in 0..self.kernels.len(){
                if let Err(x) = self.kernels[ii].enq(){
                    panic!("{}",x);
                }
            }
            self.random_offset.vector[0] = self.rgen.gen_range(0_i32,100000_i32);
            self.random_offset.write();
       }
    }
    pub fn check_accept_reject(&mut self){
        self.distresult.read();

        let mut rsum:f64 = 0.0;
        let mut wsum:f64 = 0.0;
        for aa in 0..self.distresult.vector.len(){
            wsum += self.ref_weights.vector[aa] as f64;
            rsum += (self.ref_weights.vector[aa]*((self.distresult.vector[aa] - self.ref_distresult.vector[aa]).powf(2.0)+EPSILON).sqrt()) as f64;
        }
        let mut dscore:f64 = 0.0;
        
        if wsum > 0.0{
            dscore = 1.0/(rsum/wsum);
        }

        let mut accept:bool = dscore > self.prev_score;
        if !accept{
            let dd = self.rgen.gen_range(0.0,1.0);
            accept =  1.0/(1.0+(self.prev_score/dscore).exp()*10.0) > dd;
        }

        println!("{} {} ",self.prev_score,dscore);
        if accept{
            self.prev_score = dscore;
            //if let Err(x) = self.coordinates_buffer.read(&mut self.coordinates).enq(){
            //    panic!("error when read from buff {:?}",x);
            //}
            self.coordinates_prev.vector.copy_from_slice(self.coordinates.vector.as_slice());
        }else{
            self.coordinates.vector.copy_from_slice(self.coordinates_prev.vector.as_slice());
            //if let Err(x) = self.coordinates_buffer.write(&self.coordinates).enq(){
            //    panic!("error when write to buff {:?}",x);
            //}
        }
    }
    //与えられた atom の x, y, z を cooridinates に入っている値にセットする
    pub fn map_to_atom(&self,atoms:&mut Vec<PDBAtom>){
        for ii in 0..atoms.len(){
            atoms[ii].set_xyz(self.coordinates.vector[ii*3] as f64
                ,self.coordinates.vector[ii*3+1] as f64
                ,self.coordinates.vector[ii*3+2] as f64);
        }
    }
    pub fn random_movement(&mut self,movmax:f32){
        for ii in 0..self.coordinates.vector.len(){
            self.coordinates.vector[ii] += self.rgen.gen_range(0.0_f32,movmax*2.0_f32)-movmax;
        }
        self.coordinates.write();
    }
    pub fn get_drmsd(&mut self) -> f64{
        self.coordinates.read();
        let num_atoms = self.coordinates.vector.len()/3;
        let mut ret:f64 = 0.0;
        let mut wsum:f64 = 0.0;

        for xx in 0..num_atoms{
            for yy in 0..(xx+1){
                let g0:usize = xx.max(yy);
                let g1:usize = xx.min(yy);
                if xx == yy {
                    continue;
                }
                let jj:i64 = g0 as i64*(g0 as i64+1)/2 + g1 as i64;
                let dist:f64 =(
                    ((self.coordinates.vector[g0*3+0]-self.coordinates.vector[g1*3+0]).powf(2.0)
                    +(self.coordinates.vector[g0*3+1]-self.coordinates.vector[g1*3+1]).powf(2.0)
                    +(self.coordinates.vector[g0*3+2]-self.coordinates.vector[g1*3+2]).powf(2.0)
                    +EPSILON) as f64).sqrt();
                if self.ref_distresult.vector[jj as usize] <= 0.0{
                    continue;
                }
                wsum += self.ref_weights.vector[jj as usize] as f64;
                ret += (self.ref_distresult.vector[jj as usize] as f64 - dist).powf(2.0)*self.ref_weights.vector[jj as usize] as f64;
               
            }
        }
        return (ret/(wsum+EPSILON64)+EPSILON64).sqrt();

    }
}


fn _trivial() -> ocl::Result<()> {
    let src = r#"
        __kernel void add(__global float* buffer,__global float* buffer2, float scalar) {
            buffer[get_global_id(0)] += scalar+buffer2[0];
        }
    "#;

    let mut pro_que = ocl::ProQue::builder()
        .src(src).dims(16)
        .build()?;

    let buffer = pro_que.create_buffer::<f32>()?;
    buffer.write(&vec![4.0;16]).enq().unwrap();
    pro_que.set_dims(2);
    let buffer2 = pro_que.create_buffer::<f32>()?;
    buffer2.write(&vec![1.0,2.0]).enq().unwrap();
    pro_que.set_dims(16);

    let kernel = pro_que.kernel_builder("add")
    .arg(&buffer)
    .arg(&buffer2)
        .arg(10.0f32)
        .build()?;

    unsafe { kernel.enq()?; }

    let mut vec = vec![0.0f32; buffer.len()];
    buffer.read(&mut vec).enq()?;
    
    println!("{}",vec.len());
    println!("{:?}",vec);
    Ok(())
}


#[test]
fn opencltest(){
    println!("{:?}",_trivial());
}



#[test]
fn calc_mc_test_gpu(){
    let pdbb:PDBEntry = load_pdb("D:/dummy/vbox_share/bioo/tbmtest/2r75.pdb");
    let mut ca_atoms:Vec<PDBAtom> = vec![];
    for (_cii,cc) in pdbb.chains.iter().enumerate(){
        for (_rii,rr) in cc.residues.iter().enumerate(){
            assert_eq!(_cii as i64,rr.parent_chain.unwrap());
            if rr.residue_name == "HOH"{
                continue;
            }
            for (_aii,aa) in rr.iter_atoms().enumerate(){
                if aa.atom_code == "CA"{
                    let mut tmpp:PDBAtom = PDBAtom::new();
                    tmpp.set_xyz(aa.get_x(),aa.get_y(),aa.get_z());
                    ca_atoms.push(tmpp);
                }
            }
        }
    }
    let alen:usize = ca_atoms.len();
    let mut ref_dist:Vec<Vec<f64>> = vec![vec![-1.0;alen];alen];
    let mut rgen:StdRng =  SeedableRng::seed_from_u64(100);
    let mut atoms_unfolded:Vec<PDBAtom> = vec![];
    for ii in 0..alen{
        let mut tmpp:PDBAtom = PDBAtom::new();
        tmpp.set_xyz(rgen.gen_range(0.0,100.0), rgen.gen_range(0.0,100.0), rgen.gen_range(0.0,100.0));
        atoms_unfolded.push(tmpp);
        for jj in 0..alen{
            if ca_atoms[ii].distance(&ca_atoms[jj]) < 12.0{
                ref_dist[ii][jj] = ca_atoms[ii].distance(&ca_atoms[jj]);
            }
        }
    }

    let mut mcdd = MCDistGpuBuffers::new(&(ca_atoms.iter().collect()),&ref_dist,None,None);
    for _ii in 0..80000{
        mcdd.calc();
        /*if _ii%10001 == 10000{
            println!("{}",mcdd.get_drmsd());
            //mcdd.movement_factor.vector[0] = 3.0/(_ii as f32/(1000 as f32));
            mcdd.movement_factor.vector[0] = 3.0/(_ii as f32/(10000 as f32));
            mcdd.movement_factor.write();
        }*/
    }
    println!("result: {}",mcdd.get_drmsd());
    mcdd.map_to_atom(&mut atoms_unfolded);
    let mut distsum:f64 = 0.0;
    let mut lcou = 0_usize;
    for ii in 0..alen{
        for jj in 0..(ii){
            if ca_atoms[ii].distance(&ca_atoms[jj]) < 12.0{
                lcou += 1;
                distsum += (atoms_unfolded[ii].distance(&atoms_unfolded[jj]) - ca_atoms[ii].distance(&ca_atoms[jj])).powf(2.0);
            }
        }
    }
    println!("result: {} {}",(distsum/lcou as f64).sqrt(),lcou);
}

#[allow(unused_parens)]
#[test]
fn arraymaptest(){
    let mut gall_o:i64 = -1;
    for xx in 0..100{
        for yy in 0..(xx+1){
            gall_o += 1;
            let gall_i = gall_o+1;
            let gall = gall_i as f64;
            let mut pp:i64 = ((gall*2.0).sqrt()+0.25) as i64;
            let rr = pp as i64*(1+pp as i64)/2;
            let mut  g0 = pp;
            let mut  g1 = gall_i-rr-1;
            
            println!("{} {} {} {} {} {}",xx,yy,pp,rr,g0,g1);
            if(gall_i == rr){
                g0 = pp-1;
                g1 = pp-1;
            }
            if(gall_i < rr){
                pp -= 1;
                let rr = pp as i64*(1+pp as i64)/2;
                g0 = pp;
                g1 = gall_i-rr-1;
            }
            assert_eq!(g0,xx);
            assert_eq!(g1,yy);
            
            let xx2 = xx.max(yy);
            let yy2 = xx.min(yy);
            if(xx == yy){{
                continue;
            }}

            let jj:i64 = xx2 as i64*(xx2 as i64+1)/2 + yy2 as i64;
            assert_eq!(jj ,gall_o);

        }

    }
}
*/



lazy_static! {
    static ref  OUTSIDE:Vec<u8> = vec![0,0,0,0];
    static ref  BLACK:Vec<u8> = vec![0,0,0,255];
    static ref  WHITE:Vec<u8> = vec![255,255,255,255];
    static ref  RED:Vec<u8> = vec![255,0,0,255];
    static ref  GREEN:Vec<u8> = vec![0,255,0,255];
    static ref  BLUE:Vec<u8> = vec![0,0,255,255];
}
pub struct Image2D{
    pub pixels:Vec<Vec<Vec<u8>>>
}

impl Image2D{
    pub fn new(w:usize,h:usize,c:usize)->Image2D{
        let r:Vec<Vec<Vec<u8>>> = vec![vec![vec![0;c];h];w];
        return Image2D{
            pixels:r
        };
    }
    pub fn get_width(&self)->usize{
        return self.pixels.len();
    }
    pub fn get_height(&self)->usize{
        if self.pixels.len() == 0{
            return 0;
        }
        return self.pixels[0].len();
    }
    pub fn is_inside(&self,x:i128,y:i128)->bool{
        if self.pixels.len() == 0{
            return false;
        }
        if x < 0{
            return false;
        }
        if y < 0{
            return false;
        }
        if x >= self.pixels.len() as i128{
            return false;
        }
        if y >= self.pixels[0].len() as i128{
            return false;
        }
        return true;
    }
    pub fn set_color(&mut self,x_:i128,y_:i128,c:&Vec<u8>){
        if self.is_inside(x_, y_){
            let x:usize = x_ as usize;
            let y:usize = y_ as usize;
            assert_eq!(c.len(),self.pixels[0][0].len());
            for ii in 0..c.len(){
                self.pixels[x][y][ii] = c[ii];
            }
        }
    }
    pub fn get_color(&self,x:i128,y:i128)->&Vec<u8>{
        if self.is_inside(x,y){
            return &OUTSIDE;
        }
        return &self.pixels[x as usize][y as usize];
    }
    pub fn fill_rect(&mut self,x:i128,y:i128,w:usize,h:usize,c:&Vec<u8>){
        for xx in x..(x+w as i128){
            for yy in y..(y+h as i128){
                self.set_color(xx,yy,c);
            }
        }
    }
}

use std::fs::File;
use std::io::prelude::*;
use flate2::write::*;
use flate2::Compression;



//https://www.setsuki.com/hsp/ext/png.htm
pub struct PngExporter{

}
impl PngExporter{
    pub fn export(filename:&str,pixels:&Vec<Vec<Vec<u8>>>){
        let mut file = File::create(filename).unwrap_or_else(|e|{panic!("Can not create {}. {:?}",filename,e);});
        file.write_all(&[0x89,0x50,0x4E,0x47,0x0D,0x0A,0x1A,0x0A]).unwrap();

        let u64_to_slice = |c:u64|->[u8;4]{
            [
                (c >> 24 ) as u8 &0xff,
                (c >> 16 ) as u8 &0xff,
                (c >> 8 ) as u8 &0xff,
                (c as u8) &0xff
                ]
        };

        let h:usize = pixels[0].len();
        let w:usize = pixels.len();
        let num_colors:usize = pixels[0][0].len();
        //IHDR
        file.write_all(&[0,0,0,13]).unwrap();
        let mut bbuf:Vec<u8> = vec![];
        bbuf.write_all(&[0x49,0x48,0x44,0x52]).unwrap();
        bbuf.write_all(&u64_to_slice(w as u64)).unwrap();
        bbuf.write_all(&u64_to_slice(h as u64)).unwrap();
        bbuf.write_all(&[8]).unwrap();
        if num_colors == 3{
            bbuf.write_all(&[2]).unwrap();
        }else{
            bbuf.write_all(&[6]).unwrap();
        }
        bbuf.write_all(&[0]).unwrap();
        bbuf.write_all(&[0]).unwrap();
        bbuf.write_all(&[0]).unwrap();
        
        let crc = calc_crc(&bbuf);
        file.write_all(bbuf.as_slice()).unwrap();
        file.write_all(&u64_to_slice(crc)).unwrap();



        //IDAT

        if num_colors != 3 && num_colors != 4{
            panic!("Num colors must be 3 or 4!");
        }
        //ラベル 1 の場合左との差、2 の場合上との差
        let mut varr_:Vec<u8> = vec![0;h*w*num_colors+h];
        let mut ppos:usize = 0;
        for yy in 0..h{
            if yy == 0{
                varr_[ppos] = 1;
                ppos+=1;
                for nn in 0..num_colors{
                    varr_[ppos] = pixels[0][yy][nn];
                    ppos+=1;
                }
                for xx in 1..w{
                    for nn in 0..num_colors{
                        varr_[ppos] = (pixels[xx][yy][nn] as usize +256 -pixels[xx-1][yy][nn] as usize & 0xff) as u8;
                        ppos += 1;
                    }
                }
            }else{
                varr_[ppos] = 2;
                ppos+=1;
                for xx in 1..w{
                    for nn in 0..num_colors{
                        varr_[ppos] = (pixels[xx][yy][nn] as usize +256 -pixels[xx][yy-1][nn] as usize & 0xff) as u8;
                        ppos += 1;
                    }
                }

            }
        }

        let mut e = ZlibEncoder::new(Vec::new(), Compression::new(3));
        e.write_all(varr_.as_slice()).unwrap();
        let mut bytes = e.finish().unwrap();


        //IDAT write
        
        file.write_all(&u64_to_slice(bytes.len() as u64)).unwrap();

        let mut bbuf:Vec<u8> = vec![];
        bbuf.write_all(&[0x49,0x44,0x41,0x54]).unwrap();
        bbuf.append(&mut bytes);

        let crc = calc_crc(&bbuf);
        file.write_all(bbuf.as_slice()).unwrap();
        file.write_all(&u64_to_slice(crc)).unwrap();
        
        //file.write_all(&[0x00,0x00,0x00,0x19,0x49,0x44,0x41,0x54,0x08,0x5B,0x63,0xFC,0xFF,0xFF,0x3F,0x03,0x0C,0xB0,0x30,0xC2,0x99,0x0C,0x0C,0x8C,0x08,0x71,0x06,0x06,0x00,0x7B,0x61,0x04,0x04,0xC1,0x5F,0x8B,0x88]).unwrap();

        //IEND
        file.write_all(&[0,0,0,0]).unwrap();
        let mut bbuf:Vec<u8> = vec![];
        bbuf.write_all(&[0x49,0x45,0x4E,0x44]).unwrap();
        let crc = calc_crc(&bbuf);
        file.write_all(bbuf.as_slice()).unwrap();
        file.write_all(&u64_to_slice(crc)).unwrap();
        
    }

}



//https://www.w3.org/TR/PNG/#D-CRCAppendix
pub fn calc_crc(vv:&Vec<u8>) -> u64{
    let mut crc_table:Vec<u64> = vec![0;256];
    for n in 0..256{
        let mut c:u64 = n as u64;
        for _ in 0..8{
            if (c & 1_u64) == 1{
                c = 0xedb88320 ^ (c >> 1);
            }else{
                c = c >> 1;
            }
        }
        crc_table[n] = c;
    }

    let mut c:u64 = 0xffffffff;
    for n in 0..vv.len() {
        c = crc_table[((c ^ (vv[n] as u64)) & 0xff) as usize] ^ (c >> 8);
    }
    return c ^ 0xffffffff;
}


#[test]
fn pngtest(){

    let mut writer = Vec::new();
    let mut deflater = ZlibDecoder::new(writer);
    deflater.write_all(&[
        0x78,0xDA,0x63,0xFC,0xFF,0xFF,0x3F,0x03,0x0C,0x30,0x32,0x20,0x73,0x10,0x4C,0x06,0x06,0x00,0xA9,0xCB,0x05,0xFE
        ]).unwrap();
    writer = deflater.finish().unwrap();
    println!("{:?}",writer);


    let pix:Vec<Vec<Vec<u8>>> = vec![
        vec![vec![255,255,255],vec![0,255,255],vec![255,0,0]],
        vec![vec![255,255,255],vec![0,255,255],vec![255,0,0]],
        vec![vec![255,255,255],vec![0,255,255],vec![255,0,0]],
        vec![vec![255,255,255],vec![0,255,255],vec![255,0,0]]
    ];
    PngExporter::export("example_files/example_output/test.png",&pix);
 
}
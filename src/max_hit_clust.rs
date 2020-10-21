
use super::sequence_alignment;

//cd-hit のフィルタリング弱い奴
//テストしてない
//kmer_filter = 0 を指定すると完全な SW するのでそれでいいかも

pub struct KmerChecker{
    knum:usize,
    num_types:usize,
    kmer_count1:Vec<Vec<usize>>,
    kmer_count2:Vec<Vec<usize>>,
}
impl KmerChecker{
    fn build(knum:usize,num_types:usize)->KmerChecker{
        let mut kmer_count1:Vec<Vec<usize>> = vec![];
        let mut kmer_count2:Vec<Vec<usize>> = vec![];
        for kk in 0..knum{
            let mut unum = num_types;
            for _ in 0..kk{
                unum *= num_types;
            }
            kmer_count1.push(vec![0;unum]);
            kmer_count2.push(vec![0;unum]);
        }
        return KmerChecker{
            knum,num_types,kmer_count1,kmer_count2
        };
    }

    //seq を使用して kmer の数を数える
    //set_2 が true の時は kmer_count2 の値を更新する
    pub fn calc_composition(&mut self,seq:&Vec<usize>,set_2:bool){
        for kk in 0..self.knum{
            if set_2{
                for bb in self.kmer_count2[kk].iter_mut(){
                    *bb = 0;
                }
            }else{
                for bb in self.kmer_count1[kk].iter_mut(){
                    *bb = 0;
                }
            }
            for ii in 0..(seq.len()-kk){
               
                let mut ppos:usize = 0;
                let step:usize = self.num_types;
                if kk == 0{
                    ppos = seq[ii];
                }else if kk == 1{
                    ppos = seq[ii]+seq[ii+1]*step;
                
                }else{
	                for kkk in 0..kk{
	                    ppos += seq[ii+kkk]*step.pow(kkk as u32);
	                }
                }
                if set_2{
                    self.kmer_count2[kk][ppos] += 1;
                }else{
                    self.kmer_count1[kk][ppos] += 1;
                }
            }
        }
    }
    
    //共通して見つかる kmer の数を数える。
    //kmer の数の少ない方を足す。
    pub fn calc_maximum_match_count_from_composition(buff1:&Vec<usize>,buff2:&Vec<usize>)->usize{
        let mut maxmatch:usize = 0;
        for zz in buff1.iter().zip(buff2.iter()){
            maxmatch += zz.0.min(zz.1);
        }
        return maxmatch; 
    }
    pub fn get_min_max_identity(&mut self,len_shorter:usize)->f64{
        let mut ret:f64 = 1.0;
        for kk in 0..self.knum{
            let kkf = kk as f64+1.0;
            let kmermatch = KmerChecker::calc_maximum_match_count_from_composition(
                &self.kmer_count1[kk]
                ,&self.kmer_count2[kk]);
            if kk != 0{
                let maxident = (kmermatch as f64
                +(len_shorter as f64 -(kkf-1.0)-kmermatch as f64) as f64
                *((kkf-1.0)/kkf)+kkf-1.0)/(len_shorter as f64);
                ret = ret.min(maxident);
            }else{
                ret = ret.min(kmermatch as f64 /(len_shorter as f64));
            }
        }
        return ret;
    }
}


pub fn count_cov_match(aligner:&sequence_alignment::SequenceAlignment
,seq1:&sequence_alignment::SeqData,seq2:&sequence_alignment::SeqData)->
AlignStatus{
    let mut startx_:i64 = -1;
    let mut starty_:i64 = -1;
    let mut maxscore:f32;
    let mut max_place:usize;
    let mut ret_cov1:usize = 0;
    let mut ret_cov2:usize = 0;
    let mut ret_match:usize = 0;

    if aligner.alignment_type == sequence_alignment::ALIGN_LOCAL{
        let xlen = aligner.seqlen_a +1;
        let ylen = aligner.seqlen_b +1;
        maxscore = 0.0;
        for ii in 0..xlen{
            for jj in 0..ylen{
                if aligner.cells[ii][jj].scores[sequence_alignment::CELL_MATCH] > maxscore{
                    maxscore = aligner.cells[ii][jj].scores[sequence_alignment::CELL_MATCH];
                    startx_ = ii as i64;
                    starty_ = jj as i64;
                }
            }
        }
        max_place = sequence_alignment::CELL_MATCH;
    }else if aligner.alignment_type == sequence_alignment::ALIGN_GLOBAL{
        let xlen = aligner.seqlen_a +1;
        let ylen = aligner.seqlen_b +1;
        startx_ = xlen as i64-1;
        starty_ = ylen as i64-1;
        let slen = aligner.cells[startx_ as usize][starty_ as usize].scores.len();
        maxscore = aligner.cells[startx_ as usize][starty_ as usize].scores[0];
        max_place = 0;
        for ii in 1..slen{
            if maxscore < aligner.cells[startx_ as usize][starty_ as usize].scores[ii]{
                maxscore = aligner.cells[startx_ as usize][starty_ as usize].scores[ii];
                max_place = ii;
            }
        }
        //eprintln!("max:{} {}++++++++++++++++++++",max_room,maxscore);
            
                    
    }else{
        let xlen = aligner.seqlen_a +1;
        let ylen = aligner.seqlen_b +1;
        maxscore = 0.0;
        for jj in 0..ylen{
            if aligner.cells[xlen-1][jj].scores[sequence_alignment::CELL_MATCH] > maxscore{
                maxscore = aligner.cells[xlen-1][jj].scores[sequence_alignment::CELL_MATCH];
                startx_ = xlen as i64-1;
                starty_ = jj as i64;
            }
        }
        for ii in 0..xlen{
            if aligner.cells[ii][ylen-1].scores[sequence_alignment::CELL_MATCH] > maxscore{
                maxscore = aligner.cells[ii][ylen-1].scores[sequence_alignment::CELL_MATCH];
                startx_ = ii as i64;
                starty_ = ylen as i64-1;
            }
        }
        max_place = sequence_alignment::CELL_MATCH;
    }
    if startx_ < 0{
        return AlignStatus{
            num_aligned_elements_1:0
            ,num_aligned_elements_2:0
            ,length_aligned_segment:0
            ,num_elements_match:0

        };
    }
    let mut currentx:usize = startx_ as usize;
    let mut currenty:usize = starty_ as usize;
    let mut current_direc = max_place;
    let mut aligned_len:usize = 0;
    loop{
        let prev_direc = aligner.cells[currentx][currenty].prev[current_direc];
        if aligner.alignment_type == sequence_alignment::ALIGN_LOCAL 
        && aligner.cells[currentx][currenty].scores[current_direc] == 0.0_f32{
            break;
        }
        //println!("direc:{} x:{} y:{}",current_direc,currentx,currenty);
        if current_direc == sequence_alignment::CELL_MATCH{
            
            if seq1.seq[currentx-1] == seq2.seq[currenty-1] {
                ret_match += 1;
            }
            aligned_len += 1;
            ret_cov1 += 1;
            ret_cov2 += 1;

            currentx -= 1;
            currenty -= 1;
        }else if  current_direc == sequence_alignment::CELL_GAPINX{
            currenty -= 1;
            ret_cov2 += 1;
            aligned_len += 1;
        }else if  current_direc == sequence_alignment::CELL_GAPINY{
            if currentx == 0{
                panic!("{}",currenty);
            }
            ret_cov1 += 1;
            aligned_len += 1;
            currentx -= 1;
        }else{
            panic!("???");
        }
        if currentx == 0 && currenty == 0{
            break;
        }
        current_direc = prev_direc;
        
    }
    return AlignStatus{
        num_aligned_elements_1:ret_cov1
        ,num_aligned_elements_2:ret_cov2
        ,length_aligned_segment:aligned_len
        ,num_elements_match:ret_match
    };
}

pub struct AlignStatus{
    num_aligned_elements_1:usize,
    num_aligned_elements_2:usize,
    length_aligned_segment:usize,
    num_elements_match:usize
}



pub fn cluster(seqs_:Vec<sequence_alignment::SeqData>
    ,identity:f64
    ,coverage_long:f64
    ,coverage_short:f64
    ,kmer_filt:usize)->Vec<SequenceCluster>{
    
        
    let sm = sequence_alignment::SubstitutionMatrix::get_blosum62_matrix();
    let mut kmer_checker:KmerChecker = KmerChecker::build(
        kmer_filt,sm.index_to_string.len());
    let mut len_seqs:Vec<(f64,usize,sequence_alignment::SeqData,Vec<usize>)>
    =seqs_.into_iter().enumerate().map(|m| (m.1.seq.len() as f64 + 1.0/(1.0+(m.0 as f64)/10000.0),m.1.seq.len(),m.1,vec![])).collect();
    for ss in len_seqs.iter_mut(){
        let useq:Vec<usize> = ss.2.seq.iter().map(|m| sm.get_string_index(m)).collect();
        ss.3 = useq;
    }
    len_seqs.sort_by(|a,b| b.0.partial_cmp(&a.0).unwrap());//名前の登場順であって欲しいので長い方が上に来るようにする
    //len_seqs.sort_by(|a,b| a.0.cmp(&b.0));
    //len_seqs.reverse();
    let num_seqs:usize = len_seqs.len();
    let mut is_member_of:Vec<(i64,f64)> = vec![(-1,0.0);num_seqs]; 
    let mut sw = sequence_alignment::SequenceAlignment::new(Box::new(sm),10.0,0.5,sequence_alignment::ALIGN_LOCAL);
    let mut clusters:Vec<SequenceCluster> = vec![];

    //let num_types:usize = sw.scoring_matrix.index_to_string.len();
    for ii in 0..num_seqs{
        if (ii+1)%100 == 0{
            eprintln!("{}/{}",ii+1,num_seqs);
        }
        if is_member_of[ii].0 > -1{
            continue;
        }
        let mut cluster = SequenceCluster{
            name_representative:len_seqs[ii].2.name.clone(),
            members:vec![]
        };
        if kmer_filt > 0{
            kmer_checker.calc_composition(&len_seqs[ii].3,false);
        }

        for jj in (ii+1)..num_seqs{
            if is_member_of[jj].0 > -1{
                continue;
            }
            if kmer_filt > 0{
                kmer_checker.calc_composition(&len_seqs[jj].3,true);
                let jlen = len_seqs[jj].1;
                let maxident = kmer_checker.get_min_max_identity(jlen);
                if maxident < identity*coverage_short{
                    //println!("skipped");
                    //panic!("{} {} {} ",ii,jj,maxident);
                    continue;
                }
            }

            //sw.prepare(&len_seqs[ii].1,&len_seqs[jj].1);
            sw.fill_matrix(&len_seqs[ii].3,&len_seqs[jj].3,None);
            let res1 = count_cov_match(&sw,&len_seqs[ii].2,&len_seqs[jj].2);
            

            let seq1_aligned:f64 = res1.num_aligned_elements_1 as f64;
            let seq2_aligned:f64 = res1.num_aligned_elements_2 as f64;
            let aligned_len:f64 = res1.length_aligned_segment as f64;
            let num_matches:f64 = res1.num_elements_match as f64;
            
            //println!("{} {} {}",seq1_aligned,seq2_aligned,num_matches);
            //ii の方が長い
            if  seq1_aligned/(len_seqs[ii].2.seq.len() as f64) < coverage_long {
                continue;
            }
            if  seq2_aligned/(len_seqs[jj].2.seq.len() as f64) < coverage_short{
                continue;
            }
            let ident:f64 = num_matches/aligned_len;
            if  ident < identity{
                continue;
            }
            is_member_of[jj] = (ii as i64 , ident);
            //eprintln!("{} {:?}",jj,(ii as i64 , ident));
            cluster.members.push((len_seqs[jj].2.name.clone(),ident));
        }
        clusters.push(cluster);
    }

    return clusters;
}

pub struct SequenceCluster{
    pub name_representative:String,
    pub members:Vec<(String,f64)>
}



/*
#[test]
fn clustering_test(){
    let fas = sequence_alignment::SeqData::load_fasta("test/pfam02545.cdd.fas",false);
    let res = cluster(fas,0.8,-1.0,0.8,3);
    for rr in res.iter(){
        println!(">{}",rr.name_representative);
        for mm in rr.members.iter(){
            println!("-\t{}\t{}",mm.0,mm.1);
        }
    }
}
*/
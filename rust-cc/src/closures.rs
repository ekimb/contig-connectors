use std::io::{self};
use std::error::Error;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use crate::BufReadDecompressor;
use histo::Histogram;
use std::fs::{File};
use std::sync::{Arc};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::DashMap;
use thread_id;
use super::mers;
use crate::mers::{Match, Mer, seq_to_kmers_nthash};
use std::path::PathBuf;
use super::Params;
use super::cluster;
use crate::{get_reader, hash_id};
use indicatif::ProgressBar;
use crate::Mutex;
use dashmap::DashSet;
use bv::BitVec;
use crate::SeqFileType;
use std::collections::HashMap;
use crate::ThreadIdType;
use glob::glob;

fn read_lines_lz4<P>(filename: P) -> io::Result<io::Lines<BufReadDecompressor<'static, BufReader<File>>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(BufReadDecompressor::new(BufReader::new(file)).unwrap().lines())
}

pub fn index(filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, index_path: &str, query_path: &str) {
    let mut resc : Arc<Mutex<Vec<(String, String)>>> = Arc::new(Mutex::new(Vec::new()));
    let mut r_c = 4;
    let mut clusters : Arc<DashMap<Vec<(usize, u64)>, Vec<Vec<(usize, u64)>>>> = Arc::new(DashMap::new());
    let mut ids : Arc<DashMap<Vec<(usize, u64)>, Vec<String>>> = Arc::new(DashMap::new());
    let seq_write = |file: &mut SeqFileType, s| {let _res = write!(file, "{}", s);};
    let sequences_files : Arc<DashMap<ThreadIdType, SeqFileType>> = Arc::new(DashMap::new());
    let create_sequences_file = |thread_id: ThreadIdType| -> SeqFileType {
        let seq_path = PathBuf::from(format!("{}.{}.sequences", index_path, thread_id));
        let file = match File::create(&seq_path) {
            Err(why) => panic!("Couldn't create file: {}.", why.to_string()),
            Ok(file) => file,
        };
        //let mut sequences_file = BufWriter::new(file);
        let mut sequences_file = SeqFileType::new(file); // regular lz4f
        //let mut sequences_file = WriteCompressor::new(&mut file, PreferencesBuilder::new().compression_level(CLEVEL_HIGH).build()).unwrap();  // too slow
        seq_write(&mut sequences_file, "# [node name]\t[list of minimizers]\t[origin]\n".to_string());
        sequences_file
    };
    let add_sketch =|sk: &Vec<(usize, u64)>, seq: Option<&str>, origin: &str, sequences_file: &mut SeqFileType, thread_id: usize|{
        if ids.get(sk).is_some() {
            ids.get_mut(sk).unwrap().push(origin.to_string());
        }
        else {
            ids.insert(sk.to_vec(), vec![origin.to_string()]);
        }
        cluster::insert(r_c, sk, &clusters);

       // let seq_line = format!("{}\t{:?}", origin, sk);
        //println!("{}", seq_line);
        //seq_write(sequences_file, format!("{}\n", seq_line));
    };
    let index_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<usize> {
        let thread_id :usize =  thread_id::get();
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
        let sk = seq_to_kmers_nthash(seq_str, seq_id, params);
        let origin = seq_id.to_string();
        add_sketch(&sk, Some(std::str::from_utf8(seq_str).unwrap()), &seq_id, &mut sequences_file, thread_id);
       return None;
    }; 
    let index_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = index_read_aux_mer(&seq_str, &seq_id);
    
    };
    let index_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = index_read_aux_mer(&seq_str, &seq_id);;
    };
    let mut main_thread_mer = |found: &mut Option<usize>| { // runs in main thread
       return None::<usize>;
    };
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, index_read_fasta_mer, |record, found| {
            main_thread_mer(found)
        });
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, index_read_fastq_mer, |record, found| {
            main_thread_mer(found)
        });
    } 
    for e in clusters.iter() {
        println!("Center: {:?}, member count: {}", ids.get(e.key()).unwrap().value(), e.value().len());

        for v in e.value().iter() {
            println!("Member: {:?}", ids.get(v).unwrap().value());
        }
    }

    let query_sketch =|sk: &Vec<(usize, u64)>, seq: Option<&str>, origin: &str|{
        let res = cluster::query(r_c, &sk, &clusters);
        if !res.is_empty() {
            for (big, small) in res.iter() {
                ids.get_mut(big).unwrap().sort_unstable();
                ids.get_mut(big).unwrap().dedup();
                ids.get_mut(small).unwrap().sort_unstable();
                ids.get_mut(small).unwrap().dedup(); 
                for id_big in ids.get(big).unwrap().iter() {
                    for id_small in ids.get(small).unwrap().iter() {
                        if id_big == id_small {continue;}
                        println!("{}\t{}", id_big, id_small);
                    }
                }
                //println!("{}\n{:?}\nis contained in\n{}\n{:?}", i);
                //println!("{}\tis contained in\t{}", ids.get(small).unwrap().value(), ids.get(big).unwrap().value());
            }
        }

       // let seq_line = format!("{}\t{:?}", origin, sk);
        //println!("{}", seq_line);
        //seq_write(sequences_file, format!("{}\n", seq_line));
    };

    let query_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<usize> {
        let k = params.k;
        let thread_id :usize =  thread_id::get();
        let sk = seq_to_kmers_nthash(seq_str, seq_id, params);
        if ids.get(&sk).is_some() {ids.get_mut(&sk).unwrap().push(seq_id.to_string());}
        else {ids.insert(sk.to_vec(), vec![seq_id.to_string()]);}
        if sk.is_empty() {return None;}
        query_sketch(&sk, Some(std::str::from_utf8(seq_str).unwrap()), seq_id);
        let origin = seq_id.to_string();
       return None;
    }; 
    let query_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_read_aux_mer(&seq_str, &seq_id);
    
    };
    let query_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_read_aux_mer(&seq_str, &seq_id);;
    };

    let mut main_thread_mer = |found: &mut Option<usize>| { // runs in main thread
       return None::<usize>;
    };

    let buf = get_reader(&PathBuf::from(query_path));
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_read_fasta_mer, |record, found| {
            main_thread_mer(found)
        });
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_read_fastq_mer, |record, found| {
            main_thread_mer(found)
        });
    } 
    /*let res_path = format!("{}{}", index_path, ".cc");
    let mut res_file =  match File::create(&res_path) {
        Err(why) => panic!("Couldn't create {}: {}", res_path, why.description()),
        Ok(res_file) => res_file,
    };
    let mut resc_view = resc.lock().unwrap();
    resc_view.sort();
    for (big, small) in resc_view.iter() {
        write!(res_file, "{}\t{}", big, small);
    }*/

}
/*pub fn query(filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, index_path: &str, res_file: &mut File) {
    let mut sequences : Arc<DashMap<Vec<u64>, Vec<(usize, String)>>> = Arc::new(DashMap::new());
     let mut process_sequence_line = |line: &str| {
        if line.starts_with("#") {return;}
        let v : Vec<&str> = line.split('\t').collect();
        let idx = v[0].parse::<usize>().unwrap();
        let mut v1 = v[1].to_string();
        let mut mer_raw : Vec<&str>  = v1.split('[').collect();
        let mut m1 = mer_raw[1].to_string();
        let mut mer_raw2 : Vec<&str> = m1.split(']').collect();
        let mut mer : Vec<u64> = mer_raw2[0].split(", ").collect::<Vec<&str>>().iter().map(|x| x.parse::<u64>().unwrap()).collect();
        let origin = v[2].to_string();
       let e =  sequences.get_mut(&mer);
       if e.is_some() {e.unwrap().push((idx, origin));}
       else {sequences.insert(mer, vec![(idx, origin)]);}
    };
    
    for path in glob(&format!("{}.*.sequences", &index_path)).expect("Failed to read glob pattern for .sequences files.") {
        let path = path.unwrap();
        let path = path.to_str().unwrap(); // rust really requires me to split the let statement in two..
        if let Ok(lines) = read_lines_lz4(path) {
            for line in lines {
                if let Ok(line_contents) = line {process_sequence_line(&line_contents);}
            }
        }
    }
    let query_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> Option<usize> {
        let k = params.k;
        let thread_id :usize =  thread_id::get();
        let sk = seq_to_kmers_nthash(seq_str, seq_id, params);
        let origin = seq_id.to_string();
       return None;
    }; 
    let query_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_read_aux_mer(&seq_str, &seq_id);
    
    };
    let query_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<usize> | {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_read_aux_mer(&seq_str, &seq_id);;
    };
    let mut main_thread_mer = |found: &mut Option<usize>| { // runs in main thread
       return None::<usize>;
    };
    let buf = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_read_fasta_mer, |record, found| {
            main_thread_mer(found)
        });
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_read_fastq_mer, |record, found| {
            main_thread_mer(found)
        });
    } 
}*/
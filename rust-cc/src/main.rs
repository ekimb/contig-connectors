#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
use indicatif::ProgressBar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fs::{DirEntry, File};
use std::io;
use std::io::stderr;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;
extern crate array_tool;
use core::hash::{Hash, Hasher};
use dashmap::DashMap;
use flate2::read::GzDecoder;
use lzzzz::lz4f::{BufReadDecompressor, Preferences, WriteCompressor};
use seq_io::BaseRecord;
use std::cell::UnsafeCell;
use std::collections::hash_map::DefaultHasher;
use std::fs;
use std::io::Result;
use std::mem::MaybeUninit;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use structopt::StructOpt;
use xx_bloomfilter::Bloom;
mod closures;
mod cluster;
mod mers;
const revcomp_aware: bool = true; // shouldn't be set to false except for strand-directed data or for debugging
pub struct SeqFileType(WriteCompressor<File>);
unsafe impl Sync for SeqFileType {} // same trick as below. we won't share files among threads but Rust can't know that.
impl SeqFileType {
    fn new(v: File) -> SeqFileType {
        SeqFileType(WriteCompressor::new(v, Preferences::default()).unwrap())
    }
}
impl Write for SeqFileType {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> Result<()> {
        self.0.flush()
    }
}

type ThreadIdType = usize;
pub struct Params {
    k: usize,
    density: f64,
}

pub fn hash_id<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

/// Try to get memory usage (resident set size) in bytes using the `getrusage()` function from libc.
// from https://github.com/digama0/mm0/blob/bebd670c5a77a1400913ebddec2c6248e76f90fe/mm0-rs/src/util.rs
fn get_memory_rusage() -> usize {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as usize * 1024
}

// thread helpers
fn thread_update_hashmap<U, V>(
    hashmap_all: &Arc<Mutex<HashMap<usize, HashMap<U, V>>>>,
    hashmap: HashMap<U, V>,
    thread_num: usize,
) {
    let mut hashmap_all = hashmap_all.lock().unwrap();
    let entry = hashmap_all.entry(thread_num).or_insert_with(HashMap::new);
    *entry = hashmap; // I believe hashmap is moved in this function as per https://stackoverflow.com/a/29490907
}

pub fn thread_update_vec<U>(
    vec_all: &Arc<Mutex<HashMap<usize, Vec<U>>>>,
    vec: Vec<U>,
    thread_num: usize,
) {
    let mut vec_all = vec_all.lock().unwrap();
    let entry = vec_all.entry(thread_num).or_insert_with(Vec::new);
    *entry = vec;
}

fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
        Ok(file) => file,
        Err(error) => panic!("Error opening compressed file: {:?}.", error),
    };
    if filename_str.ends_with(".gz") {
        filetype = "zip";
    }
    if filename_str.ends_with(".lz4") {
        filetype = "lz4";
    }
    let reader: Box<dyn BufRead + Send> = match filetype {
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))),
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ => Box::new(BufReader::new(file)),
    };
    reader
}

#[derive(StructOpt)]
enum cc {
    #[structopt(name = "index")]
    Index {
        #[structopt(parse(from_os_str), short, long)]
        directory: Option<PathBuf>,
        #[structopt(short, long, default_value="12")]
        k: usize,
        #[structopt(short="e", long, default_value="0.01")]
        density: f64,
        #[structopt(long, default_value="8")]
        threads: usize,
        #[structopt(parse(from_os_str), short, long)]
        prefix: Option<PathBuf>,
    },
    #[structopt(name = "query")]
    Query {
        #[structopt(parse(from_os_str), short, long)]
        index: Option<PathBuf>,
        #[structopt(parse(from_os_str), short, long)]
        query: Option<PathBuf>,
        #[structopt(short, long, default_value="12")]
        k: usize,
        #[structopt(short="e", long, default_value="0.01")]
        density: f64,
        #[structopt(long, default_value="8")]
        threads: usize,
        #[structopt(parse(from_os_str), short, long)]
        prefix: Option<PathBuf>,
    },
}

fn main() {
    let start = Instant::now();

    match cc::from_args() {
        cc::Index {
            directory,
            k,
            density,
            threads,
            prefix,
        } => {
            let mut output_prefix;
            let kval = k;
            let dval = density;
            let tval = threads;

            let mut index_dir = PathBuf::from(".");
            if directory.is_some() {
                index_dir = directory.unwrap();
            }

            output_prefix = PathBuf::from(format!("cc-k{}-d{}", kval, dval));
            if prefix.is_some() {
                output_prefix = prefix.unwrap();
            } else {
                println!(
                    "Warning: Using default output prefix ({}).",
                    output_prefix.to_str().unwrap()
                );
            }
            let params = Params {
                k: kval,
                density: dval,
            };
            let queue_len = 200;
            let index_path = format!("{}{}", output_prefix.to_str().unwrap(), ".clust");
            let mut index_file = match File::create(&index_path) {
                Err(why) => panic!("Couldn't create {}: {}", index_path, why.description()),
                Ok(index_file) => index_file,
            };
            let query_path = format!("query.fa");
            let mut fasta_reads: bool = false;
            let filename_str = index_dir.to_str().unwrap();
            if filename_str.contains(".fasta.")
                || filename_str.contains(".fa.")
                || filename_str.ends_with(".fa")
                || filename_str.ends_with(".fasta")
            {
                // not so robust but will have to do for now
                fasta_reads = true;
                println!("Input file: {}", filename_str);
                println!("Format: FASTA");
            }
            closures::index(
                &index_dir,
                &params,
                tval,
                queue_len,
                fasta_reads,
                &index_path,
                &query_path,
            );
            // idx \t mer \t origin
        }
        cc::Query {
            index,
            query,
            k,
            density,
            threads,
            prefix,
        } => {
            let mut output_prefix;
            let kval = k;
            let dval = density;
            let tval = threads;
            let mut qfile = PathBuf::new();
            let mut ip = PathBuf::new();

            if query.is_some() {
                qfile = PathBuf::from(query.unwrap());
            } else {
                panic!("Please provide a query.");
            }
            if index.is_some() {
                ip = PathBuf::from(index.unwrap());
            } else {
                panic!("Please provide an index prefix.");
            }

            output_prefix = PathBuf::from(format!("cc-k{}-d{}", kval, dval));
            if prefix.is_some() {
                output_prefix = prefix.unwrap();
            } else {
                println!(
                    "Warning: Using default output prefix ({}).",
                    output_prefix.to_str().unwrap()
                );
            }
            let params = Params {
                k: kval,
                density: dval,
            };
            let res_path = format!("{}{}", output_prefix.to_str().unwrap(), ".cc");
            let mut res_file = match File::create(&res_path) {
                Err(why) => panic!("Couldn't create {}: {}", res_path, why.description()),
                Ok(res_file) => res_file,
            };
            let mut queue_len = 200;
            let mut fasta_reads: bool = false;
            let filename_str = qfile.to_str().unwrap();
            if filename_str.contains(".fasta.")
                || filename_str.contains(".fa.")
                || filename_str.ends_with(".fa")
                || filename_str.ends_with(".fasta")
            {
                // not so robust but will have to do for now
                fasta_reads = true;
                println!("Input file: {}", filename_str);
                println!("Format: FASTA");
            }

            //closures::query(&qfile, &params, tval, queue_len, fasta_reads, &ip.to_str().unwrap(), &mut res_file);
        }
        _ => (),
    }
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!(
        "Maximum RSS: {:?}GB",
        (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0 / 1024.0
    );
}

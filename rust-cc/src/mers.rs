use crate::DashMap;
use crate::File;
use std::collections::VecDeque;

use std::io::Write;

use super::Params;
use nthash::NtHashIterator;
use std::collections::HashMap;
use std::collections::HashSet;
pub type Match = (
    String,
    String,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    bool,
);
pub type Mer = (u64, usize, usize, usize, bool);
#[derive(Clone, Debug)]
pub struct Hit {
    query_id: String,
    ref_id: String,
    query_s: usize,
    query_e: usize,
    ref_s: usize,
    ref_e: usize,
    hit_count: usize,
    is_rc: bool,
    query_span: usize,
    ref_span: usize,
    query_offset: usize,
    ref_offset: usize,
}

const seq_nt4_table: [u8; 256] = [
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];
// copy from http://www.cse.yorku.ca/~oz/hash.html:

pub fn hash(mut key: u64, mask: u64) -> u64 {
    key = (!key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

pub fn update_window(
    q: &mut VecDeque<u64>,
    q_pos: &mut VecDeque<usize>,
    q_min_val: u64,
    q_min_pos: i32,
    new_strobe_hashval: u64,
    i: usize,
    new_minimizer: bool,
) -> (u64, i32, bool) {
    q.pop_front();
    let popped_index = q_pos.pop_front();
    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    let mut min_val = q_min_val;
    let mut min_pos = q_min_pos;
    let mut new_minim = new_minimizer;
    if min_pos == popped_index.unwrap() as i32 {
        min_val = u64::max_value();
        min_pos = i as i32;
        for j in (0..q.len()).rev() {
            if q[j] < min_val {
                min_val = q[j];
                min_pos = q_pos[j] as i32;
                new_minim = true;
            }
        }
    } else if new_strobe_hashval < min_val {
        // the new value added to queue is the new minimum
        min_val = new_strobe_hashval;
        min_pos = i as i32;
        new_minim = true;
    }
    (min_val, min_pos, new_minim)
}

pub fn extract_mers(seq: &[u8], params: &Params) -> Vec<Vec<u8>> {
    let l = params.k;
    let hash_bound = ((params.density as f64) * 4_usize.pow(l as u32) as f64) as u64;
    let lmask: u64 = ((1 as u64) << 2 * l) - 1;
    let t = 3;
    let mut hash_count = 0;
    let mut seq_hashes = Vec::new();
    let seq_len = seq.len();
    let mut xl: [u64; 2] = [0; 2];
    let mut lp = 0;
    let lshift: u64 = (l as u64 - 1) * 2;
    for i in 0..seq.len() - l + 1 {
        let c = seq_nt4_table[seq[i] as usize];
        if c < 4 {
            xl[0] = (xl[0] << 2 | c as u64) & lmask; // forward strand
            xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift; // reverse strand
            lp += 1;
            //kminmer
            let yl: u64 = match xl[0] < xl[1] {
                true => xl[0],
                false => xl[1],
            };
            let hash_l = hash(yl, lmask);
            if hash_l <= hash_bound {
                seq_hashes.push(seq[i..i + l].to_vec());
                hash_count += 1;
            }
        } else {
            lp = 0;
            xl = [0; 2];
        }
    }
    return seq_hashes;
}
pub fn seq_to_kmers_nthash(seq: &[u8], id: &str, params: &Params) -> Vec<(usize, u64)> {
    let l = params.k;
    let density = params.density;
    let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
    let pos_hashes: Vec<(usize, u64)> = NtHashIterator::new(seq, l)
        .unwrap()
        .enumerate()
        .filter(|(i, x)| *x <= hash_bound)
        .collect();
    return pos_hashes;
}

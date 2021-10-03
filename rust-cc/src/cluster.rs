//use bk_tree::{BKTree, metrics};
use dashmap::DashMap;
use itertools::zip;
use rand::seq::SliceRandom;
use rand::thread_rng;

pub fn insert(
    r_c: usize,
    sk: &Vec<(usize, u64)>,
    clusters: &DashMap<Vec<(usize, u64)>, Vec<Vec<(usize, u64)>>>,
) {
    if clusters.len() == 0 {
        clusters.insert(sk.to_vec(), vec![]);
        return;
    } else {
        let mut min_center = Vec::new();
        let mut min_dist = r_c;
        for c in clusters.iter() {
            let center = c.key();
            let mut dist = hamming(center, sk, r_c);
            if dist < min_dist {
                min_dist = dist;
                min_center = center.to_vec();
            }
        }
        if min_dist == r_c {
            clusters.insert(sk.to_vec(), vec![]);
        } else if clusters.get_mut(&min_center).is_some() {
            clusters.get_mut(&min_center).unwrap().push(sk.to_vec());
        }
    }
}

pub fn hamming(a: &Vec<(usize, u64)>, b: &Vec<(usize, u64)>, r_c: usize) -> usize {
    let mut min_dist = r_c;
    if a == b {
        return min_dist;
    }
    let la = a.len();
    let lb = b.len();
    let ml = match la > lb {
        true => la,
        false => lb,
    };
    if ml == la {
        let w = lb;
        for i in 0..a.len() - w + 1 {
            let mut dist = 0;
            let sa = &a[i..i + w];
            for (i, j) in zip(b, sa) {
                if i.1 == j.1 {
                    continue;
                } else {
                    dist += 1;
                    if dist >= r_c {
                        break;
                    }
                }
            }
            if dist < min_dist {
                min_dist = dist;
            }
        }
    } else if ml == lb {
        let w = la;
        for i in 0..b.len() - w + 1 {
            let mut dist = 0;
            let sb = &b[i..i + w];
            for (i, j) in zip(a, sb) {
                if i.1 == j.1 {
                    continue;
                } else {
                    dist += 1;
                    if dist >= r_c {
                        break;
                    }
                }
            }
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }
    /*println!("{:?}", a);
    println!("{:?}", b);
    println!("{}", min_dist);*/
    return min_dist;
}

// A BK-tree using the Levenshtein distance metric.
pub fn levenshtein(a: &Vec<(usize, u64)>, b: &Vec<(usize, u64)>, r_c: usize) -> usize {
    let mut dist = 0;
    if a == b {
        return dist;
    }
    let la = a.len();
    let lb = b.len();
    if la == 0 {
        return lb;
    }
    if lb == 0 {
        return la;
    }
    /* Initialize the vector.
     *
     * This is why it’s fast, normally a matrix is used,
     * here we use a single vector. */
    let ml = match la > lb {
        true => la,
        false => lb,
    };
    let mut cache: Vec<usize> = (1..).take(ml).collect();
    let mut da;
    let mut db;

    /* Loop. */
    for (ib, tupb) in b.iter().enumerate() {
        dist = ib;
        da = ib;
        for (ia, tupa) in a.iter().enumerate() {
            db = if tupa.1 == tupb.1 { da } else { da + 1 };
            da = cache[ia];
            dist = if da > dist {
                if db > dist {
                    dist + 1
                } else {
                    db
                }
            } else if db > da {
                da + 1
            } else {
                db
            };
            if dist >= r_c {
                return r_c;
            }
            cache[ia] = dist;
        }
    }
    dist
}
pub fn query(
    r_c: usize,
    sk: &Vec<(usize, u64)>,
    clusters: &DashMap<Vec<(usize, u64)>, Vec<Vec<(usize, u64)>>>,
) -> Vec<(Vec<(usize, u64)>, Vec<(usize, u64)>)> {
    let mut res = Vec::<(Vec<(usize, u64)>, Vec<(usize, u64)>)>::new();
    if clusters.len() == 0 || sk.len() == 0 {
        return res;
    }
    for e in clusters.iter() {
        let center = e.key();
        let membs = e.value();
        if membs.len() == 0 {
            continue;
        }
        if hamming(center, sk, r_c) < r_c {
            for mem in membs.iter() {
                if hamming(mem, sk, r_c) == 0 {
                    let ml = match mem.len() > sk.len() {
                        true => (mem.to_vec(), sk.to_vec()),
                        false => (sk.to_vec(), mem.to_vec()),
                    };
                    res.push(ml.clone());
                }
            }
        }
    }
    return res.to_vec();
}

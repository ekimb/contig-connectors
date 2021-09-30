#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <chrono>  // for high_resolution_clock
#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include "sketch.hpp"

static uint64_t parse(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_id &ids, std::string fn) {
    uint64_t ref_size = 0;
    std::ifstream file(fn);
    std::string line, seq;
    int ref_idx = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (seq.length() > 0) {
                seqs.push_back(seq);
                lengths.push_back(seq.length());
                ref_size += seq.length();
            }
            ids[ref_idx] = line.substr(1, line.length() -1);
            ref_idx++;
            seq = "";
        }
        else {seq += line;}
    }
    if (seq.length() > 0) {
        seqs.push_back(seq);
        lengths.push_back(seq.length());
        ref_size += seq.length();
    }
    file.close();
    return ref_size;
}

void help() {
    std::cerr << "\n";
    std::cerr << "contig-connectors v0.0.1\n";
    std::cerr << "\n";
    std::cerr << "contig-connectors [options] <ref.fa> <query.fast[a/q.gz]>\n";
    std::cerr << "Options:\n";
    std::cerr << "\t-t INT Number of threads [4]\n";
    std::cerr << "\t-k INT Minimizer length [14]\n";
    std::cerr << "\t-o Name of output file to print results to [mapped.cc]\n";
    std::cerr << "\t-d INT Density parameter for sampling k-mers [0.01]. \n";
    std::cerr << "\t-f FLOAT Top fraction of repetitive k-mers to filter out from sampling [0.0002]\n";
}

int main (int argc, char **argv) {
    if (argc < 3) {
        help();
        return 0;
    }
    // Default parameters
    int threads = 4;
    int k = 14;
    float d = 0.01;
    float f = 0.0002;
    std::string out_name = "mapped.cc";

    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
            if (argv[opn][1] == 't') {
                threads = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'k') {
                k = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'o') {
                out_name = argv[opn + 1];
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'd') {
                d = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'f') {
                f = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            }
            else {help();}
        }
        if (!flag) {break;}
    }
    omp_set_num_threads(threads);
    std::cout << "Using" << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "d: " << d << std::endl;
    std::cout << "f: " << f << std::endl;
    assert(k > 7 && "K-mer size too small.");
    assert(k <= 32 && "Only up to k = 32 is supported.");
    // File name to reference
    std::string file = argv[opn];
    opn++;
    const char *query_file = argv[opn];
    auto start_p = std::chrono::high_resolution_clock::now();
    std::vector<std::string> ref_seqs;
    std::vector<unsigned int> ref_lengths;
    uint64_t ref_size;
    idx_to_id ids;
    ref_size = parse(ref_seqs, ref_lengths, ids, file);
    auto finish_p = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_p = finish_p - start_p;
    std::cout << "Read " << ref_seqs.size() << " reference contigs in " << elapsed_p.count() << "s." <<  std::endl;
    auto start_sk = std::chrono::high_resolution_clock::now();
    std::vector<sketch> ref_sk;
    for (size_t i = 0; i < ref_seqs.size(); ++i) {
        auto ref_seq = ref_seqs[i];
        sketch sk;
        sk = get_sketch(k, ref_seq, d, i);
        ref_sk.push_back(sk);
    }    
    auto finish_sk = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sk = finish_sk - start_sk;
    std::cout << "Built reference sketches in " << elapsed_sk.count() << " s." << std::endl;
}
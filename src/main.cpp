#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <unordered_map>
#include <math.h>
#include <chrono>  // for high_resolution_clock
#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include "sketch.hpp"
#include "bin.hpp"
#include "kseq++.hpp"
#include "query.hpp"

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
    std::string out_name = "containments.cc";

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
    std::vector<Sketch> ref_sk;
    Minimizer m;
    for (size_t i = 0; i < ref_seqs.size(); ++i) {
        auto ref_seq = ref_seqs[i];
        Sketch sk;
        sk = get_sketch(k, ref_seq, d, i, false);
        ref_sk.push_back(sk);
    }  
    
    auto finish_sk = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sk = finish_sk - start_sk;
    std::cout << "Built reference sketches in " << elapsed_sk.count() << " s." << std::endl;
    auto start_b = std::chrono::high_resolution_clock::now();
    Bin b;
    create_bins(ref_sk, b);
    print(b);
    auto finish_b = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_b = finish_b - start_b;
    std::cout << "Built bins in " << elapsed_b.count() << " s." << std::endl; 
    auto start_q = std::chrono::high_resolution_clock::now(); 
    std::vector<std::string> output_streams(threads);
    int chunk_size = 5000;
    for (int i = 0; i < threads; ++i){
        output_streams[i].reserve((chunk_size / threads + 1) * 450); //approx
    }
    gzFile fp = gzopen(query_file, "r");
    auto ks = klibpp::make_ikstream(fp, gzread);
    klibpp::KSeq record;
    std::ofstream output_file;
    output_file.open(out_name);
    std::string line, seq, seq_rc, prev_acc;
    Sketch query_mers;
    unsigned int q_id = 0;
    while (ks) {
        auto records = ks.read(chunk_size);
        int n_it = records.size();
        std::cout << "Processing chunk of " << n_it << " query sequences... ";
        #pragma omp parallel for num_threads(threads) shared(b, ref_sk, output_streams, output_file, q_id) private(record, query_mers)
        for (int i = 0; i < n_it; ++i) {
            auto record = records[i];
            query_mers = get_sketch(k, record.seq, d, q_id, true);
            std::vector<unsigned int> containments;
            containments = query_containments(query_mers, b, ref_sk);
            for (int i = 0; i < containments.size(); i++) {
                output_streams[omp_get_thread_num()].append(record.name);
                output_streams[omp_get_thread_num()].append("\t");
                output_streams[omp_get_thread_num()].append(ids[containments[i]]);
                output_streams[omp_get_thread_num()].append("\n");

            }
            
            q_id++;
        }
        std::cout << "Done!" << std::endl;
        // Output results
        for (int i = 0; i < threads; ++i) {
            output_file << output_streams[i];
            output_streams[i].clear();
        }
    }
    auto finish_q = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_q = finish_q - start_q;
    std::cout << "Finished queries in " << elapsed_q.count() << " s." << std::endl; 

}
#include "bin.hpp"

int32_t hamming(std::vector<Minimizer>::iterator v1_start, 
				std::vector<Minimizer>::iterator v1_end, 
				std::vector<Minimizer>::iterator v2_start, 
				std::vector<Minimizer>::iterator v2_end);

std::pair<int32_t, unsigned int> rolling_hamming_dist(std::vector<Minimizer> hashes_query, std::vector<Minimizer> hashes_subject);
std::vector<unsigned int> query_containments(std::string &seq, Sketch &query, Bin &bin, std::vector<Sketch> &ref_sk, std::vector<std::string> &seqs, unsigned int k);
float get_ani(unsigned int i, std::string &seq, std::vector<Minimizer> &query, std::vector<Minimizer> &subject, unsigned int subj_id, std::vector<std::string> &seqs, unsigned int k);
std::string revcomp(std::string seq);
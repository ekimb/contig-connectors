#include "bin.hpp"

int32_t hamming(std::vector<Minimizer>::iterator v1_start, 
				std::vector<Minimizer>::iterator v1_end, 
				std::vector<Minimizer>::iterator v2_start, 
				std::vector<Minimizer>::iterator v2_end);

int32_t rolling_hamming_dist(std::vector<Minimizer> hashes_query, std::vector<Minimizer> hashes_subject);
std::vector<unsigned int> query_containments(Sketch &query, Bin &bin, std::vector<Sketch> &ref_sk);
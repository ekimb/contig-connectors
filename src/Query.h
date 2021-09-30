#include "BinMin.h"

int32_t Hamming(std::vector<minimizer>::iterator v1_start, 
				std::vector<minimizer>::iterator v1_end, 
				std::vector<minimizer>::iterator v2_start, 
				std::vector<minimizer>::iterator v2_end);

int32_t Rolling_Hamming_Distance(std::vector<minimizer> hashes_query, std::vector<minimizer> hashes_subject);
std::vector<unsigned int> Query_Containments(sketch &query, bin &Bin_Map, std::vector<sketch> &Reference_Sketch);
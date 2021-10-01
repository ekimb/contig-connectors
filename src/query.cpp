#include "query.hpp"
#include <set>

int32_t hamming(std::vector<Minimizer>::iterator v1_start, 
				std::vector<Minimizer>::iterator v1_end, 
				std::vector<Minimizer>::iterator v2_start, 
				std::vector<Minimizer>::iterator v2_end){
	//Function to compute the hamming distance between two hash vectors of equal size;
	int32_t dist = -1;
    for (auto it = v1_start, it2 = v2_start; it != v1_end && it2 != v2_end; it++, it2++){
    	if ((*it).hash != (*it2).hash)
    		if (dist == -1) {
                dist += 2;
            } 
            else {
                dist++;
            }
    }
    if (dist == -1) {dist = 0;}
    return dist;
}

int32_t rolling_hamming_dist(std::vector<Minimizer> hashes_query, std::vector<Minimizer> hashes_subject){
	//Function to compute the min hamming distance between two hash vectors of unequal sizes. 
	//Assumes the query is the smaller sequence and the subject is the larger sequence.
	int32_t dist = hashes_subject.size();
	for(int i = 0; i < hashes_subject.size()-hashes_query.size(); i++){
		int32_t temp_dist = hamming(hashes_subject.begin()+i, hashes_subject.begin()+hashes_query.size(), 
							hashes_query.begin(), hashes_query.end());
		if (temp_dist < dist && temp_dist != -1)
			dist = temp_dist; 
	}
	return dist;
}

std::vector<unsigned int> query_containments(Sketch &query, Bin &b, std::vector<Sketch> &ref_sk){
	std::vector<Minimizer> query_minimizers = query.mins;
   	std::vector<Minimizer> query_minimizers_rev = query.mins_rev; 
	unsigned int query_id = query.ref_id;
	std::vector<unsigned int> containments;
	for (auto min = query_minimizers.begin(); min != query_minimizers.end(); min++){
		uint64_t min_value = min->hash; 
		if (b.find(min_value) != b.end()) {
			//checks if a hash exists in the bins. If not moves on to the next hash. 
			std::vector<std::tuple<Minimizer, unsigned int>> seqs_in_bin = b[min_value]; 
			for (auto s = seqs_in_bin.begin(); s != seqs_in_bin.end(); s++) {
                unsigned int ref_id = std::get<1>(*s);
				std::vector<unsigned int>::iterator it;
				std::vector<Minimizer> subject_minimizers = (ref_sk[ref_id].mins);
				it = std::find(query.cont_ids.begin(), query.cont_ids.end(), ref_id);
				if (it == query.cont_ids.end()){
					if (query_minimizers.size() < subject_minimizers.size()) {
						int32_t dist = rolling_hamming_dist(query_minimizers, subject_minimizers);
						if (dist == 0){
							containments.push_back(ref_id);
							//std::cout<<"Query Contained in "<<ref_id<<". Hamming distance is "<<dist<<std::endl;
						}
                        /*int32_t dist_rev = rolling_hamming_dist(query_minimizers_rev, subject_minimizers);
						if (dist_rev == 0){
							containments.push_back(ref_id);
							//std::cout<<"Query contained in "<<ref_id<<" in reverse. Hamming distance is "<<dist_rev<<std::endl;
						}	*/
					}
					else {
						int32_t dist = rolling_hamming_dist(subject_minimizers, query_minimizers);
						if (dist == 0) {//dist_thresh - parameterize
							containments.push_back(ref_id);
							//std::cout<<ref_id<<" Contained in query. Hamming distance is "<<dist<<std::endl;
						}
                        /*int32_t dist_rev = rolling_hamming_dist(subject_minimizers, query_minimizers_rev);
						if (dist_rev == 0) {//dist_thresh - parameterize
							containments.push_back(ref_id);
							//std::cout<<ref_id<<" Contained in query in reverse. Hamming distance is "<<dist_rev<<std::endl;
						}*/
					}
					query.cont_ids.push_back(ref_id);
				}
			}
		}
	}
	return containments;
}

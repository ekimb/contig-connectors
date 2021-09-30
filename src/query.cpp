#include "Query.h"
#include <set>

int32_t Hamming(std::vector<minimizer>::iterator v1_start, 
				std::vector<minimizer>::iterator v1_end, 
				std::vector<minimizer>::iterator v2_start, 
				std::vector<minimizer>::iterator v2_end){
	//Function to compute the hamming distance between two hash vectors of equal size;
	int32_t dist = 0;
    for (auto  it = v1_start, it2 = v2_start; it != v1_end && it2 != v2_end; it++, it2++){
    	if (std::get<0>(*it) != std::get<0>(*it2))
    		dist += 1;
    }
    return dist;
}

int32_t Rolling_Hamming_Distance(std::vector<minimizer> hashes_query, std::vector<minimizer> hashes_subject){
	//Function to compute the min hamming distance between two hash vectors of unequal sizes. 
	//Assumes the query is the smaller sequence and the subject is the larger sequence.
	int32_t dist = hashes_subject.size();
	for(int i = 0; i < hashes_subject.size()-hashes_query.size(); i++){
		int32_t temp_dist = Hamming(hashes_subject.begin()+i, hashes_subject.begin()+hashes_query.size(), 
							hashes_query.begin(), hashes_query.end());
		if (temp_dist < dist)
			dist = temp_dist; 
	}
	return dist;
}

std::vector<unsigned int> Query_Containments(sketch &query, bin &Bin_Map, std::vector<sketch> &Reference_Sketch){
	std::set<unsigned int> checked;
	std::vector<minimizer> query_minimizers = std::get<0>(query);
	unsigned int query_id = std::get<1>(query);

	std::vector<unsigned int> containments;
	for (auto hash = query_minimizers.begin(); hash != query_minimizers.end(); hash++){
		uint64_t min_value = std::get<0>(*hash); 
		if (Bin_Map.find(min_value) != Bin_Map.end())
		{
			//checks if a hash exists in the bins. If not moves on to the next hash. 
			std::vector<unsigned int> seqs_in_bin = Bin_Map[min_value]; 
			for(auto s = seqs_in_bin.begin(); s!=seqs_in_bin.end(); s++){

				//Checks if the subject is already checked for containments.
				if (checked.find(*s) == checked.end()){
					std::vector<minimizer> subject_minimizers = std::get<0>(Reference_Sketch[*s]);
					checked.insert(*s);
					
					if(query_minimizers.size()<subject_minimizers.size()){
						int32_t dist = Rolling_Hamming_Distance(query_minimizers, subject_minimizers);
						if (dist == 0){//dist_thresh - parameterize
							containments.push_back(*s);
							std::cout<<"Query Contained in "<<*s<<". Hamming distance is "<<dist<<std::endl;
						}	
					}
					else{
						int32_t dist = Rolling_Hamming_Distance(subject_minimizers, query_minimizers);
						if (dist == 0){//dist_thresh - parameterize
							containments.push_back(*s);
							std::cout<<*s<<" Contained in query. Hamming distance is "<<dist<<std::endl;
						}
					}
				}
			}
		}
	}
	return containments;
}

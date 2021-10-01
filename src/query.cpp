#include "query.hpp"
#include <set>

std::string revcomp(std::string seq) {   
	reverse(seq.begin(), seq.end());
	for (std::size_t i = 0; i < seq.length(); ++i) {
		switch (seq[i]){
		case 'A':
			seq[i] = 'T';
			break;    
		case 'C':
			seq[i] = 'G';
			break;
		case 'G':
			seq[i] = 'C';
			break;
		case 'T':
			seq[i] = 'A';
			break;
		}
	}
	return seq;
}  

std::string norm(std::string kmer) {
	std::string rev = revcomp(kmer);
	if (rev < kmer) {
		return rev;
	}
	else {
		return kmer;
	}
}



float get_ani(unsigned int i, std::string &seq, std::vector<Minimizer> &query, std::vector<Minimizer> &subject, unsigned int subj_id, std::vector<std::string> &seqs, unsigned int k) {
	std::string subj_seq = seqs[subj_id];
	unsigned int equals = 0;
	float percent = 0;
	std::string sub_seq;
	std::string sub_q;
	size_t pos;
	unsigned int r_offset;
	unsigned int q_offset;
	if (query.size() < subject.size()) {
		sub_q = seq.substr(query[0].pos); 
		std::string seed = seq.substr(query[0].pos, k); 
		pos = subj_seq.find(seed);
		if (pos==std::string::npos){
			subj_seq = revcomp(subj_seq);
			pos = subj_seq.find(seed);
		}
		if (pos==std::string::npos) {
			return 0.0;
		};
		sub_seq = subj_seq.substr(pos);
		r_offset = pos;
		q_offset = query[0].pos;
	}
	else {
		sub_seq = subj_seq.substr(subject[0].pos); 
		std::string seed = subj_seq.substr(subject[0].pos, k); 
		pos = seq.find(seed);
		if (pos==std::string::npos) {
			seq = revcomp(seq);
			pos = seq.find(seed, 0);
		}
		if (pos==std::string::npos){
			return 0.0;
		}
		sub_q = seq.substr(pos);	
		q_offset = pos;
		r_offset = subject[0].pos;
		
	}
	unsigned int max_len = subj_seq.length() - r_offset  < seq.length() - q_offset ? subj_seq.length() - r_offset : seq.length() - q_offset;
	for (unsigned int j = 0; j < max_len; j++){
		if (subj_seq[j+r_offset] == seq[j+q_offset]) {equals++;}
		//std::cout << j << "\t" << r_offset << "\t" << q_offset << std::endl;
	}
	percent =  (float)equals / (float)max_len;
	//std::cout << percent << std::endl;
	return percent;
}

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

std::pair<int32_t, unsigned int> rolling_hamming_dist(std::vector<Minimizer> hashes_query, std::vector<Minimizer> hashes_subject){
	//Function to compute the min hamming distance between two hash vectors of unequal sizes. 
	//Assumes the query is the smaller sequence and the subject is the larger sequence.
	int32_t dist = hashes_subject.size();
	unsigned int fin_i;
	for(int i = 0; i < hashes_subject.size()-hashes_query.size(); i++){
		int32_t temp_dist = hamming(hashes_subject.begin()+i, hashes_subject.begin()+hashes_query.size(), 
							hashes_query.begin(), hashes_query.end());
		if (temp_dist < dist && temp_dist != -1)
			dist = temp_dist; 
			fin_i = i;
	}
	std::pair<int32_t, unsigned int> p = std::make_pair(dist, fin_i);
	return p;
}

std::vector<unsigned int> query_containments(std::string &seq, Sketch &query, Bin &b, std::vector<Sketch> &ref_sk, std::vector<std::string> &seqs, unsigned int k){
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
						std::pair<int32_t, unsigned int> p = rolling_hamming_dist(query_minimizers, subject_minimizers);
						if (std::get<0>(p) == 0){
							float ani = get_ani(std::get<1>(p), seq, query_minimizers, subject_minimizers, ref_id, seqs, k);
							if (ani >= 0.95)
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
						std::pair<int32_t, unsigned int> p = rolling_hamming_dist(subject_minimizers, query_minimizers);
						if (std::get<0>(p) == 0){
							float ani = get_ani(std::get<1>(p), seq, query_minimizers, subject_minimizers, ref_id, seqs, k);
							if (ani >= 0.95)
								containments.push_back(ref_id);
						
							//std::cout<<"Query Contained in "<<ref_id<<". Hamming distance is "<<dist<<std::endl;
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

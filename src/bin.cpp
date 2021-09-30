#include "bin.hpp"

void create_bins(std::vector<Sketch> ref_sk, Bin &b){
    for (auto it = ref_sk.begin(); it != ref_sk.end(); it++) {
        unsigned int seq_id = (*it).ref_id;
        for (auto min = (*it).mins.begin(); min != (*it).mins.end(); min++){
            uint64_t hash = min->hash;
            std::tuple<Minimizer, unsigned int> t ((*min), seq_id);
            if (b.find(hash) != b.end())
                b[hash].push_back(t);
            else {
                std::vector<std::tuple<Minimizer, unsigned int>> e;
                e.push_back(t);
                b[hash] = e;
            }
        }
    }
}

void print(Bin &bin){
	//Temporary function to get stats on density of bins as the size of the references grow...
	for (auto it = bin.begin(); it != bin.end(); it++){
		std::cout<<"Bin ID:\t "<<std::get<0>(*it)<<"\t"<<std::get<1>(*it).size()<<"\n";
	}
}
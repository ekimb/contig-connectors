#include "BinMin.h"

void Create_Bins(std::vector<sketch> Reference_Sketch, bin &b){
	for(auto it = Reference_Sketch.begin(); it != Reference_Sketch.end(); it ++){
		unsigned int seq_id = std::get<1>(*it);
		std::vector<minimizer> minimizers = std::get<0>(*it);
		for (auto min = minimizers.begin(); min != minimizers.end(); min++){
			uint64_t hash = std::get<0>(*min);
			if (b.find(hash) != b.end())
				b[hash].push_back(seq_id);
			else{
				std::vector<unsigned int> ref_ids;
				ref_ids.push_back(seq_id);
				b[hash] = ref_ids;
			}
		}
	}
}

void Print_Bin_Summary(bin &Bin_Map){
	//Temporary function to get stats on density of bins as the size of the references grow...
	for (auto it = Bin_Map.begin(); it != Bin_Map.end(); it++){
		std::cout<<"Bin id "<<std::get<0>(*it)<<"\t"<<std::get<1>(*it).size()<<"\n";
	}
}
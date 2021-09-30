#include "sketch.hpp"
#include <iostream>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<uint64_t, std::vector<unsigned int>> bin;

void Create_Bins(std::vector<sketch> Reference_Sketch, bin &b);

void Print_Summary(std::unordered_map<uint64_t, std::vector<std::string>> &Bin_Map){
	//Print Bin summary(As a function of the growth in reference database)
}
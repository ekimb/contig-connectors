#include "sketch.hpp"
#include <iostream>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<uint64_t, std::vector<std::tuple<Minimizer, unsigned int>>> Bin;

void create_bins(std::vector<Sketch> ref_sk, Bin &b);
void print(Bin &bin);
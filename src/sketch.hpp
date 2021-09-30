#ifndef SKETCH_HPP
#define SKETCH_HPP
#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>

typedef std::unordered_map< unsigned int, std::string > idx_to_id;
typedef std::tuple <uint64_t, unsigned int, bool> minimizer;
typedef std::tuple <std::vector<minimizer>, unsigned int> sketch;

inline uint64_t hash64(uint64_t key, uint64_t mask);
static inline void minimizers(int k, std::string &seq, float density, std::vector<uint64_t> &hashes, std::vector<unsigned int> &coords);
sketch get_sketch(int k, std::string &seq, float density, unsigned int ref_idx, bool is_query);

#endif
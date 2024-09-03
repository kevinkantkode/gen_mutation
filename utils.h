#ifndef UTILS_H
#define UTILS_H

#include "linkedSequence.h"

/*
    Generate indel at each base with prob avg_mut_rate
    objects are always pass by reference, this method shouldn't modify any of the vectors
*/
void gen_INDEL(std::vector<LinkedSequence*>& linkedseqs,
               std::ofstream& mut_record,
               std::vector<double>& ins_prob,
               std::vector<double>& del_prob,
               double avg_mut_rate,
               std::mt19937& gen);

/*
    directly modify the base pair data within sequence
    objects are always pass by reference, this method shouldn't modify any of the vectors
*/
void gen_SNP(std::vector<Sequence*>& sequences, 
             std::ofstream& mut_record,
             std::vector<std::vector<double>>& snp_prob,
             double avg_mut_rate,
             std::mt19937& gen);

/*
    generate length n nuleotide string and store it on the stack
    probability of ATCG are passed in via base_prob, in that order on the vector
    this method shouldn't modify the vector
    
    returns a pointer to that string
    caller responsible of the data
*/
std::string* gen_n_nucleotides(std::vector<double>& base_prob, size_t n, std::mt19937& gen);

/*
    Copy length number of characters starting from pos in ls, move to .next ls if needed
    return the resulting string
    used to capture deleted segment (could potentially do more)
*/
std::string deep_copy_string(const LinkedSequence* ls, size_t pos, size_t length);

constexpr char index_to_nucleotide(int i) {
    assert((0 <= i && i <= 3) && "Invalid index");
    switch (i) {
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'C';
        case 3: return 'G';
        default: return '\0';
    }
}
#endif // UTILS_H
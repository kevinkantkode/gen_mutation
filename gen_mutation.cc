#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <vector>

// custom header files
#include "io.h"
#include "linkedSequence.h"
#include "utils.h"

// code for running gen mutation with LinekedSequence
int gen_mutation(int argc, char* argv[]) {

    /*---------------command line parsing----------*/

    // turn this to three if you take in a mutation model file 
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fasta_file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    /*---------------input files parsing----------*/
    // parse input fasta
    std::vector<Sequence*> sequences = parse_data(argv[1]);
    std::vector<LinkedSequence*> linkedseqs = init_vector_LinkedSequence(sequences);

    /*---------------Add mutations----------*/
    // mutation setup
    std::random_device rd;
    std::mt19937 gen(rd());
    // set up mutation record
    std::ofstream mut_record = init_mutation_record(argv[1]);

    // set avg mutation rate, this determines the probability of each mutation
    double avg_mut_rate_SV = 0.1;
    double avg_mut_rate_CNV = 0.1;
    double avg_mut_rate_INDEL = 0.1;
    double avg_mut_rate_SNP = 0.1;
    
    // start from large scale mutation to smaller
    // SV -> CNV -> indel -> SNP
    /*----------SV----------*/
    // TODO: Not sure how to do it yet, but involving breaking up larger LS into pieces via split and moving the next_
    //       pointer to manipulate LS between chromozones 
    /*----------CNV----------*/
    // TODO: idea is to use LS.split to break off LS we want, and use shallow copy to make copies of this segment
    // then perform similar action to insert_seq, but instead of inserting a newly created LS from sequence,
    // insert a copy of a LS

    /*----------INDEL----------*/
    // indel probability, 0 index must be 0.0 (no point ins/del 0 base pairs)
    // does NOT need to sum to 1
    std::vector<double> ins_prob = {0.0, 0.3, 0.2, 0.2, 0.05, 0.05};
    std::vector<double> del_prob = {0.0, 0.5, 0.2, 0.2, 0.05, 0.05};
    // call indel mutation
    gen_INDEL(linkedseqs, mut_record, ins_prob, del_prob, avg_mut_rate_INDEL, gen);

    /*----------SNP----------*/
    // snp probabilities
    std::vector<std::vector<double>> snp_prob(4, std::vector<double>(4, 1.0));
    for (size_t i = 0; i < 4; i++) {
        snp_prob[i][i] = 0.0;  // can't mut to itself
    }
    // call snp mutation
    gen_SNP(sequences, mut_record, snp_prob, avg_mut_rate_SNP, gen);

    /*---------------Output mutated reference----------*/
    write_mutated_ref(argv[1], linkedseqs);

    // clean up
    mut_record.close();
    // free objects
    free_vector(linkedseqs);
    free_vector(sequences);
    // ALL DONE
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    // CALL CHOICE OF MAIN HERE
    gen_mutation(argc, argv);

    /*---------------Performance----------*/
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Program elapsed time: " << duration.count() << " ms" << std::endl;
    return EXIT_SUCCESS;
}
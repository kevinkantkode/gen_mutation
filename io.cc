#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "io.h"

/*---------------FA parsing---------------*/
std::vector<Sequence*> parse_data(const char *file_path) {
    std::ifstream fastaFile(file_path);
    std::string line;
    std::ostringstream oss;
    bool inSequence = false;
    std::vector<Sequence*> sequences;
    Sequence* currentSequence;

    if (fastaFile.is_open()) {
        while (std::getline(fastaFile, line)) {
            if (line.empty()) {
                continue;  // Skip empty lines
            }

            if (line.back() == '\r') {
                line.pop_back();  // Remove the '\r' character if exist
            }

            if (line[0] == '>') {
                // Save the previous sequence before starting a new one
                if (inSequence) {
                    currentSequence->data = new std::string(std::move(oss.str()));
                    oss.str("");  // Clear the buffer
                    sequences.push_back(currentSequence);
                }
                currentSequence = new Sequence();
                currentSequence->id = line.substr(1);  // Skip the '>' character
                inSequence = true;
            } else {
                oss << line;  // Concatenate the sequence data
            }
        }
        // Don't forget to save the last sequence
        if (inSequence) {
            currentSequence->data = new std::string(std::move(oss.str()));
            oss.str("");  // Clear the buffer, not really need on last iteration
            sequences.push_back(currentSequence);
        }
        fastaFile.close();
    } else {
        std::cerr << "Unable to open fasta file" << std::endl;
    }
    return sequences;
}

std::ofstream init_mutation_record(const char *file_path) {
    std::ofstream mut_record("mutation_record");
    if (mut_record.is_open()) {
        // set up header info for mutation record
        mut_record << "##reference=" << file_path << std::endl;
        mut_record << "#CHROM\tPOS\tREF\tALT\tINFO" << std::endl;
    } else {
        std::cerr << "Unable to open mutation record file" << std::endl;
    }
    return mut_record;
}
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <vector>

#include "linkedSequence.h"
#include "utils.h"

std::string* gen_n_nucleotides(std::vector<double>& base_prob, size_t n, std::mt19937& gen) {
    assert(base_prob.size() == 4 && 
        "Base probability must be size 4 for ATCG probability, respectively");
    std::discrete_distribution<int> gen_nucleotide(base_prob.begin(), base_prob.end());
    std::ostringstream oss;

    for (size_t i = 0; i < n; i++) {
        int index = gen_nucleotide(gen);
        oss << index_to_nucleotide(index);
    }
    return new std::string(std::move(oss.str()));
}

std::string deep_copy_string(const LinkedSequence* ls, size_t pos, size_t length) {
    assert(ls->valid_pos(pos) && "Invalid input: start <= pos <= end");

    std::ostringstream oss;

    // Copy the substring, move to next LS if needed
    size_t writen = 0;
    while (writen < length && ls != NULL) {
        if (pos <= ls->get_end()) {
            oss << ls->get_seq_at(pos);
            pos++;
            writen++;
        } else {
            ls = ls->get_next();
            if (ls == NULL) {
                break;
            }
            pos = ls->get_start();
        }
    }

    return oss.str();  
}

void gen_INDEL(std::vector<LinkedSequence*>& linkedseqs,
               std::ofstream& mut_record,
               std::vector<double>& ins_prob,
               std::vector<double>& del_prob,
               double avg_mut_rate,
               std::mt19937& gen) {
    // create a bernoulli distribution to determine mutation or not
    std::bernoulli_distribution bd(avg_mut_rate);
    // split 50-50 between insert or delete, can change or pass in as variable if desired
    std::bernoulli_distribution coinflip(0.5);
    // distribution for the indel length
    std::discrete_distribution<size_t> gen_ins_len(ins_prob.begin(), ins_prob.end());
    std::discrete_distribution<size_t> gen_del_len(del_prob.begin(), del_prob.end());

    // start simulating indel
    for (LinkedSequence* cur_ls : linkedseqs) {
        std::string cur_chrom = cur_ls->get_seq_id();
        while (cur_ls != NULL) {
            // only mutate on valid positions, notice if ls is_empty it will exit forloop by default
            for (size_t pos = cur_ls->get_start(); pos <= cur_ls->get_end(); pos++) {
                if (bd(gen)) {
                    if (coinflip(gen)) {  // 50-50 for insert of del
                        // insert
                        size_t ins_len = gen_ins_len(gen);
                        // random base insertion (FOR NOW)
                        std::vector<double> atcg_prob = {0.25, 0.25, 0.25, 0.25};
                        std::string* data = gen_n_nucleotides(atcg_prob, ins_len, gen);
                        // ensure newseq.id is empty, so it gets cleaned up by ~LS and hence ~Sequence
                        Sequence* newseq = new Sequence(std::string(), data);
                        // write mutation to record
                        // chrom pos ref alt info
                        mut_record << cur_chrom << '\t'
                                   << pos << '\t'
                                   << cur_ls->get_seq_at(pos) << '\t'
                                   << *data << '\t'
                                   << "INS" << std::endl;
                        // actual mutation
                        cur_ls->insert_seq(newseq, pos);
                        // skip over the section we just added
                        cur_ls = cur_ls->get_next();
                        break;
                    } else {
                        // delete
                        size_t del_len = gen_del_len(gen);
                        std::string copy_del_seg = deep_copy_string(cur_ls, pos, del_len);
                        
                        // chrom pos ref alt info
                        mut_record << cur_chrom << '\t'
                                   << pos << '\t'
                                   << copy_del_seg << '\t'
                                   << "." << '\t'
                                   << "DEL" << std::endl;
                        // actual mutation
                        cur_ls->delete_section(pos, del_len);
                        // over delete is ok since the for loop condition will just skip it
                        break;
                    }
                }
            }
            if (cur_ls != NULL) {
                cur_ls = cur_ls->get_next();
            }
        }
    }
}

void gen_SNP(std::vector<Sequence*>& sequences, 
             std::ofstream& mut_record,
             std::vector<std::vector<double>>& snp_prob,
             double avg_mut_rate,
             std::mt19937& gen) {
    assert(snp_prob.size() == 4 && "SNP prob needs to be 4");
    // create a bernoulli distribution to determine mutation or not
    std::bernoulli_distribution bd(avg_mut_rate);
    // create individual mutation distribution for each ATCG
    std::discrete_distribution<size_t> mut_A(snp_prob[0].begin(), snp_prob[0].end());
    std::discrete_distribution<size_t> mut_T(snp_prob[1].begin(), snp_prob[1].end());
    std::discrete_distribution<size_t> mut_C(snp_prob[2].begin(), snp_prob[2].end());
    std::discrete_distribution<size_t> mut_G(snp_prob[3].begin(), snp_prob[3].end());


    for (Sequence* cur_seq : sequences) {
        std::string cur_chrom = cur_seq->id;
        for (size_t pos = 0; pos < cur_seq->data->size(); pos++) {
            if (bd(gen)) {
                // a snp mutation occur at cur_seq[pos]
                char ref_base = cur_seq->data->at(pos);
                char new_base = '\0';
                switch (ref_base) {
                    case 'A': new_base = index_to_nucleotide(mut_A(gen)); break;
                    case 'T': new_base = index_to_nucleotide(mut_T(gen)); break;
                    case 'C': new_base = index_to_nucleotide(mut_C(gen)); break;
                    case 'G': new_base = index_to_nucleotide(mut_G(gen)); break;
                }
                assert(ref_base != new_base && new_base != '\0');

                // actual mutation
                cur_seq->data->at(pos) = new_base;

                // chrom pos ref alt info
                mut_record << cur_chrom << '\t'
                           << pos << '\t'
                           << ref_base << '\t'
                           << new_base << '\t'
                           << "SNP" << std::endl;
            }
        }
    }
}
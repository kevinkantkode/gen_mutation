#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "linkedSequence.h"
/*----------Class member functions----------*/
// ctor
LinkedSequence::LinkedSequence(Sequence* seq, const size_t start, const size_t end, LinkedSequence* next)
    : seq_(seq), start_(start), end_(end), next_(next) {
    assert((start <= seq->data->size()-1 && end <= seq->data->size()-1) && "Invalid start and/or end value");
}

LinkedSequence::LinkedSequence(Sequence* seq) 
    : seq_(seq), start_(0), end_(seq->data->size()-1), next_(nullptr){}

// cctor
LinkedSequence::LinkedSequence(const LinkedSequence& ls)
    : seq_(ls.seq_), start_(ls.start_), end_(ls.end_), next_(ls.next_) {
    // Note: shallow copy 
}

// dtor
LinkedSequence::~LinkedSequence(){
    if (seq_->id.empty()) {
        delete seq_;
    }
    delete next_;
}

// assignment
LinkedSequence& LinkedSequence::operator=(const LinkedSequence& ls) {
    if (this != &ls) {
        seq_ = ls.seq_;
        start_ = ls.start_;
        end_ = ls.end_;
        next_ = ls.next_;
    }
    return *this;
}

bool LinkedSequence::contain_cycle() const {
    const LinkedSequence* slow = this;
    const LinkedSequence* fast = next_;
    while (slow != fast) {
        if (fast == nullptr || fast->next_ == nullptr) {
            return false;
        }
        slow = slow->next_;
        fast = fast->next_->next_;
    }
    return true;
}

LinkedSequence* LinkedSequence::split(size_t new_start) {
    assert(valid_pos(new_start) && "Invalid input: start <= new_start <= end");

    if (start_ == new_start) {
        // already split, nothing to do
        return this;
    }
    // call copy constructor
    LinkedSequence* newls = new LinkedSequence(*this);

    end_ = new_start-1;
    newls->start_ = new_start;

    newls->next_ = next_;
    next_ = newls;

    return newls;
}

LinkedSequence* LinkedSequence::insert_seq(Sequence* seq, size_t insert_pos) {
    assert(valid_pos(insert_pos) && "Invalid input: start <= insert_pos <= end");
    LinkedSequence* newls = new LinkedSequence(seq);
    LinkedSequence* nextls = split(insert_pos);

    if (insert_pos == this->start_) {
        // meaning split did not create a new LS
        assert(nextls == this);
        LinkedSequence* copyls = new LinkedSequence(*this);
        copyls->next_ = next_;
        // mark current ls as empty
        start_ = end_ + 1;
        next_ = newls;
        newls->next_ = copyls;
        return copyls;
    } else {
        next_ = newls;
        newls->next_ = nextls;
        return nextls;
    } 
}

LinkedSequence* LinkedSequence::delete_section(size_t delete_start, size_t size) {
    assert(valid_pos(delete_start) && "Invalid input: start <= delete_start <= end");

    LinkedSequence* nextls = split(delete_start);
    
    if (delete_start + size > nextls->end_) {
        // nextls is fully deleted, mark it as empty
        nextls->start_ = nextls->end_ + 1;
        assert(nextls->is_empty() && "Nextls isn't empty");

        size_t over = delete_start + size - nextls->end_ - 1;  // -1 because end is inclusive
        LinkedSequence* nnext = nextls->next_;
        if (nnext == nullptr) {
            return nullptr;
        }
        return nnext->delete_section(nnext->start_, over);
    }
    // effectively skip over the (delete_start,delete_start+size-1) segment
    nextls->start_ = delete_start + size;
    return nextls;
}

void LinkedSequence::reverse() {
    // reference binding, do NOT copy sequence
    Sequence* seq = seq_;
    size_t left = start_;
    size_t right = end_;
    while (left < right) {
        // Swap characters at left and right indices
        std::swap(seq->data[left], seq->data[right]);
        left++;
        right--;
    }
}

std::string LinkedSequence::to_string_all(bool debug) const {
    std::ostringstream oss;
    const LinkedSequence* runner = this;
    assert(contain_cycle() == false && "LinkedSequence contains cycle");
    while (runner != nullptr) {
        if (runner->is_empty() == false) {
            oss << runner->to_string();
            if (debug) {
                oss << "->";
            }
        }
        runner = runner->next_;
    }
    return std::move(oss.str());
}

/*----------Functions that uses the class----------*/

std::vector<LinkedSequence*> init_vector_LinkedSequence(std::vector<Sequence*>& sequences) {
    std::vector<LinkedSequence*> linkedseqs;
    for (Sequence* seq : sequences) {
        linkedseqs.push_back(new LinkedSequence(seq));
    }
    return linkedseqs;
}

void write_mutated_ref(const char *ref_path, std::vector<LinkedSequence*>& linkedseqs) {
    // create mut_ ref filename
    std::string original(ref_path);
    size_t dot_pos = original.rfind('.');
    //                           mut_ prefix    filename      everything after ".", should just be fa
    std::string mutated_filename = "mut_" + original.substr(0, dot_pos) + original.substr(dot_pos);

    std::ofstream mut_file(mutated_filename);
    if (mut_file.is_open()) {
        for (LinkedSequence* ls : linkedseqs) {
            // write id
            mut_file << ls->get_seq_id() << std::endl;
            // write data, to_string_all(true) to include "->"
            mut_file << ls->to_string_all(true) << std::endl; 
        }
        mut_file.close();
    } else {
        std::cerr << "Unable to open mutated output file" << std::endl;
    }
}
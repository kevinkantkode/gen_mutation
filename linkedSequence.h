#ifndef LINKEDSEQUENCE
#define LINKEDSEQUENCE

#include "io.h"

/*
    Class representing an LinkedSequence object, who points to a Sequence on the stack
    a start and end position of that Sequence (inclusive), and a .next pointer to another LinkedSequence

    for simplicity, the class will often be refered to as LS instead
*/
class LinkedSequence {
    public:
        /* 
            Construct LS with all object field specified
        */
        LinkedSequence(Sequence* seq, const size_t start, const size_t end, LinkedSequence* next);
        /*
            Construct LS with just Sequence
        */
        LinkedSequence(Sequence* seq);

        /*
            Copy constructor for LS
        */
        LinkedSequence(const LinkedSequence& ls);

        /*
            Destructor
            Only delete Sequence if LS owns the Sequence (ex. created via inserts)
            Created Sequence have nullptr for Sequence.id
            Recursively delete LS that follows after next_
        */
        ~LinkedSequence();

        /*
            Assignment Operator
        */
        // Note: do i even need this?
        LinkedSequence& operator=(const LinkedSequence& ls);

        /*
            return size of LinkedSequence
        */
        size_t size() const {
            return end_ - start_ + 1;
        }

        /*
            Check if the LS is empty (contains no data)
            return true if empty, false otherwise
        */
        bool is_empty() const {
            return size() <= 0;
        }

        /*
            split "this" into two linked LS with new LS starting at new_start
            before: prev LS -> this LS(start,end) -> next LS
            after:  prev LS -> this LS(start, nextStart-1) -> new LS(nextStart,end) -> next LS
            return the LS with new nextStart------------------^

            be VERY careful, as mishandling pointers can cause data to be lost and not freed at program exit
            or causing cycles in LS
        */
        LinkedSequence* split(size_t new_start);

        /*
            Insert given LS at insert_pos of "this"
            before: prev LS -> this LS (start, end) -> next LS
            after:  prev LS -> this LS (start, pos-1) -> new LS (pos, end) -> next LS 
            return LS after insertion
        */
        LinkedSequence* insert_seq(Sequence* seq, size_t insert_pos);

        /*
            "Delete" size bases starting at delete_start
                before: (start -> end)
                after: (start -> delete_start-1) <deleted> (delete_start+size -> end)
                                               |--.next_---^
                if delete_start + size > end, delete remaining from future LinkedSequence, 
                until reach end
            Note: no actual data is deleted, just won't be visible to any LinkedSequence
            return the LS after the deleted segment
        */
        LinkedSequence* delete_section(size_t delete_start, size_t size);

        /*
            reverse the data in Sequence in-place
        */
        void reverse();

        /*
            return the string representation all the way to end
            debug = true sets "->" delimiters between different LS
        */
        std::string to_string_all(bool debug = false) const;

        // return if the pos is valid within the context of this LS
        // notice it must be valid to THIS ls, not just valid to the sequence
        bool valid_pos(size_t pos) const {
            return start_ <= pos && pos <= end_;
        }
        // getters
        size_t get_start() const {
           return start_;
        }

        size_t get_end() const {
            return end_;
        }

        // be careful of calling non const method on this pointer
        // mainly should be used to advances a local LS*
        // ex. LS* ls = ls.get_next()
        LinkedSequence* get_next() const {
            return next_;
        }

        std::string get_seq_at(size_t pos) const {
            assert(valid_pos(pos));
            return seq_->data->substr(pos,1);
        }

        std::string get_seq_id() const {
            return seq_->id;
        }
        
    private:
        // object fields
        Sequence* seq_;
        size_t start_;
        size_t end_;
        LinkedSequence* next_;

        /*
            return true if LS contains a cycle, false if not
        */
        bool contain_cycle() const;

        /*
            return the string representation of this LinkedSequence
        */
        std::string to_string() const {
            return seq_->data->substr(start_, size());
        }
};


/*
    Create a vector<LinkedSequence*> stored on the stack using the read in data of vector<Sequence*>
*/
std::vector<LinkedSequence*> init_vector_LinkedSequence(std::vector<Sequence*>& sequences);

/*
    Write output FA using the heads of each LinkedSequence
*/
void write_mutated_ref(const char *ref_path, std::vector<LinkedSequence*>& linkedseqs);

#endif // LINKEDSEQUENCE
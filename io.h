#ifndef IO_H
#define IO_H

/*
    Class representing a sequence of DNA, which holds a string id (could be empty if Sequence is created from INS)
    and a string* to memory on stack
*/
struct Sequence{
    std::string id;
    std::string* data;

    Sequence() : Sequence("", nullptr) {}

    Sequence(std::string id_, std::string* data_):id(id_), data(data_) {}

    ~Sequence(){
        delete data;
    }
};

/*
        
    <const char*> filepath represent the path to FA/FQ file

    vector<Sequence*> represents an array of Sequences*
    Caller is responsible for freeing the resource
*/
std::vector<Sequence*> parse_data(const char *file_path);

/*
    Template function to free all vector stored resources
    ~LS and ~Sequence will handle freeing resources they own
*/
template <typename T>
void free_vector(std::vector<T*>& vec) {
    for (auto* e : vec) {
        delete e;
    }
}

/*
    create mutation record file and write the header line
    TODO: right now the output file name is "mutation record", can change in future to
          some variation of the input .FA filename
          take a look at first few lines in write_mutated_ref in LS.cc
*/
std::ofstream init_mutation_record(const char *file_path);

#endif // IO_H
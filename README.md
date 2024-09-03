Not much is required for this project, almost all are standard library

just make sure C++ v11 or higher should be enough

This project is still considered incomplete, while it can successfully run and output results,
it is still lacking 2 major types of mutation, CNV and SV.

to run, first install all the files onto your local linux enviorment, then enter

$ make

$ ./gen_mutation <fasta>

this version will output the mutated fast file to <mut_fasta> with '->' representing the connection between segments, this can be turned off via removing the 'true' argument for the to_string_all method of write_mutated_ref in linkedSequence.cc

all mutations will be documented in "mutation_record" file, note the pos value is 0-index based, adjust if desire 1-index based

Please direct any questions towards: kevinshi1118@gmail.com or create an issue under this repo

/*

    projectA:
    abPOA.hpp
    This file holds the definitions for the connector between projectA and abPOA.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <string>

#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"
#include "abPOA/abpoa.h"

using namespace std;

#ifndef PROJECTA_ABPOA_HPP
#define PROJECTA_ABPOA_HPP

// Struct to hold parameters for abPOA
struct projectA_abpoa_parameters_t {
    abpoa_t *ab;
    abpoa_para_t *abpt;
    int n_seqs;
    char **seq_names;
    int *seq_lens;
    uint8_t **seqs;
    int ** qual_weights;
    FILE *out_fp;

};

// Struct to hold the inputs and outputs for abPOA
struct projectA_abpoa_io_t {
    // TODO
};


// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gssw structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by abPOA including reserved space for the results.
void* projectA_abpoa_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by abPOA including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the abPOA alignment as well as the abPOA structs needed for alignment.
void* projectA_abpoa_calculate_batch(void* ptr, int32_t thread_index);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the abPOA alignment as well as the abPOA structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_abpoa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for abPOA.
projectA_algorithm_t* projectA_get_abpoa();


// PRE:     abpoa
//      abpoa:          Pointer to an existing gssw project.
// POST:    abpoa
//      abpoa:          Nullpointer and gssw has been destroyed.
void projectA_abpoa_destroy(projectA_algorithm_t* abpoa);


#endif
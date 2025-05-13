/*

    projectA:
    s_gwfa.hpp
    This file holds the definitions for the connector between projectA and s_gwfa.
    Author: Frederic zur Bonsen <fzurbonsen@ethz.ch>

*/
#include <vector>
#include <string>
#include <unordered_map>

#include "s_gwfa/s_gwfa.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"

using namespace std;

#ifndef PROJECTA_S_GWFA_HPP
#define PROJECTA_S_GWFA_HPP

// PRE:     alignment
//      alignment:      Pointer to a projectA_alignment_t.
// POST:    alignment
//      alignment:      Pointer to a projectA_alignment with the full CIGAR of the given path.
void projectA_s_gwfa(projectA_alignment_t* alignment);


// Struct to hold inputs for gwfa
struct projectA_s_gwfa_parameters_t {
    s_gwfa_graph_t* graph;
    const char* seq;
    int32_t len;

    projectA_hash_graph_t* projectA_hash_graph;

    // Constructor for projectA_s_gwfa_parameters_t
    projectA_s_gwfa_parameters_t(s_gwfa_graph_t* graph, const char* seq, int32_t len, projectA_hash_graph_t* projectA_hash_graph);
};


// Struct to hold complete path information
struct projectA_s_gwfa_path_t {
    s_gwfa_path_t* path;
    int32_t score;
};

// Struct to hold the inputs and outputs for s_gwfa
struct projectA_s_gwfa_io_t {
    vector<vector<projectA_s_gwfa_parameters_t>> parameters;
    vector<vector<projectA_s_gwfa_path_t>> paths;
    int32_t size;
};


// create a s_gwfa graph
s_gwfa_graph_t* projectA_s_gwfa_graph_build(projectA_hash_graph_t* in_graph);

// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gssw structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gwfa including reserved space for the results.
void* projectA_s_gwfa_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by s_gwfa including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the s_gwfa alignment as well as the gwfa structs needed for alignment.
void* projectA_s_gwfa_calculate_batch(void* ptr, int32_t thread_index);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the s_gwfa alignment as well as the gwfa structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_s_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for s_gwfa.
projectA_algorithm_t* projectA_get_s_gwfa();


// PRE:     s_gwfa
//      s_gwfa:       Pointer to an existing s_gwfa project.
// POST:    s_gwfa
//      s_gwfa:       Nullpointer and gssw has been destroyed.
void projectA_s_gwfa_destroy(projectA_algorithm_t* s_gwfa);

#endif
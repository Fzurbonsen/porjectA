/*

    projectA:
    gssw.hpp
    This file holds the definitions for the connector between projectA and gssw.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/
#include <vector>
#include <string>

#include "gwfa/gwfa.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"

using namespace std;


#ifndef PROJECTA_GWFA_HPP
#define PROJECTA_GWFA_HPP

// Struct to hold inputs for gwfa
struct projectA_gwfa_parameters_t {
    void* km;
    gwf_graph_t* graph;
    int32_t ql;
    const char *q;
    int32_t v0;
    int32_t v1;
    uint32_t max_lag;
    int32_t traceback;

    projectA_hash_graph_t* projectA_hash_graph;

    // Constructor for projectA_gwfa_parameters_t
    projectA_gwfa_parameters_t(void* km, gwf_graph_t* graph, int32_t ql, const char* q, int32_t v0, int32_t v1,
                                uint32_t max_lag, int32_t traceback, projectA_hash_graph_t* projectA_hash_graph);
};


// Struct to hold complete path information
struct projectA_gwfa_path_t {
    gwf_path_t path;
    int32_t score;
};

// Struct to hold the inputs and outputs for gwfa
struct projectA_gwfa_io_t {
    vector<vector<projectA_gwfa_parameters_t>> parameters;
    vector<vector<projectA_gwfa_path_t>> paths;
    int32_t size;
};


// PRE:     g
//      g:              Pointer to a gwf_graph_t.
// POST:    g
//      g:              Nullpointer and all the memory g pointed to is freed.
void gwf_free(gwf_graph_t *g);


// 
gwf_graph_t* projectA_hash_graph_to_gwf_graph(projectA_hash_graph_t* in_graph);


// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gssw structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gwfa including reserved space for the results.
void* projectA_gwfa_init(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by gwfa including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the gwfa alignment as well as the gwfa structs needed for alignment.
void* projectA_gwfa_calculate_batch(void* ptr, int32_t thread_index);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the gwfa alignment as well as the gwfa structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for gssw.
projectA_algorithm_t* projectA_get_gwfa();


// PRE:     gssw
//      gssw:       Pointer to an existing gssw project.
// POST:    gssw
//      gssw:       Nullpointer and gssw has been destroyed.
void projectA_gwfa_destroy(projectA_algorithm_t* gssw);


#endif
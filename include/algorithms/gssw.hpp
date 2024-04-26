/*

    projectA:
    gssw.hpp
    This file holds the definitions for the connector between projectA and gssw.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/
#include <vector>
#include <string>

#include "gssw/gssw.h"
#include "algorithm.hpp"
#include "graph.hpp"

using namespace std;

#ifndef PROJECTA_GSSW_HPP
#define PROJECTA_GSSW_HPP

// Sturct to hold inputs for gssw
struct projectA_gssw_parameters_t {
    gssw_graph* graph;
    const char* read;
    int8_t* nt_table;
    int8_t* mat;
    uint8_t gap_open;
    uint8_t gap_extension;

    // Constructor for projectA_gssw_parameters_t
    projectA_gssw_parameters_t(gssw_graph* graph, const char* read, int8_t* nt_table,
                                 int8_t* mat, uint8_t gap_open, uint8_t gap_extension);
};


// Struct to hold the inputs and outputs for gssw
struct projectA_gssw_io_t {
    vector<projectA_gssw_parameters_t> parameters;
    vector<gssw_graph_mapping*> gms;
};



gssw_graph* projectA_hash_graph_to_gssw_graph(projectA_hash_graph_t* in_graph, int8_t* nt_table, 
                                                int8_t* mat, uint8_t gap_open, uint8_t gap_extension);


void* projectA_gssw_init(vector<projectA_algorithm_input_t>& graphs);


void* projectA_gssw_calculate_batch(void* ptr);


void projectA_gssw_post(void* ptr);



// PRE:     
// POST:    return
//      return:     Pointer to a projectA_algorithm_t struct that holds the function pointers for gssw.
projectA_algorithm_t* projectA_get_gssw();


// PRE:     gssw
//      gssw:       Pointer to an existing gssw project.
// POST:    gssw
//      gssw:       Nullpointer and gssw has been destroyed.
void projectA_gssw_destroy(projectA_algorithm_t* gssw);


#endif // PROJECTA_GSSW_HPP
/*

    projectA:
    gt_gwfa.hpp
    This file holds the definitions for the connector between projectA and gt_gwfa.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/
#include <vector>
#include <string>
#include <utility>
#include <tuple>

#include "gt_gwfa/graphs.h"
#include "algorithm.hpp"
#include "graph.hpp"

using namespace std;

#ifndef PROJECTA_GT_GWFA_HPP
#define PROJECTA_GT_GWFA_HPP

// Struct to hold inputs for gt_gwfa
struct projectA_gt_gwfa_parameters_t {
    gssw_graph* gssw;
    gwf_graph_t* gwf;
    const char* read;

    // Constructor for projectA_gt_gwfa_t
    projectA_gt_gwfa_parameters_t(gssw_graph* gssw, gwf_graph_t* gwf, const char* r);
};


// Struct to hold the inputs and outputs for gt_gwfa
struct projectA_gt_gwfa_io_t {
    vector<projectA_gt_gwfa_parameters_t> parameters;
    vector<gssw_graph_mapping*> gms;
};


gssw_graph* projectA_hash_graph_to_gt_gssw_graph(projectA_hash_graph_t* in_graph);


void* projectA_gt_gwfa_init(vector<projectA_algorithm_input_t>& graphs);


void* projectA_gt_gwfa_calcualte_batch(void* ptr);


void projectA_gt_gwfa_post(void* ptr);



// PRE:     
// POST:    return 
//      return:     projectA_algorithm struct that holds the function pointers to use
//                  gt_gwfa init/calculate_batch/post.
projectA_algorithm_t* projectA_get_gt_gwfa();


// PRE:     gt_gwfa
//      gt_gwfa:    Pointer to an existing gssw project.
// POST:    gt_gwfa
//      gt_gwfa:    Nullpointer and gssw has been destroyed.
void projectA_gt_gwfa_destroy(projectA_algorithm_t* gt_gwfa);

#endif // PROJECTA_GT_GWFA_HPP
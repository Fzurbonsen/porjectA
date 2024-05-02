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


// PRE:     in_graph, nt_table, mat, gap_open, gap_extension
//      in_graph:       Pointer to a projectA_hash_graph_t.
// POST:    return
//      return:         Pointer to gssw graph with the same information as the in_graph and the relevant paramaters given as inputs.
gssw_graph* projectA_hash_graph_to_gt_gssw_graph(projectA_hash_graph_t* in_graph);


// PRE:     cigar
//      cigar:          Pointer to a valid gssw cigar.
// POST:    return
//      return:         ProjectA cigar struct that holds an equivalent cigar to the cigar in the input struct.
projectA_cigar_t projectA_gt_gwfa_get_cigar(gssw_cigar* cigar);


// PRE:     graph, gm
//      graph:          Pointer to a valid projectA hash graph struct.
//      gm:             Pointer to a valid gssw graph mapping on the input graph or a derivative of it graph.
// POST:    retrun:
//      return:         Pointer to a projectA alignment struct that holds the information of the gssw graph mapping.
projectA_alignment_t* projectA_gt_gwfa_graph_mapping_to_alignment(projectA_hash_graph_t* graph, gssw_graph_mapping* gm);


// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gt_gwfa structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gt_gwfa including reserved space for the results.
void* projectA_gt_gwfa_init(vector<projectA_algorithm_input_t>& graphs);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by gt_gwfa including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the gt_gwfa alignment as well as the gssw structs needed for alignment.
void* projectA_gt_gwfa_calcualte_batch(void* ptr);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the gt_gwfa alignment as well as the gssw structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.
void projectA_gt_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments);


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
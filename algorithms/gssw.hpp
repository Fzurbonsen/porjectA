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
#include "alignment.hpp"


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

    projectA_hash_graph_t* projectA_hash_graph;

    // Constructor for projectA_gssw_parameters_t
    projectA_gssw_parameters_t(gssw_graph* graph, const char* read, int8_t* nt_table,
                                 int8_t* mat, uint8_t gap_open, uint8_t gap_extension,
                                 projectA_hash_graph_t* projectA_hash_graph);
};


// Struct to hold the inputs and outputs for gssw
struct projectA_gssw_io_t {
    vector<projectA_gssw_parameters_t> parameters;
    vector<gssw_graph_mapping*> gms;
};


// PRE:     in_graph, nt_table, mat, gap_open, gap_extension
//      in_graph:       Pointer to a projectA_hash_graph_t.
//      nt_table:       Pointer to int8_t holding a valid nt table created by a gssw helper function for this graph.
//      mat:            Pointer to int8_t jolding a valid mat created by a gssw helper function for this graph.
//      gap_open:       Uint8_t that holds the gap open parameter for this graph alignment.
//      gap-extension:  Uint8_t that holds the gap extension paramter for this graph alignment.
// POST:    return
//      return:         Pointer to gssw graph with the same information as the in_graph and the relevant paramaters given as inputs.
gssw_graph* projectA_hash_graph_to_gssw_graph(projectA_hash_graph_t* in_graph, int8_t* nt_table, 
                                                int8_t* mat, uint8_t gap_open, uint8_t gap_extension);


// PRE:     cigar
//      cigar:          Pointer to a valid gssw cigar.
// POST:    return
//      return:         ProjectA cigar struct that holds an equivalent cigar to the cigar in the input struct.
projectA_cigar_t projectA_gssw_get_cigar(gssw_cigar* cigar);


// PRE:     graph, gm
//      graph:          Pointer to a valid projectA hash graph struct.
//      gm:             Pointer to a valid gssw graph mapping on the input graph or a derivative of it graph.
// POST:    retrun:
//      return:         Pointer to a projectA alignment struct that holds the information of the gssw graph mapping.
projectA_alignment_t* projectA_gssw_graph_mapping_to_alignment(projectA_hash_graph_t* graph, gssw_graph_mapping* gm);



// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gssw structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gssw including reserved space for the results.
void* projectA_gssw_init(vector<projectA_algorithm_input_t>& graphs);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by gssw including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the gssw alignment as well as the gssw structs needed for alignment.
void* projectA_gssw_calculate_batch(void* ptr);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the gssw alignment as well as the gssw structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_gssw_post(void* ptr, vector<projectA_alignment_t*>& alignments);


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for gssw.
projectA_algorithm_t* projectA_get_gssw();


// PRE:     gssw
//      gssw:       Pointer to an existing gssw project.
// POST:    gssw
//      gssw:       Nullpointer and gssw has been destroyed.
void projectA_gssw_destroy(projectA_algorithm_t* gssw);


#endif // PROJECTA_GSSW_HPP
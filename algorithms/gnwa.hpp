/*

    projectA:
    gssw.hpp
    This file holds the definitions for the connector between projectA and gssw.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/
#include <vector>
#include <string>

#include "GNWA/gnwa.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"


using namespace std;

#ifndef PROJECTA_GNWA_HPP
#define PROJECTA_GNWA_HPP

// Sturct to hold inputs for gnwa
struct projectA_gnwa_parameters_t {
    gnwa_graph_t* graph;
    const char* read;
    int8_t* nt_table;
    int8_t* mat;
    uint8_t gap_open;
    uint8_t gap_extension;

    projectA_hash_graph_t* projectA_hash_graph;

    // Constructor for projectA_gnwa_parameters_t
    projectA_gnwa_parameters_t(gnwa_graph_t* graph, const char* read, int8_t* nt_table,
                                 int8_t* mat, uint8_t gap_open, uint8_t gap_extension,
                                 projectA_hash_graph_t* projectA_hash_graph);
};


// Struct to hold the inputs and outputs for gnwa
struct projectA_gnwa_io_t {
    vector<vector<projectA_gnwa_parameters_t>> parameters;
    vector<vector<gnwa_alignment_t*>> gms;
    int32_t size;
};


// PRE:     in_graph, nt_table, mat, gap_open, gap_extension
//      in_graph:       Pointer to a projectA_hash_graph_t.
//      nt_table:       Pointer to int8_t holding a valid nt table created by a gnwa helper function for this graph.
//      mat:            Pointer to int8_t jolding a valid mat created by a gnwa helper function for this graph.
//      gap_open:       Uint8_t that holds the gap open parameter for this graph alignment.
//      gap-extension:  Uint8_t that holds the gap extension paramter for this graph alignment.
// POST:    return
//      return:         Pointer to gnwa graph with the same information as the in_graph and the relevant paramaters given as inputs.
gnwa_graph_t* projectA_hash_graph_to_gnwa_graph(projectA_hash_graph_t* in_graph, int8_t* nt_table, 
                                                int8_t* mat, uint8_t gap_open, uint8_t gap_extension);


// PRE:     cigar
//      cigar:          Pointer to a valid gnwa cigar.
// POST:    return
//      return:         ProjectA cigar struct that holds an equivalent cigar to the cigar in the input struct.
projectA_cigar_t projectA_gnwa_get_cigar(gnwa_cigar_t* cigar);



// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gnwa structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gnwa including reserved space for the results.
void* projectA_gnwa_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by gnwa including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the gnwa alignment as well as the gnwa structs needed for alignment.
void* projectA_gnwa_calculate_batch(void* ptr, int32_t thread_index);


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the gnwa alignment as well as the gnwa structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_gnwa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads);


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for gnwa.
projectA_algorithm_t* projectA_get_gnwa();


// PRE:     gnwa
//      gnwa:       Pointer to an existing gnwa project.
// POST:    gnwa
//      gnwa:       Nullpointer and gnwa has been destroyed.
void projectA_gnwa_destroy(projectA_algorithm_t* gnwa);


#endif // PROJECTA_GNWA_HPP
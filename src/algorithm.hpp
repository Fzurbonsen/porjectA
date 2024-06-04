/*

    projectA:
    algorithm.hpp
    This file holds the definition for the algorithm struct.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <vector>
#include <utility>
#include <string>

#include "graph.hpp"
#include "alignment.hpp"

using namespace std;

#ifndef PROJECTA_ALGORITHM_HPP
#define PROJECTA_ALGORITHM_HPP


// Struct to hold the input into an algorithm
struct projectA_algorithm_input_t {

    string read; // Read to be aligned
    projectA_hash_graph_t* graph; // Graph on which to align the read

};


// Strutct to hold the funciton calls to use an alignment algorithm.
struct projectA_algorithm_t {

    // Void pointer that holds the function pointer to initialize an algorithm. 
    // This function gets passed a vector of reads with corresponding graphs that should be aligned.
    // This function gets a void* return value that can be used to store any relevant information.
    void* (*init)(vector<projectA_algorithm_input_t>&, int32_t); 

    // Void pointer that holds the function pointer to start the calculation of the initialized batch.
    // This function gets passed the void* return form the init function.
    // This funciotn gets a void* return to hold the results from the alignment.
    void* (*calculate_batch)(void*, int32_t);

    // Void pointer that holds the function pointer to execute the post of the alignment.
    // This function gets passed the void* return from the calculate_batch function.
    // This funciont returns a vector that holds the results from the alignment.
    void (*post)(void*, vector<projectA_alignment_t*>&, int32_t);

};


#endif // PROJECTA_ALGORITHM_HPP
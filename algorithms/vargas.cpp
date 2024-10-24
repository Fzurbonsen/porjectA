/*

    projectA:
    vargas.cpp
    This file holds the implementation of the connector between projectA and vargas.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <string>

// #include "algorithms/vargas.hpp"
#include "algorithm.hpp"
#include "graph.hpp"
#include "file_io.hpp"


void* projectA_vargas_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    // TODO
    return nullptr;
}

void* projectA_vargas_calculate_batch(void* ptr, int32_t thread_index) {
    // TODO
    return nullptr;
}

void projectA_vargas_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    // TODO
}

// Function to get the algorithm sturct for abPOA
projectA_algorithm_t* projectA_get_vargas() {
    // Create new object
    projectA_algorithm_t* vargas = new projectA_algorithm_t;

    // Assign function pointers
    vargas->init = projectA_vargas_init;
    vargas->calculate_batch = projectA_vargas_calculate_batch;
    vargas->post = projectA_vargas_post;

    return vargas;
}

void projectA_vargas_destroy(projectA_algorithm_t* vargas) {
    // TODO
}
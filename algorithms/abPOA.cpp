/*

    projectA:
    abPOA.cpp
    This file holds the implementation of the connector between projectA and abPOA.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <string>

#include "algorithms/abPOA.hpp"
#include "algorithm.hpp"
#include "graph.hpp"
#include "file_io.hpp"

void* projectA_abpoa_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    // TODO
    return nullptr;
}

void* projectA_abpoa_calculate_batch(void* ptr, int32_t thread_index) {
    // TODO
    return nullptr;
}

void projectA_abpoa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    // TODO
}

// Function to get the algorithm sturct for abPOA
projectA_algorithm_t* projectA_get_abpoa() {
    // Create new object
    projectA_algorithm_t* abpoa = new projectA_algorithm_t;

    // Assign function pointers
    abpoa->init = projectA_abpoa_init;
    abpoa->calculate_batch = projectA_abpoa_calculate_batch;
    abpoa->post = projectA_abpoa_post;

    return abpoa;
}

void projectA_abpoa_destroy(projectA_algorithm_t* abpoa) {
    // TODO
}
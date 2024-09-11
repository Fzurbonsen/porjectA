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

projectA_algorithm_t* projectA_get_abpoa() {
    // TODO
    return nullptr;
}

void projectA_abpoa_destroy(projectA_algorithm_t* abpoa) {
    // TODO
}
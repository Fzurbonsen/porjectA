/*

    projectA:
    edlib.hpp
    This file holds the definitions for the connector between proejctA and edlib.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <vector>
#include <string>

#include "edlib/edlib.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"

using namespace std;

#ifndef PROJECTA_EDLIB_HPP
#define PROJECTA_EDLIB_HPP

// PRE:     alignment
//      alignment:      Pointer to a projectA_alignment_t.
// POST:    alignment
//      alignment:      Pointer to a projectA_alignment with the full CIGAR of the given path.
void projectA_edlib(projectA_alignment_t* alignment);

#endif
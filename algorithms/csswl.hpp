/*

    projectA:
    csswl.hpp
    This file holds the definitions for the connector between projectA and csswl.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <vector>
#include <string>

#include "csswl/ssw.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "alignment.hpp"

using namespace std;

#ifndef PROJECTA_CSSWL_HPP
#define PROJECTA_CSSWL_HPP

// As this is only a S2S algorithm it is not housed in the algorithm struct.

// PRE:     alignment
//      alignment:      Pointer to a projectA_alignment_t.
// POST:    alignment
//      alignment:      Pointer to a projectA_alignment with the full CIGAR of the given path.
void projectA_csswl(projectA_alignment_t* alignment);

#endif
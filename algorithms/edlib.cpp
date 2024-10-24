/*

    projectA:
    edlib.cpp
    This file holds the implementation of the connector between projectA and edlib.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <iostream>
#include <string.h>
#include <stdio.h>

#include "algorithms/edlib.hpp"
#include "file_io.hpp"

using namespace std;

void projectA_edlib(projectA_alignment_t* alignment) {

    // Define local variables
    auto& read = alignment->read;
    auto& reference = alignment->reference;
    auto& cigar = alignment->cigar_string;
    cigar.elements.clear();
    cigar.len = 0;

    EdlibAlignResult result = edlibAlign(read.c_str(), strlen(read.c_str()), 
                                            reference.c_str(), strlen(reference.c_str()), 
                                            edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    char* edlib_cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    string cigar_str = edlib_cigar;
    int32_t score = result.editDistance;
    free(edlib_cigar);
    edlibFreeAlignResult(result);
    cigar = projectA_parse_cigar_string(cigar_str);
    alignment->offset = 0;
    alignment->score = score;
}
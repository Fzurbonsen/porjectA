/*

    projectA:
    csswl.cpp
    This file holds the implementation of the connector between projectA and csswl.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <iostream>
#include <string.h>
#include <stdio.h>

#include "algorithms/csswl.hpp"
#include "csswl/kseq.h"
#include "file_io.hpp"

using namespace std;


// Function to use csswl with an alignment from projectA
void projectA_csswl(projectA_alignment_t* alignment) {

    // Define local variables
    auto& read = alignment->read;
    auto& reference = alignment->reference;
    auto& cigar = alignment->cigar_string;
    cigar.elements.clear();
    cigar.len = 0;

    // Define parameters
    int32_t l, m, k;
    // int32_t match = 1, mismatch = 1000, gap_open = 0, gap_extension = 1000;
    int32_t match = 1, mismatch = 127, gap_open = 1, gap_extension = 127;
    // int32_t match = 1, mismatch = -2, gap_open = 2, gap_extension = 1;

    // Reference and read sequences
    const char* ref_seq = reference.c_str();
    const char* read_seq = read.c_str();

    s_profile* profile;
    int read_len = strlen(read_seq);
    int ref_len = strlen(ref_seq);

    // Allocate memory for numeric sequences
    int8_t* num = (int8_t*)malloc(read_len * sizeof(int8_t));   // the read sequence represented in numbers
    int8_t* ref_num = (int8_t*)malloc(ref_len * sizeof(int8_t)); // the reference sequence represented in numbers
    s_align* result;

    /* This table is used to transform nucleotide letters into numbers. */
    static const int8_t nt_table[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    // Initialize scoring matrix for genome sequences
    //  A  C  G  T  N (or other ambiguous code)
    //  2 -2 -2 -2  0  A
    // -2  2 -2 -2  0  C
    // -2 -2  2 -2  0  G
    // -2 -2 -2  2  0  T
    //  0  0  0  0  0  N (or other ambiguous code)
    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : -mismatch;  /* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    // Convert the read sequence to numerical values
    for (m = 0; m < read_len; ++m) num[m] = nt_table[(int)read_seq[m]];

    // Initialize the profile with the numerical read sequence
    profile = ssw_init(num, read_len, mat, 5, 2);

    // Convert the reference sequence to numerical values
    for (m = 0; m < ref_len; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

    // Perform the alignment
    result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension, 1, 0, 0, 15);

    // Iterate over the cigar
    cigar.len = result->cigarLen;
    for (int i = 0; i < result->cigarLen; ++i) {
        // Create element
        projectA_cigar_element_t element;
        element.len = cigar_int_to_len(result->cigar[i]);
        element.type = cigar_int_to_op(result->cigar[i]);
        cigar.elements.push_back(element);
    }

    // Copy offset
    alignment->offset = result->ref_begin1;
    alignment->read_start = result->read_begin1;

    // Cleanup
    align_destroy(result);
    init_destroy(profile);
    free(mat);
    free(ref_num);
    free(num);
}
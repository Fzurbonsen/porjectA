/*

    projectA:
    ksw2.cpp
    This file holds the implementation of the connector between projectA and ksw2.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <iostream>
#include <string.h>
#include <stdio.h>

#include "algorithms/ksw2.hpp"
#include "file_io.hpp"

using namespace std;


// Function to use ksw2 with an alignment from projectA
void projectA_ksw2(projectA_alignment_t* alignment) {

    // Define local variables
    auto& read = alignment->read;
    auto& reference = alignment->reference;
    auto& cigar = alignment->cigar_string;
    cigar.elements.clear();
    cigar.len = 0;

    // Prepare parameters for ksw2
    const char *tseq = read.c_str();
    const char *qseq = reference.c_str();
    int sc_mch = 1;
    int sc_mis = -2;
    int gapo = 2;
    int gape = 1;
    // int sc_mch = 1;
    // int sc_mis = 1000;
    // int gapo = 0;
    // int gape = 1000;


    int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	// for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
	// 	fprintf(stderr, "%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
    // fprintf(stderr, "\n");

    // Store CIGAR
    cigar.len = ez.n_cigar;
    for (i = 0; i < ez.n_cigar; ++i) {
        // Create element
        projectA_cigar_element_t element;
        element.len = ez.cigar[i]>>4;
        element.type = "MID"[ez.cigar[i]&0xf];
        cigar.elements.push_back(element);
    }

	free(ez.cigar); free(ts); free(qs);

}
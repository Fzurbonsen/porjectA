/*

    projectA:
    alignment.cpp
    This file holds the implementation for the projectA alignment struts helperfunctions.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include "alignment.hpp"
#include "file_io.hpp"


// Function to compare two CIGARs
bool projectA_compare_cigar(projectA_cigar_t& cigar1, projectA_cigar_t& cigar2) {

    // Check whether the two CIGARs have the same size
    if (cigar1.len != cigar2.len) return false;

    // Iterate over all the CIGAR elements
    for (int i = 0; i < cigar1.len; ++i) {
        
        // Compare element length and type
        if (!(cigar1.elements[i].len == cigar2.elements[i].len
            && cigar1.elements[i].type == cigar2.elements[i].type)) return false;
    }

    return true;
}


// Function to compare two alignments
int projectA_compare_alignments(bool print, FILE* file, projectA_alignment_t* alignment1, 
                                                            projectA_alignment_t* alignment2) {

    // If print is enabled we check wether the file pointer is valid
    if (print && file == nullptr) {
        cerr << "Error: invalid file pointer!\n";
        exit(1);
    }

    // // Compare score, offset and size
    // if (!(alignment1->offset == alignment2->offset
    //     && alignment1->score == alignment2->score
    //     && alignment1->size == alignment2->size)) {
        
    //     // If print is enabled we print the mismatch to file
    //     if (print) {
    //         fprintf(file, "Mismatched alignment:\n");
    //         fprintf(file, "\t");
    //         projectA_print_alignment(file, alignment1);
    //         fprintf(file, "\n\t");
    //         projectA_print_alignment(file, alignment2);
    //         fprintf(file, "\n");
    //     }

    //     // If the alignment parameters don't match we return zero to indicate a mismatch in the alignments
    //     return 0;
    // }

    // Compare offset and size
    if (!(alignment1->size == alignment2->size)) {
        
        // If print is enabled we print the mismatch to file
        if (print) {
            fprintf(file, "Mismatched alignment:\n");
            fprintf(file, "size\n");
            fprintf(file, "\t");
            projectA_print_alignment(file, alignment1);
            fprintf(file, "\n\t");
            projectA_print_alignment(file, alignment2);
            fprintf(file, "\n");
        }

        // If the alignment parameters don't match we return zero to indicate a mismatch in the alignments
        return 0;
    }

    // Iterate over all CIGARS and NODES
    for (int i = 0; i < alignment1->size; ++i) {

        // Check if the two node ids are identical
        if (alignment1->nodes[i] != alignment2->nodes[i]) {

            // If print is enabled we print the mismatch to file
            if (print) {
                fprintf(file, "Mismatched alignment:\n");
                fprintf(file, "nodes\n");
                fprintf(file, "\t");
                projectA_print_alignment(file, alignment1);
                fprintf(file, "\n\t");
                projectA_print_alignment(file, alignment2);
                fprintf(file, "\n");
            }

            // If the referenced nodes don't match we retunr zero to indicate a mismatch in the alignments
            return 0;
        }


        // // Check if the two CIGARs are identical
        // if (!projectA_compare_cigar(alignment1->cigar[i], alignment2->cigar[i])) {

        //     // If print is enabled we print the mismatch to file
        //     if (print) {
        //         fprintf(file, "Mismatched alignment:\n");
        //         fprintf(file, "cigar\n");
        //         fprintf(file, "\t");
        //         projectA_print_alignment(file, alignment1);
        //         fprintf(file, "\n\t");
        //         projectA_print_alignment(file, alignment2);
        //         fprintf(file, "\n");
        //     }

        //     // If the two CIGARs don't match we return a zero to indicate a mismatch in the alignments
        //     return 0;
        // }
    }

    return 1;
}


// Function to compare the path of two alignments
int projectA_compare_alignments_path(bool print, FILE* file, projectA_alignment_t* alignment1, 
                                                            projectA_alignment_t* alignment2) {
    
    // If print is enabled we check wether the file pointer is valid
    if (print && file == nullptr) {
        cerr << "Error: invalid file pointer!\n";
        exit(1);
    }

    int size = (alignment1->size < alignment2->size) ? alignment1->size : alignment2->size;
    
    // Iterate over the nodes and check whether they are the same
    for (int i = 0; i < size; ++i) {
        if (alignment1->nodes[i] != alignment2->nodes[i]) {

            // If print is enabled we print the mismatch to file
            if (print) {
                fprintf(file, "Mismatched alignment:\n");
                fprintf(file, "nodes\n");
                fprintf(file, "\t");
                projectA_print_alignment(file, alignment1);
                fprintf(file, "\n\t");
                projectA_print_alignment(file, alignment2);
                fprintf(file, "\n");
            }

            cerr << alignment1->nodes[i] << "\t" << alignment2->nodes[i] << endl;

            return 0;
        }
    }

    return 1;
}


/*

    projectA:
    alignment.cpp
    This file holds the implementation for the projectA alignment struts helperfunctions.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include "alignment.hpp"
#include "file_io.hpp"


// Function to concatenate two CIGARs
void projectA_concat_cigar(projectA_cigar_t* cigar1, projectA_cigar_t* cigar2) {

    // Increase length of CIGAR
    cigar1->len = cigar1->len + cigar2->len;

    // Add CIGAR2 to CIGAR1
    for (int32_t i = 0; i < cigar2->len; ++i) {
        cigar1->elements.push_back(cigar2->elements[i]);
    }
}


// Function to generate aligned pairs
set<int32_t> projectA_generate_aligned_positions(projectA_cigar_t* cigar) {
    // Variables to keep track of positions
    int32_t ref_pos = 0;
    int32_t query_pos = 0;
    set<int32_t> aligned_positions;

    // Iterate over CIGAR
    for (auto& element : cigar->elements) {
        
        // Differentiate between element types
        switch (element.type) {

            // Match operation is both in the reference and the query
            case 'M' :
                // Iterate over element length to add the positions to the set
                for (int32_t i = 0; i < element.len; ++i) {
                    aligned_positions.insert(ref_pos + i);
                }
                ref_pos += element.len;
                query_pos += element.len;
                break;

            // Insertion is only in the query
            case 'I' :
                // Iterate over element length to add the positions to the set
                for (int32_t i = 0; i < element.len; ++i) {
                    aligned_positions.insert(query_pos + i);
                }
                query_pos += element.len;
                break;
            
            // Deletion is only in the reference
            case 'D' :
                // Iterate over element length to add the positions to the set
                for (int32_t i = 0; i < element.len; ++i) {
                    aligned_positions.insert(ref_pos + i);
                }
                ref_pos += element.len;
                break;

            // Skipped region is only in the reference
            case 'N' :
                // Iterate over element length to add the positions to the set
                for (int32_t i = 0; i < element.len; ++i) {
                    aligned_positions.insert(ref_pos + i);
                }
                ref_pos += element.len;
                break;

            default :
                cerr << "Error: element unknown!\n";
                cerr << "Element: " << element.type << " is not known!\n";
                exit(1);
        }
    }

    return aligned_positions;
}


// Function to calculate the Jaccard index between two sets
double projectA_jaccard_index (set<int32_t>& set1, set<int32_t>& set2) {
    set<int32_t> intersection;
    set<int32_t> union_set;

    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                        inserter(intersection, intersection.begin()));
    set_union(set1.begin(), set1.end(), set2.begin(), set2.end(),
                inserter(union_set, union_set.begin()));

    return !union_set.empty() ? static_cast<double>(intersection.size()) / union_set.size() : 0.0;
}


// Function to calculate the overlap of two CIGARs
double projectA_cigar_accuracy(projectA_cigar_t* cigar1, projectA_cigar_t* cigar2) {
    set<int32_t> aligned_positions1 = projectA_generate_aligned_positions(cigar1);
    set<int32_t> aligned_positions2 = projectA_generate_aligned_positions(cigar2);

    return projectA_jaccard_index(aligned_positions1, aligned_positions2);
}


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


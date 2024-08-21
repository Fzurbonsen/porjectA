/*

    projectA:
    gssw.cpp
    This file holds the implementation of the connector between projectA and gssw.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cstring>
#include <cstdint>

#include "algorithms/gssw.hpp"
#include "gssw/gssw.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "file_io.hpp"

using namespace std;


// Constructor for parameter struct
projectA_gssw_parameters_t::projectA_gssw_parameters_t(gssw_graph* graph, const char* read, int8_t* nt_table,
                                                        int8_t* mat, uint8_t gap_open, uint8_t gap_extension,
                                                        projectA_hash_graph_t* projectA_hash_graph) :
                                                        graph(graph), read(read), nt_table(nt_table), mat(mat),
                                                        gap_open(gap_open), gap_extension(gap_extension),
                                                        projectA_hash_graph(projectA_hash_graph) {}


// Function to convert projectA_hash_graph_t into gssw_graph
gssw_graph* projectA_hash_graph_to_gssw_graph(projectA_hash_graph_t* in_graph, int8_t* nt_table, 
                                                int8_t* mat, uint8_t gap_open, uint8_t gap_extension) {

    unordered_map<projectA_node_t*, gssw_node*> node_map;

    gssw_sse2_disable();

    gssw_node* nodes[in_graph->n_nodes];

    // Iterate over all nodes in the graph to create the corresponding gssw nodes
    int i = 0;
    for (auto& curr : in_graph->nodes_in_order) {

        // Fill node
        nodes[i] = (gssw_node*)gssw_node_create(NULL, i, curr->seq.c_str(), nt_table, mat);

        // Add to node map
        node_map[curr] = nodes[i];
        i++;
    }

    // Iterate over all nodes
    for (auto& curr : in_graph->nodes_in_order) {
        
        // Iterate over all outgoing edges
        for (auto& next : curr->next) {
            gssw_nodes_add_edge(node_map[curr], node_map[next]);
        }
    }

    // Create gssw graph
    gssw_graph* out_graph = gssw_graph_create(in_graph->n_nodes);

    // Iterate over all nodes
    for (i = 0; i < in_graph->n_nodes; ++i) {
        // Add node to graph
        gssw_graph_add_node(out_graph, nodes[i]);
    }

    // Check if graph sizes match
    if (out_graph->size != in_graph->n_nodes) {
        cerr << "Error: graph size does not match!\n";
        exit(1);
    }

    return out_graph;
}


// Function to count number of matches, mismatches, insertions, deletions, non-aligned
void procetA_gssw_count_cigar_op(gssw_cigar* cigar, int32_t& matches, int32_t& mismatches,
                                                    int32_t& insertions, int32_t& deletions,
                                                    int32_t& length) {

    // Iterate over cigar elements
    for (int32_t i = 0; i < cigar->length; ++i) {
        
        // Switch case to distinguish the operations
        switch (cigar->elements[i].type) {
            case 'I' :
                insertions += cigar->elements[i].length;
                break;
            case 'D' :
                deletions += cigar->elements[i].length;
                break;
            case 'M' : 
                matches += cigar->elements[i].length;
                break;
            case 'X' :
                mismatches += cigar->elements[i].length;
                break;
            case 'S' :
                break;
            case 'N' :
                break;
            default :
                cerr << "Error: Invalid CIGAR element!\n";
                cerr << "Element: " << cigar->elements[i].type << " not known!\n";
                exit(1);
        }
        length += cigar->elements[i].length;
    }


}

// Function to convert gssw cigar to projectA cigar
projectA_cigar_t projectA_gssw_get_cigar(gssw_cigar* cigar) {

    projectA_cigar_t projectA_cigar;

    // Iterate over cigar elements
    for (int i = 0; i < cigar->length; ++i) {
        // Local variables
        projectA_cigar_element_t projectA_cigar_element;

        // Convert CIGAR elements
        switch (cigar->elements[i].type) {
            case 'I' :
                projectA_cigar_element.type = 'I';
                break;
            case 'D' :
                projectA_cigar_element.type = 'D';
                break;
            case 'M' : 
                projectA_cigar_element.type = 'M';
                break;
            case 'X' :
                projectA_cigar_element.type = 'M';
                break;
            case 'S' :
                projectA_cigar_element.type = 'M';
                break;
            case 'N' :
                projectA_cigar_element.type = 'N';
                break;
            default :
                cerr << "Error: Invalid CIGAR element!\n";
                cerr << "Element: " << cigar->elements[i].type << " not known!\n";
                exit(1);
        }

        // Copy cigar length
        projectA_cigar_element.len = cigar->elements[i].length;

        if (projectA_cigar.elements.size() ==  0) {

            // Push element to CIGAR vector
            projectA_cigar.elements.push_back(projectA_cigar_element);

        } else {

            // Check if last element and current element are the same
            if (projectA_cigar_element.type == projectA_cigar.elements.back().type) {

                // Increase last element length
                projectA_cigar.elements.back().len = projectA_cigar_element.len + projectA_cigar.elements.back().len;

            } else {
                
                // Push element to CIGAR vector
                projectA_cigar.elements.push_back(projectA_cigar_element);

            }
        }
    }

    projectA_cigar.len = projectA_cigar.elements.size();

    return projectA_cigar;
}


// Function to convert the gssw graph mapping to the projectA_alignment_t struct
void projectA_gssw_graph_mapping_to_alignment(projectA_hash_graph_t* graph, gssw_graph_mapping* gm, projectA_alignment_t* alignment) {

    // Prepare CIGAR string field
    auto& cigar = alignment->cigar_string;
    cigar.len = 0;
    cigar.operations_length = 0;
    alignment->n_matches = 0;
    alignment->n_mismatches = 0;

    // Copy offset and score
    alignment->offset = gm->position;
    alignment->score = gm->score;
    alignment->size = gm->cigar.length;

    // Reserve memory for alignment vectors
    alignment->nodes.reserve(alignment->size);
    alignment->cigar.reserve(alignment->size);

    // Iterate over the graph mapping graph CIGAR
    for (int i = 0; i < alignment->size; ++i) {
        auto& node_cigar = gm->cigar.elements[i];

        // Add node to the nodes vector
        alignment->nodes.push_back(graph->nodes_in_order[node_cigar.node->id]->id);

        // Add CIGAR to the CIGARs vector
        projectA_cigar_t new_cigar = projectA_gssw_get_cigar(node_cigar.cigar);
        alignment->cigar.push_back(new_cigar);

        // Update reference string
        alignment->reference = alignment->reference + node_cigar.node->seq;

        // Add CIGAR to total CIGAR
        projectA_concat_cigar(&cigar, &new_cigar);

        // Count matches and mismatches
        int32_t n_insertions;
        int32_t n_deletions;
        procetA_gssw_count_cigar_op(node_cigar.cigar, alignment->n_matches, alignment->n_mismatches, 
                                                    n_insertions, n_deletions, cigar.operations_length);
    }
}


// Function to initialize gssw
void* projectA_gssw_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    // Create our io vectors for the gssw algorithm
    projectA_gssw_io_t* out = new projectA_gssw_io_t;
    for (int i = 0; i < numThreads; ++i) {
        vector<projectA_gssw_parameters_t> params;
        vector<gssw_graph_mapping*> gms;
        params.reserve(alignments.size()/numThreads);
        gms.reserve(alignments.size()/numThreads);
        out->parameters.push_back(params);
        out->gms.push_back(gms);
    }

    // Assign size
    out->size = alignments.size();

    // Create the specifications for gssw.
    // int8_t match = 1;
    // int8_t mismatch = 1000;
    // uint8_t gap_open = 0;
    // uint8_t gap_extension = 1000;

    // int8_t match = 1;
    // int8_t mismatch = 4;
    // uint8_t gap_open = 4;
    // uint8_t gap_extension = 4;

    int8_t match = 1;
    int8_t mismatch = 1;
    uint8_t gap_open = 1;
    uint8_t gap_extension = 1;

    gssw_sse2_disable();
    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);

    // Iterate over the input graphs
    int32_t thread_index = 0;
    for (auto& itr : alignments) {
        // Construct gssw graph
        gssw_graph* new_gssw = projectA_hash_graph_to_gssw_graph(itr->graph, nt_table, mat, gap_open, gap_extension);

        // Construct parameter entry
        projectA_gssw_parameters_t entry(new_gssw, itr->read.c_str(), nt_table, mat, gap_open, gap_extension, itr->graph);

        // Append to parameter vector
        out->parameters[thread_index].push_back(entry);
        
        // Update thread index
        thread_index = (thread_index == numThreads - 1) ? 0 : thread_index + 1;
    }

    // Cast return value into void pointer
    return static_cast<void*>(out);
}


// Function to calculate a batch in gssw
void* projectA_gssw_calculate_batch(void* ptr, int32_t thread_index) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gssw_io_t*>(ptr);
    auto& parameters = input->parameters[thread_index];
    auto& gms = input->gms[thread_index];

    // Run gssw on every set of parameters
    for (int i = 0; i < parameters.size(); ++i) {
        auto& parameter = parameters[i];
        gssw_graph_fill(parameter.graph, parameter.read, parameter.nt_table, parameter.mat, 
                        parameter.gap_open, parameter.gap_extension, 0, 0, 15, 2, true);
        gms.push_back(gssw_graph_trace_back(parameter.graph,
                                            parameter.read,
                                            strlen(parameter.read),
                                            parameter.nt_table,
                                            parameter.mat,
                                            parameter.gap_open,
                                            parameter.gap_extension,
                                            0, 0));
    }

    // Cast return value into void pointer
    return static_cast<void*>(ptr);
}


// Function to handle post of gssw
void projectA_gssw_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gssw_io_t*>(ptr);
    auto& parameters_vec = input->parameters;
    auto& gms_vec = input->gms;

    // Iterate over all entries in the graph mapping vectors
    int32_t thread_index = 0;
    int32_t j = 0;
    for (int i = 0; i < input->size; ++i) {
        auto& gms = gms_vec[thread_index];
        auto& parameters = parameters_vec[thread_index];

        // Check for out of bounds error
        if (j >= gms.size() || j >= parameters.size()) {
            cerr << "Error: The index is out of bound!\n";
            cerr << "index: " << j << "\tin vector of the thread: " << thread_index << endl;
            cerr << "at value: " << i << endl;
            exit(1);
        }

        // Update alignments
        projectA_gssw_graph_mapping_to_alignment(parameters[j].projectA_hash_graph, gms[j], alignments[i]);

        // Update thread index
        if (thread_index == numThreads - 1) {
            thread_index = 0;
            ++j;
        } else {
            ++thread_index;
        }
    }

    // Free memory allocated for mappings
    for (auto& gms : gms_vec) {
        for (auto& gm : gms) {
            gssw_graph_mapping_destroy(gm);
        }
        gms.clear();
    }
    gms_vec.clear();

    // Free memory allocated for parameters
    free(parameters_vec[0][0].nt_table);
    free(parameters_vec[0][0].mat);
    for (auto& parameters : parameters_vec) {
        for (auto& parameter : parameters) {
            gssw_graph_destroy(parameter.graph);
        }
        parameters.clear();
    }
    parameters_vec.clear();

    // Delete the object that holds the inputs and outputs.
    delete input;
    ptr = nullptr;
}


// Function to get the algorithm struct for gssw
projectA_algorithm_t* projectA_get_gssw() {
    // Create new object
    projectA_algorithm_t* gssw = new projectA_algorithm_t;

    // Assign function pointers
    gssw->init = projectA_gssw_init;
    gssw->calculate_batch = projectA_gssw_calculate_batch;
    gssw->post = projectA_gssw_post;

    return gssw;
}


// Function to delete projectA_algorithm_t struct
void projectA_gssw_destroy(projectA_algorithm_t* gssw) {

    // Check if the struct still exists
    if (gssw == nullptr) {
        cerr << "Warning: gssw struct is null pointer!\n";
        return;
    }

    // Delete struct
    delete gssw;
}
/*

    projectA:
    gt_gwfa.cpp
    This file holds the implementation of the connector between projectA and gt_gwfa.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <tuple>
#include <unordered_map>

#include "algorithms/gt_gwfa.hpp"
#include "gt_gwfa/graphs.h"
#include "algorithm.hpp"
#include "graph.hpp"
#include "file_io.hpp"

using namespace std;

// Constructor for parameter struct
projectA_gt_gwfa_parameters_t::projectA_gt_gwfa_parameters_t(gssw_graph* gssw , gwf_graph_t* gwf, const char* r,
                                                                projectA_hash_graph_t* projectA_hash_graph) : 
                                                                gssw(gssw), gwf(gwf), read(r), 
                                                                projectA_hash_graph(projectA_hash_graph) {}


// Function to convert projectA_hash_graph_t into gssw_graph
gssw_graph* projectA_hash_graph_to_gt_gssw_graph(projectA_hash_graph_t* in_graph) {

    // Define parameters for gssw
    int8_t match = 1, mismatch = 4;
    uint8_t gap_open = 6, gap_extension = 1;
    gssw_sse2_disable();
    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);
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

    free(nt_table);
    free(mat);
    return out_graph;
}


// Function to convert gssw cigar to projectA cigar
projectA_cigar_t projectA_gt_gwfa_get_cigar(gssw_cigar* cigar) {

    projectA_cigar_t projectA_cigar;

    // Copy cigar length
    projectA_cigar.len = cigar->length;

    // Reserve cigar size
    projectA_cigar.elements.reserve(projectA_cigar.len);

    // Iterate over cigar elements
    for (int i = 0; i < projectA_cigar.len; ++i) {

        projectA_cigar_element_t projectA_cigar_element;

        // Copy cigar element
        projectA_cigar_element.len = cigar->elements[i].length;
        projectA_cigar_element.type = cigar->elements[i].type;

        projectA_cigar.elements.push_back(projectA_cigar_element);
    }

    return projectA_cigar;
}


// Function to convert the gssw graph mapping to the projectA_alignment_t struct
projectA_alignment_t* projectA_gt_gwfa_graph_mapping_to_alignment(projectA_hash_graph_t* graph, gssw_graph_mapping* gm) {

    // Create new projectA_alignment_t
    projectA_alignment_t* alignment = new projectA_alignment_t;

    // Prepare CIGAR string field
    auto& cigar = alignment->cigar_string;
    cigar.len = 0;

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
        projectA_cigar_t new_cigar = projectA_gt_gwfa_get_cigar(node_cigar.cigar);
        alignment->cigar.push_back(new_cigar);

        // Add CIGAR to total CIGAR
        projectA_concat_cigar(&cigar, &new_cigar);
    }

    return alignment;
};


// Function to initialize gt_gwfa
void* projectA_gt_gwfa_init(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {

    // Create the io vectors for the gt_gwfa algorithm
    projectA_gt_gwfa_io_t* out = new projectA_gt_gwfa_io_t;
    for (int i = 0; i < numThreads; ++i) {
        vector<projectA_gt_gwfa_parameters_t> params;
        vector<gssw_graph_mapping*> gms;
        params.reserve(graphs.size()/numThreads);
        gms.reserve(graphs.size()/numThreads);
        out->parameters.push_back(params);
        out->gms.push_back(gms);
    }

    // Assign size
    out->size = graphs.size();

    // Iterate over the input graphs
    int32_t thread_index = 0;
    for (auto& itr : graphs) {

        // Construct gssw graph
        gssw_graph* new_gssw = projectA_hash_graph_to_gt_gssw_graph(itr.graph);

        // Create gwf graph
        gwf_graph_t* new_gwf = gt_gssw_to_gwf_direct(new_gssw);

        // Construct parameter entry
        projectA_gt_gwfa_parameters_t entry(new_gssw, new_gwf, itr.read.c_str(), itr.graph);

        // Append entry to parameter vector
        out->parameters[thread_index].push_back(entry);

        // Update thread index
        thread_index = (thread_index == numThreads - 1) ? 0 : thread_index + 1;
    }

    // Cast return value into void pointer
    return static_cast<void*>(out);
}


// Function to calculate a batch in gt_gwfa
void* projectA_gt_gwfa_calculate_batch(void* ptr, int32_t thread_index) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gt_gwfa_io_t*>(ptr);
    auto& parameters = input->parameters[thread_index];
    auto& gms = input->gms[thread_index];

    // Run gt_gwfa on every set of parameters
    for (int i = 0; i < parameters.size(); ++i) {
        auto& parameter = parameters[i];
        gms.push_back(gt_read_align_gwf(parameter.gwf, parameter.gssw, parameter.read));
    }
    
    // Cast return value into void pointer
    return static_cast<void*>(ptr);
}


// Function to handle post of gt_gwfa
void projectA_gt_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gt_gwfa_io_t*>(ptr);
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

        alignments.push_back(projectA_gt_gwfa_graph_mapping_to_alignment(parameters[j].projectA_hash_graph, gms[j]));

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
    for (auto& parameters : parameters_vec) {
        for (auto& parameter : parameters) {
            gssw_graph_destroy(parameter.gssw);
            gt_gwf_free(parameter.gwf);
        }
        parameters.clear();
    }
    parameters_vec.clear();

    // Free the input pointer itself
    delete input;
    ptr = nullptr;
}

// Function to get the algorithm struct for gt_gwfa
projectA_algorithm_t* projectA_get_gt_gwfa() {

    // Create new algorithm object
    projectA_algorithm_t* gt_gwfa = new projectA_algorithm_t;

    // Assign function pointers 
    gt_gwfa->init = projectA_gt_gwfa_init;
    gt_gwfa->calculate_batch = projectA_gt_gwfa_calculate_batch;
    gt_gwfa->post = projectA_gt_gwfa_post;

    return gt_gwfa;
}


// Function to delete projectA_algorithm_t struct
void projectA_gt_gwfa_destroy(projectA_algorithm_t* gt_gwfa) {

    // Check if the struct still exists
    if (gt_gwfa == nullptr) {
        cerr << "Warning: gssw struct is null pointer!\n";
        return;
    }

    // Delete struct
    delete gt_gwfa;
}
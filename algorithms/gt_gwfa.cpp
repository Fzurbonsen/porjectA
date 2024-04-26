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

using namespace std;

// Constructor for parameter struct
projectA_gt_gwfa_parameters_t::projectA_gt_gwfa_parameters_t(gssw_graph* gssw , gwf_graph_t* gwf, const char* r) : 
                                                                gssw(gssw), gwf(gwf), read(r) {}


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
        nodes[i] = (gssw_node*)gssw_node_create(NULL, curr->index, curr->seq.c_str(), nt_table, mat);

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


// Function to initialize gt_gwfa
void* projectA_gt_gwfa_init(vector<projectA_algorithm_input_t>& graphs) {

    // Create the io vectors for the gt_gwfa algorithm
    projectA_gt_gwfa_io_t* out = new projectA_gt_gwfa_io_t;

    // Iterate over the input graphs
    for (auto& itr : graphs) {

        // Construct gssw graph
        gssw_graph* new_gssw = projectA_hash_graph_to_gt_gssw_graph(itr.graph);

        // Create gwf graph
        gwf_graph_t* new_gwf = gt_gssw_to_gwf_direct(new_gssw);

        // Construct parameter entry
        projectA_gt_gwfa_parameters_t entry(new_gssw, new_gwf, itr.read.c_str());

        // Append entry to parameter vector
        out->parameters.push_back(entry);
    }

    // Reserve space for results
    out->gms.reserve(out->parameters.size());

    // Cast return value into void pointer
    return static_cast<void*>(out);
}


// Function to calculate a batch in gt_gwfa
void* projectA_gt_gwfa_calculate_batch(void* ptr) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gt_gwfa_io_t*>(ptr);
    auto& parameters = input->parameters;
    auto& gms = input->gms;

    // Run gt_gwfa on every set of parameters
    for (int i = 0; i < parameters.size(); ++i) {
        auto& parameter = parameters[i];
        gms.push_back(gt_read_align_gwf(parameter.gwf, parameter.gssw, parameter.read));
    }
    
    // Cast return value into void pointer
    return static_cast<void*>(ptr);
}


// Function to handle post of gt_gwfa
void projectA_gt_gwfa_post(void* ptr) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gt_gwfa_io_t*>(ptr);
    auto& parameters = input->parameters;
    auto& gms = input->gms;


    // TODO: Create usable results
    for (auto& gm : gms) {
        gssw_print_graph_mapping(gm, stderr);
    }



    // Free memory allocated for mappings
    for (auto& gm : gms) {
        gssw_graph_mapping_destroy(gm);
    }
    gms.clear();

    // Free memory allocated for parameters
    for (auto& parameter : parameters) {
        gssw_graph_destroy(parameter.gssw);
        gt_gwf_free(parameter.gwf);
    }
    parameters.clear();

    // Free the input pointer itself
    delete input;
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
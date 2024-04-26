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

using namespace std;


// Constructor for parameter struct
projectA_gssw_parameters_t::projectA_gssw_parameters_t(gssw_graph* graph, const char* read, int8_t* nt_table,
                                                        int8_t* mat, uint8_t gap_open, uint8_t gap_extension) :
                                                        graph(graph), read(read), nt_table(nt_table), mat(mat),
                                                        gap_open(gap_open), gap_extension(gap_extension) {}


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

    return out_graph;
}


// Function to initialize gssw
void* projectA_gssw_init(vector<projectA_algorithm_input_t>& graphs) {

    // Create our io vectors for the gssw algorithm
    projectA_gssw_io_t* out = new projectA_gssw_io_t;

    // Create the specifications for gssw.
    int8_t match = 1;
    int8_t mismatch = 4;
    uint8_t gap_open = 6;
    uint8_t gap_extension = 1;
    gssw_sse2_disable();
    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);


    // Iterate over the input graphs
    for (auto& itr : graphs) {
        // Construct gssw graph
        gssw_graph* new_gssw = projectA_hash_graph_to_gssw_graph(itr.graph, nt_table, mat, gap_open, gap_extension);

        // Construct parameter entry
        projectA_gssw_parameters_t entry(new_gssw, itr.read.c_str(), nt_table, mat, gap_open, gap_extension);

        // Append to parameter vector
        out->parameters.push_back(entry);
    }

    // Reserve space for results
    out->gms.reserve(out->parameters.size());

    // Cast return value into void pointer
    return static_cast<void*>(out);
}


// Function to calculate a batch in gssw
void* projectA_gssw_calculate_batch(void* ptr) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gssw_io_t*>(ptr);
    auto& parameters = input->parameters;
    auto& gms = input->gms;

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
void projectA_gssw_post(void* ptr) {

    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gssw_io_t*>(ptr);
    auto& parameters = input->parameters;
    auto& gms = input->gms;


    // TODO: Postprocessing data
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
        gssw_graph_destroy(parameter.graph);
    }
    free(parameters[0].nt_table);
    free(parameters[0].mat);
    parameters.clear();

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
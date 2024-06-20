/*

    projectA:
    gwfa.cpp
    This file holds the implementation of the connector between projectA and gwfa.
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

#include "algorithms/gwfa.hpp"
#include "gwfa/gwfa.h"
#include "gwfa/gfa.h"
#include "gwfa/gfa-priv.h"
#include "gwfa/ketopt.h"
#include "gwfa/kalloc.h"
#include "gwfa/kseq.h"
#include "algorithm.hpp"
#include "graph.hpp"

using namespace std;


// Constructor for parameter struct
projectA_gwfa_parameters_t::projectA_gwfa_parameters_t(void* km, gwf_graph_t* graph, int32_t ql, const char* q, int32_t v0,
                                                        int32_t v1, uint32_t max_lag, int32_t traceback, 
                                                        projectA_hash_graph_t* projectA_hash_graph) :
                                                        km(km), graph(graph), ql(ql), q(q), v0(v0), v1(v1), 
                                                        max_lag(max_lag), traceback(traceback),
                                                        projectA_hash_graph(projectA_hash_graph) {}


// Function to gree gwf graph
void gwf_free(gwf_graph_t *g) {
	int32_t i;
	for (i = 0; i < g->n_vtx; ++i) free(g->seq[i]);
	free(g->len); free(g->seq); free(g->arc); free(g->src); free(g);
    g = nullptr;
}


// Function to transform hash graph to gwf graph
gwf_graph_t* projectA_hash_graph_to_gwf_graph(projectA_hash_graph_t* in_graph) {

    // Map to get node id's for edges
    unordered_map<projectA_node_t*, int32_t> node_map;

    // Initialize graph
    gwf_graph_t* graph
    GFA_CALLOC(graph, 1);
    graph->n_vtx = in_graph->n_nodes;
    graph->n_arc = in_graph->n_edges;

    // Allocate memory for graph struct fields
    GFA_MALLOC(graph->len, graph->n_vtx);
    GFA_MALLOC(graph->src, graph->n_vtx);
    GFA_MALLOC(graph->seq, graph->n_vtx);
    GFA_MALLOC(graph->arc, graph->n_arc);

    // Iterate over all nodes to create gwf nodes
    int32_t i, k;
    for (i = k = 0; i < in_graph->nodes_in_order.size(); ++i) {
        // Assign local variables
        auto& node = in_graph->nodes_in_order[i];
        uint32_t v;
        uint32_t len = node->len; // Sequence length
        uint32_t count_prev = node->prev.size(); // Number of outgoing edges
        uint32_t edge_origin_id;

        // Add node to node_map
        node_map[node] = i;

        // Assign node length
        graph->len[i] = len;

        // Allocate memory for sequence
        GFA_MALLOC(graph->seq[i], len + 1);

        // Copy sequence
        for (int32_t j = 0; j < len; ++j) {
            graph->seq[i][j] = node->seq[j];
        }
        graph->seq[i][len] = 0; // 0 terminated string

        // Iterate over all incoming edges
        for (int32_t j = 0; j < count_prev; ++j, ++k) {
            edge_origin_id = node_map[node->prev[j]];
            graph->arc[k].a = (uint64_t)edge_origin_id<<32 | i;
            graph->arc[k].o = 0;
        }
    }

    return graph;
}


// Function to convert path to alignment
projectA_alignment_t* projectA_gwfa_path_to_alignment(projectA_hash_graph_t* graph, gwf_path_t* path) {

    // Create new projectA_alignment_t
    projectA_alignment_t* alignment = new projectA_alignment_t;

    // Copy alignment size
    alignment->size = path->nv;

    // Reserve memory for alignmnet vectors
    alignment->nodes.reserve(alignment->size);
    alignment->cigar.reserve(alignment->size);

    // Iterate over the path
    for (int i = 0; i < alignment->size; ++i) {
        auto& node = path->v[i];

        // Add node to the nodes vector
        alignment->nodes.push_back(graph->nodes_in_order[node]->id);
    }

    return alignment;
}


// Function to initialize gwfa
void* projectA_gwfa_init(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {

    // Create io vectors for the gwfa algorithm
    projectA_gwfa_io_t* out = new projectA_gwfa_io_t;
    for (int i = 0; i < numThreads; ++i) {
        // Prepare vectors
        vector<projectA_gwfa_parameters_t> params;
        vector<projectA_gwfa_path_t> paths;

        // Reserve memory
        params.reserve(graphs.size()/numThreads);
        paths.reserve(graphs.size()/numThreads);

        // Push vectors to I/O struct
        out->parameters.push_back(params);
        out->paths.push_back(paths);
    }

    // Assign size
    out->size =  graphs.size();

    // Create the specifications for gwfa.
    int32_t v0 = 0;
    int32_t v1 = -1;
    uint32_t max_lag = 0;
    int32_t traceback = 1;

    // Iterate over the input graphs
    int32_t thread_index = 0;
    for (auto& itr : graphs) {
        // Construct gwfa graph
        gwf_graph_t* new_gwf_graph = projectA_hash_graph_to_gwf_graph(itr.graph);

        // Construct parameter entry
        projectA_gwfa_parameters_t entry(km_init(), new_gwf_graph, itr.read.size(), itr.read.c_str(), v0, v1, max_lag, traceback, itr.graph);

        // Append entry to parameter vector
        out->parameters[thread_index].push_back(entry);

        // Update thread index
        thread_index = (thread_index == numThreads - 1) ? 0 : thread_index + 1;
    }

    // Cast return value int void pointer
    return static_cast<void*>(out);
}


// Function to calculate a batch in gwfa
void* projectA_gwfa_calculate_batch(void* ptr, int32_t thread_index) {

    // Check if ptr is valid
    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gwfa_io_t*>(ptr);
    auto& parameters = input->parameters[thread_index];
    auto& paths = input->paths[thread_index];

    // Run gwfa on every set of parameters
    for (int32_t i = 0; i < parameters.size(); ++i) {
        // Create local variables
        auto& parameter = parameters[i];
        projectA_gwfa_path_t path;
        int32_t score;
        gwf_path_t gwf_path;

        // Run gwfa
        gwf_ed_index(parameter.km, parameter.graph);
        score = gwf_ed(parameter.km, parameter.graph, parameter.ql, parameter.q, parameter.v0, parameter.v1,
                        parameter.max_lag, parameter.traceback, &gwf_path);

        path.score = score;
        path.path = gwf_path;

        paths.push_back(path);
    }
    
    // Cast return value into void pointer
    return static_cast<void*>(ptr);
}


// Function to handle post of gwfa
void projectA_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    // Check if ptr is valid
    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_gwfa_io_t*>(ptr);
    auto& parameters_vec = input->parameters;
    auto& paths_vec = input->paths;

    // Iterate over all entries in the paths
    int32_t thread_index = 0;
    int32_t j = 0;
    for (int32_t i = 0; i < input->size; ++i) {
        auto& paths = paths_vec[thread_index];
        auto& parameters = parameters_vec[thread_index];

        // Check for out of bounds error 
        if (j >= paths.size() || j >= parameters.size()) {
            cerr << "Error: The index is out of bound!\n";
            cerr << "index: " << j << "\tin vector of the thread: " << thread_index << endl;
            cerr << "at value: " << i << endl;
            exit(1);
        }

        // Push alignment to alignments vector
        alignments.push_back(projectA_gwfa_path_to_alignment(parameters[j].projectA_hash_graph, &(paths[j].path)));

        // Update thread index
        if (thread_index == numThreads - 1) {
            thread_index = 0;
            ++j;
        } else {
            ++thread_index;
        }
    }


    for (auto& parameters : parameters_vec) {
        for (auto& parameter : parameters) {
            gwf_cleanup(parameter.km, parameter.graph);
            gwf_free(parameter.graph);
            km_destroy(parameter.km);
        }
    }

    delete input;
    ptr = nullptr;
}


// Function to get the algorithm struct for gwfa
projectA_algorithm_t* projectA_get_gwfa() {
    // Create new object
    projectA_algorithm_t* gwfa = new projectA_algorithm_t;

    // Assign function pointers
    gwfa->init = projectA_gwfa_init;
    gwfa->calculate_batch = projectA_gwfa_calculate_batch;
    gwfa->post = projectA_gwfa_post;

    return gwfa;
}


// Function to delete projectA_algorithm_t struct
void projectA_gwfa_destroy(projectA_algorithm_t* gwfa) {

    // Check if the struct still exists
    if (gwfa == nullptr) {
        cerr << "Warning: gwfa struct is null pointer!\n";
        return;
    }

    // Delete struct
    delete gwfa;
}
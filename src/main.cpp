/*

    projectA:
    main.cpp
    This file holds the main function for projectA.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <iostream>
#include <string>

#include "file_io.hpp"
#include "graph.hpp"
#include "algorithms/gssw.hpp"
#include "alignment.hpp"


int main(int argc, char* argv[]) {

    // Check format of the input
    if (argc != 3) {
        cerr << "Error: input must have the following format:\n";
        cerr << "projectA reference_graph.gfa node_list.txt\n";
        exit(1);
    }

    std::string ref_graph_file = argv[1];
    std::string node_list_file = argv[2];


    // Check validity of input file names
    if (ref_graph_file.size() < 5) {
        cerr << "Error: " << ref_graph_file << " is not a valid .gfa file!\n";
        exit(1);
    }
    if (ref_graph_file.substr(ref_graph_file.length() - 4) != ".gfa") {
        cerr << "Error: " << ref_graph_file << " is not a valid .gfa file!\n";
        exit(1);
    }


    // Read files

    // Read reference graph file
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa(ref_graph_file);
    // Index reference graph
    projectA_index_hash_graph(ref_graph);

    // Read node list file
    vector<projectA_node_list_t> clusters;
    projectA_read_node_list(clusters, node_list_file);


    // Build graph form cluster information
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);


    // Run gssw
    vector<projectA_alignment_t*> alignments_gssw;
    projectA_algorithm_t* gssw = projectA_get_gssw();
    void* ptr = gssw->init(graphs);
    gssw->calculate_batch(ptr);
    gssw->post(ptr, alignments_gssw);
    projectA_gssw_destroy(gssw);
    for (auto& alignment : alignments_gssw) {
        projectA_print_alignment(stderr, alignment);
    }


    // Cleanup
    for (auto& alignment : alignments_gssw) {
        delete alignment;
    }
    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
    cerr << "run succesfull!\n";
    return 0;
}
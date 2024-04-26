/*

    projectA:
    file_io.hpp
    This file holds the definitions for the projectA I/O functions used to read from and write to files.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/


#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "graph.hpp"

using namespace std;

// Basic IO:
// PRE:     fileName, testCases{TODO: struct/class that can hold the test cases}
//      fileName:       String containing the path to test case graph file.
//      testCases:      Instance of {testCase struct/class}.
// POST:    testCases{TODO: struct/class that can hold the test cases}
//      testCases:      Instance of {testCase struct/class} that holds the loaded test cases. 
void projectA_read_graph_file(string& fileName);


// PRE:     fileName, positions
//      fileName:       String containing the path to test case cluster file.
//      positions:      Vector holding projectA_positions_chain_t.
// POST:    positions
//      positions:      Vector of projectA_position_chain_t that holds the loaded test cases.
void projectA_read_cluster_file(string& fileName, vector<projectA_position_t>& positions);


// IO used for file comparison
// PRE:     fileName, graphs_map
//      fileName:       String containing the path to test case graph file.
//      graphs_map:     Empty unordered map with a string key and a vector holding a vector holding a int as value.
// POST:    graphs_map
//      graphs_map:     Unordered map holding the different reads as keys and the corresponing 
//                      node ids stored in a vector.
void projectA_read_graph_file(string& fileName, unordered_map<string, vector<vector<string>>>& graphs_map);


// PRE:     fileName, cluster_map
//      fileName:       String containing the path to test case cluster file.
//      cluster_map:    Empty unordered map with a strig key and int vector as value.
// POST:    cluster_map
//      cluster_map:    Unordered map holding the different reads as keys and the correspongind
//                      node ids stored in a vector.
void projectA_read_cluster_file(string& fileName, unordered_map<string, vector<string>>& cluster_map);


// PRE:     fileName
//      fileName:       String containing the path to a .gfa file that contains a reference graph
// POST:    return
//      return:         Pointer to a projectA_graph_t that holds the graph from the .gfa file.
projectA_graph_t* projectA_read_gfa(const string& fileName);


// PRE:     fileName
//      fileName:       String containing the path to a .gfa file that contains a reference graph
// POST:    return
//      return:         Pointer to a projectA_hash_graph_t that holds the graph from the .gfa file.
projectA_hash_graph_t* projectA_hash_read_gfa(const string& fileName);


// PRE:     node_lists, fileName
//      node_lists:     Reference to a vector holding projectA_node_list_t.
//      fileName:       String containing the path to a file in the node_list format.
// POST:    node_lists
//      node_lists:     Reference to a vector holding projectA_node_list_t where the node_lists
//                      and corresponding reads stored in the file at fileName have been added.
void projectA_read_node_list(vector<projectA_node_list_t>& node_lists, const string& fileName);


// PRE:     file, graph
//      file:           Valid pointer to an output file.
//      graph:          Pointer to a projectA_graph_t.
// POST:    file
//      file:           The graph in graph is printed to file
void projectA_print_graph(FILE* file, projectA_graph_t* graph);


// PRE:     file, graph
//      file:           Valid pointer to an output file.
//      graph:          Pointer to a projectA_hash_graph_t.
// POST:    file
//      file:           The graph in graph is printed to file.
void projectA_print_graph(FILE* file, projectA_hash_graph_t* graph);
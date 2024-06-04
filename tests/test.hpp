/*

    projectA:
    test.hpp
    This file holds the definitions of the functions used in the projectA test suit.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <unordered_map>

#include "file_io.hpp"
#include "graph.hpp"

using namespace std;

// PRE:     id, vec
//      id:             Unisgned 32 bit integer
//      vec:            Vector holding a vector of unsigned 32 bit integers
// POST:    return
//      return:         1 if id is contained in vec 0 if it is not.
int find_id_match (uint32_t id, vector<vector<uint32_t>> vec);


// PRE:     graph_file, cluster_file
//      graph_file:     String holding the path to a graph file.
//      cluster_file:   String holding the path to a cluster file.
// POST:    return 
//      return:         1 if all nodes in the cluster file can also be found
//                      in the graph file. 0 if not all nodes could be found.
int file_node_id_check (string graph_file, string cluster_file);
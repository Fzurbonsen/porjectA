/*

    projectA:
    extract_graph.hpp
    This file holds the definitions for the extract graph algorithm used to construct a graph from a position
    and is part of projectA.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <string>
#include <vector>
#include <utility>
#include <set>

#include "graph.hpp"


// PRE:     reference_graph, positions
//      reference_graph:    Instance of projectA_graph_t holding the reference graph.
//      positions:          Vector holding filled projectA_position_t.
// POST:    return
//      return:             Returns a vector holding pair of set holding unsigned integers 
//                          of size 32 bit holding all the node ids corresponding to a position 
//                          and the corresponding reads.
vector<pair<string, set<projectA_node_t*>>> projectA_extract_graph(projectA_hash_graph_t* ref_graph, vector<projectA_position_t>& positions);


// PRE:     reference_graph, positions
//      reference_graph:    Instance of projectA_graph_t holding the reference graph.
//      positions:          Filled projectA_position_t struct.
// POST:    return
//      return:             Returns a set holding unsigned integers of size 32 bit holding
//                          all the node ids corresponding to the position.
set<projectA_node_t*> projectA_extract_nodes(projectA_hash_graph_t* ref_graph, projectA_position_t& pos);
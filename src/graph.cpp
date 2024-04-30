/*

    projectA:
    graph.cpp
    This file holds the implemntations of the projectA graph helper functions.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/


#include <cstring>
#include <iostream>
#include <utility>
#include <set>
#include <queue>
#include <stack>

#include "graph.hpp"
#include "algorithm.hpp"

using namespace std;


// Function to delete a node
void projectA_delete_node(projectA_node_t* node) {
    if (node != nullptr) {
        delete node;
    }
}


// Function to delete a hash graph
void projectA_delete_hash_graph(projectA_hash_graph_t* graph) {
    if (graph != nullptr) {
        for (auto& node : graph->nodes) {
            projectA_delete_node(node.second);
        }
        delete graph;
    }
}


// Function to append a node to a projectA_hash_graph_t
void projectA_hash_graph_append_node(projectA_hash_graph_t* graph, string id, uint32_t len, string& seq, uint32_t index) {

    // Check if input seq matches the supposed length
    if (len != seq.size()) {
        cerr << "Error: sequence does not match sequence length!" << endl;
        exit(1);
    }
    // Check if a node with this id already exists
    if (graph->nodes.count(id)) {
        cerr << "Error: node with id: " << id << " already exists!\n";
        exit(1);
    }

    projectA_node_t* node = new projectA_node_t;

    // Fill node
    node->id = id;
    node->len = len;
    node->seq = seq;
    node->index = index;

    // Add node information to node map
    graph->nodes[id] = node;

    // Increase node counter
    graph->n_nodes++;
}


// Function to append an edge to a projectA_hash_graph_t
void projectA_hash_graph_append_edge(projectA_hash_graph_t* graph, string start, string end) {

    // Add edges to the nodes
    graph->nodes[start]->next.push_back(graph->nodes[end]);
    graph->nodes[end]->prev.push_back(graph->nodes[start]);

    // Increase edge counter
    graph->n_edges++;
}


// Function to store traversed nodes in the forward direction 
void projectA_graph_traverse_forward_distance(projectA_node_t* node , set<projectA_node_t*>& nodes, int64_t len) {

    // Check if we still need to traverse further
    if (len > 0) {

        len = len - node->len;
        nodes.insert(node);

        for (int i = 0; i < node->next.size(); ++i) {
            projectA_graph_traverse_forward_distance(node->next[i], nodes, len);
        }
    }
}


// Function to store traversed nodes in the backward direction 
void projectA_graph_traverse_backward_distance(projectA_node_t* node , set<projectA_node_t*>& nodes, int64_t len) {

    // Check if we still need to traverse further
    if (len > 0) {

        nodes.insert(node);
        len = len - node->len;

        for (int i = 0; i < node->prev.size(); ++i) {
            projectA_graph_traverse_forward_distance(node->prev[i], nodes, len);
        }
    }
}


// Function to build cluster information
void projectA_build_graph_from_cluster(projectA_hash_graph_t* graph, projectA_hash_graph_t* ref_graph, 
                                        projectA_node_list_t& node_list) {
    
    // cerr << node_list.n_nodes << endl;

    // Append nodes to graph
    for (int i = 0; i < node_list.n_nodes; ++i) {
        string& id = ref_graph->nodes[node_list.nodes[i]]->id;
        uint32_t& len = ref_graph->nodes[node_list.nodes[i]]->len;
        string& seq = ref_graph->nodes[node_list.nodes[i]]->seq;
        uint32_t& index = ref_graph->nodes[node_list.nodes[i]]->index;
        projectA_hash_graph_append_node(graph, id, len, seq, index);
    }

    // Iterate over all nodes to extract the contained edges.
    for (int i = 0; i < node_list.n_nodes; ++i) {
        const auto& node = ref_graph->nodes[node_list.nodes[i]];

        // Iterate over all edges linked to the node in the reference graph
        for (const auto& next : node->next) {
            
            // If the next node is included in the graph we add the edge.
            if (graph->nodes.find(next->id) != graph->nodes.end()) {
                projectA_hash_graph_append_edge(graph, node->id, next->id);
            }
        }
    }

    projectA_hash_graph_in_order_nodes(graph);
    // cerr << graph->n_edges << "\t" << graph->n_nodes << endl;
}


// Function to build a vector of cluster information
void projectA_build_graph_from_cluster(vector<projectA_algorithm_input_t>& graphs, projectA_hash_graph_t* ref_graph, 
                                        vector<projectA_node_list_t>& node_lists) {
    
    // Iterate over all node lists
    for (auto& node_list : node_lists) {
        projectA_algorithm_input_t algorithm_inputs;
        projectA_hash_graph_t* graph = new projectA_hash_graph_t;
        graph->n_edges = 0;
        graph->n_nodes = 0;

        // Build new graph for each node list
        projectA_build_graph_from_cluster(graph, ref_graph, node_list);

        // Fill algorithm input struct
        algorithm_inputs.read = node_list.read;
        algorithm_inputs.graph = graph;
        
        // Append the graph to the graph vector
        graphs.push_back(algorithm_inputs);
    }
}


// Function to index a hash graph
void projectA_index_hash_graph(projectA_hash_graph_t* graph) {

    uint32_t i = 0;
    for (auto& curr_node : graph->nodes) {
        curr_node.second->index = i;
        i++;
    }
}


// Function implementing DFS for topo sort
void projectA_hash_graph_topo_sort_DFS(projectA_node_t* node, stack<projectA_node_t*>& stack) {

    // Mark the node as visited
    node->visited = true;

    // Call the DFS on all following nodes
    for (auto& next : node->next) {

        // Check if the next nodes have already been visited
        if (!next->visited) {
            projectA_hash_graph_topo_sort_DFS(next, stack);
        } 
    }

    // Push the node to the stack
    stack.push(node);
}


// Function to create the in order vector of all the nodes in a hash graph
void projectA_hash_graph_in_order_nodes(projectA_hash_graph_t* graph) {
    
    // Reste the traversed bool in all nodes.
    for (auto& itr : graph->nodes) {
        auto& node = itr.second;
        node->visited = false;
    }

    // Stack to hold sorted nodes
    stack<projectA_node_t*> stack;

    // Perform topological sort:
    // Itterate over all nodes
    for (auto& itr : graph->nodes) {
        auto& node = itr.second;
        // Check if the node has already been visited
        if (!node->visited) {
            projectA_hash_graph_topo_sort_DFS(node, stack);
        }
    }


    // Check if all nodes have been visited
    for (auto& itr : graph->nodes) {
        auto& node = itr.second;
        if (!node->visited) {
            cerr << "Error: not all nodes could be reached from the top node!\n";
            exit(1);
        }
        // Reset node->visited
        node->visited = false;
    }

    // Flip order for in order vector
    while(!stack.empty()) {
        graph->nodes_in_order.push_back(stack.top());
        stack.pop();
    }
}
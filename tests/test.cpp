/*

    projectA:
    test.cpp
    This file holds the test suit for projectA.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/

#include <iostream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <unordered_map>
#include <utility>
#include <set>

#include "test.hpp"
#include "extract_graph.hpp"
#include "algorithm.hpp"
#include "algorithms/gt_gwfa.hpp"
#include "algorithms/gssw.hpp"

using namespace std;



// Function to find a match in the id vectors
int find_id_match (string id, vector<vector<string>> vec) {

    // Itterate over outer vector
    for (auto it1 = vec.begin(); it1 != vec.end(); ++it1) {

        // Itterate over inner vector
        for (auto it2 = (*it1).begin(); it2 != (*it1).end(); ++it2) {

            // If we find a match we return 1
            if (id == *it2) {
                return 1;
            }
        }
    }

    // If we don't find a match we return 0
    return 0;
}



// Function to compare input files to check whether they hold the same nodes.
int file_node_id_check (string graph_file, string cluster_file) {

    unordered_map<string, vector<vector<string>>> graph_map;
    unordered_map<string, vector<string>> cluster_map;

    int i = 0;

    // Read graph and cluster file
    projectA_read_graph_file(graph_file, graph_map);
    projectA_read_cluster_file(cluster_file, cluster_map);

    // Itterate over all reads in the cluster map
    for (auto cluster_it = cluster_map.begin(); cluster_it != cluster_map.end(); ++cluster_it) {
        // cluster_it->first;   key
        // cluster_it->second;  value

        // Check if the current read is in the graphs map
        if (graph_map.count(cluster_it->first) > 0) {

            // Itterate over all ids that correspond to the current read
            for (auto id_it = cluster_it->second.begin(); id_it != cluster_it->second.end(); ++id_it) {

                // Check if the node id can be found in the vectors coresponding to the read in the graphs map
                if (find_id_match(*id_it, graph_map.at(cluster_it->first)) == 0) {
                    cerr << "id not found:\t" << cluster_it->first << endl;
                } else {
                    i++;
                }
            }
        } else {
            cerr << "read not found:\t"<<cluster_it->first << endl;
        }
    }

    return i;
}



// Function to import tests form an external file
int import_tests(string fileName, vector<string>& read_vector) {
    // TODO
    return 0;
}



// Function to get gold results
void gold_results(vector<string>& gold_vector, vector<string>& read_vector)  {
    // TODO
}



// Function to get test results
void test_results(vector<string>& result_vector, vector<string>& read_vector) {
    // TODO
}



// Function to compare the test and gold results
void compare_results(int& n_tests, int& correct_tests, 
                        vector<string>& gold_vector, vector<string>& result_vector) {
    
    // Ensure that the gold vector has the size of the number of tests
    if (gold_vector.size() != n_tests) {
        cerr << "Error: Gold size does not match n_tests\n";
        return;
    }

    // Ensure that the result vector has the size of the number of tests
    if (result_vector.size() != n_tests) {
        cerr << "Error: Result size does not match n_tests\n";
        return;
    }
    
    // Iterate through the vectors and compare the elements
    for (size_t i = 0; i < n_tests; ++i) {
        if (gold_vector[i] == result_vector[i]) {
            // If the elements match we increment the counter
            correct_tests++;
        } else {
            // Print the mismatched test case
            cerr << "Mismatch at index " << i << ":\n";
            cerr << "Gold result: " << gold_vector[i] << "\n";
            cerr << "Result:      " << result_vector[i] << "\n";
        }
    }
}



// Function that handles all the test cases
int run_tests(string fileName) {
    cerr << "running tests from file: " << fileName << "\n";
    
    // Variables used to track performance
    double runtime;

    // record the starting time
    clock_t start = clock();


    
    // record the ending time
    clock_t end = clock();
    // calculate the runtime in seconds
    runtime = double(end - start) / CLOCKS_PER_SEC;

    // output the results of the tests
    // cerr << "time: " << runtime << " seconds\n" << correct_tests << " out of " << n_tests << " are correct\n";

    return 0;
}

int main() {


    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    projectA_read_node_list(clusters, "./test_cases/tests.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);

    // for (auto& graph : graphs) {
    //     projectA_print_graph(stderr, graph.graph);
    // }



    // FILE* file = fopen("test.txt", "w");
    // for (auto& graph : graphs) {
    //     projectA_print_graph(file, graph.graph);
    // }
    // fclose(file);

    void* ptr;
    vector<projectA_alignment_t*> alignments_gt_gwfa;

    // Tests for gt_gwfa:
    cerr << "testing gt_gwfa!\n";
    projectA_algorithm_t* gt_gwfa = projectA_get_gt_gwfa();
    cerr << "loading inputs\n";
    ptr = gt_gwfa->init(graphs);
    cerr << "calculate batch\n";
    ptr = gt_gwfa->calculate_batch(ptr);
    cerr << "entering post\n";
    gt_gwfa->post(ptr, alignments_gt_gwfa);
    cerr << "destroying algorithm struct\n";
    projectA_gt_gwfa_destroy(gt_gwfa);
    







    vector<projectA_alignment_t*> alignments_gssw;

    // Tests for gssw:
    cerr << "testing gssw!" << endl;
    projectA_algorithm_t* gssw = projectA_get_gssw();
    cerr << "loading inputs\n";
    ptr = gssw->init(graphs);
    cerr << "calculating batch\n";
    gssw->calculate_batch(ptr);
    cerr << "entering post\n";
    gssw->post(ptr, alignments_gssw);
    cerr << "destroying algorithm struct\n";
    projectA_gssw_destroy(gssw);



    for (auto& alignment : alignments_gt_gwfa) {
        projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
    for (auto& alignment : alignments_gssw) {
        projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
    cerr << "run succesfull!\n";
    return 0;
}
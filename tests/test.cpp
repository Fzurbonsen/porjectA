/*

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
    int n_tests, correct_tests = 0;

    // record the starting time
    clock_t start = clock();


    // Vector filled with graphs
    // vector<"graph_class"> graph_vector;
    // Vector filled with the reads
    vector<string> read_vector;
    // // Vector for gold and normal results
    vector<string> result_vector, gold_vector;

    // Load file with test graphs and alignments
    n_tests = import_tests(fileName, read_vector);

    // Get gold results
    gold_results(gold_vector, read_vector);

    // Get test results
    test_results(result_vector, read_vector);

    // Compare the tests and gold standard.
    compare_results(n_tests, correct_tests, gold_vector, result_vector);
    
    // record the ending time
    clock_t end = clock();
    // calculate the runtime in seconds
    runtime = double(end - start) / CLOCKS_PER_SEC;

    // output the results of the tests
    cerr << "time: " << runtime << " seconds\n" << correct_tests << " out of " << n_tests << " are correct\n";

    return 0;
}

int main() {

    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    projectA_index_hash_graph(ref_graph);

    projectA_read_node_list(clusters, "./test_cases/node_list.txt");


    vector<pair<const string, projectA_hash_graph_t*>> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);

    // FILE* file = fopen("test.txt", "w");
    // for (auto& graph : graphs) {
    //     projectA_print_graph(file, graph.second);
    // }
    // fclose(file);



    // Tests for gt_gwfa:
    // projectA_algorithm gt_gwfa;
    // projectA_get_gt_gwfa(gt_gwfa);
    // void* pointer = nullptr;
    // pointer = gt_gwfa.init(graphs);

    // Cast the void pointer back to its original type
    // auto ptr = static_cast<projectA_gt_gwfa_io_t*>(pointer);

    // pointer = gt_gwfa.calculate_batch(pointer);
    // gt_gwfa.post(pointer);



    // Tests for gssw:
    projectA_algorithm_t* gssw = projectA_get_gssw();
    void* ptr = gssw->init(graphs);
    gssw->calculate_batch(ptr);
    gssw->post(ptr);
    projectA_gssw_destroy(gssw);



    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.second);
    }
    projectA_delete_hash_graph(ref_graph);
    cerr << "run succesfull!\n";
    return 0;
}
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
#include <unordered_map>
#include <utility>
#include <set>
#include <chrono>
#include <thread>

#include "test.hpp"
#include "extract_graph.hpp"
#include "algorithm.hpp"
#include "algorithms/gt_gwfa.hpp"
#include "algorithms/gssw.hpp"
#include "alignment.hpp"

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


// Function to run the gssw algorithm
void projectA_run_gssw(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {
    void* ptr2;
    vector<projectA_alignment_t*> alignments_gssw;
    vector<thread> threads;
    
    projectA_algorithm_t* gssw = projectA_get_gssw();
    ptr2 = gssw->init(graphs, numThreads);
    for (int i = 0; i < numThreads; ++i) {
        // gssw->calculate_batch(ptr2, i);
        threads.push_back(thread(gssw->calculate_batch, ptr2, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    threads.clear();
    gssw->post(ptr2, alignments_gssw, numThreads);
    projectA_gssw_destroy(gssw);

    for (auto& alignment : alignments_gssw) {
        // projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
}

// Function to run gt_gwfa algoritm
void projectA_run_gt_gwfa(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {
    void* ptr1;
    vector<projectA_alignment_t*> alignments_gt_gwfa;
    vector<thread> threads;

    projectA_algorithm_t* gt_gwfa = projectA_get_gt_gwfa();
    ptr1 = gt_gwfa->init(graphs, numThreads);
    // gt_gwfa->calculate_batch(ptr1, 0);
    for (int i = 0; i < numThreads; ++i) {
        // gssw->calculate_batch(ptr2, i);
        threads.push_back(thread(gt_gwfa->calculate_batch, ptr1, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    gt_gwfa->post(ptr1, alignments_gt_gwfa, numThreads);
    projectA_gt_gwfa_destroy(gt_gwfa);

    for (auto& alignment : alignments_gt_gwfa) {
        // projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
}


// Function to get a timed run of gssw
void timed_run_gssw(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {
    typedef std::chrono::high_resolution_clock Clock;
    void* ptr2;
    vector<projectA_alignment_t*> alignments_gssw;
    vector<thread> threads;
    
    projectA_algorithm_t* gssw = projectA_get_gssw();
    ptr2 = gssw->init(graphs, numThreads);
    auto t1 = Clock::now();
    for (int i = 0; i < numThreads; ++i) {
        // gssw->calculate_batch(ptr2, i);
        threads.push_back(thread(gssw->calculate_batch, ptr2, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    threads.clear();
    auto t2 = Clock::now();
    gssw->post(ptr2, alignments_gssw, numThreads);
    projectA_gssw_destroy(gssw);

    for (auto& alignment : alignments_gssw) {
        // projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
    // auto t1 = Clock::now();
    // projectA_run_gssw(graphs, numThreads);
    // auto t2 = Clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "gssw:\t" << "threads: " << numThreads << "\t" << duration1.count() << " ms\t" << endl;
}


// Function to get a timed run of gt_gwfa
void timed_run_gt_gwfa(vector<projectA_algorithm_input_t>& graphs, int32_t numThreads) {
    typedef std::chrono::high_resolution_clock Clock;
    void* ptr1;
    vector<projectA_alignment_t*> alignments_gt_gwfa;
    vector<thread> threads;

    projectA_algorithm_t* gt_gwfa = projectA_get_gt_gwfa();
    ptr1 = gt_gwfa->init(graphs, numThreads);
    auto t3 = Clock::now();
    for (int i = 0; i < numThreads; ++i) {
        // gssw->calculate_batch(ptr2, i);
        threads.push_back(thread(gt_gwfa->calculate_batch, ptr1, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    auto t4 = Clock::now();
    gt_gwfa->post(ptr1, alignments_gt_gwfa, numThreads);
    projectA_gt_gwfa_destroy(gt_gwfa);

    for (auto& alignment : alignments_gt_gwfa) {
        // projectA_print_alignment(stderr, alignment);
        delete alignment;
    }
  
    // auto t3 = Clock::now();
    // projectA_run_gt_gwfa(graphs, 1);
    // auto t4 = Clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
    std::cerr << "gt_gwfa:\t" << "threads: " << numThreads << "\t" << duration2.count() << " ms\t" << endl;
}



int main() {

    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    projectA_read_node_list(clusters, "./test_cases/node_list.txt");
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

    timed_run_gssw(graphs, 1);
    timed_run_gssw(graphs, 2);
    timed_run_gssw(graphs, 4);
    timed_run_gssw(graphs, 8);
    timed_run_gssw(graphs, 16);
    timed_run_gssw(graphs, 32);
    timed_run_gssw(graphs, 64);
    timed_run_gssw(graphs, 128);
    timed_run_gssw(graphs, 256);
    // timed_run_gt_gwfa(graphs, 1);
    // timed_run_gt_gwfa(graphs, 2);
    // timed_run_gt_gwfa(graphs, 4);
    // timed_run_gt_gwfa(graphs, 8);
    // timed_run_gt_gwfa(graphs, 16);
    // timed_run_gt_gwfa(graphs, 32);
    // timed_run_gt_gwfa(graphs, 64);
    // timed_run_gt_gwfa(graphs, 128);
    // timed_run_gt_gwfa(graphs, 256);
    // timed_run_gt_gwfa(graphs, 512);
    // timed_run_gt_gwfa(graphs, 1024);
    // timed_run_gt_gwfa(graphs, 2048);


    // if (alignments_gssw.size() != graphs.size()) {
    //     cerr << "Error: gssw size does not match grpah size.\t" << graphs.size() << " " << alignments_gssw.size() << endl;
    //     exit(1);
    // }

    // int match = 0;
    // for (int i = 0; i < graphs.size(); ++i) {

    //     match += projectA_compare_alignments(true, stderr, alignments_gssw[i], 
    //                                                         alignments_gt_gwfa[i]);
    // }

    // cerr << match << " matches of " << graphs.size() << " alignments.\n";

    // for (auto& alignment : alignments_gt_gwfa) {
    //     // projectA_print_alignment(stderr, alignment);
    //     delete alignment;
    // }
    // for (auto& alignment : alignments_gssw) {
    //     // projectA_print_alignment(stderr, alignment);
    //     delete alignment;
    // }
    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
    cerr << "run succesfull!\n";
    return 0;
}
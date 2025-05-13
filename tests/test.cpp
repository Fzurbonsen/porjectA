/*

    projectA:
    test.cpp
    This file holds the test suite for projectA.
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
#include <stack>
#include <queue>

#include "test.hpp"
#include "extract_graph.hpp"
#include "algorithm.hpp"
// #include "algorithms/gt_gwfa.hpp"
#include "algorithms/gssw.hpp"
#include "algorithms/gwfa.hpp"
#include "alignment.hpp"
// #include "algorithms/ksw2.hpp"
#include "algorithms/csswl.hpp"
// #include "algorithms/abPOA.hpp"
// #include "algorithms/vargas.hpp"
#include "algorithms/gnwa.hpp"

// #include "gt_gwfa/edlib.h"
#include "algorithms/edlib.hpp"

#include "algorithms/s_gwfa.hpp"
#include "s_gwfa/s_gwfa.h"
#include "edlib/edlib.h"


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


// Function to get the aligned reference from gssw
void projectA_get_gssw_reference(projectA_alignment_t* alignment, string& reference) {
    // Define local variables
    auto& cigar = alignment->cigar_string;
    int32_t pos = alignment->offset;

    // Iterate over the CIGAR
    for (int32_t i = 0; i < cigar.len; i++) {
        auto& element = cigar.elements[i];

        // Swith case to differentiate between the different operations
        switch(element.type) {

            // Match
            case('M') :
                pos += element.len;
                break;

            // Deletion
            case('D') :
                pos += element.len;
                break;

            // Default
            default :
                break;
        }
    }

    // Cut the reference
    reference = alignment->reference.substr(alignment->offset, pos - alignment->offset);
    // alignments2[i]->reference = alignments2[i]->reference.substr(0, element_length);
    // cerr << reference << endl;
    // cerr << alignment->offset << "\t" << pos << endl;
    // cerr << pos - alignment->offset << "\t" << strlen(reference.c_str()) << endl;
}


// Function to get the aligned read from gssw
void projectA_get_gssw_read(projectA_alignment_t* alignment, string& read) {
    // Define local variables
    auto& cigar = alignment->cigar_string;
    int32_t pos = 0;

    // Iterate over the CIGAR
    for (int32_t i = 0; i < cigar.len; i++) {
        auto& element = cigar.elements[i];

        // Swith case to differentiate between the different operations
        switch(element.type) {

            // Match
            case('M') :
                pos += element.len;
                break;

            // Insertion
            case('I') :
                pos += element.len;
                break;

            // Default
            default :
                break;
        }
    }

    // Cut the read
    read = alignment->read.substr(0,pos);
    // cerr << pos << endl;
}


// Function to get the aligned read from csswl
void projectA_get_csswl_read(projectA_alignment_t* alignment, string& read) {
    // Define local variables
    auto& cigar = alignment->cigar_string;
    int32_t pos = alignment->read_start;

    // Iterate over the CIGAR
    for (int32_t i = 0; i < cigar.len; i++) {
        auto& element = cigar.elements[i];

        // Swith case to differentiate between the different operations
        switch(element.type) {

            // Match
            case('M') :
                pos += element.len;
                break;

            // Insertion
            case('I') :
                pos += element.len;
                break;

            // Default
            default :
                break;
        }
    }

    // Cut the read
    read = alignment->read.substr(alignment->read_start, pos - alignment->read_start);
    // cerr << alignment->read_start << "\t" << pos << endl;
}


// Function to run the gssw algorithm
void projectA_get_alignment_gssw(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr2;
    vector<thread> threads;

    projectA_algorithm_t* gssw = projectA_get_gssw();
    ptr2 = gssw->init(alignments, numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(gssw->calculate_batch, ptr2, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    threads.clear();
    gssw->post(ptr2, alignments, numThreads);
    projectA_gssw_destroy(gssw);
}

// Function to run and time the gssw algorithm
int projectA_get_timed_alignment_gssw(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr2;
    vector<thread> threads;
    typedef std::chrono::high_resolution_clock Clock;

    projectA_algorithm_t* gssw = projectA_get_gssw();
    ptr2 = gssw->init(alignments, numThreads);

    auto t0 = Clock::now();
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(gssw->calculate_batch, ptr2, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    threads.clear();
    auto t1 = Clock::now();

    gssw->post(ptr2, alignments, numThreads);
    projectA_gssw_destroy(gssw);
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    return duration;
}


// Function to run gwfa algorithm
void projectA_get_alignment_gwfa(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr;
    vector<thread> threads;
    
    // Get algorithm struct
    projectA_algorithm_t* gwfa = projectA_get_gwfa();
    
    // Initialize data
    ptr = gwfa->init(alignments, numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(gwfa->calculate_batch, ptr, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }

    // gwfa post
    gwfa->post(ptr, alignments, numThreads);

    // Destroy algorithm struct
    projectA_gwfa_destroy(gwfa);
}


// Function to run and time the gwfa algorithm
int projectA_get_timed_alignment_gwfa(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr;
    vector<thread> threads;
    typedef std::chrono::high_resolution_clock Clock;
    
    // Get algorithm struct
    projectA_algorithm_t* gwfa = projectA_get_gwfa();
    
    // Initialize data
    ptr = gwfa->init(alignments, numThreads);

    auto t0 = Clock::now();
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(gwfa->calculate_batch, ptr, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    auto t1 = Clock::now();

    // gwfa post
    gwfa->post(ptr, alignments, numThreads);

    // Destroy algorithm struct
    projectA_gwfa_destroy(gwfa);
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    return duration;
}

// Function to run gwfa algorithm
void projectA_get_alignment_s_gwfa(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr;
    vector<thread> threads;
    
    // Get algorithm struct
    projectA_algorithm_t* s_gwfa = projectA_get_s_gwfa();
    
    // Initialize data
    ptr = s_gwfa->init(alignments, numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(s_gwfa->calculate_batch, ptr, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }

    // gwfa post
    s_gwfa->post(ptr, alignments, numThreads);

    // Destroy algorithm struct
    projectA_gwfa_destroy(s_gwfa);
}


// Function to run gwfa algorithm
int projectA_get_timed_alignment_s_gwfa(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr;
    vector<thread> threads;
    typedef std::chrono::high_resolution_clock Clock;
    
    // Get algorithm struct
    projectA_algorithm_t* s_gwfa = projectA_get_s_gwfa();
    
    // Initialize data
    ptr = s_gwfa->init(alignments, numThreads);
    auto t0 = Clock::now();
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(s_gwfa->calculate_batch, ptr, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    auto t1 = Clock::now();

    // gwfa post
    s_gwfa->post(ptr, alignments, numThreads);

    // Destroy algorithm struct
    projectA_gwfa_destroy(s_gwfa);
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    return duration;
}


// Function to run the GNWA algorithm
void projectA_get_alignment_gnwa(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {
    void* ptr;
    vector<thread> threads;

    projectA_algorithm_t* gnwa = projectA_get_gnwa();
    ptr = gnwa->init(alignments, numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(thread(gnwa->calculate_batch, ptr, i));
    }
    for (auto& th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }
    threads.clear();
    gnwa->post(ptr, alignments, numThreads);

    projectA_gnwa_destroy(gnwa);
}


// Function to write the alignment information to a test file
void projectA_write_to_test_file(FILE* file, vector<projectA_alignment_t*>& alignments1, vector<projectA_alignment_t*>& alignments2, int32_t i) {

    fprintf(file, "gssw:\t");
    fprintf(file, "%i\t", alignments1[i]->offset);
    projectA_print_cigar(file, &alignments1[i]->cigar_string);
    // fprintf(file, "score:\t%d\t\tmatch/lenght:\t%f\n", alignments1[i]->score, match_per_length);
    // fprintf(file, "matches:\t%d\t\tmismatches:\t%d\t\tlength:\t%d\n", alignments1[i]->n_matches, alignments1[i]->n_mismatches, 
    //                                                                     alignments1[i]->cigar_string.operations_length);
    fprintf(file, "read:\t%s\n", alignments1[i]->read.c_str());
    fprintf(file, "ref:\t%s\n", alignments1[i]->reference.c_str());
    projectA_print_path(file, alignments1[i]->nodes);
    

    fprintf(file, "gwfa:\t");
    fprintf(file, "%i\t", alignments2[i]->offset);
    projectA_print_cigar(file, &alignments2[i]->cigar_string);
    // fprintf(file, "score:\t%d\n", alignments2[i]->score);
    fprintf(file, "read:\t%s\n", alignments2[i]->read.c_str());
    fprintf(file, "ref:\t%s\n", alignments2[i]->reference.c_str());
    projectA_print_path(file, alignments2[i]->nodes);
}


// Function to create read sets
void projectA_create_read_sets(unordered_map<string, vector<projectA_alignment_t*>>& map, projectA_alignment_t* alignment) {
    
    // Check if the pointers are valid
    if (alignment == nullptr) {
        cerr << "Error: invalid alignment pointer!\n";
        exit(1);
    }

    // Check if the read is already in the map
    if (map.count(alignment->read) > 0) {

        // If the read is int the map then we push the alignment to the vector in the map
        map[alignment->read].push_back(alignment);

    } else {

        // If the read is not in the map yet, we add it together with a vector of alignment pointers
        vector<projectA_alignment_t*> alignment_vec;
        alignment_vec.push_back(alignment);
        map[alignment->read] = alignment_vec;
    }
}


// Function to determine max score read
void projectA_determine_max_score_alignment(unordered_map<string, projectA_alignment_t*>& max_score_alignments_map,
                                            unordered_map<string, vector<projectA_alignment_t*>>& aligned_reads_map) {
    
    //Iterate over all the reads in the map
    for (auto& read_pair : aligned_reads_map) {

        projectA_alignment_t* max_score_alignment = read_pair.second[0];
        // Iterate over all the alignments for a specific read
        for (auto& alignment : read_pair.second) {
            // Check if the alignment has a higher score than the previous max score
            if (alignment->score > max_score_alignment->score) {
                max_score_alignment = alignment;
            }
        }
        // Add max score alignment to the map
        max_score_alignments_map[read_pair.first] = max_score_alignment;
    }
}


// Function to sort a vector of alignments
vector<projectA_alignment_t*> projectA_sort_alignment_vector(vector<projectA_alignment_t*>& alignments) {
    if (alignments.empty()) {
        cerr << "[projectA] Warning: Vector is empty!\n";
        return {};
    }

    // Use sort with a custom comparator
    sort(alignments.begin(), alignments.end(), [](const projectA_alignment_t* a, const projectA_alignment_t* b) {
        if (a == nullptr || b == nullptr) return a != nullptr; // Nulls are treated as lowest priority
        return a->score > b->score; // Sort descending by score
    });

    return alignments; // Return the sorted vector
}


// Function to determine the top five max score read
void projectA_determine_n_max_score_alignment(unordered_map<string, vector<projectA_alignment_t*>>& n_max_score_alignments_map,
                                            unordered_map<string, vector<projectA_alignment_t*>>& aligned_reads_map,
                                            int32_t n) {
    
    // Iterate over all the reads in the map
    for (auto& read_pair : aligned_reads_map) {

        read_pair.second = projectA_sort_alignment_vector(read_pair.second);
        
        vector<projectA_alignment_t*> n_max_score_alignments;

        // Check if there are at least n entries
        if (read_pair.second.size() < n) {
            // Copy all the entries
            for (auto& entry : read_pair.second) {
                n_max_score_alignments.push_back(entry);
            }
        } else {
            // Copy the first n entries
            for (int32_t i = 0; i < n; ++i) {
                n_max_score_alignments.push_back(read_pair.second[i]);
            }
        }

        n_max_score_alignments_map[read_pair.first] = n_max_score_alignments;
    }
}


// Function to extract graph from sim_position
projectA_hash_graph_t* projectA_extract_graph_from_sim_pos(projectA_hash_graph_t* ref_graph, vector<string> nodes) {

    // Create new graph
    projectA_hash_graph_t* graph = new projectA_hash_graph_t;
    graph->n_edges = 0;
    graph->n_nodes = 0;
    for (auto& node : nodes) {
        string id = ref_graph->nodes[node]->id;
        uint32_t len = ref_graph->nodes[node]->len;
        string seq = ref_graph->nodes[node]->seq;
        uint32_t index = ref_graph->nodes[node]->index;
        projectA_hash_graph_append_node(graph, id, len, seq, index);
    }

    // Iterate over all nodes to extract the contained edges.
    for (int i = 0; i < nodes.size(); ++i) {
        const auto& node = ref_graph->nodes[nodes[i]];

        // Iterate over all edges linked to the node in the reference graph
        for (const auto& next : node->next) {
            
            // If the next node is included in the graph we add the edge.
            if (graph->nodes.find(next->id) != graph->nodes.end()) {
                projectA_hash_graph_append_edge(graph, node->id, next->id);
            }
        }
    }
    projectA_hash_graph_in_order_nodes(graph);
    return graph;
}


void run_standard_tests(string graphFile, string positionFile, string simPositionFile, int match, int mismatch, uint gap_open, uint gap_extend, FILE* outputFile) {
    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa(graphFile);
    // projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/linear_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    // projectA_read_node_list(clusters, "./test_cases/node_list_2.txt");
    // projectA_read_node_list(clusters, "./test_cases/node_list_2_small.txt");
    projectA_read_node_list(clusters, positionFile);
    // projectA_read_node_list(clusters, "./test_cases/tests.txt");
    // projectA_read_node_list(clusters, "./test_cases/linear_node_list.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);


    // Read simulated positions
    unordered_map<string, vector<string>> sim_positions;
    projectA_read_sim_positions(sim_positions, simPositionFile);
    // projectA_read_sim_positions_from_two_files(sim_positions, "./test_cases/sim_reads_new.txt", "./test_cases/sim_paths.txt");
 
    // Maps to hold the different read sets
    unordered_map<string, vector<projectA_alignment_t*>> read_set_map1;
    unordered_map<string, vector<projectA_alignment_t*>> read_set_map2;
    unordered_map<string, vector<projectA_alignment_t*>> read_set_map3;
    // Maps to hold the max scoring alignments for each read
    unordered_map<string, projectA_alignment_t*> max_scoring_alignments1;
    unordered_map<string, projectA_alignment_t*> max_scoring_alignments2;
    unordered_map<string, projectA_alignment_t*> max_scoring_alignments3;
    // Maps to hold the top 5 max scoring alignments for each read
    unordered_map<string, vector<projectA_alignment_t*>> five_max_scoring_alignments1;
    unordered_map<string, vector<projectA_alignment_t*>> five_max_scoring_alignments2;
    unordered_map<string, vector<projectA_alignment_t*>> five_max_scoring_alignments3;


    vector<projectA_alignment_t*> alignments1;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments1.push_back(alignment);
        projectA_create_read_sets(read_set_map1, alignment);
    }
    vector<projectA_alignment_t*> alignments2;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments2.push_back(alignment);
        projectA_create_read_sets(read_set_map2, alignment);
    }
    vector<projectA_alignment_t*> alignments3;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments3.push_back(alignment);
        projectA_create_read_sets(read_set_map3, alignment);
    }

    vector<projectA_algorithm_input_t> sim_graphs;
    for (auto& pos : sim_positions) {
        // projectA_algorithm_input_t pair;
        // pair.graph = projectA_extract_graph_from_sim_pos(ref_graph, pos.second);
        // pair.read = pos.first;
        // sim_graphs.push_back(pair);
    }

    int32_t base_size = alignments1.size();
    for (int32_t i = 0; i < sim_graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = base_size + i;
        alignment->graph = sim_graphs[i].graph;
        alignment->read = sim_graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments1.push_back(alignment);
        projectA_create_read_sets(read_set_map1, alignment);
    }

    base_size = alignments2.size();
    for (int32_t i = 0; i < sim_graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = base_size + i;
        alignment->graph = sim_graphs[i].graph;
        alignment->read = sim_graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments2.push_back(alignment);
        projectA_create_read_sets(read_set_map2, alignment);
    }

    base_size = alignments3.size();
    for (int32_t i = 0; i < sim_graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = base_size + i;
        alignment->graph = sim_graphs[i].graph;
        alignment->read = sim_graphs[i].read;
        alignment->match = match;
        alignment->mismatch = mismatch;
        alignment->gap_open = gap_open;
        alignment->gap_extend = gap_extend;
        alignments3.push_back(alignment);
        projectA_create_read_sets(read_set_map3, alignment);
    }

    typedef std::chrono::high_resolution_clock Clock;

    auto t0 = Clock::now();
    // // projectA_get_alignment_gssw(alignments1, 48);
    // projectA_get_alignment_gssw(alignments1, 32);
    // // projectA_get_alignment_gnwa(alignments1, 16);
    // // projectA_get_alignment_gwfa(alignments1, 32);
    auto t1 = Clock::now();
    // cerr << "gssw time: " << projectA_get_timed_alignment_gssw(alignments2, 1) << endl;
    // cerr << "gssw time: " << projectA_get_timed_alignment_gssw(alignments1, 32) << endl;
    // cerr << "s_gwfa time: " << projectA_get_timed_alignment_s_gwfa(alignments1, 1) << endl;
    // cerr << "gwfa time: " << projectA_get_timed_alignment_gwfa(alignments1, i) << endl;

    auto t2 = Clock::now();
    // projectA_get_alignment_gwfa(alignments3, 16);
    projectA_get_alignment_gwfa(alignments1, 32);
    projectA_get_alignment_gssw(alignments3, 32);
    // projectA_get_alignment_s_gwfa(alignments1, 48);
    // projectA_get_alignment_s_gwfa(alignments1, 32);

    // for (int i = 1; i <= 30; ++i) {
    //     cerr << i << "\t";
    //     auto time_1 = Clock::now();
    //     projectA_get_alignment_gwfa(alignments1, i);
    //     for (auto alignment : alignments1) {
    //         // projectA_edlib(alignment);
    //         projectA_csswl(alignment);
    //     }
    //     auto time_2 = Clock::now();
    //     auto duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(time_2 - time_1);
    //     cerr << duration_.count() << endl;
    // }
    auto t3 = Clock::now();


    // Truncate references for gwfa
    // for (int i = 0; i < alignments1.size(); ++i) {

        // Strings to hold the cut reference and read
        // string new_reference;
        // string new_read;

        // projectA_get_gssw_reference(alignments1[i], new_reference);
        // projectA_get_gssw_read(alignments1[i], new_read);

        // alignments2[i]->reference = new_reference;
        // alignments2[i]->read = new_read;
    // }

    for (int32_t i = 0; i < alignments2.size(); ++i) {
        auto& alignment1 = alignments1[i];
        auto& alignment2 = alignments2[i];
        auto& alignment3 = alignments3[i];

        // csswl
        // alignment2->reference = alignment1->reference;
        // alignment2->read = alignment1->read;
        // projectA_csswl(alignment1);
        // projectA_csswl(alignment3);

        // // recut read
        // string new_read;
        // projectA_get_csswl_read(alignment2, new_read);
        // alignment2->read = new_read;

        // string new_reference;
        // string new_read;

        // projectA_get_gssw_reference(alignment3, new_reference);
        // projectA_get_gssw_read(alignment3, new_read);

        // alignment2->reference = new_reference;
        // alignment2->read = new_read;

        // edLib
        projectA_edlib(alignment3);
        projectA_edlib(alignment1);

        // // ksw2
        // alignment2->reference = alignment1->reference;
        // alignment2->read = alignment1->read;
        // projectA_ksw2(alignment1);
        // alignment2->offset = 0;

    }
    auto t4 = Clock::now();


    // Find the max scoring alignment for each read
    projectA_determine_max_score_alignment(max_scoring_alignments1, read_set_map1);
    projectA_determine_max_score_alignment(max_scoring_alignments2, read_set_map2);
    projectA_determine_max_score_alignment(max_scoring_alignments3, read_set_map3);

    projectA_determine_n_max_score_alignment(five_max_scoring_alignments1, read_set_map1, 5);
    projectA_determine_n_max_score_alignment(five_max_scoring_alignments2, read_set_map2, 5);
    projectA_determine_n_max_score_alignment(five_max_scoring_alignments3, read_set_map3, 5);



    FILE* file = fopen("./files/CIGAR.txt", "w");
    double cigar_accuracy = 0;
    double node_accuracy = 0;
    int count = 0;
    int n_high_score = 0;
    for (int i = 0; i < alignments1.size(); ++i) {

        double match_per_length = (float)(alignments1[i]->n_matches)/(float)(alignments1[i]->cigar_string.operations_length);

        // if (alignments3[i]->score > 200) {
        // if (match_per_length > 0.6) {
        if (true) {
            n_high_score += 1;


            // Alignment 1/2
            // Calculate accuracy of the two CIGARs and print the CIGARs
            // double acc = 0;
            // acc = projectA_cigar_accuracy(&alignments1[i]->cigar_string, &alignments2[i]->cigar_string);
            // cigar_accuracy += acc;
            // projectA_write_to_test_file(file, alignments1, alignments2, i);
            // fprintf(file, "#accuracy:\t%f\n", acc);
            // fprintf(file, "\n\n");



            // Alignment 2/3
            // Calculate accuracy of the two CIGARs and print the CIGARs
            // double acc = 0;
            // acc = projectA_cigar_accuracy(&alignments3[i]->cigar_string, &alignments2[i]->cigar_string);
            // cigar_accuracy += acc;
            // projectA_write_to_test_file(file, alignments3, alignments2, i);
            // fprintf(file, "#accuracy:\t%f\n", acc);
            // fprintf(file, "\n\n");



            // Alignment 1/3
            // Calculate accuracy of the two CIGARs and print the CIGARs
            double acc = 0;
            acc = projectA_cigar_accuracy(&alignments1[i]->cigar_string, &alignments3[i]->cigar_string);
            cigar_accuracy += acc;
            projectA_write_to_test_file(file, alignments1, alignments2, i);
            fprintf(file, "#accuracy:\t%f\n", acc);
            // fprintf(file, "\n\n");
            // projectA_print_cigar(stderr, &alignments1[i]->cigar);


            // Calculate accuracy of the two paths and print paths
            double n_acc = 0;
            // node_accuracy += projectA_path_accuracy(alignments1[i]->nodes, alignments2[i]->nodes);
            // node_accuracy += projectA_path_accuracy(alignments1[i]->nodes, alignments3[i]->nodes);
            n_acc = project_weighted_path_accuracy(alignments1[i], alignments3[i]);
            node_accuracy += n_acc;
            fprintf(file, "#node accuracy:\t%f\n", n_acc);
            fprintf(file, "\n\n");

            // projectA_print_path(stderr, alignments1[i]->nodes);
            // projectA_print_path(stderr, alignments3[i]->nodes);
            // fprintf(stderr, "\n");
            // fprintf(stderr, "[%s]\n\n", alignments1[i]->graph->nodes_in_order[alignments1[i]->graph->nodes_in_order.size()-1]->id.c_str());
        }
    }
    fclose(file);



    // Compare the max scoring alignments of alignments1 and alignments3
    double pos_accuracy = 0;
    int32_t pos_count = 0;
    for (auto& read_pair : max_scoring_alignments1) {
        string read = read_pair.first;

        // Compare the max alignments for this read
        if (max_scoring_alignments1[read]->id == max_scoring_alignments3[read]->id) {
            pos_count++;
            node_accuracy += project_weighted_path_accuracy(max_scoring_alignments1[read], max_scoring_alignments3[read]);
            cigar_accuracy += projectA_cigar_accuracy(&max_scoring_alignments1[read]->cigar_string, &max_scoring_alignments3[read]->cigar_string);
        }
   }
    pos_accuracy = (double)pos_count / max_scoring_alignments1.size();
    node_accuracy = (double)node_accuracy / max_scoring_alignments1.size();
    cigar_accuracy = (double)cigar_accuracy / max_scoring_alignments1.size();


    // for (const auto& nodes_pair : sim_positions) {
    //     const auto& nodes = nodes_pair.second;
    //     for (const auto& node : nodes) {
    //         fprintf(stderr, "%s\t", node.c_str());
    //     }
    //     fprintf(stderr, "\n");
    // }
    // fprintf(stderr, "\n");




    // Compare the max scoring alignments to the simulated positions
    // double pos_accuracy = 0;
    // int32_t pos_count = 0;
    // for (auto& read_pair : max_scoring_alignments1) {
    //     string read = read_pair.first;
    //     // fprintf(stderr ,"%s\n", read.c_str());
    //     // for (int i = 0; i < read_pair.second->nodes.size(); ++i) {
    //     //     fprintf(stderr, "%s\t", read_pair.second->nodes[i].c_str());
    //     // }
    //     // fprintf(stderr, "\n");
    //     // for (int i = 0; i < sim_positions[read].size(); ++i) {
    //     //     fprintf(stderr, "%s\t", sim_positions[read][i].c_str());
    //     // }
    //     // fprintf(stderr, "\n");

    //     // Compare the max alignment with the simulated alignment
    //     pos_accuracy += projectA_node_sub_set(sim_positions[read], max_scoring_alignments1[read]->nodes);
    // }
    // pos_accuracy = (double)pos_accuracy / max_scoring_alignments1.size();



    // // Compare the top 5 max scoring alignments to the simulated positions
    // double pos_accuracy = 0;
    // int32_t pos_count = 0;
    // for (auto& read_pair : five_max_scoring_alignments1) {
    //     string read = read_pair.first;

    //     // Compare the five max alignments with the simulated alignment
    //     for (auto& alignment : five_max_scoring_alignments1[read]) {
    //         int32_t is_sub_set = projectA_node_sub_set(sim_positions[read], alignment->nodes);
    //         if (is_sub_set) {
    //             pos_accuracy += is_sub_set;
    //             break;
    //         }
    //         // cerr << five_max_scoring_alignments1[read].size() << endl;
    //     }
    // }
    // pos_accuracy = (double)pos_accuracy / five_max_scoring_alignments1.size();



    auto duration_gssw = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    auto duration_gwfa = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    auto duration_s2s = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
    auto duration_gwfa_s2s = std::chrono::duration_cast<std::chrono::milliseconds>( t4 - t2);

    // Print result values to console
    // cigar_accuracy = cigar_accuracy / n_high_score;
    // node_accuracy = node_accuracy / n_high_score;
    fprintf(outputFile, "%f\t", cigar_accuracy);
    fprintf(outputFile, "%f\t", node_accuracy);
    fprintf(outputFile, "%f\n", pos_accuracy);

    fprintf(outputFile, "n alignments:\t%d\t", n_high_score);
    fprintf(outputFile, "all alignments:\t%zu\n", alignments1.size());

    fprintf(outputFile, "gssw: %lld\tgwfa+s2s: %lld\tgwfa: %lld\ts2s: %lld\n",
            duration_gssw.count(),
            duration_gwfa_s2s.count(),
            duration_gwfa.count(),
            duration_s2s.count());
    
    // fprintf(outputFile, "E:\t%i\t%i\t%u\t%u\t%f\n",
    //         match,
    //         mismatch,
    //         gap_open,
    //         gap_extend,
    //         pos_accuracy);
    // cerr << count << "\t" << alignments1.size() << endl;

    // // Iterate over the max_scoring_alignments to print the best scores
    // for (auto& read_pair : max_scoring_alignments1) {
    //     cerr << read_pair.first << "\n" << read_pair.second->score << endl;
    // }
    // for (auto& read_pair : max_scoring_alignments3) {
    //     cerr << read_pair.first << "\n" << read_pair.second->score << endl;
    // }

    for (auto& alignment : alignments1) {
        delete alignment;
    }
    for (auto& alignment : alignments2) {
        delete alignment;
    }
    for (auto& alignment : alignments3) {
        delete alignment;
    }




    // projectA_cigar_t test_cigar1 = projectA_parse_cigar_string("4D3M2I2D");
    // projectA_cigar_t test_cigar2 = projectA_parse_cigar_string("4D2M2I2D");

    // projectA_print_cigar(stderr, &test_cigar1);
    // projectA_print_cigar(stderr, &test_cigar2);

    // cerr << projectA_cigar_accuracy(&test_cigar1, &test_cigar2) << endl;


    for (auto& graph : graphs) {
        // projectA_print_graph(stderr, graph.graph);
        projectA_delete_hash_graph(graph.graph);
    }

    for (auto& graph : sim_graphs) {
        // projectA_print_graph(stderr, graph.graph);
        projectA_delete_hash_graph(graph.graph);
    }

    projectA_delete_hash_graph(ref_graph);
}

void run_tests_gnwa() {
    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    // projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/linear_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    // projectA_read_node_list(clusters, "./test_cases/node_list.txt");
    projectA_read_node_list(clusters, "./test_cases/node_list_small.txt");
    // projectA_read_node_list(clusters, "./test_cases/tests.txt");
    // projectA_read_node_list(clusters, "./test_cases/linear_node_list.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);



    vector<projectA_alignment_t*> alignments1;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignments1.push_back(alignment);
    }

    projectA_get_alignment_gnwa(alignments1, 46);



    for (auto& alignment : alignments1) {
        delete alignment;
    }


    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
}



void run_benchmark(string graphFile, string positionFile, string simPositionFile, int match, int mismatch, uint gap_open, uint gap_extend, FILE* outputFile) {
    
    int32_t runtime = 0;
    int32_t setup_runtime = 0;
    typedef std::chrono::high_resolution_clock Clock;
    
    
    
    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa(graphFile);
    // projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/linear_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    // projectA_read_node_list(clusters, "./test_cases/node_list_2.txt");
    // projectA_read_node_list(clusters, "./test_cases/node_list_2_small.txt");
    projectA_read_node_list(clusters, positionFile);
    // projectA_read_node_list(clusters, "./test_cases/tests.txt");
    // projectA_read_node_list(clusters, "./test_cases/linear_node_list.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);


    // Read simulated positions
    unordered_map<string, vector<string>> sim_positions;
    projectA_read_sim_positions(sim_positions, simPositionFile);
    // projectA_read_sim_positions_from_two_files(sim_positions, "./test_cases/sim_reads_new.txt", "./test_cases/sim_paths.txt");
 
    // Maps to hold the different read sets
    unordered_map<string, vector<projectA_alignment_t*>> read_set_map;
    // Maps to hold the max scoring alignments for each read
    unordered_map<string, projectA_alignment_t*> max_scoring_alignments;
    // Maps to hold the top 5 max scoring alignments for each read
    unordered_map<string, vector<projectA_alignment_t*>> five_max_scoring_alignments;

    fprintf(outputFile, "n_threads,runtime,setup_runtime\n");

    for (int i = 0; i <= 1; ++i) {
        vector<projectA_alignment_t*> alignments;
        for (int32_t i = 0; i < graphs.size(); ++i) {
            projectA_alignment_t* alignment = new projectA_alignment_t;
            alignment->id = i;
            alignment->graph = graphs[i].graph;
            alignment->read = graphs[i].read;
            alignment->match = match;
            alignment->mismatch = mismatch;
            alignment->gap_open = gap_open;
            alignment->gap_extend = gap_extend;
            alignments.push_back(alignment);
            projectA_create_read_sets(read_set_map, alignment);
        }

        runtime = projectA_get_timed_alignment_s_gwfa(alignments, 2^i);
        auto t0 = Clock::now();
        // for (auto alignment : alignments) {
        //     projectA_edlib(alignment);
        // }
        auto t1 = Clock::now();
        projectA_get_alignment_s_gwfa(alignments, 2^i);
        // for (auto alignment : alignments) {
        //     projectA_edlib(alignment);
        // }
        auto t2 = Clock::now();
        runtime += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        setup_runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

        fprintf(outputFile, "%i,%i,%i\n", 2^i, runtime, setup_runtime);

        for (auto& alignment : alignments) {
            delete alignment;
        }
    }
    

    for (auto& graph : graphs) {
        // projectA_print_graph(stderr, graph.graph);
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
}

void sim_reads_and_path(string inputFile, string outputFile1, string outputFile2, int32_t n_reads, int32_t read_length) {

    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa(inputFile);
    projectA_index_hash_graph(ref_graph);
    vector<vector<string>> paths;
    vector<string> reads;

    srand(time(0));

    for (int i = 0; i < n_reads; ++i) {
        int32_t r = rand() % ref_graph->n_nodes;
        int32_t len = 0;
        vector<string> path;
        auto it = ref_graph->nodes.begin();
        advance(it, r);
        projectA_node_t* node = it->second;

        path.push_back(node->id);
        len += node->len;
        reads.push_back(node->seq);
        
        while (len < read_length) {
            projectA_node_t* next_node = node->next[r % (node->next.size())];
            path.push_back(next_node->id);
            reads[i] += next_node->seq;
            len += next_node->len;
            node = next_node;
        }

        paths.push_back(path);
    }

    cerr << paths.size() << endl;

    // Write paths to the first output file
    std::ofstream pathFile(outputFile1);
    if (!pathFile.is_open()) {
        std::cerr << "Error: Unable to open output file for paths: " << outputFile1 << std::endl;
        return;
    }
    for (const auto& path : paths) {
        for (size_t j = 0; j < path.size(); ++j) {
            pathFile << path[j];
            if (j < path.size() - 1) {
                pathFile << "\t";
            }
        }
        pathFile << "\n";
    }
    pathFile.close();

    // Write reads to the second output file
    std::ofstream readsFile(outputFile2);
    if (!readsFile.is_open()) {
        std::cerr << "Error: Unable to open output file for reads: " << outputFile2 << std::endl;
        return;
    }
    for (const auto& read : reads) {
        readsFile << read << "\n";
    }
    readsFile.close();
    
    projectA_delete_hash_graph(ref_graph);
}

char gen_random_base() {
    switch (rand() % 4) {
        case 0: {
            return 'A';
        }
        case 1: {
            return 'C';
        }
        case 2: {
            return 'G';
        }
        case 3: {
            return 'T';
        }
        default: {
            cerr << "error: rand is not in scope!\n";
            exit(1);
        }
    }
}

void gen_ref_seq_pair(string& ref, string& seq, int32_t len, int32_t n_changes) {
    for (int i = 0; i < len; ++i) {
        int32_t r_char = rand() % 4;
        switch (r_char) {
            case 0:
                ref += 'A';
                break;
            case 1:
                ref += 'C';
                break;
            case 2:
                ref += 'G';
                break;
            case 3:
                ref += 'T';
                break;
            default:
                cerr << "error: rand is not in scope!\n";
                exit(1);
        }
    }
    seq = ref;
    for (int i = 0; i < n_changes; ++i) {
        int32_t type = rand() % 3;
        switch (type) {
            case 0: {
                // Mismatch
                int32_t pos = rand() % seq.size();
                char r_base = gen_random_base();
                while (r_base == seq[pos]) {
                    r_base = gen_random_base();
                }
                seq[pos] = r_base;
                break;
            }
            case 1: {
                // Insertion
                int32_t pos = rand() % seq.size();
                seq.insert(pos, string(1, gen_random_base()));
                break;
            }
            case 2: {
                // Deletion
                int32_t pos = rand() % seq.size();
                seq.erase(pos, 1);
                break;
            }
            default: {
                cerr << "error: rand is not in scope!\n";
                exit(1);
            }
        }
    }
}

void test_s_gwfa_edlib() {

    typedef std::chrono::high_resolution_clock Clock;

    srand(time(0));
    int32_t r = rand();

    int32_t len = 10000;
    int32_t n_changes = 100;
    string ref;
    string seq;
    int32_t n_equal = 0;


    // s2s tests
    // ed tests
    for (int i = 0; i < 100; ++i) {
        ref = "";
        seq = "";
        gen_ref_seq_pair(ref, seq, len, n_changes);
        s_gwfa_node_t* node = s_gwfa_node_create(0, strlen(ref.c_str()), ref.c_str());
        s_gwfa_graph_t* graph = s_gwfa_graph_create();
        s_gwfa_graph_add_node(graph, node);
        EdlibAlignResult result = edlibAlign(seq.c_str(), strlen(seq.c_str()), 
                                            ref.c_str(), strlen(ref.c_str()), 
                                            edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
        int32_t edlib_ed = result.editDistance;

        // int32_t wfa_editdistance = wfa_ed(seq.c_str(), strlen(seq.c_str()), ref.c_str(), strlen(ref.c_str()));

        void* km = km_init();
        k_gwfa_path_t path_k_gwfa;
        int32_t wfa_editdistance = k_gwfa_ed(km, graph, strlen(seq.c_str()), seq.c_str(), 1, &path_k_gwfa);
        km_destroy(km);

        cerr << "Edlib: " << edlib_ed << "\tWFA_ed: " << wfa_editdistance << endl;

        if (edlib_ed == wfa_editdistance) {
            n_equal++;
        } else {
            cerr << "\t\t\tref: " << ref << endl << "\t\t\tseq: " << seq << endl;
        }
    }

    // infix tests
    // for (int i = 0; i < 100; ++i) {
    //     ref = "";
    //     seq = "";
    //     gen_ref_seq_pair(ref, seq, len, n_changes);

    //     ref = "XXXXXXXXXXXXXXXXXXXXX" + ref + "XXXXXXXXXXXXXXXXXXXXXX";

    //     EdlibAlignResult result = edlibAlign(seq.c_str(), strlen(seq.c_str()), 
    //                                         ref.c_str(), strlen(ref.c_str()), 
    //                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    //     int32_t edlib_ed = result.editDistance;

    //     int32_t wfa_editdistance = wfa_ed_infix(seq.c_str(), strlen(seq.c_str()), ref.c_str(), strlen(ref.c_str()));

    //     cerr << "Edlib: " << edlib_ed << "\tWFA_ed: " << wfa_editdistance << endl;

    //     if (edlib_ed == wfa_editdistance) {
    //         n_equal++;
    //     } else {
    //         cerr << "\t\t\tref: " << ref << endl << "\t\t\tseq: " << seq << endl;
    //     }
    // }

    // // prefix tests
    // for (int i = 0; i < 100; ++i) {
    //     ref = "CTGAAGATTTTCG";
    //     seq = "TGCTACTTGGATGGGCATGTGC";
    //     // gen_ref_seq_pair(ref, seq, len, n_changes);
    //     s_gwfa_node_t* node = s_gwfa_node_create(0, strlen(ref.c_str()), ref.c_str());
    //     s_gwfa_graph_t* graph = s_gwfa_graph_create();
    //     s_gwfa_graph_add_node(graph, node);

    //     ref = ref;

    //     EdlibAlignResult result = edlibAlign(seq.c_str(), strlen(seq.c_str()), 
    //                                         ref.c_str(), strlen(ref.c_str()), 
    //                                         edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    //     int32_t edlib_ed = result.editDistance;

    //     // int32_t wfa_editdistance = wfa_ed_prefix(seq.c_str(), strlen(seq.c_str()), ref.c_str(), strlen(ref.c_str()));
    //     void* km = km_init();
    //     int32_t wfa_editdistance = k_gwfa_ed(km, graph, strlen(seq.c_str()), seq.c_str());
    //     km_destroy(km);

    //     s_gwfa_graph_destroy(graph);

    //     cerr << "Edlib: " << edlib_ed << "\tWFA_ed: " << wfa_editdistance << endl;

    //     if (edlib_ed == wfa_editdistance) {
    //         n_equal++;
    //     } else {
    //         cerr << "\t\t\tref: " << ref << endl << "\t\t\tseq: " << seq << endl;
    //     }
    // }
    cerr  << endl << "n_equal:\t" << n_equal << endl;
}

void print_node(FILE* file, s_gwfa_node_t* node) {
    fprintf(file, "S\t%i\t%s\n", node->id, node->seq);
}


void print_edge(FILE* file, s_gwfa_edge_t* edge) {
    fprintf(file, "L\t%i\t+\t%i\t+\t0M\n", edge->start->id, edge->end->id);
}


void print_graph(FILE* file, s_gwfa_graph_t* graph) {

    s_gwfa_edge_t** edges = (s_gwfa_edge_t**)malloc(0);
    int n_edges = 0;

    // Print all the nodes and gather all the edges
    for (int i = 0; i < graph->size; ++i) {
        print_node(file, graph->nodes[i]);
        edges = (s_gwfa_edge_t**)realloc(edges, sizeof(s_gwfa_edge_t*) * (n_edges + graph->nodes[i]->n_edges));
        for (int j = 0; j < graph->nodes[i]->n_edges; ++j) {
            edges[n_edges + j] = graph->nodes[i]->edges[j];
        }
        n_edges += graph->nodes[i]->n_edges;
    }

    // Print all the edges
    for (int i = 0; i < n_edges; ++i) {
        print_edge(file, edges[i]);
    }

    free(edges);
}


void s_gwfa_test_gwfa() {

    typedef std::chrono::high_resolution_clock Clock;

    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    // projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/linear_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    // projectA_read_node_list(clusters, "./test_cases/node_list_2.txt");
    // projectA_read_node_list(clusters, "./test_cases/node_list_2_small.txt");
    projectA_read_node_list(clusters, "./test_cases/node_list_single_node.txt");
    // projectA_read_node_list(clusters, "./test_cases/tests.txt");
    // projectA_read_node_list(clusters, "./test_cases/linear_node_list.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);

    projectA_hash_graph_t* test_graph = projectA_hash_read_gfa("./test_cases/test_graph_2.gfa");
    projectA_index_hash_graph(test_graph);
    projectA_hash_graph_in_order_nodes(test_graph);
    projectA_algorithm_input_t inp;
    inp.graph = test_graph;
    inp.read = "AGCAATTCTGTAGCCGTACA";
    // graphs.push_back(inp);

    vector<projectA_alignment_t*> alignments1;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignments1.push_back(alignment);
    }

    // void* ptr = projectA_s_gwfa_init(alignments1, 1);
    // ptr = projectA_s_gwfa_calculate_batch(ptr, 0);
    // projectA_s_gwfa_post(ptr, alignments1, 1);
    // projectA_get_alignment_s_gwfa(alignments1, 1);
    int32_t time_gwfa = projectA_get_timed_alignment_gwfa(alignments1, 1);
    // projectA_get_alignment_gwfa(alignments1, 1);
    int duration_k = 0;
    int duration_s = 0;

    for (auto& alignment: alignments1) {
        s_gwfa_graph_t* graph = projectA_s_gwfa_graph_build(alignment->graph);
        // s_gwfa_path_t** path = (s_gwfa_path_t**)malloc(sizeof(s_gwfa_path_t*));
        print_graph(stderr, graph);
        cerr << alignment->read << endl;
        // int32_t score = g_wfa_ed_infix(alignment->read.c_str(), alignment->read.length(), graph, path);
        // projectA_print_graph(stderr, alignment->graph);
        // void* km = km_init();
        // k_gwfa_path_t path_k_gwfa;
        // auto t0 = Clock::now();
        // int32_t score2 = k_gwfa_ed(km, graph, alignment->read.length(), alignment->read.c_str(), 0, &path_k_gwfa);
        // auto t1 = Clock::now();
        // km_destroy(km);

        auto t2 = Clock::now();
        // int32_t score = g_wfa_ed_infix(alignment->read.c_str(), alignment->read.length(), graph, path);
        auto t3 = Clock::now();

        // cerr << alignment->score << "\t" << score << "\t" << score2 << endl;
        cerr << alignment->score << endl;
        // string g_wfa_ref = "";
        // for (int i = 0; i < (*path)->size; ++i) {
        //     fprintf(stderr, "[%i]->", (*path)->nodes[i]->id);
        //     g_wfa_ref += (*path)->nodes[i]->seq;
        // }
        // fprintf(stderr, "\n");
        // s_gwfa_path_destroy(*path);
        // free(path);
        s_gwfa_graph_destroy(graph);
        for (auto& node : alignment->nodes) {
            fprintf(stderr, "[%s]->", node.c_str());
        }
        fprintf(stderr, "\n");
        // EdlibAlignResult result1 = edlibAlign(alignment->read.c_str(), strlen(alignment->read.c_str()), 
        //                                     alignment->reference.c_str(), strlen(alignment->reference.c_str()), 
        //                                     edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        // int32_t edlib_ed = result1.editDistance;
        // EdlibAlignResult result2 = edlibAlign(alignment->read.c_str(), strlen(alignment->read.c_str()), 
        //                                     g_wfa_ref.c_str(), strlen(g_wfa_ref.c_str()), 
        //                                     edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        // int32_t edlib_g_wfa = result2.editDistance;
        // cerr << edlib_ed << "\t" << edlib_g_wfa << endl;
        // cerr << edlib_g_wfa << endl;
        // cerr << endl;
        // duration_k += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        // duration_s += std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    }

    cerr << duration_k << "\t" << duration_s << "\t" << time_gwfa << endl;

    for (auto& alignment : alignments1) {
        delete alignment;
    }


    for (auto& graph : graphs) {
        projectA_delete_hash_graph(graph.graph);
    }
    projectA_delete_hash_graph(ref_graph);
}


void test_gwf_ed_infix() {

    vector<projectA_node_list_t> clusters;
    projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/reference_graph.gfa");
    // projectA_hash_graph_t* ref_graph = projectA_hash_read_gfa("./test_cases/linear_graph.gfa");
    projectA_index_hash_graph(ref_graph);


    projectA_read_node_list(clusters, "./test_cases/node_list_2.txt");
    // projectA_read_node_list(clusters, "./test_cases/node_list_2_small.txt");
    // projectA_read_node_list(clusters, "./test_cases/node_list_single_node.txt");
    // projectA_read_node_list(clusters, "./test_cases/tests.txt");
    // projectA_read_node_list(clusters, "./test_cases/linear_node_list.txt");
    vector<projectA_algorithm_input_t> graphs;
    projectA_build_graph_from_cluster(graphs, ref_graph, clusters);

    projectA_hash_graph_t* test_graph = projectA_hash_read_gfa("./test_cases/test_graph_2.gfa");
    projectA_index_hash_graph(test_graph);
    projectA_hash_graph_in_order_nodes(test_graph);
    projectA_algorithm_input_t inp;
    inp.graph = test_graph;
    inp.read = "AGCAATTCTGTAGCCGTACA";
    // graphs.push_back(inp);

    vector<projectA_alignment_t*> alignments1;
    for (int32_t i = 0; i < graphs.size(); ++i) {
        projectA_alignment_t* alignment = new projectA_alignment_t;
        alignment->id = i;
        alignment->graph = graphs[i].graph;
        alignment->read = graphs[i].read;
        alignments1.push_back(alignment);
    }

    projectA_get_alignment_gwfa(alignments1, 1);

    int n = 0;

    for (auto alignment : alignments1) {
        EdlibAlignResult result2 = edlibAlign(alignment->read.c_str(), strlen(alignment->read.c_str()), 
                                            alignment->reference.c_str(), strlen(alignment->reference.c_str()), 
                                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        int32_t ed = result2.editDistance;
        if (ed != alignment->score) {
            cerr << ed << "\t" << alignment->score << endl;
            cerr << endl;
            n++;
        }
    }
    
    cerr << (double)n/alignments1.size() << endl;
    
}



int main() {

    int match = 1;
    int mismatch = 1;
    uint gap_open = 1;
    uint gap_extend = 1;

    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2_small.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2_medium.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/tests.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/linear_node_list.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend);
    // run_tests_gnwa();
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_3.txt", "", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_4.txt", "", match, mismatch, gap_open, gap_extend, stderr);

    // FILE* outputFile = fopen("./files/gssw_benchmark.csv", "w");
    // run_benchmark("./test_cases/reference_graph.gfa", "./test_cases/node_list_2.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // fclose(outputFile);


    // projectA_get_alignment_abpoa(alignments2, 1);
    // FILE* outputFile = fopen("./files/gssw_results.txt", "w");
    // fprintf(outputFile, "H:match//mismatch//gap_open//gap_extend//pos_accuracy\n");

    // for (int i = 1; i <= 40; i += 10) {
    //     for (int j = 1; j <= 50; j += 10) {
    //         for (int k = 1; k <= 50; k += 10) {
    //             fprintf(stderr, "Running test with i=%d, j=%d, k=%d\n", i, j, k);

    //             try {
    //                 // Call the function and handle exceptions if applicable
    //                 run_standard_tests("./test_cases/reference_graph.gfa", 
    //                                 "./test_cases/node_list_2.txt", 
    //                                 "./test_cases/1000.new.sim.txt", 
    //                                 i, j, k, k, outputFile);
    //             } catch (const std::exception& e) {
    //                 // Catch and log standard exceptions
    //                 fprintf(stderr, "Error: Exception during test run: %s\n", e.what());
    //             } catch (...) {
    //                 // Catch any other unknown errors
    //                 fprintf(stderr, "Error: Unknown exception during test run.\n");
    //             }
    //         }
    //     }
    // }

    // fclose(outputFile);

    // sim_reads_and_path("./test_cases/reference_graph.gfa", "./files/sim_paths.txt", "./files/sim_reads.txt", 100, 1000);

    // test_s_gwfa_edlib();
    // s_gwfa_test_gwfa();
    test_gwf_ed_infix();

    // string seq = "AAGTAAAACCCTTTACCATTAAGCAGCATTTCCCCTCCCCTAAGCCCTTCCCCTCAGCCCCACCAACCTGCATCTGACTCCATGGGCTTATCTACTCTGGATATCTCATAAATATGGAACCATACAATATGTGACCTTTTGTGTCTGGTTTCTTTCATTTAGCATGATGTTTTCAAGGATATCAGGATTTTTAAAACTTTCCAAGTGATTCTGATATGGAATCAAGATATCCACAAAGCAGGTGTGAGGCAGGGTAGACAATCAGATTCCATTTTAACTAAGCACTTGTGTGCAAGGAAACCACAAATATCTCATTTAATCCTTACAACAACTTGTGTGTATATCATTCCCATTTTACAGGCCGGGAAACTGAAGCTCAGAGAGGTGAAATCACCCACCTGACATCACATAGGTATTTAAAAGCAAATTGGTCTGACTCTCCAGCACCATGTCTCCTCTGCAACCCAGGCCCCTCTCCCAAATTAGCACTTGGTAACCATCAGAAAGGAAGAGTGGGCTGTGCCTACTGGCTGGCAAGCCTGGAATCTATAGGCAGAGTTGATCTACTGGCTAGAGAACATAGGACTTAAGGGTATTTGAGCTCCTTAACTTTTAAAATACCCTGCAACTAGAAAGAAAGTCAACAAACATGAAAGCTGGAGAGAACTTAGACCAATTTGTAAAAGCACTCATTTTTAAAACCAAAAACCCAAGGCCCAAAGACAGCAGACTGGCAAAACATGGTCAGATGCCTTCTGAATTTCAAAAAAGCTATCAGTCATTTTTACTTTTCATTGATCCTGAGAACAGGGTCTTCCCCAAACATAAAAACAAAACTCATTTTTACACAGTTTGACTACCACCGTATTATTAGCCATGGTTTCAATTCATTTAGGCAACAGATACCATGCATATACTATGTTCCAGGCACATTCTAGGGAGGAGGGATATTGCAATGAACAAAGACAGATAACACCCCCTGTCCCTCACAGAGCTACAT";
    // string ref = "TTGTCACCTTGTCTACTCCTGTTCCTCAGGGACTCACCTGCCTTAGCTCAGTCCTTGCTGGGCAATTCTCTGTGTCTTTTCTAGGGCCTGCAAGTTCTGCCACTGGAGTTTACAAATGACCACTATCTTTAGCCAGATAACCAGGTATGCCCTAGCATTGGGCAGTTTTCCTCTACATCTTGCCCTCTCTGTCTCAAGGATTTAGGACCTAGATTATGTTGCTGCTTCAGGCAGTTCTGCTGCTGACAAGATAAATTCTGTTAGCATTTCACTAAGCTGATGGGTATAATTACTTTAGTTGGAGGTCACTAAAACCTGAATTTACAAATAAAAACTGGATCAGGGAGCAGTCTGTGTGTGTGTTGTGGAGGAGTGTGGGTGAAGGTACAACAGATAATATCCCTCCTTTGCTGAGATCCTTCCAGTGGCTCCCTAAATCACTCATAGTAAAAGCCAACTTCCTCACAGGGCCCTTCATCATGTTCTGCTGCTAGCCTCCATCTGCTTCTTTCCCTTACTCAACCTAATTCTAATTTCAGTTCCATTGGAGCCATACAGATCTCCTTGCTATTTCTTTATTTCAAATTCCAGTTATTCATTGCTGGCATATAGAAAAGCAATAGATTGTTGTATATTAACTTTGTGTCTTGTAACCTTGCTGTAAGTGCTGTTCTGGCAGTCTTTGTTTTGATCAATTCTTTGGGATATTCTACATAGACAGATCATCTGCAAACAAAGATAGTGTTATTTCTTCCATTTCAATCTGTAGAGCTTTTATTTCCTTTTCTGGCTTATTCCACTCGCTAAGATTTCCAGTAAAATGTTGAATAGGAGTAGAGAGAGAGGCACCTTGCTTTGTTCCTGATCTTAAGGGGAAAGCATCCAATTTCTTACCATTAAGTATGAAGTTGGCTGCAAAGAAGTTCTTTATCCAGTTGAGGAAATTCCCTCCTATTCATAGTTTGCTGCCCCATCAGCTTTTTAGGTTGGTGCAAAAGTAATTGTGGTTTTTGCTGTGACTTT";

    // EdlibAlignResult result = edlibAlign(seq.c_str(), strlen(seq.c_str()), 
    //                                         ref.c_str(), strlen(ref.c_str()), 
    //                                         edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    // int32_t edlib_ed = result.editDistance;

    // cerr << edlib_ed << endl;
    

    cerr << "run succesfull!\n";
    return 0;
}
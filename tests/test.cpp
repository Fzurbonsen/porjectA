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

#include "test.hpp"
#include "extract_graph.hpp"
#include "algorithm.hpp"
// #include "algorithms/gt_gwfa.hpp"
#include "algorithms/gssw.hpp"
#include "algorithms/gwfa.hpp"
#include "alignment.hpp"
#include "algorithms/ksw2.hpp"
#include "algorithms/csswl.hpp"
// #include "algorithms/abPOA.hpp"
// #include "algorithms/vargas.hpp"
#include "algorithms/gnwa.hpp"

// #include "gt_gwfa/edlib.h"
#include "algorithms/edlib.hpp"


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

    typedef std::chrono::high_resolution_clock Clock;

    auto t0 = Clock::now();
    // projectA_get_alignment_gssw(alignments1, 16);
    // projectA_get_alignment_gnwa(alignments1, 16);
    auto t1 = Clock::now();

    auto t2 = Clock::now();
    // projectA_get_alignment_gwfa(alignments3, 16);
    // projectA_get_alignment_gwfa(alignments1, 32);
    projectA_get_alignment_gssw(alignments1, 16);
    auto t3 = Clock::now();


    // Truncate references for gwfa
    for (int i = 0; i < alignments1.size(); ++i) {

        // Strings to hold the cut reference and read
        // string new_reference;
        // string new_read;

        // projectA_get_gssw_reference(alignments1[i], new_reference);
        // projectA_get_gssw_read(alignments1[i], new_read);

        // alignments2[i]->reference = new_reference;
        // alignments2[i]->read = new_read;
    }

    for (int32_t i = 0; i < alignments2.size(); ++i) {
        auto& alignment1 = alignments1[i];
        auto& alignment2 = alignments2[i];
        auto& alignment3 = alignments3[i];

        // csswl
        // alignment2->reference = alignment1->reference;
        // alignment2->read = alignment1->read;
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
        // projectA_ksw2(alignment2);
        // alignment2->offset = 0;

    }
    auto t4 = Clock::now();


    // Find the max scoring alignment for each read
    projectA_determine_max_score_alignment(max_scoring_alignments1, read_set_map1);
    projectA_determine_max_score_alignment(max_scoring_alignments2, read_set_map2);
    projectA_determine_max_score_alignment(max_scoring_alignments3, read_set_map3);



    FILE* file = fopen("./files/CIGAR.txt", "w");
    double cigar_accuracy = 0;
    double node_accuracy = 0;
    int count = 0;
    int n_high_score = 0;
    for (int i = 0; i < alignments1.size(); ++i) {

        double match_per_length = (float)(alignments1[i]->n_matches)/(float)(alignments1[i]->cigar_string.operations_length);

        // if (alignments1[i]->score > 200) {
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
            projectA_write_to_test_file(file, alignments1, alignments3, i);
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
            // fprintf(stderr, "[%s]\n\n", alignments1[i]->graph->nodes_in_order[alignments1[i]->graph->nodes_in_order.size()-1]->id.c_str());
        }
    }
    fclose(file);



    // // Compare the max scoring alignments of alignments1 and alignments2
    // double pos_accuracy = 0;
    // int32_t pos_count = 0;
    // for (auto& read_pair : max_scoring_alignments1) {
    //     string read = read_pair.first;

    //     // Compare the max alignments for this read
    //     if (max_scoring_alignments1[read]->id == max_scoring_alignments3[read]->id) {
    //         pos_count++;
    //     }
    // }
    // pos_accuracy = (double)pos_count / max_scoring_alignments1.size();

    // for (const auto& nodes_pair : sim_positions) {
    //     const auto& nodes = nodes_pair.second;
    //     for (const auto& node : nodes) {
    //         fprintf(stderr, "%s\t", node.c_str());
    //     }
    //     fprintf(stderr, "\n");
    // }
    // fprintf(stderr, "\n");

    // Compare the max scoring alignments to the simulated positions
    double pos_accuracy = 0;
    int32_t pos_count = 0;
    for (auto& read_pair : max_scoring_alignments1) {
        string read = read_pair.first;
        fprintf(stderr ,"%s\n", read.c_str());
        for (int i = 0; i < read_pair.second->nodes.size(); ++i) {
            fprintf(stderr, "%s\t", read_pair.second->nodes[i].c_str());
        }
        fprintf(stderr, "\n");
        for (int i = 0; i < sim_positions[read].size(); ++i) {
            fprintf(stderr, "%s\t", sim_positions[read][i].c_str());
        }
        fprintf(stderr, "\n");

        // Compare the max alignment with the simulated alignment
        pos_accuracy += projectA_node_sub_set(sim_positions[read], max_scoring_alignments1[read]->nodes);
    }
    pos_accuracy = (double)pos_accuracy / max_scoring_alignments1.size();

    auto duration_gssw = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    auto duration_gwfa = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
    auto duration_s2s = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
    auto duration_gwfa_s2s = std::chrono::duration_cast<std::chrono::milliseconds>( t4 - t2);

    // Print result values to console
    cigar_accuracy = cigar_accuracy / n_high_score;
    node_accuracy = node_accuracy / n_high_score;
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


int main() {

    int match = 1;
    int mismatch = 100;
    uint gap_open = 100;
    uint gap_extend = 100;


    run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2_small.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_2_medium.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/tests.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend);s
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/linear_node_list.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend);
    // run_tests_gnwa();
    // run_standard_tests("./test_cases/reference_graph.gfa", "./test_cases/node_list_3.txt", "./test_cases/1000.new.sim.txt", match, mismatch, gap_open, gap_extend, stderr);


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

    cerr << "run succesfull!\n";
    return 0;
}
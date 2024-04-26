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


// PRE:     fileName, read_vector, graph_vector
//      fileName:       String containing the path to a file that holds the neccesary information
//                      to build the test graphs and reads.
//      read_vector:    Empty vector that has the string type.
//      graph_vector:   Empty vector that has the graph type.
// POST:    read_vector, graph_vector
//      read_vector:    Vector containing the reads for each test alignment.
//      graph_vector:   Vector containing the graphs for each test alignment.
//      return:         Integer value of the size of the test vectors.
int import_tests(string fileName, vector<string>& read_vector);


// PRE:     gold_vector, read_vector, graph_vector
//      gold_vector:    Empty vector that has the string type.
//      read_vector:    Vector holding the reads for each test alignment.
//      graph_vector:   Vector holding the graphs for each test alignment.
// POST:    gold_vector
//      gold_vector:    Vector holing the gold results for each test alignment.
void gold_results(vector<string>& gold_vector, vector<string>& read_vector);


// PRE:     test_vector, read_vector, graph_vector
//      test_vector:    Empty vector that has the string type.
//      read_vector:    Vector holding the reads for each test alignment.
//      graph_vector:   Vector holding the graphs for each test alignment.
// POST:    test_vector
//      test_vector:    Vector holing the test results for each test alignment.
void test_results(vector<string>& result_vector, vector<string>& read_vector);


// PRE:     n_tests, correct_tests, gold_vector, result_vector
//      n_tests:        Integer holding the numer of test cases in the test vector.
//      correct_tests:  Integer of value 0.
//      gold_vector:    Vector holing the gold results for each test alignment.
//      test_vector:    Vector holing the test results for each test alignment.
// POST:    correct_tests
//      correct_tests:  Integer holding the amount of test cases where the gold and test vector
//                      hold the same result.
void compare_results(int& n_tests, int& correct_tests, 
                        vector<string>& gold_vector, vector<string>& result_vector);

// PRE:     fileName
//      fileName:       String containing the path to a file that holds the neccesary information
//                      to build the test graphs and reads. 
// POST:
int run_tests(string fileName);
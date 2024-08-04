/*

    projectA:
    file_io.cpp
    This file holds the implementation of the projectA I/O functions to read from and write to files.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>
    
*/


#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>

#include "file_io.hpp"

using namespace std;

// Function to read contents of a graph file
void projectA_read_graph_file(string& fileName) {


    // Initialise variables
    bool skip;
    int8_t state = 0; // 0 = default, 1 = graph, 2 = read
    ifstream file(fileName);
    string line, ignore, read;
    char type, dir;
    uint32_t seq_len;
    string sequence, id, id1, id2;

    // Check if file could be opened
    if (!file.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        exit(0);
    }

    // Read file
    while (getline(file, line)) {

        // Skip empty lines
        if (line.empty() || line[0] == '#') continue;

        // Check if we have a graph next
        if (line[0] == 'N' && line[1] == 'G' && state == 0) {
            state++;
        }
            
        // Check if line is node
        if (line[0] == 'S' && state == 1) {
            // Initialise string stream
            istringstream iss(line);
            
            // Read values from stream
            iss >> type >> id >> sequence;
            // Extract sequence length
            size_t length_index = line.find("LN:i:");
            if (length_index != string::npos) {
                istringstream length_stream(line.substr(length_index + 5));
                length_stream >> seq_len;
            }

            // Compare sequence length if they don't then we skip to the next graph
            if (sequence.length() != seq_len) {
                cerr << "Error: sequnce length does not match, skipping to next graph!\n";
                skip = true;
                break;
            }
            // TODO:
            // Create corresponding node
            // cerr << type << "\t" << id << "\t" << sequence << "\t" << seq_len << endl;
        }

        // Check if line is edge
        if (line[0] == 'L' && state == 1) {
            // Initialise string stream
            istringstream iss(line);

            // Read values from stream
            iss >> type >> id1 >> dir >> id2 >> ignore;
            // TODO:
            // Create corresponding edges
            cerr << type << "\t" << id1 << "\t" << dir << "\t" << id2 << endl;
        }


        // Check if we have a read next
        if (line[0] == 'N' && line[1] == 'R' && state == 1) {
            state++;
        }

        // Check if line is read
        if (line[0] == 'R' && state == 2) {
            istringstream iss(line);
            iss >> type >> read;
            state = 0;
            // TODO:
            // Create corresponding reads
            cerr << type << "\t" << read << endl;
        }

        // Check for errors
        if ((state > 2) || (line[0] != 'N' && line[0] != 'S' && line[0] != 'L' && line[0] != 'R')) {
            cerr << "Error: Graph was not formated propperly, skipping to next graph!\n";
            state = 0;
        }
    }
}


// Function to read a cluster file
void projectA_read_cluster_file(string& fileName, vector<projectA_position_t>& positions) {
    
    // Initialise variables
    bool skip;
    int8_t state = 0; // 0 = default, 1 = graph, 2 = read
    ifstream file(fileName);
    string line, read;
    char type;
    uint32_t is_reverse, offset, forward_max_dist, backward_max_dist;
    string id, sequence;
    projectA_position_t pos;

    // Check if file could be opened
    if (!file.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        exit(0);
    }

    // Read file
    while (getline(file, line)) {

        // Skip empty lines
        if (line.empty() || line[0] == '#') continue;

        // Check if we have a position next
        if (line[0] == 'N' && line[1] == 'C' && state == 0) {
            state++;
        }
            
        // Check if line is position
        if (line[0] == 'P' && state == 1) {
            // Initialise string stream
            istringstream iss(line);
            
            // Read values from stream
            iss >> type >> id >> is_reverse >> offset >> forward_max_dist >> backward_max_dist;

            // cerr << type << "\t" << id << "\t" << is_reverse << "\t" << offset << "\t" << forward_max_dist << "\t" << backward_max_dist << endl;
            pos.id.push_back(id);
            pos.is_reverse.push_back(is_reverse);
            pos.offset.push_back(offset);
            pos.forward_search_lengths.push_back(forward_max_dist);
            pos.backward_search_lengths.push_back(backward_max_dist);
        }


        // Check if we have a read next
        if (line[0] == 'N' && line[1] == 'R' && state == 1) {
            state++;
        }

        // Check if line is read
        if (line[0] == 'R' && state == 2) {
            istringstream iss(line);
            iss >> type >> read;
            state = 0;
            // cerr << type << "\t" << read << endl;

            // Assign read to position chain
            pos.read = read;

            // Push new position chain to vector
            positions.push_back(pos);

            // Clear vectors and string to reuse
            pos.read.clear();
            pos.id.clear();
            pos.is_reverse.clear();
            pos.offset.clear();
            pos.forward_search_lengths.clear();
            pos.backward_search_lengths.clear();
        }

        // Check for errors
        if ((state > 2) || (line[0] != 'N' && line[0] != 'P' && line[0] != 'R')) {
            cerr << "Error: Position was not formated propperly, skipping to next Position!\n";
            state = 0;
        }
    }
}


// Function that reads a grpah file and stores the reads and the corresponding node ids in a umap
void projectA_read_graph_file(string& fileName, unordered_map<string, vector<vector<string>>>& graphs_map) {

    // Initialise variables
    bool skip;
    int8_t state = 0; // 0 = default, 1 = graph, 2 = read
    ifstream file(fileName);
    string line, ignore, read;
    char type, dir;
    uint32_t seq_len;
    string id, id1, id2, gssw_id, sequence;
    vector<string> id_vec;
    vector<vector<string>> graph_vec;

    // Check if file could be opened
    if (!file.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        exit(0);
    }

    // Read file
    while (getline(file, line)) {

        // Skip empty lines
        if (line.empty() || line[0] == '#') continue;

        // Check if we have a graph next
        if (line[0] == 'N' && line[1] == 'G' && state == 0) {
            state++;
        }
            
        // Check if line is node
        if (line[0] == 'S' && state == 1) {
            // Initialise string stream
            istringstream iss(line);
            
            // Read values from stream
            iss >> type >> id >> gssw_id >> sequence;
            // Extract sequence length
            size_t length_index = line.find("LN:i:");
            if (length_index != string::npos) {
                istringstream length_stream(line.substr(length_index + 5));
                length_stream >> seq_len;
            }

            // Compare sequence length if they don't then we skip to the next graph
            if (sequence.length() != seq_len) {
                cerr << "Error: sequnce length does not match, skipping to next graph!\n";
                skip = true;
                break;
            }

            // append current node to our id vector
            id_vec.push_back(gssw_id);

            // cerr << type << "\t" << id << "\t" << sequence << "\t" << seq_len << endl;
        }

        // Check if line is edge
        if (line[0] == 'L' && state == 1) {
            // Initialise string stream
            istringstream iss(line);

            // Read values from stream
            iss >> type >> id1 >> dir >> id2 >> ignore;
            // TODO:
            // Create corresponding edges
            // cerr << type << "\t" << id1 << "\t" << dir << "\t" << id2 << endl;
        }


        // Check if we have a read next
        if (line[0] == 'N' && line[1] == 'R' && state == 1) {
            state++;
        }

        // Check if line is read
        if (line[0] == 'R' && state == 2) {
            istringstream iss(line);
            iss >> type >> read;
            state = 0;

            // check if the read is already in the map
            if (graphs_map.count(read) > 0) {

                // if the key already exists then we append our current vector to the existing one
                graphs_map.at(read).push_back(id_vec);
                id_vec.clear();
            } else {

                // if the key does not already exist we create a new key value pair
                graph_vec.push_back(id_vec);
                graphs_map.insert(make_pair(read, graph_vec));
                graph_vec.clear();
                id_vec.clear();
            }
            
            // cerr << type << "\t" << read << endl;
        }

        // Check for errors
        if ((state > 2) || (line[0] != 'N' && line[0] != 'S' && line[0] != 'L' && line[0] != 'R')) {
            cerr << "Error: Graph was not formated propperly, skipping to next graph!\n";
            state = 0;
        }
    }
}




// Function that reads a cluster file and stores the reads and corresponding node ids in a umap
void projectA_read_cluster_file(string& fileName, unordered_map<string, vector<string>>& cluster_map) {
    
    // Initialise variables
    bool skip;
    int8_t state = 0; // 0 = default, 1 = graph, 2 = read
    ifstream file(fileName);
    string line, read;
    char type;
    uint32_t is_reverse, offset, forward_max_dist, backward_max_dist;
    string sequence, id;
    vector<string> cluster_vec;

    // Check if file could be opened
    if (!file.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        exit(0);
    }

    // Read file
    while (getline(file, line)) {

        // Skip empty lines
        if (line.empty() || line[0] == '#') continue;

        // Check if we have a position next
        if (line[0] == 'N' && line[1] == 'C' && state == 0) {
            state++;
        }
            
        // Check if line is position
        if (line[0] == 'P' && state == 1) {
            // Initialise string stream
            istringstream iss(line);
            
            // Read values from stream
            iss >> type >> id >> is_reverse >> offset >> forward_max_dist >> backward_max_dist;

            // cerr << type << "\t" << id << "\t" << is_reverse << "\t" << offset << "\t" << forward_max_dist << "\t" << backward_max_dist << endl;
        }


        // Check if we have a read next
        if (line[0] == 'N' && line[1] == 'R' && state == 1) {
            state++;
        }

        // Check if line is read
        if (line[0] == 'R' && state == 2) {
            istringstream iss(line);
            iss >> type >> read;
            state = 0;
            
            // Check if the read already exists in the map
            if (cluster_map.count(read) > 0) {

                // If it already exists then we append the corresponding id vector
                cluster_map.at(read).push_back(id);

            } else {

                // If it doesn't exist already we add a new key value pair
                cluster_vec.push_back(id);
                cluster_map.insert(make_pair(read, cluster_vec));
                cluster_vec.clear();
            }
            // cerr << type << "\t" << read << endl;
        }

        // Check for errors
        if ((state > 2) || (line[0] != 'N' && line[0] != 'P' && line[0] != 'R')) {
            cerr << "Error: Position was not formated propperly, skipping to next Position!\n";
            state = 0;
        }
    }
}


// Function to read graph from gfa and store as hash graph
projectA_hash_graph_t* projectA_hash_read_gfa(const string& fileName) {

    // Create a new projectA_hash_graph_t object
    projectA_hash_graph_t* graph = new projectA_hash_graph_t();
    // Initialize number of nodes
    graph->n_nodes = 0;

    // Open the file
    ifstream file(fileName);
    // Check if the file could be opened
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << fileName << endl;
        return nullptr;
    }

    string line;
    // Read each line of the file
    while (getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Read new line
        istringstream iss(line);
        string keyword;
        iss >> keyword;

        // Check if the line is a sequence/node
        if (keyword == "S") {
            uint32_t len;
            string id, seq, skip;
            iss >> id >> seq >> skip;
            len = seq.size();
            projectA_hash_graph_append_node(graph, id, len, seq, 0);
        }

        // Check if the line is an edge
        if (keyword == "L") {
            string start, end;
            string skip;
            iss >> start >> skip >> end >> skip;
            projectA_hash_graph_append_edge(graph, start, end);
        }
    }

    file.close();
    return graph;
}


// Function to read a node list file
void projectA_read_node_list(vector<projectA_node_list_t>& node_lists, const string& fileName) {

    projectA_node_list_t node_list;
    node_list.n_nodes = 0;

    // Open the file
    ifstream file(fileName);
    // Check if the file could be opened
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << fileName << endl;
        return;
    }

    string line;
    uint32_t line_counter = 0;
    bool n_line = false;

    // Read each line of the file
    while (getline(file, line)) {
        line_counter++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Read new line
        istringstream iss(line);
        string keyword;
        iss >> keyword;

        // Check if the line is a sequence/node
        if (keyword == "N") {


            // Check for valid preceding line
            if (n_line) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "Previous N-line was not followed by an R-line!\n";
                exit(1);
            }
            n_line = true;

            if (iss.peek() == EOF) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "N-line is empty!\n";
                exit(1);
            }

            string id;

            // Read ids until we reach the end of the line
            while(iss >> id) {
                if (!id.empty()) {
                    node_list.nodes.push_back(id); // Add the node id to the node_list struct
                    node_list.n_nodes++; // Increase the node counter
                }
            }

            // Check if there was at least one node in the N-line
            if (node_list.n_nodes == 0) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "N-line is empty!\n";
                exit(1);
            }
        }

        // Check if the line is an edge
        if (keyword == "R") {

            // Check for valid preceding line
            if (!n_line) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "R-line was not preceeded by an N-line!\n";
                exit(1);
            }
            n_line = false;

            // Check if there is an instance in the line
            if (iss.peek() == EOF) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "R-line is empty!\n";
                exit(1);
            }

            // Add read to the node_list struct
            iss >> node_list.read;
            node_list.read_len = node_list.read.size();

            // Check if the line actually holds a value
            if (node_list.read.empty()) {
                cerr << "Error: file not formatted correctly!\n";
                cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
                cerr << "R-line is empty!\n";
                exit(1);
            }

            // If we read a node we can add the struct to our node_lists vector
            node_lists.push_back(node_list);
            // Clear the node_list node vector for next input
            node_list.nodes.clear();
            node_list.read.clear();
            node_list.n_nodes = 0;
        }

        // Check for invalid line indicators
        if (keyword != "N" && keyword != "R") {
            cerr << "Error: file not formatted correctly!\n";
            cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
            cerr << keyword << " is not a valid line inidcator!\n";
            exit(1);
        }
    }

    // Check if we ended on an N-line
    if (n_line) {
        cerr << "Error: file not formatted correctly!\n";
        cerr << "Error in file: " << fileName << " in line: " << line_counter << endl;
        cerr << "Ended on an N-line!\n";
        exit(1);
    }

    file.close();
    return;
}


// Function to print a graph to a file in GFA format
void projectA_print_graph(FILE* file, projectA_hash_graph_t* graph) {
    // Check if the file is valid
    if (!file) {
        cerr << "Error: Invalid file pointer!" << endl;
        return;
    }

    // Write header
    fprintf(file, "H\tVN:Z:1.0\n");

    // Write segments
    for (auto& it : graph->nodes) {
        fprintf(file, "S\t%s\t%s\tLN:i:%d\n", it.second->id.c_str(), it.second->seq.c_str(), it.second->len);
    }

    // Write links (edges)
    for (auto& it : graph->nodes) {
        for (auto& next : it.second->next) {
            fprintf(file, "L\t%s\t+\t%s\t+\t0M\n", it.second->id.c_str(), next->id.c_str());
        }
    }
}


// Function to print alignment structs
void projectA_print_alignment(FILE* file, projectA_alignment_t* alignment) {

    // Print the score and offset in the first node
    fprintf(file, "%i@%u:", alignment->score, alignment->offset);

    // Iterate over all nodes
    for (int i = 0; i < alignment->size; ++i) {

        // Print the node id
        fprintf(file, "%s[", alignment->nodes[i].c_str());

        // Iterate over the cigar in the node
        for (auto& cigar_element : alignment->cigar[i].elements) {

            // Print cigar element
            fprintf(file, "%i%c", cigar_element.len, cigar_element.type);
        }

        fprintf(file, "]");
    }

    fprintf(file, "\n");
}


// Function to print CIGAR
void projectA_print_cigar(FILE* file, projectA_cigar_t* cigar) {

    // Print lenght
    // fprintf(file, "%i\t", cigar->len);
    
    // Iterate over all elements of the CIGAR
    for (auto& cigar_element : cigar->elements) {
        
        // Print the CIGAR element
        fprintf(file, "%i%c", cigar_element.len, cigar_element.type);
    }

    fprintf(file, "\n");
}


// Function to print alignment path
void projectA_print_path(FILE* file, vector<string>& path) {

    // Print lenght
    fprintf(file, "%i\t", path.size());

    // Iterate over all nodes
    for (auto& node : path) {

        // Print the path id
        fprintf(file, "[%s]->", node.c_str());
    }

    fprintf(file, "\n");
}
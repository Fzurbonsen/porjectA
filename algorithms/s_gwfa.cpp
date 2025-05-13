#include "algorithms/s_gwfa.hpp"

void s_gwfa_print_node(FILE* file, s_gwfa_node_t* node) {
    fprintf(file, "S\t%i\t%s\n", node->id, node->seq);
}


void s_gwfa_print_edge(FILE* file, s_gwfa_edge_t* edge) {
    fprintf(file, "L\t%i\t+\t%i\t+\t0M\n", edge->start->id, edge->end->id);
}


void s_gwfa_print_graph(FILE* file, s_gwfa_graph_t* graph) {

    s_gwfa_edge_t** edges = (s_gwfa_edge_t**)malloc(0);
    int n_edges = 0;

    // Print all the nodes and gather all the edges
    for (int i = 0; i < graph->size; ++i) {
        s_gwfa_print_node(file, graph->nodes[i]);
        edges = (s_gwfa_edge_t**)realloc(edges, sizeof(s_gwfa_edge_t*) * (n_edges + graph->nodes[i]->n_edges));
        for (int j = 0; j < graph->nodes[i]->n_edges; ++j) {
            edges[n_edges + j] = graph->nodes[i]->edges[j];
        }
        n_edges += graph->nodes[i]->n_edges;
    }

    // Print all the edges
    for (int i = 0; i < n_edges; ++i) {
        s_gwfa_print_edge(file, edges[i]);
    }

    free(edges);
}

// Constructor for parameter struct
projectA_s_gwfa_parameters_t::projectA_s_gwfa_parameters_t(s_gwfa_graph_t* graph, const char* seq, int32_t len,
                                                            projectA_hash_graph_t* projectA_hash_graph) :
                                                            graph(graph), seq(seq), len(len), 
                                                            projectA_hash_graph(projectA_hash_graph) {}



s_gwfa_graph_t* projectA_s_gwfa_graph_build(projectA_hash_graph_t* in_graph) {
    s_gwfa_graph_t* graph = s_gwfa_graph_create();
    std::unordered_map<projectA_node_t*, s_gwfa_node_t*> node_map;
    // Iterate over all nodes and add them to the graph
    int i = 0;
    for (auto& n : in_graph->nodes_in_order) {
        s_gwfa_node_t* node = s_gwfa_node_create(i, n->len, n->seq.c_str());
        s_gwfa_graph_add_node(graph, node);
        node_map[n] = node;
        for (auto& prev : n->prev) {
            s_gwfa_edge_create(node_map[prev], node_map[n], 0);
        }
        ++i;
    }
    graph->top_nodes = s_gwfa_graph_find_top_nodes(graph);
    return graph;
}



// Function to convert path to alignment
void projectA_s_gwfa_path_to_alignment(projectA_hash_graph_t* graph, s_gwfa_path_t* path, projectA_alignment_t* alignment, int32_t score) {

    // Copy alignment size and score
    alignment->size = path->size;
    alignment->score = score;

    // Reserve memory for alignmnet vectors
    alignment->nodes.reserve(alignment->size);
    alignment->cigar.reserve(alignment->size);

    // Iterate over the path
    for (int i = 0; i < path->size; ++i) {
        auto& node = path->nodes[i];

        // Add node to the nodes vector
        alignment->nodes.push_back(graph->nodes_in_order[node->id]->id);

        // Add sequence to referece sequence
        alignment->reference = alignment->reference + graph->nodes_in_order[node->id]->seq;
    }
}


// PRE:     graphs
//      graphs:         Reference of a vector of projectA_algorithm_input_t that holds the relveant information
//                      to create the gssw structs needed for alignment.
// POST:    return
//      return:         Void pointer that holds the populated structs needed for alignment by gwfa including reserved space for the results.
void* projectA_s_gwfa_init(vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    // Create io vectors for the s_gwfa algorithm
    projectA_s_gwfa_io_t* out = new projectA_s_gwfa_io_t;
    for (int i = 0; i < numThreads; ++i) {
        // Prepare vectors
        vector<projectA_s_gwfa_parameters_t> params;
        vector<projectA_s_gwfa_path_t> paths;

        // Reserve memory
        params.reserve(alignments.size()/numThreads);
        paths.reserve(alignments.size()/numThreads);

        // Push vectors to I/O struct
        out->parameters.push_back(params);
        out->paths.push_back(paths);
    }

    // Assign size
    out->size =  alignments.size();

    // Iterate over the input graphs
    int32_t thread_index = 0;
    for (auto& itr : alignments) {
        // Construct gwfa graph
        s_gwfa_graph_t* new_s_gwfa_graph = projectA_s_gwfa_graph_build(itr->graph);

        // Construct parameter entry
        projectA_s_gwfa_parameters_t entry(new_s_gwfa_graph, itr->read.c_str(), itr->read.size(), itr->graph);

        // Append entry to parameter vector
        out->parameters[thread_index].push_back(entry);

        // Update thread index
        thread_index = (thread_index == numThreads - 1) ? 0 : thread_index + 1;
    }

    // Cast return value int void pointer
    return static_cast<void*>(out);
}


// PRE:     ptr
//      ptr:            Void pointer that holds the populated structs needed for alignment by s_gwfa including reserved space for the results.
// POST:    return
//      return:         Void pointer that holds the results of the s_gwfa alignment as well as the gwfa structs needed for alignment.
void* projectA_s_gwfa_calculate_batch(void* ptr, int32_t thread_index){

    // Check if ptr is valid
    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_s_gwfa_io_t*>(ptr);
    auto& parameters = input->parameters[thread_index];
    auto& paths = input->paths[thread_index];

    // Run gwfa on every set of parameters
    for (int32_t i = 0; i < parameters.size(); ++i) {
        // Create local variables
        auto& parameter = parameters[i];
        projectA_s_gwfa_path_t path;
        int32_t score;
        s_gwfa_path_t** final_path = (s_gwfa_path_t**)malloc(sizeof(s_gwfa_path_t*));

        // Run gwfa
        // score = g_wfa_ed_infix(parameter.seq, parameter.len, parameter.graph, final_path);
        // score = g_wfa_ed(parameter.seq, parameter.len, parameter.graph, final_path);
        void* km = km_init();
        k_gwfa_path_t path_k_gwfa;
        // score = k_gwfa_ed(km, parameter.graph, parameter.len, parameter.seq, 0, &path_k_gwfa);
        score = k_gwfa_infix_ed(km, parameter.graph, parameter.len, parameter.seq, 0, &path_k_gwfa);
        km_destroy(km);

        path.score = score;
        path.path = *final_path;

        free(final_path);
        paths.push_back(path);
    }
    
    // Cast return value into void pointer
    return static_cast<void*>(ptr);
}


// PRE:     ptr, alignments
//      ptr:            Void pointer that holds the results of the s_gwfa alignment as well as the gwfa structs needed for alignment.
//      alignments:     Reference to a vector of pointers to projectA alignment structs.
// POST:    alignments
//      alignments:     Reference to a vector of pointers to projectA alignment structs that hold the information about the performed alignment.   
void projectA_s_gwfa_post(void* ptr, vector<projectA_alignment_t*>& alignments, int32_t numThreads) {

    // Check if ptr is valid
    if (ptr == nullptr) {
        cerr << "Error: input is nullptr!\n";
        exit(1);
    }

    // Define and cast local variables
    auto input = static_cast<projectA_s_gwfa_io_t*>(ptr);
    auto& parameters_vec = input->parameters;
    auto& paths_vec = input->paths;

    // Iterate over all entries in the paths
    // int32_t thread_index = 0;
    // int32_t j = 0;
    // for (int32_t i = 0; i < input->size; ++i) {
    //     auto& paths = paths_vec[thread_index];
    //     auto& parameters = parameters_vec[thread_index];

    //     // Check for out of bounds error 
    //     if (j >= paths.size() || j >= parameters.size()) {
    //         cerr << "Error: The index is out of bound!\n";
    //         cerr << "index: " << j << "\tin vector of the thread: " << thread_index << endl;
    //         cerr << "at value: " << i << endl;
    //         exit(1);
    //     }

    //     // Push alignment to alignments vector
    //     projectA_s_gwfa_path_to_alignment(parameters[j].projectA_hash_graph, paths[j].path, alignments[i], paths[j].score);

    //     // Update thread index
    //     if (thread_index == numThreads - 1) {
    //         thread_index = 0;
    //         ++j;
    //     } else {
    //         ++thread_index;
    //     }

    // }


    for (auto& parameters : parameters_vec) {
        for (auto& parameter : parameters) {
            s_gwfa_graph_destroy(parameter.graph);
        }
    }
    // for (auto& paths : paths_vec) {
    //     for (auto& path : paths) {
    //         s_gwfa_path_destroy(path.path);
    //     }
    // }

    delete input;
    ptr = nullptr;
}


// PRE:     
// POST:    return
//      return:         Pointer to a projectA_algorithm_t struct that holds the function pointers for s_gwfa.
projectA_algorithm_t* projectA_get_s_gwfa() {
    // Create new object
    projectA_algorithm_t* s_gwfa = new projectA_algorithm_t;

    // Assign function pointers
    s_gwfa->init = projectA_s_gwfa_init;
    s_gwfa->calculate_batch = projectA_s_gwfa_calculate_batch;
    s_gwfa->post = projectA_s_gwfa_post;

    return s_gwfa;
}


// PRE:     s_gwfa
//      s_gwfa:       Pointer to an existing s_gwfa project.
// POST:    s_gwfa
//      s_gwfa:       Nullpointer and gssw has been destroyed.
void projectA_s_gwfa_destroy(projectA_algorithm_t* s_gwfa) {

    // Check if the struct still exists
    if (s_gwfa == nullptr) {
        cerr << "Warning: s_gwfa struct is null pointer!\n";
        return;
    }

    // Delete struct
    delete s_gwfa;
}
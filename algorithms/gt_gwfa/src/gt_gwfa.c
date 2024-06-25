#include "graphs.h"

// gwf utility functions:
void gt_gwf_free(gwf_graph_t *g)
{
	int32_t i;
	for (i = 0; i < g->n_vtx; ++i) free(g->seq[i]);
	free(g->len); free(g->seq); free(g->arc); free(g->src);
    free(g->gssw_id);
    free(g);
}

void gt_gwf_graph_print(FILE *fp, const gwf_graph_t *g)
{
	int32_t i;
	for (i = 0; i < g->n_vtx; ++i)
		fprintf(fp, "S\t%d\t%s\tLN:i:%d\n", i, g->seq[i], g->len[i]);
	for (i = 0; i < g->n_arc; ++i)
		fprintf(fp, "L\t%d\t+\t%d\t+\t%dM\n", (uint32_t)(g->arc[i].a>>32), (uint32_t)g->arc[i].a, g->arc[i].o);
}





// create gssw grpah for testing
gssw_graph* gt_import_gssw_graph(int argc, char * const argv[]) {


    if (argc != 6) {
        fprintf(stderr, "usage: gssw_example nodeseq1 nodeseq2 nodeseq3 nodeseq4 readseq\n");
        exit(1);
    }
    
    // default parameters for genome sequence alignment
    int8_t match = 1, mismatch = 4;
    uint8_t gap_open = 6, gap_extension = 1;
    // from Mengyao's example about the importance of using all three matrices in traceback.
    // int32_t l, m, k, match = 2, mismatch = 1, gap_open = 2, gap_extension = 1;

    char *ref_seq_1 = argv[1];
    char *ref_seq_2 = argv[2];
    char *ref_seq_3 = argv[3];
    char *ref_seq_4 = argv[4];
    char *read_seq = argv[5];

    gssw_sse2_disable();
	/* This table is used to transform nucleotide letters into numbers. */
    int8_t* nt_table = gssw_create_nt_table();
    
	// initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code)
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)
    int8_t* mat = gssw_create_score_matrix(match, mismatch);

    gssw_node* nodes[4];
    nodes[0] = (gssw_node*)gssw_node_create(NULL, 1, ref_seq_1, nt_table, mat);
    nodes[1] = (gssw_node*)gssw_node_create(NULL, 2, ref_seq_2, nt_table, mat);
    nodes[2] = (gssw_node*)gssw_node_create(NULL, 3, ref_seq_3, nt_table, mat);
    nodes[3] = (gssw_node*)gssw_node_create(NULL, 4, ref_seq_4, nt_table, mat);
    
    // makes a diamond
    gssw_nodes_add_edge(nodes[0], nodes[1]);
    gssw_nodes_add_edge(nodes[0], nodes[2]);
    gssw_nodes_add_edge(nodes[1], nodes[3]);
    gssw_nodes_add_edge(nodes[2], nodes[3]);
    
    gssw_graph* graph = gssw_graph_create(4);

    gssw_graph_add_node(graph, nodes[0]);
    gssw_graph_add_node(graph, nodes[1]);
    gssw_graph_add_node(graph, nodes[2]);
    gssw_graph_add_node(graph, nodes[3]);

    gssw_graph_fill(graph, read_seq, nt_table, mat, gap_open, gap_extension, 0, 0, 15, 2, true);

    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    read_seq,
                                                    strlen(read_seq),
                                                    nt_table,
                                                    mat,
                                                    gap_open,
                                                    gap_extension,
                                                    0, 0);

    printf("Optimal local mapping:\n");
    gssw_print_graph_mapping(gm, stdout);
    // gt_print_graph_mapping(gm);
    gssw_graph_mapping_destroy(gm);

    free(nt_table);
	free(mat);
    return graph;
}


void gt_print_gssw_graph(gssw_graph* graph) {
    int size = graph->size;
    gssw_node* max_node = graph->max_node;
    gssw_node** nodes = graph->nodes;
    gssw_node* curr;
    gssw_node* next;
    
    fprintf(stderr, "gssw_graph:\n");

    for (int i = 0; i < size; i++) {
        curr = nodes[i];
        fprintf(stderr, "S\t%i\t%s\tLN:i:%i\n", curr->gwfa_index, curr->seq, curr->len);
    }

    for (int i = 0; i < size; i++) {
        curr = nodes[i];
        for (int j = 0; j < curr->count_next; j++) {
            next = curr->next[j];
            fprintf(stderr, "L\t%i\t+\t%i\t+\t0M\n", curr->gwfa_index, next->gwfa_index);
        }
    }
    return;
}


void gt_check_graphs(gssw_graph* gssw, gwf_graph_t* gwf) {
    if (gssw->size != gwf->n_vtx) {
        fprintf(stderr, "size not equal\n");
        exit(1);
    }
    int size = gssw->size;
    int len;
    gssw_node** gssw_nodes = gssw->nodes;
    gssw_node* gssw_curr;

    for (int i = 0; i < size; i++) {
        gssw_curr = gssw_nodes[i];
        if (gssw_curr->len != gwf->len[i]) {
            fprintf(stderr, "seq len not equal in seq%i:\n", i);
            gt_gwf_graph_print(stderr, gwf);
            gt_print_gssw_graph(gssw);
            exit(1);
        }
        len = gssw_curr->len;
        for (int j = 0; j < len; j++) {
            if (gssw_curr->seq[j] != gwf->seq[i][j]) {
                fprintf(stderr, "seq not equal\tnode: %i\tposition: %i\n", i, j);
                gt_gwf_graph_print(stderr, gwf);
                gt_print_gssw_graph(gssw);
                exit(1);
            }
        }
    }
    return;
}

// utility function to compare two cigar elements
void gt_check_cigar_element(gssw_cigar_element ce1, gssw_cigar_element ce2) {

    if (ce1.length != ce2.length) {
        fprintf(stderr, "cigar element length is not equal c1:\t%i\tc2:\t%i\n", ce1.length, ce2.length);
    }

    if (ce1.type != ce2.type) {
        fprintf(stderr, "cigar element is not equal ce1:\t%c\tce2:\t%c\n", ce1.type, ce2.type);
    }
}

// utility function to compare two cigars
void gt_check_cigar(gssw_cigar* c1, gssw_cigar* c2) {

    if (c1->length != c2->length) {
        fprintf(stderr, "cigar length is not equal c1:\t%i\tc2:\t%i\n", c1->length, c2->length);
    }

    int32_t length = c1->length;

    for (int i = 0; i < length; i++) {
        gt_check_cigar_element(c1->elements[i], c2->elements[i]);
    }
}

// utility function to compare two node cigars
void gt_check_node_cigar(gssw_node_cigar nc1, gssw_node_cigar nc2) {

    if (nc1.node != nc2.node) {
        fprintf(stderr, "node cigar nodes not equal nc1:\t%p\tnc2:\t%p\n", nc1.node, nc2.node);
    }

    gt_check_cigar(nc1.cigar, nc2.cigar);
}

// utility function to compare two graph cigars
void gt_check_graph_cigar (gssw_graph_cigar gc1, gssw_graph_cigar gc2) {

    if (gc1.length != gc2.length) {
        fprintf(stderr, "graph cigar length is not equal gm1:\t%u\tgm2:\t%i\n", gc1.length, gc2.length);
    }

    uint32_t length = gc1.length;

    for (int i = 0; i < length; i++) {
        fprintf(stderr, "%i\t", i);
        gt_check_node_cigar(gc1.elements[i], gc2.elements[i]);
    }
}

// utility function to compare two graph mappings
void gt_check_alignment(gssw_graph_mapping gm1, gssw_graph_mapping gm2) {

    if (gm1.position != gm2.position) {
        fprintf(stderr, "mapping position is not equal gm1:\t%i\tgm2:\t%i\n", gm1.position, gm2.position);
    }

    // if (gm1.score != gm2.score) {
    //     fprintf(stderr, "mapping score is not equal gm1:\t%i\tgm2:\t%i\n", gm1.score, gm2.score);
    //     exit(1);
    // }

    gt_check_graph_cigar(gm1.cigar, gm2.cigar);
    fprintf(stderr, "\n\n");
}


void gt_check_node_length(gssw_graph_mapping* gm) {
    gssw_node_cigar* nc = gm->cigar.elements;
    gssw_cigar* c;
    int node_length, cigar_length;
    for (int i = 0; i < gm->cigar.length; i++) {
        cigar_length = 0;
        node_length = strlen(nc[i].node->seq);
        c = nc[i].cigar;
        for (int j = 0; j < c->length; j++) {
            if (c->elements[j].type == 'M' || c->elements[j].type == 'I') {
                cigar_length += c->elements[j].length;
            }
        }
        if (i == 0) {
            node_length -= gm->position;
        }
        if (cigar_length > node_length) {
            fprintf(stderr, "cigar and node are not the same length\nid:\t%i\nnode:\t%i\tcigar:\t%i\n",
                            nc[i].node->id, node_length, cigar_length);
        }
    }
}

void gt_check_alignment_length(gssw_graph_mapping* gm, int32_t expected_length) {
    gssw_node_cigar* nc = gm->cigar.elements;
    gssw_cigar* c;
    int alignment_length = 0;
    for (int i = 0; i < gm->cigar.length; i++) {
        c = nc[i].cigar;
        for (int j = 0; j < c->length; j++) {
            if (c->elements[j].type == 'M' || c->elements[j].type == 'I') {
                alignment_length += c->elements[j].length;
            }
        }
    }
    if (alignment_length != expected_length) {
        fprintf(stderr, "alignment is not as long as expected\naligned bases:\t%i\texpected:\t%i\n",
                        alignment_length, expected_length);
        gssw_print_graph_mapping(gm, stderr);
    }
}

void gt_gwfa_print_to_file(FILE *fp, const gwf_graph_t *g)
{
	int32_t i;
    fprintf(fp, "NG:\n");
	for (i = 0; i < g->n_vtx; ++i)
		fprintf(fp, "S\t%d\t%lu\t%s\tLN:i:%d\n", i, g->gssw_id[i], g->seq[i], g->len[i]);
	for (i = 0; i < g->n_arc; ++i)
		fprintf(fp, "L\t%d\t+\t%d\t+\t%dM\n", (uint32_t)(g->arc[i].a>>32), (uint32_t)g->arc[i].a, g->arc[i].o);
}


gssw_graph* gt_read_gfa(const char* filename) {
    int8_t match = 1, mismatch = 4;
    uint8_t gap_open = 6, gap_extension = 1;
    gssw_sse2_disable();
    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(match, mismatch);
    uint32_t size = 0;

    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file!");
        return;
    }

    char line[2048];
    while (fgets(line, sizeof(line), file) != NULL) {
        char* newline = strchr(line, '\n');
        if (newline != NULL) *newline = '\0';
        char* token = strtok(line, " \t\n");
        while (token != NULL) {
            if (token[0] == 'S' && strlen(token) == 1) {
                size++;
            }
            token = strtok(NULL, " \t\n");
        }
    }
    fclose(file);
    fprintf(stderr, "%i\n", size);


    gssw_node* nodes[size];

    FILE* file2 = fopen(filename, "r");
    if (file2 == NULL) {
        fprintf(stderr, "Error opening file!");
        return;
    }

    uint32_t n_nodes = 0;
    while (fgets(line, sizeof(line), file2) != NULL) {
        char* newline = strchr(line, '\n');
        if (newline != NULL) *newline = '\0';

        char* token = strtok(line, " \t\n");
        while (token != NULL) {
            // printf("Token: %s\n", token);

            if (token[0] == 'S' && strlen(token) == 1) {

                char* id_str = strtok(NULL, " \t\n");
                int32_t id = atoi(id_str);
                char* seq = strtok(NULL, " \t\n");

                if (n_nodes >= size) {
                    fprintf(stderr, "Error: to many nodes!\n");
                    return NULL;
                }

                nodes[n_nodes] = (gssw_node*)gssw_node_create(NULL, id, seq, nt_table, mat);
                nodes[n_nodes]->gwfa_index = id;
                n_nodes++;
            }

            if (token[0] == 'L' && strlen(token) == 1) {


                char* start_str = strtok(NULL, " \t\n");
                int32_t start = atoi(start_str);
                char* dir_start = strtok(NULL, " \t\n");
                char* end_str = strtok(NULL, " \t\n");
                int32_t end = atoi(end_str);
                char* dir_end = strtok(NULL, " \t\n");

                // fprintf(stderr, "%i\t%i\n", start, end);
                gssw_nodes_add_edge(nodes[start], nodes[end]);
            } 

            token = strtok(NULL, " \t\n");
        }
    }
    fclose(file2);
    

    gssw_graph* graph = gssw_graph_create(size);

    for (int i = 0; i < size; ++i) {
        gssw_graph_add_node(graph, nodes[i]);
    }

    free(mat);
    free(nt_table);
    return graph;
}






// change gssw graph to gwf graph
gwf_graph_t* gt_gssw_to_gwf(gssw_graph* input_graph) {
    gwf_graph_t* out;
    GFA_CALLOC(out, 1);

    uint32_t i, k, j;

    uint32_t size;
    gssw_node* current_node;
    uint64_t current_id;
    char* current_seq;
    int32_t current_len;
    int32_t current_count_next;
    uint64_t current_edge_source_id;
    uint64_t current_edge_destination_id;

    size = input_graph->size;
    int32_t total_edges = 0;
    for (uint32_t i = 0; i < size; i++) {
        current_node = input_graph->nodes[i];
        total_edges += current_node->count_next;
        current_node->gwfa_index = i;
    }
    out->n_vtx = size;
    out->n_arc = total_edges;
    // out->gssw_node_id = current_node->id;
    GFA_MALLOC(out->len, out->n_vtx);
	GFA_MALLOC(out->src, out->n_vtx);
	GFA_MALLOC(out->seq, out->n_vtx);
	GFA_MALLOC(out->arc, out->n_arc);
    GFA_MALLOC(out->gssw_id, out->n_vtx);

    //iterating over all nodes to store the node information
    for (i = k = 0; i < out->n_vtx; i++) {
        current_node = input_graph->nodes[i];
        uint32_t v;
        uint32_t current_len =  current_node->len;
        uint32_t current_count_next;
        uint32_t current_edge_destination_id;
        // fprintf(stderr, "%u\n", current_node->gwfa_index);
        out->len[i] = current_len;
        current_count_next = current_node->count_next;

        out->gssw_id[i] = current_node->id;

        GFA_MALLOC(out->seq[i], current_len + 1);

        for (j = 0; j < current_len; j++) {
            out->seq[i][j] = current_node->seq[j];
        }
        out->seq[i][current_len] = 0;

        for (j = 0; j < current_count_next; j++, k++) {
            current_edge_destination_id = current_node->next[j]->gwfa_index;
            out->arc[k].a = (uint64_t)i<<32 | current_edge_destination_id;
            out->arc[k].o = 0;
            //fprintf(stderr, "\n arcs: \n%lu\t%u\n", out->arc[k].a, out->arc[k].o);
        }
    }

    // gt_check_graphs(input_graph, out);
    return out;
}



// propper handling of graph transition
gwf_graph_t* gt_gssw_to_gwf_direct(gssw_graph* in_graph) {
    gwf_graph_t* out_graph;

    out_graph = gt_gssw_to_gwf(in_graph);

    return out_graph;
}




// flatten CIGAR-string
char* gt_flatten_cigar(const char *cigarString) {
    char* flattened_cigar = (char*)malloc(1); // Allocate memory for an empty string
    flattened_cigar[0] = '\0'; // Ensure the string is properly null-terminated
    int32_t number, length, current_length;
    int32_t i, j;
    char operation;
    length = strlen(cigarString);
    i = 0;
    number = 0;
    
    while (i < length) {
        if (isdigit(cigarString[i])) {
            while (isdigit(cigarString[i])) {
                number = number * 10 + (cigarString[i] - '0');
                i++;
            }
        } else {
            operation = cigarString[i];
            current_length = strlen(flattened_cigar);
            
            // Expand memory to accommodate the flattened_cigar string
            flattened_cigar = (char*)realloc(flattened_cigar, strlen(flattened_cigar) + number + 1);
            
            // Append the current operation number times to the flattened_cigar string
            for (j = 0; j < number; j++) {
                flattened_cigar[current_length + j] = operation;
            }
            flattened_cigar[current_length + number] = '\0'; // Add null terminator
            
            // Reset number for the next pair
            number = 0;
            i++;
        }
    }
    
    flattened_cigar[strlen(flattened_cigar)] = '\0';
    return flattened_cigar;
}





// Append CIGAR operation while rebuilding
char* gt_append_re_cigar(int32_t count, char last_operation, char* cigar) {
    int32_t num_len;

    // Reallocate memory to account for longer string
    char countStr[20];
    sprintf(countStr, "%d", count);
    cigar = (char*)realloc(cigar, strlen(cigar) + strlen(countStr) + 1);

    // Concatenate the new additons to the CIGAR
    strcat(cigar, countStr);
    char lastOpStr[2]; // Temporary string to hold the last_operation
    lastOpStr[0] = last_operation;
    lastOpStr[1] = '\0';
    strcat(cigar, lastOpStr);


    return cigar;
}





// Rebuild CIGAR from flattened
char* gt_rebuild_cigar(char* flattened_cigar) {

    // Initialies variables
    int32_t length, count;
    int32_t i;
    char operation, last_operation;
    char* cigar;

    // Allocate memeory and assign default values
    cigar = (char*)malloc(1);
    cigar[0] = '\0';
    length = strlen(flattened_cigar);
    last_operation = '\0';
    count = 0;

    // Iterate over the flattened CIGAR to rebuild the CIGAR
    for (i = 0; i < length; i++) {
        operation = flattened_cigar[i];

        // If the last operation is the same as the current operation we can increase the count by one
        if (operation == last_operation) {
            count++;

        // If the last operation isn't the same as the current operation we need to save the current operation in the CIGAR
        // and start a new count for the next operation
        } else {

            // If this is the first iteration we don't need to safe the last operation but only need to start a new count
            if (last_operation == '\0') {
                count = 1;
                last_operation = operation;

            // Save last operation and start new count
            } else {

                // Append operation to CIGAR
                cigar = gt_append_re_cigar(count, last_operation, cigar);

                // Reset count and operation
                count = 1;
                last_operation = operation;
            }
        }
    }

    // Append last operation of the CIGAR
    if (count > 0) {
        cigar = gt_append_re_cigar(count, last_operation, cigar);
    }

    cigar[strlen(cigar)] = '\0';
    return cigar;
}


// Function to create a CIGAR string from gssw_cigar
char* gt_gssw_create_cigar_string(gssw_cigar* c) {
    char* cigar = (char*)malloc(1);  // Allocate memory for an empty string
    cigar[0] = '\0';  // Initialize the string

    int i;
    int l = c->length;
    gssw_cigar_element* e = c->elements;

    for (i = 0; i < l; ++i, ++e) {
        char length_str[20];  // Assuming a reasonable maximum length for an integer
        sprintf(length_str, "%d", e->length);
        
        // Allocate memory for the concatenated string
        cigar = (char*)realloc(cigar, strlen(cigar) + strlen(length_str) + 2);  // Explicit cast to char*

        strcat(cigar, length_str);
        strncat(cigar, &e->type, 1);
    }

    return cigar;
}


// calculate offset in a cigar string
int gt_calculate_offset(const char *cigar_string) {
    int offset = 0;
    int current_length = 0;

    while (*cigar_string != '\0') {
        if (isdigit(*cigar_string)) {
            current_length = current_length * 10 + (*cigar_string - '0');
        } else {
            if (*cigar_string == 'M' || *cigar_string == 'D' || *cigar_string == 'N') {
                offset += current_length;
            } else if (*cigar_string == 'I' || *cigar_string == 'S') {
                // Insertions and soft clippings in the read do not contribute to the offset
            }

            // Reset the current length for the next operation
            current_length = 0;
        }

        cigar_string++;
    }

    return offset;
}



// Rebuild CIGAR from flattened and add to graph mapping
int gt_flattened_cigar_to_gm(char* flattened_cigar, gssw_node_cigar* nc, int32_t current_index) {

    // Initialies variables
    int32_t length, count;
    int32_t i, j;
    char operation, last_operation;

    // Allocate memeory and assign default values
    length = strlen(flattened_cigar);
    last_operation = '\0';
    count = 0;
    j = 0;

    // Iterate over the flattened CIGAR to rebuild the CIGAR
    for (i = 0; i < length; i++) {
        operation = flattened_cigar[i];

        // If the last operation is the same as the current operation we can increase the count by one
        if (operation == last_operation) {
            count++;

        // If the last operation isn't the same as the current operation we need to save the current operation in the CIGAR
        // and start a new count for the next operation
        } else {

            // If this is the first iteration we don't need to safe the last operation but only need to start a new count
            if (last_operation == '\0') {
                count = 1;
                last_operation = operation;

            // Save last operation and start new count
            } else {
                
                // Append operation to CIGAR
                gssw_cigar_push_back(nc->cigar, last_operation, count);
                // fprintf(stderr, "%i%c", count, last_operation);

                // Reset count and operation
                count = 1;
                last_operation = operation;
                j++;
            }
        }
    }

    // Append last operation of the CIGAR
    if (count > 0) {
        gssw_cigar_push_back(nc->cigar, last_operation, count);
        // fprintf(stderr, "%i%c", count, last_operation);
    }
    // fprintf(stderr, "\n");

    return j;
}

void gt_print_graph_mapping(gssw_graph_mapping* gm) {
    fprintf(stderr, "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
    fprintf(stderr, "position: %i\tscore: %i\tlength: %u\n", gm->position, gm->score, gm->cigar.length);
    gssw_node_cigar* e = gm->cigar.elements;
    for (int i = 0; i < gm->cigar.length; i++) {
        gssw_cigar* c = e[i].cigar;
        fprintf(stderr, "node id: %i\t node index: %i\t cigar: ", e[i].node->id, e[i].node->gwfa_index);
        gssw_print_cigar(c, stderr);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
}

// Create transfer CIGAR from gwfa to gssw
gssw_graph_mapping* gt_gwf_traceback_to_gm(gwf_path_t *path, gwf_graph_t *graph, gssw_graph *gssw_g, const char* read) {
    int32_t i, j, k, l;
    int32_t count, index_shift, position, position_buffer;
    char* flattened_cigar = NULL;
    char* aligned_seq = (char*)malloc(1);
    char* curr_cigar_flattened = (char*)malloc(1);
    char* substring = NULL;
    aligned_seq[0] = '\0';
    curr_cigar_flattened[0] = '\0';

    gssw_graph_mapping* gm = gssw_graph_mapping_create();
    gssw_graph_cigar* gc = &gm->cigar;
    gm->cigar.length = path->nv;
    gssw_node_cigar* nc = gc->elements;

    gm->score = path->s;

    for (i = 0; i < path->nv; i++) {
        int32_t curr_node = path->v[i];
        char* curr_seq = graph->seq[curr_node];
        aligned_seq = (char*)realloc(aligned_seq, strlen(aligned_seq) + strlen(curr_seq) + 1);
        strcat(aligned_seq, curr_seq);
    }

    EdlibAlignResult result = edlibAlign(read, strlen(read), aligned_seq, strlen(aligned_seq), edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
    char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    gm->position = result.startLocations[0];
    position = result.startLocations[0];

    flattened_cigar = gt_flatten_cigar(cigar);

    uint32_t graph_cigar_bufsize = 16;
    gc->elements = (gssw_node_cigar*)malloc(graph_cigar_bufsize * sizeof(gssw_node_cigar));
    gc->length = 0;

    for (i = 0, index_shift = 0; i < path->nv; i++) {
        count = 0;
        j = 0;
        position_buffer = 0;
        l = path->v[i];

        if (position > 0) {
            if (position <= graph->len[l]) {
                position_buffer = position;
                position = 0;
            } else {
                position -= graph->len[l];
                gm->position -= graph->len[l];
                continue;
            }
        }

        if (gc->length == graph_cigar_bufsize) {
            graph_cigar_bufsize *= 2;
            gc->elements = (gssw_node_cigar*)realloc(gc->elements, graph_cigar_bufsize * sizeof(gssw_node_cigar));
        }

        nc = gc->elements + gc->length;

        while (count < graph->len[l] - position_buffer && flattened_cigar[index_shift + j] != '\0') {
            if (flattened_cigar[index_shift + j] == 'M' || flattened_cigar[index_shift + j] == 'I') {
                count++;
            }
            j++;
        }

        while (flattened_cigar[index_shift + j] != 'M' && flattened_cigar[index_shift + j] != 'I' && flattened_cigar[index_shift + j] != '\0') {
            j++;
        }

        if (index_shift < strlen(flattened_cigar)) {
            if (j < 0) {
                fprintf(stderr, "error: count per node < 0\n");
                exit(1);
            }

            if (j == 0 || position > 0) {
                gm->position -= graph->len[l];
            } else {
                substring = (char*)malloc(j + 1);

                if (index_shift + j > strlen(flattened_cigar)) {
                    fprintf(stderr, "error: out of bounds access\n");
                    exit(1);
                }

                for (k = 0; k < j; k++) {
                    substring[k] = flattened_cigar[index_shift + k];
                }
                substring[k] = '\0';

                nc->cigar = (gssw_cigar*)calloc(1, sizeof(gssw_cigar));
                gt_flattened_cigar_to_gm(substring, nc, l);
                nc->node = gssw_g->nodes[l];

                index_shift += j;
                gc->length++;
                free(substring);
                substring = NULL;
            }
        } else {
            break;
        }
    }

    free(flattened_cigar);
    free(cigar);
    free(curr_cigar_flattened);
    edlibFreeAlignResult(result);
    free(aligned_seq);
    return gm;
}




// hand gwf graph to gwfa for alignment
gssw_graph_mapping* gt_read_align_gwf (gwf_graph_t *graph, gssw_graph *gssw_g, const char* read) {
    void *km = 0;
    int32_t s;
    size_t read_len;
    uint32_t max_lag = 0;
    int traceback = 1;
    gwf_path_t path;
    km = km_init();

    read_len = strlen(read);
    gwf_ed_index(km, graph);

    s = gwf_ed(km, graph, read_len, read, 0, -1, max_lag, traceback, &path);

    gssw_graph_mapping* gm = gt_gwf_traceback_to_gm(&path, graph, gssw_g, read);
	gwf_cleanup(km, graph);
	km_destroy(km);

    // fprintf(stderr, "\n INTERNAL CIGAR:\n");
    // for (int i = 0; i < gm->cigar.length; i++) {
    //     fprintf(stderr, "%i:\t", i);
    //     gssw_print_cigar(gm->cigar.elements[i].cxcigar, stderr);
    //     fprintf(stderr, "\n");
    // }
    // gssw_print_graph_mapping(gm, stderr);
    // gt_print_graph_mapping(gm);
    return gm;
}
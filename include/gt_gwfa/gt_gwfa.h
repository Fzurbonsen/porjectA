#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <ctype.h>

#include "gwfa/gwfa.h"
#include "gwfa/gfa.h"
#include "gwfa/gfa-priv.h"
#include "gssw.h"
#include "edlib.h"

#ifdef __cplusplus
extern "C" {
#endif

void gt_gwf_free(gwf_graph_t *g);
void gt_gwf_graph_print(FILE *fp, const gwf_graph_t *g);

gssw_graph* gt_import_gssw_graph(int argc, char * const argv[]);
// void gt_print_gssw_graph(gssw_graph* graph, char* read_seq);
void gt_print_gssw_graph(gssw_graph* graph);
void gt_check_graphs(gssw_graph* gssw, gwf_graph_t* gwf);
void gt_check_alignment(gssw_graph_mapping gm1, gssw_graph_mapping gm2);
void gt_check_node_length(gssw_graph_mapping* gm);
void gt_check_alignment_length(gssw_graph_mapping* gm, int32_t expected_length);
gssw_graph* gt_read_gfa(const char* filename);

void gt_gwfa_print_to_file(FILE *fp, const gwf_graph_t *g);

void gt_print_graph_mapping (gssw_graph_mapping* gm);

gwf_graph_t* gt_gssw_to_gwf(gssw_graph* input_graph);
gwf_graph_t* gt_gssw_to_gwf_direct(gssw_graph* in_graph);

gssw_graph_mapping* gwf_traceback_to_gm(gwf_path_t *path, gwf_graph_t *graph, gssw_graph *gssw_g, const char* gfa_read);
gssw_graph_mapping* gt_read_align_gwf(gwf_graph_t *graph, gssw_graph *gssw_g, const char* read);

#ifdef __cplusplus
}
#endif
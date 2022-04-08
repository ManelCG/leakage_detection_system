#include <graph.h>
#include <stdlib.h>

#include <stdio.h>
#include <stdbool.h>

int main(int argc, char *argv[]){
  int   sorig[] = {1, 2, 2,  2, 12, 3, 3, 6, 6,  9,  9, 12};
  int   torig[] = {2, 6, 9, 12,  3, 4, 5, 7, 8, 10, 11, 13};
  float dorig[] = {0.15, 0.15, 0.15, 0.15, 0.16, 0.075, 0.075, 0.075, 0.1, 0.1, 0.075, 0.075};
  float rough[] = {0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015};
  float outfw[] = {3.6, 2.7, 0, 23.1, 0.5, 1.8, 3.6};
  // float outfw[] = {1, 1, 1, 1, 1, 1, 1};
  int num_pipes = 12;

  Graph *g = graph_new(NULL, num_pipes, sorig, torig);
  graph_set_diameters(g, dorig);
  graph_set_roughness(g, rough);

  int num = graph_get_n_output_nodes(g);
  printf("%d demand nodes\n", num);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_output_node(g, i);
    node_set_flowrate(n, outfw[i]);
    node_set_flowrate_measured(n, true);
  }

  graph_set_inflow_evenly(g);

  graph_print(g);
}

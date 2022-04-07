#include <graph.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
  int   sorig[] = {1, 2, 2,  2, 12, 3, 3, 6, 6,  9,  9, 12};
  int   torig[] = {2, 6, 9, 12,  3, 4, 5, 7, 8, 10, 11, 13};
  float dorig[] = {0.15, 0.15, 0.15, 0.15, 0.16, 0.075, 0.075, 0.075, 0.1, 0.1, 0.075, 0.075};
  float rough[] = {0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015, 0.0000015};
  int num_pipes = 12;

  Graph *g = graph_new(NULL, num_pipes, sorig, torig);
  graph_set_diameters(g, dorig);
  graph_set_roughness(g, rough);


  graph_print(g);
}

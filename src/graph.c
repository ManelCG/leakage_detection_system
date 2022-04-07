#include <graph.h>

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define PI 3.1415926536

typedef struct Pipe{
  Node *orig;
  Node *dest;

  float diam;
  float area;
  float rough;

  int ID;

} Pipe;

typedef struct Node{
  int n_pipes_in;
  int n_pipes_out;

  Pipe **pipes_in;
  Pipe **pipes_out;

  _Bool is_junction;
  _Bool is_connected;
  _Bool is_input;
  _Bool is_output;

  int ID;
} Node;

typedef struct Graph{
  Node **nodes;
  Pipe **pipes;

  int *inc_matrix;

  int n_pipes;
  int n_nodes;
} Graph;


//Constructors
Pipe *pipe_new(Pipe **ret, Node *orig, Node *dest){
  Pipe *new = malloc(sizeof(Pipe));

  new->orig = orig;
  new->dest = dest;

  pipe_set_diam(new, 0);
  pipe_set_rough(new, 0);

  if (ret != NULL){
    *ret = new;
  }
  return new;
}
Node *node_new(Node **ret){
  Node *new = malloc(sizeof(Node));

  new->n_pipes_in = 0;
  new->n_pipes_out = 0;
  new->pipes_in = NULL;
  new->pipes_out = NULL;

  new->is_junction = false;
  new->is_connected = false;
  new->is_input = false;
  new->is_output = false;

  if (ret != NULL){
    *ret = new;
  }
  return new;
}
Graph *graph_new(Graph **ret, int n_pipes, int *sorig, int *torig){
  Graph *g = malloc(sizeof(Graph));

  //Get num of nodes
  int n_nodes = 0;
  for (int i = 0; i < n_pipes; i++){
    n_nodes = fmax(n_nodes, sorig[i]);
    n_nodes = fmax(n_nodes, torig[i]);
  }
  n_nodes++;  //If biggest node is 13 -> There are 14 nodes (counting 0)

  g->n_nodes = n_nodes;
  g->n_pipes = n_pipes;

  //Alloc memory for arrays
  g->nodes = malloc(sizeof(Node *) * n_nodes);
  g->pipes = malloc(sizeof(Pipe *) * n_pipes);

  g->inc_matrix = calloc(sizeof(int) * n_nodes * n_pipes, 1);

  //Create nodes
  for (int i = 0; i < n_nodes; i++){
    g->nodes[i] = node_new(NULL);
    node_set_id(g->nodes[i], i);
  }

  for (int i = 0; i < n_pipes; i++){  //Create pipes.
    g->pipes[i] = pipe_new(NULL, g->nodes[sorig[i]], g->nodes[torig[i]]);
    pipe_set_id(g->pipes[i], i);

    node_add_pipe_out(g->nodes[sorig[i]], g->pipes[i]);
    node_add_pipe_in(g->nodes[torig[i]], g->pipes[i]);

    g->inc_matrix[i*n_nodes + sorig[i]] = -1;
    g->inc_matrix[i*n_nodes + torig[i]] = 1;
  }

  if (ret != NULL){
    *ret = g;
  }
  return g;
}

//Node functions
void node_add_pipe_in(Node *n, Pipe *p){
  if (n->pipes_in == NULL){
    n->pipes_in = malloc(sizeof(Pipe*));
  } else {
    n->pipes_in = realloc(n->pipes_in, sizeof(Pipe*) * n->n_pipes_in + 1);
  }

  n->pipes_in[n->n_pipes_in] = p;
  n->n_pipes_in += 1;

  //Update node status
  n->is_connected = true;
  if (n->pipes_out > 0){
    n->is_junction = true;
    n->is_input = false;
    n->is_output = false;
  } else {
    n->is_junction = false;
    n->is_input = false;
    n->is_output = true;
  }
}
void node_add_pipe_out(Node *n, Pipe *p){
  if (n->pipes_out == NULL){
    n->pipes_out = malloc(sizeof(Pipe*));
  } else {
    n->pipes_out = realloc(n->pipes_out, sizeof(Pipe*) * n->n_pipes_out + 1);
  }

  n->pipes_out[n->n_pipes_out] = p;
  n->n_pipes_out += 1;

  //Update node status
  n->is_connected = true;
  if (n->pipes_in > 0){
    n->is_junction = true;
    n->is_input = false;
    n->is_output = false;
  } else {
    n->is_junction = false;
    n->is_input = true;
    n->is_output = false;
  }
}
void node_set_id(Node *n, int id){
  n->ID = id;
}

//Pipe functions
void pipe_set_diam(Pipe *p, float d){
  p->diam = d;
  p->area = 1.0/4 * pow(PI, 2) * d;
}
void pipe_set_id(Pipe *p, int id){
  p->ID = id;
}
void pipe_set_rough(Pipe *p, float r){
  p->rough = r;
}


//Graph functions
void graph_set_diameters(Graph *g, float *d){
  for (int i = 0; i < g->n_pipes; i++){
    pipe_set_diam(g->pipes[i], d[i]);
  }
}
void graph_set_roughness(Graph *g, float *r){
  for (int i = 0; i < g->n_pipes; i++){
    pipe_set_rough(g->pipes[i], r[i]);
  }
}
void graph_print(Graph *g){
  printf("Printing graph: \n");
  printf("PIPES: \n");
  for (int i = 0; i < g->n_pipes; i++){
    Pipe *p = g->pipes[i];

    printf("Pipe nº %d with ID %d:\n", i, p->ID);
    printf("Diam %g, area %g\n", p->diam, p->area);
    printf("Roughness: %.7f\n", p->rough);
    printf("%3d -> %3d\n", p->orig->ID, p->dest->ID);
    printf("\n");
  }

  printf("NODES: \n");
  for (int i = 0; i < g->n_nodes; i++){
    Node *n = g->nodes[i];

    if (n->is_connected){
      printf("Node nº %d with ID %d\n", i, n->ID);

      if (n->is_junction){
        printf("Junction node\n");
      }

      if (! n->is_input){
        printf("Fed by: ");
        for (int j = 0; j < n->n_pipes_in; j++){
          printf("%d ", n->pipes_in[j]->ID);
        }
      } else {
        printf("Input node");
      }
      printf("\n");

      if (! n->is_output){
        printf("Out to: ");
        for (int j = 0; j < n->n_pipes_out; j++){
          printf("%d ", n->pipes_out[j]->ID);
        }
      } else {
        printf("Output node");
      }
      printf("\n");
      printf("\n");
    }
  }

  printf("INCIDENCE MATRIX:\n");
  for (int j = 0; j < g->n_nodes; j++){
  for (int i = 0; i < g->n_pipes; i++){
    if (g->inc_matrix[i*g->n_nodes + j] < 0){
      printf("%d ", g->inc_matrix[i*g->n_nodes + j]);
    } else {
      printf(" %d ", g->inc_matrix[i*g->n_nodes + j]);
    }
  }
  printf("\n");
  }
  printf("\n");
}

#include <graph.h>

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define PI 3.1415926536

#define GEOMETRY_CIRCULAR 0
#define GEOMETRY_SQUARE 1
#define GEOMETRY_RECTANGULAR 2
#define GEOMETRY_CIRCULAR_ANNULUS 3
#define GEOMETRY_CUSTOM 4

union dimensions{
  float circ_diam;
  float squa_side;
  float rect_sides[2];
  float *custom_sides;
};

typedef struct Pipe{
  Node *orig;
  Node *dest;

  int geometry;
  union dimensions dimensions;
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

  float flowrate;             //m³/s. if -1 -> unknown
  _Bool flowrate_is_measured; //false -> calculated or unknown

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

  pipe_set_geometry(new, GEOMETRY_CIRCULAR);

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

  new->flowrate = -1;
  new->flowrate_is_measured = false;

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
void node_print(Node *n){
  printf("Node nº with ID %d\n", n->ID);

  if (n->is_connected){
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

    if (n->flowrate != -1){
      printf("Flowrate ");
      if (!n->flowrate_is_measured){
        printf("NOT ");
      }
      printf("measured: %g m³/s\n", n->flowrate);
    } else {
      printf("Flowrate unknown\n");
    }


  } else {    //Not connected
    printf("Disconnected node\n");
  }
}
void node_set_flowrate(Node *n, float f){
  n->flowrate = f;
}
float node_get_flowrate(Node *n){
  return n->flowrate;
}
void node_set_flowrate_measured(Node *n, _Bool m){
  n->flowrate_is_measured = true;
}

//Pipe functions
int pipe_set_diam(Pipe *p, float d){
  if (p->geometry == GEOMETRY_CIRCULAR ||
      p->geometry == GEOMETRY_CIRCULAR_ANNULUS){
    p->dimensions.circ_diam = d;
    p->area = 1.0/4 * pow(PI, 2) * d;
    return 0;
  }
  return -1;
}
int pipe_set_side(Pipe *p, float s){
  if (p->geometry == GEOMETRY_SQUARE){
    p->dimensions.squa_side = s;
    p->area = s*s;
    return 0;
  }
  return -1;
}
int pipe_set_sides(Pipe *p, float s1, float s2){
  if (p->geometry == GEOMETRY_RECTANGULAR){
    p->dimensions.rect_sides[0] = s1;
    p->dimensions.rect_sides[1] = s2;
    p->area = s1 * s2;
    return 0;
  }
  return -1;
}
void pipe_set_id(Pipe *p, int id){
  p->ID = id;
}
void pipe_set_rough(Pipe *p, float r){
  p->rough = r;
}
void pipe_set_geometry(Pipe *p, int g){
  p->geometry = g;
}
void pipe_print(Pipe *p){
  printf("Pipe nº with ID %d:\n", p->ID);
  switch(p->geometry){
    case GEOMETRY_CIRCULAR:
      printf("Diam %g, area %g\n", p->dimensions.circ_diam, p->area);
      break;
  }
  printf("Roughness: %.7f\n", p->rough);
  printf("%3d -> %3d\n", p->orig->ID, p->dest->ID);
  printf("\n");
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
void graph_print_incidence_matrix(Graph *g){
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
}
void graph_print(Graph *g){
  printf("---------------------------------------\n");
  printf("Printing graph: \n");
  printf("PIPES: \n");
  for (int i = 0; i < g->n_pipes; i++){
    pipe_print(g->pipes[i]);
  }

  printf("\n\nNODES: \n");
  for (int i = 0; i < g->n_nodes; i++){
    node_print(g->nodes[i]);
    printf("\n");
  }

  printf("\n\nINCIDENCE MATRIX:\n");
  graph_print_incidence_matrix(g);
  printf("---------------------------------------\n");
}
void graph_print_disconnected_nodes(Graph *g){
  for (int i = 0; i < g->n_nodes; i++){
    if (!g->nodes[i]->is_connected){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_connected_nodes(Graph *g){
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_connected){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_junction_nodes(Graph *g){
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_junction){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_input_nodes(Graph *g){
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_input){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_output_nodes(Graph *g){
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_output){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
int graph_get_n_disconnected_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (!g->nodes[i]->is_connected){
      sum++;
    }
  }
  return sum;
}
int graph_get_n_connected_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_connected){
      sum++;
    }
  }
  return sum;
}
int graph_get_n_junction_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_junction){
      sum++;
    }
  }
  return sum;
}
int graph_get_n_input_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_input){
      sum++;
    }
  }
  return sum;
}
int graph_get_n_output_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_output){
      sum++;
    }
  }
  return sum;
}
Node *graph_get_nth_disconnected_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (!g->nodes[i]->is_connected){
      if (index == sum){
        return g->nodes[i];
      }
      sum++;
    }
  }
  return NULL;
}
Node *graph_get_nth_connected_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_connected){
      if (index == sum){
        return g->nodes[i];
      }
      sum++;
    }
  }
  return NULL;
}
Node *graph_get_nth_junction_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_junction){
      if (index == sum){
        return g->nodes[i];
      }
      sum++;
    }
  }
  return NULL;
}
Node *graph_get_nth_input_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_input){
      if (index == sum){
        return g->nodes[i];
      }
      sum++;
    }
  }
  return NULL;
}
Node *graph_get_nth_output_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_output){
      if (index == sum){
        return g->nodes[i];
      }
      sum++;
    }
  }
  return NULL;
}

float graph_get_total_outflow(Graph *g){
  int num_outputs = graph_get_n_output_nodes(g);
  float sum = 0;
  for (int i = 0; i < num_outputs; i++){
    Node *n = graph_get_nth_output_node(g, i);
    sum += node_get_flowrate(n);
  }
  return sum;
}

void graph_set_inflow_evenly(Graph *g){
  int num_inputs = graph_get_n_input_nodes(g);
  float total_outflow = graph_get_total_outflow(g);

  float even_outflow = total_outflow / num_inputs;

  for (int i = 0; i < num_inputs; i++){
    Node *n = graph_get_nth_input_node(g, i);
    node_set_flowrate(n, even_outflow);
  }
}

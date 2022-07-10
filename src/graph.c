#include <graph.h>

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include <unistd.h>
#include <string.h>

#define PI 3.1415926536

#define GEOMETRY_CIRCULAR 0
#define GEOMETRY_SQUARE 1
#define GEOMETRY_RECTANGULAR 2
#define GEOMETRY_CIRCULAR_ANNULUS 3
#define GEOMETRY_CUSTOM 4

#define MAX_LEAK_OUTFLOW 0.01

// #define __GRAPH_C_DEBUG_
// #define __GRAPH_C_DETECTION_DEBUG_

union dimensions{
  float circ_diam;
  float squa_side;
  float rect_sides[2];
  float *custom_sides;
};

typedef struct Leaks{
  int n;
  float *outfw;
  Node **nodes;
} Leaks;

typedef struct Pipe{
  Node *orig;
  Node *dest;

  int geometry;
  union dimensions dimensions;
  float area;
  float rough;
  float length;

  float fluid_viscosity;
  float fluid_density;

  float pressure_in;
  float pressure_out;

  float flowrate;

  float fluid_velocity;

  float friction;

  float flowrate_ideal;
  float flowrate_real;

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

  _Bool has_leak;
  float leak_flowrate;

  float height;

  float fluid_viscosity;
  float fluid_density;

  float fluid_velocity;

  float pressure_measured;    //Pa.   if -1 -> not measured
  float pressure_calculated;  //Pa.   if -1 -> not yet calculated
  float flowrate_measured;    //m³/s. if -1 -> not measured
  float flowrate_calculated;  //m³/s. if -1 -> not yet calculated

  // float pressure_real;
  // float flowrate_real;

  _Bool is_measured;

  int ID;
} Node;

typedef struct Graph{
  Node **nodes;
  Pipe **pipes;

  int depth;
  int width;

  Leaks *leaks;

  int *inc_matrix;
  float *mass_conservation_matrix;

  FrictionModel friction_model;

  float fluid_viscosity;
  float fluid_density;

  int n_pipes;
  int n_nodes;
} Graph;

//Constructors
Pipe *pipe_new(Pipe **ret, Node *orig, Node *dest){
  Pipe *new = malloc(sizeof(Pipe));

  new->orig = orig;
  new->dest = dest;

  new->length = 1;
  new->fluid_viscosity = -1;
  new->fluid_density = -1;
  new->friction = -1;

  new->pressure_in = -1;
  new->pressure_out = -1;

  new->flowrate = -1;
  new->flowrate_ideal = -1;

  new->fluid_velocity = -1;

  pipe_set_geometry(new, GEOMETRY_CIRCULAR);

  pipe_set_diam(new, 0);
  pipe_set_rough(new, 0);

  if (ret != NULL){
    *ret = new;
  }
  return new;
}
Pipe *pipe_copy(Pipe **r, Pipe *s){
  Pipe *n = malloc(sizeof(Pipe));

  n->geometry = s->geometry;
  n->dimensions = s->dimensions;
  n->area = s->area;
  n->rough = s->rough;
  n->length = s->length;

  n->fluid_viscosity = s->fluid_viscosity;
  n->fluid_density = s->fluid_density;

  n->pressure_in = s->pressure_in;
  n->pressure_out = s->pressure_out;

  n->flowrate = s->flowrate;
  n->fluid_velocity = s->fluid_velocity;
  n->friction = s->friction;

  n->flowrate_ideal = s->flowrate_ideal;
  n->flowrate_real = s->flowrate_real;

  n->ID = s->ID;

  if (r != NULL){
    *r = n;
  }
  return n;
}
Node *node_new(Node **ret){
  Node *new = malloc(sizeof(Node));

  new->n_pipes_in = 0;
  new->n_pipes_out = 0;
  new->pipes_in = NULL;
  new->pipes_out = NULL;

  new->fluid_viscosity = -1;
  new->fluid_density = -1;

  new->fluid_velocity = -1;

  new->height = 0;

  new->is_junction = false;
  new->is_connected = false;
  new->is_input = false;
  new->is_output = false;

  new->has_leak = false;
  new->leak_flowrate = 0;

  new->pressure_measured = -1;
  new->pressure_calculated = -1;
  new->flowrate_measured = -1;
  new->flowrate_calculated = -1;

  // new->pressure_real = -1;
  // new->flowrate_real = -1;

  new->is_measured = false;

  if (ret != NULL){
    *ret = new;
  }
  return new;
}
Node *node_copy(Node **r, Node *s){
  Node *n = malloc(sizeof(Node));

  n->n_pipes_in = s->n_pipes_in;
  n->n_pipes_out = s->n_pipes_out;

  n->pipes_in = malloc(sizeof(Pipe *) * n->n_pipes_in);
  n->pipes_out = malloc(sizeof(Pipe *) * n->n_pipes_out);

  n->is_junction = s->is_junction;
  n->is_connected = s->is_connected;
  n->is_input = s->is_input;
  n->is_output = s->is_output;

  n->has_leak = s->has_leak;
  n->leak_flowrate = s->leak_flowrate;

  n->height = s->height;

  n->fluid_velocity = s->fluid_velocity;

  n->fluid_viscosity = s->fluid_viscosity;
  n->fluid_density = s->fluid_density;

  n->pressure_measured = s->pressure_measured;
  n->pressure_calculated = s->pressure_calculated;
  n->flowrate_measured = s->flowrate_measured;
  n->flowrate_calculated = s->flowrate_calculated;

  n->is_measured = s->is_measured;

  n->ID = s->ID;

  if (r != NULL){
    *r = n;
  }
  return n;
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
  g->mass_conservation_matrix = calloc(sizeof(float) * n_nodes * n_pipes, 1);

  g->leaks = NULL;


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

    g->inc_matrix[sorig[i]*g->n_pipes + i] = -1;
    g->inc_matrix[torig[i]*g->n_pipes + i] = 1;
  }

  if (ret != NULL){
    *ret = g;
  }
  g->depth = -1;
  g->width = -1;
  return g;
}
Leaks *leaks_new(Leaks **r, int n){
  Leaks *l = malloc(sizeof(Leaks));
  l->n = n;
  l->outfw = malloc(sizeof(float) * n);
  l->nodes = malloc(sizeof(Node *) * n);

  if (r != NULL){
    *r = l;
  }
  return l;
}
Leaks *leaks_copy(Leaks **r, Leaks *s){
  Leaks *n = malloc(sizeof(Leaks));

  n->n = s->n;
  n->outfw = malloc(sizeof(float) * s->n);
  n->nodes = malloc(sizeof(Node *) * s->n);

  if (r != NULL){
    *r = n;
  }
  return n;
}
Graph *graph_copy(Graph **r, Graph *s){
  Graph *n = malloc(sizeof(Graph));
  n->n_pipes = s->n_pipes;
  n->n_nodes = s->n_nodes;
  n->nodes = malloc(sizeof(Node *) * n->n_nodes);
  n->pipes = malloc(sizeof(Pipe *) * n->n_pipes);

  n->width = s->width;
  n->depth = s->depth;

  //Copy node and pipe values
  for (int i = 0; i < n->n_pipes; i++){
    pipe_copy(&n->pipes[i], s->pipes[i]);
  }
  for (int i = 0; i < n->n_nodes; i++){
    node_copy(&n->nodes[i], s->nodes[i]);
  }

  //Set node and pipe handles
  for (int i = 0; i < n->n_pipes; i++){
    n->pipes[i]->orig = n->nodes[s->pipes[i]->orig->ID];
    n->pipes[i]->dest = n->nodes[s->pipes[i]->dest->ID];
  }
  for (int i = 0; i < n->n_nodes; i++){
    Node *node = n->nodes[i];
    for (int j = 0; j < node->n_pipes_in; j++){
      node->pipes_in[j] = n->pipes[s->nodes[i]->pipes_in[j]->ID];
    }
    for (int j = 0; j < node->n_pipes_out; j++){
      node->pipes_out[j] = n->pipes[s->nodes[i]->pipes_out[j]->ID];
    }
  }

  n->inc_matrix = calloc(sizeof(int) * n->n_nodes * n->n_pipes, 1);
  n->mass_conservation_matrix = calloc(sizeof(float) * n->n_nodes * n->n_pipes, 1);

  for (int j = 0; j < n->n_nodes; j++){
  for (int i = 0; i < n->n_pipes; i++){
    n->mass_conservation_matrix[j*n->n_pipes + i] = s->mass_conservation_matrix[j*n->n_pipes + i];
    n->inc_matrix[j*n->n_pipes + i] = s->inc_matrix[j*n->n_pipes + i];
  }
  }

  n->friction_model = s->friction_model;

  n->fluid_viscosity = s->fluid_viscosity;
  n->fluid_density = s->fluid_density;

  //Copy leak data
  leaks_copy(&n->leaks, s->leaks);
  Leaks *l = n->leaks;

  //Set leak handlers
  for (int i = 0; i < l->n; i++){
    l->nodes[i] = n->nodes[s->leaks->nodes[i]->ID];
  }


  if (r != NULL){
    *r = n;
  }
  return n;
}


//Destructors
void pipe_destroy(Pipe *p){
  if (p != NULL){
    free(p);
  }
}
void node_destroy(Node *n){
  if (n == NULL){
    return;
  }
  free(n->pipes_in);
  free(n->pipes_out);

  free(n);
  return;
}
void leaks_destroy(Leaks *l){
  if (l == NULL){
    return;
  }

  free(l->outfw);
  free(l->nodes);
  free(l);
  return;
}
void graph_destroy(Graph *g){
  if (g == NULL){
    return;
  }

  for (int i = 0; i < g->n_pipes; i++){
    pipe_destroy(g->pipes[i]);
  }
  for (int i = 0; i < g->n_nodes; i++){
    node_destroy(g->nodes[i]);
  }
  free(g->pipes);
  free(g->nodes);
  leaks_destroy(g->leaks);

  free(g->inc_matrix);
  free(g->mass_conservation_matrix);

  free(g);
  return;
}


//Node functions
void node_add_pipe_in(Node *n, Pipe *p){
  if (n->pipes_in == NULL){
    n->pipes_in = malloc(sizeof(Pipe*));
  } else {
    n->pipes_in = realloc(n->pipes_in, sizeof(Pipe*) * (n->n_pipes_in + 1));
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
    n->pipes_out = realloc(n->pipes_out, sizeof(Pipe*) * (n->n_pipes_out + 1));
  }

  #ifdef __GRAPH_C_DEBUG_
  printf("ID: %d, pipes: %d, pipe ID: %d\n", n->ID, n->n_pipes_out, p->ID);
  #endif
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
int node_get_id(Node *n){
  return n->ID;
}
void node_set_is_measured(Node *n, _Bool m){
  n->is_measured = m;
}
_Bool node_get_is_measured(Node *n){
  return n->is_measured;
}
void node_print(Node *n){
  if (n == NULL){
    printf("(null)\n");
    return;
  }
  printf("Node with ID %d\n", n->ID);

  if (n->is_measured){
    printf("MEASURED NODE\n");
  }

  if (n->is_connected){
    if (n->is_junction){
      printf("Junction node\n");
    }

    if (! n->is_input){
      printf("Fed by pipes: ");
      for (int j = 0; j < n->n_pipes_in; j++){
        printf("%d ", n->pipes_in[j]->ID);
      }
    } else {
      printf("Input node");
    }
    printf("\n");

    if (! n->is_output){
      printf("Outputs to:   ");
      for (int j = 0; j < n->n_pipes_out; j++){
        printf("%d ", n->pipes_out[j]->ID);
      }
    } else {
      printf("Output node");
    }
    printf("\n");

    if (n->has_leak){
      printf("Leaks at %f m³/s\n", n->leak_flowrate);
    } else {
      printf("Doesn't leak\n");
    }

    if (n->pressure_measured != -1){  //Pressure measured
      printf("Pressure measured:   %g Pa\n", n->pressure_measured);
    } else {
      printf("Pressure unknown\n");
    }

    if (n->pressure_calculated != -1){  //Pressure calculated
      printf("Pressure calculated: %g Pa\n", n->pressure_calculated);
    } else {
      printf("Pressure not calculated\n");
    }

    if (n->fluid_velocity != -1){  //Velocity
      printf("Fluid velocity:      %g m/s\n", n->fluid_velocity);
    } else {
      printf("Velocity unknown\n");
    }


    if (n->flowrate_measured != -1){  //Flowrate measured
      printf("Flowrate measured:   %g m³/s\n", n->flowrate_measured);
    } else {
      printf("Flowrate unknown\n");
    }

    if (n->flowrate_calculated != -1){  //Flowrate calculated
      printf("Flowrate calculated: %g m³/s\n", n->flowrate_calculated);
    } else {
      printf("Flowrate not calculated\n");
    }

    if (n->flowrate_calculated != -1 && n->flowrate_measured != -1){
      //We known ideal and real flowrates
      if (n->flowrate_calculated != n->flowrate_measured){
        printf("Leakage detected!!!\n");
      }
    }
  } else {    //Not connected
    printf("Disconnected node\n");
  }
}
void graph_set_output_flowrates(Graph *g, float *outflows){
  int num = graph_get_n_output_nodes(g);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_output_node(g, i);
    node_set_flowrate_measured(n, outflows[i]);
  }
}
void node_set_height(Node *n, float h){
  n->height = h;
}
float node_get_height(Node *n){
  return n->height;
}
void node_set_flowrate_measured(Node *n, float f){
  n->is_measured = true;
  n->flowrate_measured = f;
}
float node_get_flowrate_measured(Node *n){
  return n->flowrate_measured;
}
void node_set_flowrate_calculated(Node *n, float f){
  n->flowrate_calculated = f;
}
float node_get_flowrate_calculated(Node *n){
  return n->flowrate_calculated;
}
float node_input_compute_pressure(Node *n){
  float p = 101325 + n->height*9.81*n->fluid_density;
  // float p = n->height*9.806*n->fluid_density;
  node_set_pressure_measured(n, p);
  return p;
}
void node_set_pressure_measured(Node *n, float p){
  n->is_measured = true;
  n->pressure_measured = p;
}
float node_get_pressure_measured(Node *n){
  return n->pressure_measured;
}
void node_set_pressure_calculated(Node *n, float p){
  n->pressure_calculated = p;
}
float node_get_pressure_calculated(Node *n){
  return n->pressure_calculated;
}
void node_set_fluid_viscosity(Node *n, float v){
  n->fluid_viscosity = v;
}
float node_get_fluid_viscosity(Node *n){
  return n->fluid_viscosity;
}
void node_set_fluid_density(Node *n, float d){
  n->fluid_density = d;
}
float node_get_fluid_density(Node *n){
  return n->fluid_density;
}
void node_set_fluid_velocity(Node *n, float v){
  n->fluid_velocity = v;
}
float node_get_fluid_velocity(Node *n){
  return n->fluid_velocity;
}
void node_set_leak_flowrate(Node *n, float f){
  n->has_leak = true;
  n->leak_flowrate = f;
}
float node_get_leak_flowrate(Node *n){
  return n->leak_flowrate;
}

//Pipe functions
int pipe_set_diam(Pipe *p, float d){
  if (p->geometry == GEOMETRY_CIRCULAR ||
      p->geometry == GEOMETRY_CIRCULAR_ANNULUS){
    p->dimensions.circ_diam = d;
    p->area = 1.0/4 * pow(PI, 2) * d;
    // p->area = PI/4 * pow(d, 2);
    // p->area = PI * pow(d/2, 2);
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
void pipe_set_length(Pipe *p, float l){
  p->length = l;
}
void pipe_set_geometry(Pipe *p, int g){
  p->geometry = g;
}
void pipe_print(Pipe *p){
  if (p == NULL){
    printf("(null)\n");
    return;
  }
  printf("Pipe with ID %d:\n", p->ID);
  printf("Length: %gm\n", p->length);
  switch(p->geometry){
    case GEOMETRY_CIRCULAR:
      printf("Geometry: Circular\n");
      printf("Diam %g, area %g\n", p->dimensions.circ_diam, p->area);
      break;
  }

  printf("Roughness: %.7f\n", p->rough);
  printf("Fluid density:   %g\n", p->fluid_density);
  printf("Fluid viscosity: %g\n", p->fluid_viscosity);

  if (p->flowrate != -1){
    printf("Flowrate:        %g m³/s\n", p->flowrate);
  } else {
    printf("Flowrate unknown\n");
  }
  if (p->fluid_velocity != -1){
    printf("Velocity:        %g m/s\n", p->fluid_velocity);
  } else {
    printf("Velocity unknown\n");
  }
  if (p->friction != -1){
    printf("Friction:        %g\n", p->friction);
  } else {
    printf("Friction unknown\n");
  }
  if (p->pressure_in != -1){
    printf("Pressure in:     %g\n", p->pressure_in);
  } else {
    printf("Pressure in unkown\n");
  }
  if (p->pressure_out != -1){
    printf("Pressure out:    %g\n", p->pressure_out);
  } else {
    printf("Pressure out unkown\n");
  }
  printf("%3d -> %3d\n", p->orig->ID, p->dest->ID);
  printf("\n");
}
void pipe_set_fluid_viscosity(Pipe *p, float v){
  p->fluid_viscosity = v;
}
float pipe_get_fluid_viscosity(Pipe *p){
  return p->fluid_viscosity;
}
void pipe_set_fluid_density(Pipe *p, float d){
  p->fluid_density = d;
}
float pipe_get_fluid_density(Pipe *p){
  return p->fluid_density;
}
void pipe_set_friction(Pipe *p, float fd){
  p->friction = fd;
}
float pipe_get_friction(Pipe *p){
  return p->friction;
}
float pipe_compute_friction(Pipe *p, FrictionModel fm){
  float fd;
  fd = fm(p->dimensions.circ_diam,
          p->rough,
          p->fluid_density,
          p->fluid_viscosity,
          p->fluid_velocity);

  p->friction = fd;

  return fd;
}


//Graph functions
void graph_compute_mass_conservation_matrix(Graph *g){
  for (int j = 0; j < g->n_nodes; j++){
  for (int i = 0; i < g->n_pipes; i++){
    g->mass_conservation_matrix[j*g->n_pipes + i] = g->inc_matrix[j*g->n_pipes + i] * g->pipes[i]->area;
  }
  }
}
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
void graph_set_lengths(Graph *g, float *l){
  for (int i = 0; i < g->n_pipes; i++){
    pipe_set_length(g->pipes[i], l[i]);
  }
}
Node **graph_get_nodes(Graph *g){
  return g->nodes;
}
Pipe **graph_get_pipes(Graph *g){
  return g->pipes;
}
void graph_print_incidence_matrix(Graph *g){
  if (g->inc_matrix == NULL){
    printf("(null)\n");
    return;
  }
  for (int j = 0; j < g->n_nodes; j++){
  for (int i = 0; i < g->n_pipes; i++){
    if (g->inc_matrix[j*g->n_pipes + i] < 0){
      printf("%d ", g->inc_matrix[j*g->n_pipes + i]);
    } else {
      printf(" %d ", g->inc_matrix[j*g->n_pipes + i]);
    }
  }
  printf("\n");
  }
}
void graph_print_mass_conservation_matrix(Graph *g){
  if (g->mass_conservation_matrix == NULL){
    printf("(null)\n");
    return;
  }
  for (int j = 0; j < g->n_nodes; j++){
  for (int i = 0; i < g->n_pipes; i++){
    if (g->mass_conservation_matrix[j*g->n_pipes + i] < 0){
      printf("%.4f ", g->mass_conservation_matrix[j*g->n_pipes + i]);
    } else {
      printf(" %.4f ", g->mass_conservation_matrix[j*g->n_pipes + i]);
    }
  }
  printf("\n");
  }
}
void graph_print(Graph *g){
  printf("---------------------------------------\n");
  printf("Printing graph: \n");
  if (g == NULL){
    printf("(null)\n");
    return;
  }
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

  // printf("\n\nMASS CONSERVATION MATRIX:\n");
  // graph_print_mass_conservation_matrix(g);

  printf("\n\nLEAK DATA: \n");
  graph_print_leaks_data(g);
  printf("---------------------------------------\n");
}
void graph_print_disconnected_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (!g->nodes[i]->is_connected){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_connected_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_connected){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_junction_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_junction){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_input_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_input){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_output_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->is_output){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
void graph_print_leak_nodes(Graph *g){
  if (g == NULL){
    printf("(null)\n");
    return;
  }
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i]->has_leak){
      node_print(g->nodes[i]);
      printf("\n");
    }
  }
}
int graph_get_n_nodes(Graph *g){
  return g->n_nodes;
}
int graph_get_n_disconnected_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (!g->nodes[i]->is_connected){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_connected_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_connected){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_junction_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_junction){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_input_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_input){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_measurement_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_measured){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_output_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_output){
        sum++;
      }
    }
  }
  return sum;
}
int graph_get_n_leak_nodes(Graph *g){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->has_leak){
        sum++;
      }
    }
  }
  return sum;
}
Node *graph_get_nth_node(Graph *g, int i){
  return g->nodes[i];
}
Node *graph_get_nth_disconnected_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (!g->nodes[i]->is_connected){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_connected_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_connected){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_measurement_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_measured){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_junction_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_junction){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_leak_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->has_leak){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_input_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_input){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}
Node *graph_get_nth_output_node(Graph *g, int index){
  int sum = 0;
  for (int i = 0; i < g->n_nodes; i++){
    if (g->nodes[i] != NULL){
      if (g->nodes[i]->is_output){
        if (index == sum){
          return g->nodes[i];
        }
        sum++;
      }
    }
  }
  return NULL;
}

float graph_get_total_outflow(Graph *g){
  int num_outputs = graph_get_n_output_nodes(g);
  float sum = 0;
  for (int i = 0; i < num_outputs; i++){
    Node *n = graph_get_nth_output_node(g, i);
    float flow = node_get_flowrate_measured(n);
    if (flow != -1){
      sum += flow;
    }
  }
  int num_leaks = graph_get_n_leak_nodes(g);
  for (int i = 0; i < num_leaks; i++){
    Node *n = graph_get_nth_leak_node(g, i);
    float flow = node_get_leak_flowrate(n);
    sum += flow;
  }
  return sum;
}
float graph_get_total_leak_outflow(Graph *g){
  int num_leaks = graph_get_n_leak_nodes(g);
  float sum = 0;
  for (int i = 0; i < num_leaks; i++){
    Node *n = graph_get_nth_leak_node(g, i);
    float flow = node_get_leak_flowrate(n);
    sum += flow;
  }
  printf("%f\n", sum);
  return sum;
}
float graph_get_total_calculated_outflow(Graph *g){
  int num_outputs = graph_get_n_output_nodes(g);
  float sum = 0;
  for (int i = 0; i < num_outputs; i++){
    Node *n = graph_get_nth_output_node(g, i);
    float flow = node_get_flowrate_calculated(n);
    if (flow != -1){
      sum += flow;
    }
  }
  return sum;
}

void graph_set_inflow_evenly(Graph *g){
  int num_inputs = graph_get_n_input_nodes(g);
  float total_outflow = graph_get_total_outflow(g);

  float even_outflow = total_outflow / num_inputs;

  for (int i = 0; i < num_inputs; i++){
    Node *n = graph_get_nth_input_node(g, i);
    node_set_flowrate_calculated(n, even_outflow);
  }
}
void graph_add_leaks_to_measured_nodes(Graph *g){
  int num_measured = graph_get_n_measurement_nodes(g);
  for (int i = 0; i < num_measured; i++){
    Node *nm = graph_get_nth_measurement_node(g, i);
    node_set_flowrate_measured(nm, node_get_flowrate_calculated(nm));
  }


  int num_leaks = graph_get_n_leak_nodes(g);

  for (int index_leak = 0; index_leak < num_leaks; index_leak++){
    Node *nl = graph_get_nth_leak_node(g, index_leak);

    if (nl->is_measured){
      if (nl->flowrate_measured != -1){
        node_set_flowrate_measured(nl, node_get_flowrate_measured(nl) + nl->leak_flowrate);
      } else {
        node_set_flowrate_measured(nl, node_get_flowrate_calculated(nl) + nl->leak_flowrate);
      }
    }

    float leak_flowrate = nl->leak_flowrate;
    float flowrate_per_area;

    Node **node_vector;
    Node **aux = NULL;
    float *area_vector;
    float *aux_area;

    int vector_len = nl->n_pipes_in;
    int aux_len = 0;
    node_vector = malloc(sizeof(Node *) * vector_len);
    area_vector = malloc(sizeof(float) * vector_len);

    float total_area_in = 0;
    for (int i = 0; i < nl->n_pipes_in; i++){
      Pipe *p = nl->pipes_in[i];
      node_vector[i] = p->orig;
      area_vector[i] = p->area;
      total_area_in += p->area;
    }

    flowrate_per_area = leak_flowrate / total_area_in;

    while (vector_len != 0){
      for (int i = 0; i < vector_len; i++){
        Node *n = node_vector[i];
        int n_pipes = n->n_pipes_in;

        printf("Checking node %d\n", n->ID);

        float node_leaked_flowrate = flowrate_per_area * area_vector[i];

        if (n->is_measured){
          node_set_flowrate_measured(n, node_get_flowrate_calculated(n) + node_leaked_flowrate);
        }

        total_area_in = 0;
        for (int j = 0; j < n_pipes; j++){
          total_area_in += n->pipes_in[j]->area;
        }
        flowrate_per_area = node_leaked_flowrate / total_area_in;

        for (int j = 0; j < n_pipes; j++){
          Pipe *p = n->pipes_in[j];

          aux_len++;
          if (aux == NULL){
            aux = malloc(sizeof(Node *) * aux_len);
            aux_area = malloc(sizeof(float) * aux_len);
          } else {
            aux = realloc(aux, sizeof(Node *) * aux_len);
            aux_area = realloc(aux_area, sizeof(float) * aux_len);
          }
          aux[aux_len - 1] = p->orig;
          aux_area[aux_len -1] = p->area;
        }
      }

      free(node_vector);
      free(area_vector);
      node_vector = aux;
      area_vector = aux_area;

      aux = NULL;
      aux_area = NULL;
      vector_len = aux_len;
      aux_len = 0;
    }

    if (aux != NULL){
      free(aux);
    }
    free(node_vector);

  }
}
void graph_add_leaks_to_inflow(Graph *g){
  int num_inputs = graph_get_n_input_nodes(g);
  float total_outflow = graph_get_total_outflow(g);

  float even_outflow = total_outflow / num_inputs;

  for (int i = 0; i < num_inputs; i++){
    Node *n = graph_get_nth_input_node(g, i);
    node_set_flowrate_measured(n, even_outflow);
  }
}

float graph_get_total_calculated_inflow(Graph *g){
  int num = graph_get_n_input_nodes(g);
  float sum = 0;
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_input_node(g, i);
    float flow = node_get_flowrate_calculated(n);
    if (flow != -1){
      sum += flow;
    }
  }
  return sum;
}
float graph_get_total_inflow(Graph *g){
  int num = graph_get_n_input_nodes(g);
  float sum = 0;
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_input_node(g, i);
    float flow = node_get_flowrate_measured(n);
    if (flow != -1){
      sum += flow;
    }
  }
  return sum;
}

void graph_inflow_real_to_calc(Graph *g){
  int num = graph_get_n_input_nodes(g);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_input_node(g, i);
    node_set_flowrate_calculated(n, node_get_flowrate_measured(n));
  }
}
void graph_inflow_calc_to_real(Graph *g){
  int num = graph_get_n_input_nodes(g);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_input_node(g, i);
    node_set_flowrate_measured(n, node_get_flowrate_calculated(n));
  }
}
void graph_outflow_real_to_calc(Graph *g){
  int num = graph_get_n_output_nodes(g);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_output_node(g, i);
    node_set_flowrate_calculated(n, node_get_flowrate_measured(n));
  }
}
void graph_outflow_calc_to_real(Graph *g){
  int num = graph_get_n_output_nodes(g);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_output_node(g, i);
    node_set_flowrate_measured(n, node_get_flowrate_calculated(n));
  }
}
void graph_set_fluid_viscosity(Graph *g, float v){
  g->fluid_viscosity  = v;
  for (int i = 0; i < g->n_pipes; i++){
    pipe_set_fluid_viscosity(g->pipes[i], v);
  }
  for (int i = 0; i < g->n_nodes; i++){
    node_set_fluid_viscosity(g->nodes[i], v);
  }
}
void graph_set_fluid_density(Graph *g, float d){
  g->fluid_density = d;
  for (int i = 0; i < g->n_pipes; i++){
    pipe_set_fluid_density(g->pipes[i], d);
  }
  for (int i = 0; i < g->n_nodes; i++){
    node_set_fluid_density(g->nodes[i], d);
  }
}
float graph_get_fluid_viscosity(Graph *g){
  return g->fluid_viscosity;
}
float graph_get_fluid_density(Graph *g){
  return g->fluid_density;
}
void graph_set_friction_model(Graph *g, FrictionModel fm){
  g->friction_model = fm;
}
void graph_backpropagate_flowrate(Graph *g){
  #ifdef __GRAPH_C_DEBUG_
  printf("BACK PROPAGATING FLOWRATE\n");
  #endif

  Node **vector;
  Node **aux = NULL;
  int n_nodes = graph_get_n_output_nodes(g);
  int n_aux = 0;
  vector = malloc(sizeof(Node *) * n_nodes);
  for (int i = 0; i < n_nodes; i++){
    vector[i] = graph_get_nth_output_node(g, i);
  }


  while (n_nodes != 0){
    for (int i = 0; i < n_nodes; i++){
      Node *n = vector[i];
      int n_pipes = n->n_pipes_in;

      //Sum flowrate demanded from outgoing pipes:
      if (! n->is_output){
        float sum = 0;
        for (int k = 0; k < n->n_pipes_out; k++){
          sum += n->pipes_out[k]->flowrate;
        }

        n->flowrate_calculated = sum;
      }

      float sum_area_in = 0;
      for (int j = 0; j < n_pipes; j++){
        sum_area_in += n->pipes_in[j]->area;
      }
      float flowrate_divided = n->flowrate_calculated / sum_area_in;


      for (int j = 0; j < n_pipes; j++){
        Pipe *p = n->pipes_in[j];
        p->flowrate = flowrate_divided * p->area;

        p->fluid_velocity = p->flowrate / p->area;

        if (p->orig->is_input == false){
          int is_added = false;
          for (int k = 0; k < n_aux; k++){
            if (aux[k] == p->orig){
              is_added = true;
            }
          }
          if (! is_added){
            n_aux ++;
            if (aux == NULL){
              aux = malloc(sizeof(Node *) * n_aux);
            } else {
              aux = realloc(aux, sizeof(Node *) *n_aux);
            }
            aux[n_aux - 1] = p->orig;
          }
        }
      }
    }

    free(vector);
    vector = aux;
    aux = NULL;
    n_nodes = n_aux;
    n_aux = 0;
  }
  if (aux != NULL){
    free(aux);
  }
  free(vector);
}
void graph_propagate_pressure(Graph *g){
  #ifdef __GRAPH_C_DEBUG_
  printf("PROPAGATING PRESSURES:\n");
  #endif

  Node **vector;
  Node **aux = NULL;
  int n_nodes = graph_get_n_input_nodes(g);
  int n_aux = 0;
  vector = malloc(sizeof(Node *) * n_nodes);

  for (int i = 0; i < n_nodes; i++){
    vector[i] = graph_get_nth_input_node(g, i);
  }

  while (n_nodes != 0){
    for (int i = 0; i < n_nodes; i++){
      Node *n = vector[i];
      int n_pipes = n->n_pipes_out;

      float sum_area_out = 0;
      for (int j = 0; j < n_pipes; j++){
        sum_area_out += n->pipes_out[j]->area;
      }

      float sum_area_in = 0;
      for (int j = 0; j < n->n_pipes_in; j++){
        sum_area_in += n->pipes_in[j]->area;
      }

      //Get pressure for every out-pipe
      float pressure_divided = n->pressure_calculated;

      // pressure_divided = n->pressure_calculated / sum_area_out;
      // if (! n->is_input){
      //   float pout = n->pressure_calculated/sum_area_in;
      //   float pout = (n->pressure_calculated*sum_area_in)/sum_area_out;
      //   pressure_divided = pout / sum_area_out;
      // } else {
      //   pressure_divided = n->pressure_calculated / sum_area_out;
      // }

      for (int j = 0; j < n_pipes; j++){
        Pipe *p = n->pipes_out[j];

        // printf("Doing pipe %d\n", p->ID);

        //Set pressure in
        p->pressure_in = pressure_divided;
        // p->pressure_in = pressure_divided * p->area;

        // printf("Pressure in  = %f\n", p->pressure_in);


        //Calculate friction
        pipe_compute_friction(p, g->friction_model);

        //Set pressure in next node
        float pressure_drop = calculate_pressure_drop(p->fluid_velocity,
                                                      p->dimensions.circ_diam,
                                                      p->length,
                                                      p->friction,
                                                      p->fluid_density);

        #ifdef __GRAPH_C_DEBUG_
        printf("Pressure drop: %f\n", pressure_drop);
        #endif

        p->pressure_out = p->pressure_in - pressure_drop;
        // p->pressure_out = p->friction / (pow(p->fluid_velocity, 2) * p->area);

        // printf("Pressure out = %f\n", p->pressure_out);

        p->dest->pressure_calculated = p->pressure_out;



        // printf("\n");
        if (p->dest->is_output == false){
          n_aux++;
          if (aux == NULL){
            aux = malloc(sizeof(Node *) * n_aux);
          } else {
            aux = realloc(aux, sizeof(Node *) * n_aux);
          }
          aux[n_aux - 1] = p->dest;
        }
      }
    }
    free(vector);
    vector = aux;
    aux = NULL;
    n_nodes = n_aux;
    n_aux = 0;
  }
  if (aux != NULL){
    free(aux);
  }
  free(vector);
}

//LEAKS FUNCTIONS
Leaks *graph_generate_random_leaks(Graph *g, int num){
  //Check if num is valid
  int n_junction = graph_get_n_junction_nodes(g);
  if (num > n_junction){
    return NULL;
  }

  //Free previous leaks memory
  if (g->leaks != NULL){
    free(g->leaks);
  }
  Leaks *leaks = malloc(sizeof(Leaks));

  //Generate leaks structure
  leaks->n = num;
  leaks->outfw = malloc(sizeof(float) * num);
  leaks->nodes = malloc(sizeof(Node *) * num);

  //Auxiliary variables
  _Bool is_valid;
  int node;
  int k = 0;
  int aux[num];
  for (int i = 0; i < num; i++){
    aux[i] = -1;
  }

  //Generate the leaks...

  //For each leak to generate...
  for (int i = 0; i < num; i++){
    is_valid = false;

    //Generate randoms until valid
    while (! is_valid){
      node = random() % n_junction;

      //check if valid
      is_valid = true;
      for (int j = 0; j < num; j++){
        if (aux[j] == node){
          is_valid = false;
        }
      }

      //Generate the leak
      if (is_valid){
        float outflow = (float)rand()/(float)(RAND_MAX/MAX_LEAK_OUTFLOW);
        Node *n = graph_get_nth_junction_node(g, node);
        node_set_leak_flowrate(n, outflow);
        leaks->nodes[k] = n;

        aux[k] = node;
        k++;
      }

    }
  }

  g->leaks = leaks;
  return leaks;
}
void graph_print_leaks_data(Graph *g){
  Leaks *l = g->leaks;
  if (l == NULL){
    printf("Leaks not yet initialized\n");
    return;
  }
  for (int i = 0; i < l->n; i++){
    Node *n = l->nodes[i];
    printf("Node %d: %f m³/s\n", node_get_id(n), node_get_leak_flowrate(n));
  }
}
_Bool graph_has_leaks(Graph *g){
  float graph_real_inflow = graph_get_total_inflow(g);
  float graph_calc_inflow = graph_get_total_calculated_inflow(g);
  float graph_real_outflow = graph_get_total_outflow(g);
  float graph_calc_outflow = graph_get_total_calculated_outflow(g);

  #ifdef __GRAPH_C_DETECTION_DEBUG_
  printf("Leak detection: \n");
  printf("Real inflow:  %f\n", graph_real_inflow);
  printf("Calc inflow:  %f\n", graph_calc_inflow);
  printf("Real outflow: %f\n", graph_real_outflow);
  printf("Calc outflow: %f\n", graph_calc_outflow);
  #endif

  if (graph_real_inflow == graph_calc_inflow &&
      graph_real_inflow == graph_real_outflow &&
      graph_real_inflow == graph_calc_outflow){
    return false;
  } else {
    return true;
  }
}
Leaks *graph_find_leaks(Graph *g){
  Leaks *l = leaks_new(NULL, 0);

  return l;
}
Graph *graph_optimize_naive(Graph *g){
  return NULL;
}
void graph_calculate_geometry(Graph *G){
  Graph *g = graph_copy(NULL, G);
  int n_input = graph_get_n_input_nodes(g);
  int width = n_input;
  int depth = 0;
  _Bool finished;

  do {
    finished = true;
    depth++;
    n_input = graph_get_n_input_nodes(g);
    int level_width = 0;
    Node **nodev = malloc(sizeof(Node *) * n_input);
    for (int i = 0; i < n_input; i++){
      nodev[i] = graph_get_nth_input_node(g, i);
    }
    for (int i = 0; i < n_input; i++){
      Node *n = nodev[i];
      n->is_input = false;
      level_width += n->n_pipes_out;
      for (int p = 0; p < n->n_pipes_out; p++){
        n->pipes_out[p]->dest->is_input = true;
        finished = false;
      }
    }
    if (level_width >= width){
      width = level_width;
    }

    free(nodev);
  } while (! finished);

  graph_destroy(g);
  G->width = width;
  G->depth = depth;
}
int graph_get_depth(Graph *g){
  if (g->depth == -1){
    graph_calculate_geometry(g);
  }
  return g->depth;
}
int graph_get_width(Graph *g){
  if (g->width == -1){
    graph_calculate_geometry(g);
  }
  return g->width;
}
void graph_cut_node(Graph *g, int node_i){
  Node *ni = g->nodes[node_i];
  if (ni == NULL){
    return;
  }

  Node **node_vector;
  Node **aux = NULL;

  int vector_len = ni->n_pipes_out;
  int aux_len = 0;

  node_vector = malloc(sizeof(Node *) * vector_len);
  for (int i = 0; i < ni->n_pipes_out; i++){
    node_vector[i] = ni->pipes_out[i]->dest;
  }

  graph_del_node(g, node_i);

  while (vector_len != 0){
    for (int i = 0; i < vector_len; i++){
      Node *n = node_vector[i];
      int n_pipes = n->n_pipes_out;

      for (int j = 0; j < n_pipes; j++){
        Pipe *p = n->pipes_out[j];

        if (g->nodes[p->dest->ID] != NULL){
          aux_len++;

          if (aux == NULL){
            aux = malloc(sizeof(Node *) * aux_len);
          } else {
            aux = realloc(aux, sizeof(Node *) * aux_len);
          }
          aux[aux_len -1] = p->dest;
        }
      }

      graph_del_node(g, n->ID);
    }

    free(node_vector);
    node_vector = aux;
    aux = NULL;
    vector_len = aux_len;
    aux_len = 0;
  }
}
Node *graph_del_node(Graph *g, int node_i){
  Node *n = g->nodes[node_i];
  g->nodes[node_i] = NULL;
  return n;
}
float node_measurement_get_diff(Node *n){
  if (n->is_measured){
    if (n->flowrate_measured != -1){
      return n->flowrate_measured - n->flowrate_calculated;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}
float node_measurement_get_successors_diff(Node *n){
  float result = 0;
  int vector_len = 1;
  int aux_len = 0;

  Node **node_vector = malloc(sizeof(Node *) * vector_len);
  Node **aux = NULL;

  node_vector[0] = n;

  while (vector_len != 0){
    for (int i = 0; i < vector_len; i++){
      Node *n = node_vector[i];
      int n_pipes = n->n_pipes_out;

      for (int j = 0; j < n_pipes; j++){
        Pipe *p = n->pipes_out[j];
        Node *dest = p->dest;

        if (dest->is_measured){
          if (dest->flowrate_measured != -1){
            result += node_measurement_get_diff(dest);
          }
        } else {
          aux_len++;
          if (aux == NULL){
            aux = malloc(sizeof(Node *) * aux_len);
          } else {
            aux = realloc(aux, sizeof(Node *) * aux_len);
          }

          aux[aux_len -1] = dest;
        }
      }
    }
    free(node_vector);
    node_vector = aux;
    aux = NULL;
    vector_len = aux_len;
    aux_len = 0;
  }

  if (aux != NULL){
    free(aux);
  }
  free(node_vector);

  return result;
}
void graph_plot(Graph *g){
  //Make copy of graph to not destroy original
  g = graph_copy(NULL, g);

  const float float_tolerance = 0.000001;

  //Variables for running graphviz
  int p[2], pid;
  if (pipe(p) < 0){
    perror("Could not pipe:");
  }
  if ((pid = fork()) < 0){
    perror("Could not fork:");
  }

  //Fork for graphviz
  if (pid == 0){
    close(0);
    dup(p[0]);
    close(p[0]);
    close(p[1]);

    close(1);
    fopen("test.png", "w");

    execlp("dot", "dot", "-Tpng", NULL);
    exit(1);
  }
  close(p[0]);

  //String for piping to graphviz
  int str_size = 1024;
  int buffer_size = 1024;

  //Graphviz header
  char *str = malloc(sizeof(char) * str_size);
  char buffer[buffer_size];
  strcpy(str, "digraph G{fontname=\"Helvetica,Arial,sans-serif\"\nnode [fontname=\"Helvetica,Arial,sans-serif\"]\nedge [fontname=\"Helvetica,Arial,sans-serif\"]\n");

  //Clone again to g2 and cut nodes after leak
  Graph *g2 = graph_copy(NULL, g);
  int *nodes_to_del = NULL;
  int nodes_to_del_l = 0;
  int *nodes_to_cut = NULL;
  int nodes_to_cut_l = 0;
  for (int i = 0; i < graph_get_n_measurement_nodes(g2); i++){
    Node *nm = graph_get_nth_measurement_node(g2, i);
    printf("Node %d has %f diff and %f succ diff\n", nm->ID, node_measurement_get_diff(nm), node_measurement_get_successors_diff(nm));
    if (nm != NULL){
      _Bool del = false;
      _Bool cut = false;
      if (fabs(node_measurement_get_diff(nm) - node_measurement_get_successors_diff(nm)) < float_tolerance && node_measurement_get_diff(nm) > 0){
        printf("Deleting node %d\n", nm->ID);
        del = true;
      }
      if (nm->flowrate_measured == nm->flowrate_calculated){
        printf("Cutting node %d\n", nm->ID);
        cut = true;
      }
      if (del){
        if (nodes_to_del == NULL){
          nodes_to_del = malloc(sizeof(int));
          nodes_to_del[nodes_to_del_l] = nm->ID;
        } else {
          nodes_to_del = realloc(nodes_to_del, sizeof(int) * (nodes_to_del_l+1));
          nodes_to_del[nodes_to_del_l] = nm->ID;
        }
        nodes_to_del_l++;
      }
      if (cut){
        if (nodes_to_cut == NULL){
          nodes_to_cut = malloc(sizeof(int));
          nodes_to_cut[nodes_to_cut_l] = nm->ID;
        } else {
          nodes_to_cut = realloc(nodes_to_cut, sizeof(int) * (nodes_to_cut_l+1));
          nodes_to_cut[nodes_to_cut_l] = nm->ID;
        }
        nodes_to_cut_l++;
      }
    }
  }
  for (int i = 0; i < nodes_to_del_l; i++){
    graph_del_node(g2, nodes_to_del[i]);
  }
  for (int i = 0; i < nodes_to_cut_l; i++){
    graph_cut_node(g2, nodes_to_cut[i]);
  }
  if (nodes_to_cut != NULL){
    free(nodes_to_cut);
  }
  if (nodes_to_del != NULL){
    free(nodes_to_del);
  }

  // g = g2;

  //FAULTY NODE FORMAT
  sprintf(buffer, "node [shape=diamond color=red]; ");
  if (strlen(str) + strlen(buffer) + 10 > str_size){
    str_size *= 2;
    str = realloc(str, str_size);
  }
  strcat(str, buffer);
  for (int i = 0; i < graph_get_n_nodes(g2); i++){
    Node *node = graph_get_nth_node(g2, i);

    if (node != NULL && g->nodes[i]->is_connected){
      sprintf(buffer, "%d; ", node->ID);
      if (strlen(str) + strlen(buffer) + 10 < str_size){
        str_size*=2;
        str = realloc(str, str_size);
      }
      strcat(str, buffer);
    }
  }

  //INPUT NODE FORMAT
  sprintf(buffer, "node [shape=ellipse color=blue]; ");
  if (strlen(str) + strlen(buffer) + 10 > str_size){
    str_size *= 2;
    str = realloc(str, str_size);
  }
  strcat(str, buffer);
  for (int i = 0; i < graph_get_n_input_nodes(g); i++){
    Node *node = graph_get_nth_input_node(g, i);
    if (node != NULL){
      sprintf(buffer, "%d; ", node->ID);
      if (strlen(str) + strlen(buffer) + 10 < str_size){
        str_size*=2;
        str = realloc(str, str_size);
      }
      strcat(str, buffer);
    }
  }

  //JUNCTION NODE FORMAT
  sprintf(buffer, "node [shape=diamond color=black]; ");
  if (strlen(str) + strlen(buffer) + 10 > str_size){
    str_size *= 2;
    str = realloc(str, str_size);
  }
  strcat(str, buffer);
  for (int i = 0; i < graph_get_n_junction_nodes(g); i++){
    Node *node = graph_get_nth_junction_node(g, i);

    if (node != NULL){
      sprintf(buffer, "%d; ", node->ID);
      if (strlen(str) + strlen(buffer) + 10 < str_size){
        str_size*=2;
        str = realloc(str, str_size);
      }
      strcat(str, buffer);
    }
  }

  //OUTPUT NODE FORMAT
  sprintf(buffer, "node [shape=box color=green]; ");
  if (strlen(str) + strlen(buffer) + 10 > str_size){
    str_size*=2;
    str = realloc(str, str_size);
  }
  strcat(str, buffer);
  for (int i = 0; i < graph_get_n_output_nodes(g); i++){
    Node *node = graph_get_nth_output_node(g, i);

    if (node != NULL){
      sprintf(buffer, "%d; ", node->ID);
      if (strlen(str) + strlen(buffer) + 10 > str_size){
        str_size*=2;
        str = realloc(str, str_size);
      }
      strcat(str, buffer);
    }
  }

  printf("Draw pupes\n");

  //DRAW PIPES
  Node **node_vector;
  Node **aux = NULL;

  int vector_len = graph_get_n_input_nodes(g);
  int aux_len = 0;
  node_vector = malloc(sizeof(Node *) * vector_len);

  for (int i = 0; i < vector_len; i++){
    node_vector[i] = graph_get_nth_input_node(g, i);
  }

  while (vector_len != 0){
    for (int i = 0; i < vector_len; i++){
      Node *n = node_vector[i];
      int n_pipes = n->n_pipes_out;

      float flowrate_inpipes = 0;
      float flowrate_outpipes = 0;
      for (int j = 0; j < n->n_pipes_out; j++){
        flowrate_outpipes += n->pipes_out[j]->flowrate;
      }
      for (int j = 0; j < n->n_pipes_in; j++){
        flowrate_inpipes += n->pipes_in[j]->flowrate;
      }

      printf("\n");
      for (int j = 0; j < n_pipes; j++){
        Pipe *p = n->pipes_out[j];
        if (g->nodes[p->dest->ID] != NULL){

          printf("Drawing %d->%d\n", p->orig->ID, p->dest->ID);

          sprintf(buffer, "%d->%d ", p->orig->ID, p->dest->ID);
          if (strlen(str) + strlen(buffer) + 10 > str_size){
            str_size*=2;
            str = realloc(str, str_size);
          }
          strcat(str, buffer);

          if (g2->nodes[p->orig->ID] != NULL && g2->nodes[p->dest->ID] != NULL){
            sprintf(buffer, "[color=red] ");
            if (strlen(str) + strlen(buffer) + 10 > str_size){
              str_size*=2;
              str = realloc(str, str_size);
            }
            strcat(str, buffer);
          }

          ////ADD PIPE INFORMATION
          //sprintf(buffer, "[label=\"Len: %.1fm\nPrssIn: %.1f Pa\nPrssOut: %.1f Pa\nFlowrate: %.4f m³/s\n\"];", p->length, p->pressure_in, p->pressure_out, p->flowrate);
          //if (strlen(str) + strlen(buffer) + 10 > str_size){
          //  str_size*=2;
          //  str = realloc(str, str_size);
          //}
          //strcat(str, buffer);


          if (p->dest->is_output == false){
            aux_len++;
            if (aux == NULL){
              aux = malloc(sizeof(Node *) * aux_len);
            } else {
              aux = realloc(aux, sizeof(Node *) * aux_len);
            }
            aux[aux_len -1] = p->dest;
          }
        }
      }
    }
    free(node_vector);
    node_vector = aux;
    aux = NULL;
    vector_len = aux_len;
    aux_len = 0;
  }

  if (aux != NULL){
    free(aux);
  }
  free(node_vector);


  // Node **nodes = graph_get_nodes(g);
  // for (int i = 0; i < graph_get_n_nodes(g); i++){
  //   Node *node = nodes[i];
  //   for (int j = 0; j < node->n_pipes_out; j++){
  //     Pipe *pipe = node->pipes_out[j];
  //     Node *dest = pipe->dest;
  //     sprintf(buffer, "%d->%d [label=\"Len: %.1fm\nPrssIn: %.1f Pa\nPrssOut: %.1f Pa\nFlowrate: %.4f m³/s\n\"];", node->ID, dest->ID, pipe->length, pipe->pressure_in, pipe->pressure_out, pipe->flowrate);
  //     // sprintf(buffer, "%d->%d ", node->ID, dest->ID);
  //     if (strlen(str) + strlen(buffer) + 10 > str_size){
  //       str_size*=2;
  //       str = realloc(str, str_size);
  //     }
  //     strcat(str, buffer);
  //   }
  // }


  //TERMINATION
  strcat(str, "}");
  write(p[1], str, strlen(str)+1);
  close(p[1]);

  graph_destroy(g);
  graph_destroy(g2);
}

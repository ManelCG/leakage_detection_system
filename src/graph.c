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

  float height;

  float fluid_viscosity;
  float fluid_density;

  float fluid_velocity;

  float pressure_measured;    //Pa.   if -1 -> not measured
  float pressure_calculated;  //Pa.   if -1 -> not yet calculated
  float flowrate_measured;    //m³/s. if -1 -> not measured
  float flowrate_calculated;  //m³/s. if -1 -> not yet calculated

  int ID;
} Node;

typedef struct Graph{
  Node **nodes;
  Pipe **pipes;

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

  new->fluid_velocity = -1;

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

  new->fluid_viscosity = -1;
  new->fluid_density = -1;

  new->fluid_velocity = -1;

  new->height = 0;

  new->is_junction = false;
  new->is_connected = false;
  new->is_input = false;
  new->is_output = false;

  new->pressure_measured = -1;
  new->pressure_calculated = -1;
  new->flowrate_measured = -1;
  new->flowrate_calculated = -1;

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
  g->mass_conservation_matrix = calloc(sizeof(float) * n_nodes * n_pipes, 1);


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
int node_get_id(Node *n){
  return n->ID;
}
void node_print(Node *n){
  printf("Node with ID %d\n", n->ID);

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
void node_set_height(Node *n, float h){
  n->height = h;
}
float node_get_height(Node *n){
  return n->height;
}
void node_set_flowrate_measured(Node *n, float f){
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
  n->pressure_measured = p;
  return p;
}
void node_set_pressure_measured(Node *n, float p){
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
void pipe_set_length(Pipe *p, float l){
  p->length = l;
}
void pipe_set_geometry(Pipe *p, int g){
  p->geometry = g;
}
void pipe_print(Pipe *p){
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
int graph_get_n_nodes(Graph *g){
  return g->n_nodes;
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
Node *graph_get_nth_node(Graph *g, int i){
  return g->nodes[i];
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
    float flow = node_get_flowrate_measured(n);
    if (flow != -1){
      sum += flow;
    }
  }
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
  printf("BACK PROPAGATING FLOWRATE\n");

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
  printf("PROPAGATING PRESSURES:\n");

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

        printf("Pressure drop: %f\n", pressure_drop);

        p->pressure_out = p->pressure_in - pressure_drop;
        // p->pressure_out = p->friction / (pow(p->fluid_velocity, 2) * p->area);

//         printf("Pressure out = %f\n", p->pressure_out);

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
_Bool graph_detect_leaks(Graph *g){
  float graph_real_inflow = graph_get_total_inflow(g);
  float graph_calc_inflow = graph_get_total_calculated_inflow(g);
  float graph_real_outflow = graph_get_total_outflow(g);
  float graph_calc_outflow = graph_get_total_calculated_outflow(g);

  if (graph_real_inflow == graph_calc_inflow &&
      graph_real_inflow == graph_real_outflow &&
      graph_real_inflow == graph_calc_outflow){
    return false;
  } else {
    return true;
  }
}

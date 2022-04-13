#include <fluid_mechanics.h>

typedef struct Pipe Pipe;
typedef struct Node Node;

typedef struct Graph Graph;

Pipe *pipe_new(Pipe **, Node *, Node *);
Node *node_new(Node **);

Graph *graph_new(Graph **ret, int n_pipes, int *sorig, int *torig);

//Node functions
void node_add_pipe_in(Node *n, Pipe *p);
void node_add_pipe_out(Node *n, Pipe *p);
void node_set_id(Node *n, int);
int node_get_id(Node *n);
void node_print(Node *n);

void node_set_height(Node *n, float h);
float node_get_height(Node *n);

void node_set_flowrate_measured(Node *n, float f);
void node_set_flowrate_calculated(Node *n, float f);
float node_get_flowrate_measured(Node *n);
float node_get_flowrate_calculated(Node *n);

void node_set_fluid_velocity(Node *n, float v);
float node_get_fluid_velocity(Node *n);

void node_set_pressure_measured(Node *n, float p);
float node_get_pressure_measured(Node *n);
void node_set_pressure_calculated(Node *n, float p);
float node_get_pressure_calculated(Node *n);
float node_input_compute_pressure(Node *n);

//Pipe functions
void pipe_set_geometry(Pipe *p, int g);
int pipe_set_diam(Pipe *p, float d);
int pipe_set_side(Pipe *p, float s);
int pipe_set_sides(Pipe *p, float s1, float s2);
int pipe_set_custom_sides(Pipe *p, float *s);
void pipe_set_rough(Pipe *p, float r);
void pipe_set_length(Pipe *p, float l);
void pipe_set_id(Pipe *p, int);
void pipe_print(Pipe *p);

void pipe_set_friction(Pipe *p, float fd);
float pipe_get_friction(Pipe *p);
float pipe_compute_friction(Pipe *p, FrictionModel fm);

//Graph functions
void graph_set_diameters(Graph *g, float *d);
void graph_set_roughness(Graph *g, float *r);
void graph_set_lengths(Graph *g, float *l);
void graph_print(Graph *g);
void graph_print_incidence_matrix(Graph *g);
void graph_print_mass_conservation_matrix(Graph *g);
void graph_print_disconnected_nodes(Graph *g);
void graph_print_connected_nodes(Graph *g);
void graph_print_junction_nodes(Graph *g);
void graph_print_input_nodes(Graph *g);
void graph_print_output_nodes(Graph *g);

void graph_compute_mass_conservation_matrix(Graph *g);

Node **graph_get_nodes(Graph *g);
Pipe **graph_get_pipes(Graph *g);

int graph_get_n_nodes(Graph *g);
int graph_get_n_disconnected_nodes(Graph *g);
int graph_get_n_connected_nodes(Graph *g);
int graph_get_n_junction_nodes(Graph *g);
int graph_get_n_input_nodes(Graph *g);
int graph_get_n_output_nodes(Graph *g);

Node *graph_get_nth_node(Graph *g, int);
Node *graph_get_nth_disconnected_node(Graph *g, int);
Node *graph_get_nth_connected_node(Graph *g, int);
Node *graph_get_nth_junction_node(Graph *g, int);
Node *graph_get_nth_input_node(Graph *g, int);
Node *graph_get_nth_output_node(Graph *g, int);

float graph_get_total_outflow(Graph *g);
float graph_get_total_calculated_outflow(Graph *g);
void graph_set_inflow_evenly(Graph *g);
float graph_get_total_inflow(Graph *g);
float graph_get_total_calculated_inflow(Graph *g);

void graph_inflow_real_to_calc(Graph *g);
void graph_inflow_calc_to_real(Graph *g);
void graph_outflow_real_to_calc(Graph *g);
void graph_outflow_calc_to_real(Graph *g);

void graph_set_fluid_viscosity(Graph *g, float v);
void graph_set_fluid_density(Graph *g, float d);
float graph_get_fluid_viscosity(Graph *g);
float graph_get_fluid_density(Graph *g);
void graph_set_friction_model(Graph *g, FrictionModel fm);
void graph_backpropagate_flowrate(Graph *g);
void graph_propagate_pressure(Graph *g);

_Bool graph_detect_leaks(Graph *g);

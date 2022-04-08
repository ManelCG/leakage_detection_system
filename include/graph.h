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
void node_print(Node *n);

void node_set_flowrate(Node *n, float f);
float node_get_flowrate(Node *n);
void node_set_flowrate_measured(Node *n, _Bool m);


//Pipe functions
void pipe_set_geometry(Pipe *p, int g);
int pipe_set_diam(Pipe *p, float d);
int pipe_set_side(Pipe *p, float s);
int pipe_set_sides(Pipe *p, float s1, float s2);
int pipe_set_custom_sides(Pipe *p, float *s);
void pipe_set_rough(Pipe *p, float r);
void pipe_set_id(Pipe *p, int);
void pipe_print(Pipe *p);

//Graph functions
void graph_set_diameters(Graph *g, float *d);
void graph_set_roughness(Graph *g, float *r);
void graph_print(Graph *g);
void graph_print_incidence_matrix(Graph *g);
void graph_print_disconnected_nodes(Graph *g);
void graph_print_connected_nodes(Graph *g);
void graph_print_junction_nodes(Graph *g);
void graph_print_input_nodes(Graph *g);
void graph_print_output_nodes(Graph *g);

int graph_get_n_disconnected_nodes(Graph *g);
int graph_get_n_connected_nodes(Graph *g);
int graph_get_n_junction_nodes(Graph *g);
int graph_get_n_input_nodes(Graph *g);
int graph_get_n_output_nodes(Graph *g);

Node *graph_get_nth_disconnected_node(Graph *g, int);
Node *graph_get_nth_connected_node(Graph *g, int);
Node *graph_get_nth_junction_node(Graph *g, int);
Node *graph_get_nth_input_node(Graph *g, int);
Node *graph_get_nth_output_node(Graph *g, int);

float graph_get_total_outflow(Graph *g);
void graph_set_inflow_evenly(Graph *g);

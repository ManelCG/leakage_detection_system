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

//Graph functions
void graph_set_diameters(Graph *g, float *d);
void graph_set_roughness(Graph *g, float *r);
void graph_print(Graph *g);

//Pipe functions
void pipe_set_diam(Pipe *p, float d);
void pipe_set_rough(Pipe *p, float r);
void pipe_set_id(Pipe *p, int);

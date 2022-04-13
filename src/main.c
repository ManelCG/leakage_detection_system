#include <graph.h>
#include <stdlib.h>

#include <stdio.h>
#include <stdbool.h>

#include <fluid_mechanics.h>

int main(int argc, char *argv[]){
  //INITIAL DATA:
  float water_viscosity = 0.08903; float water_density = 997.08;   //Asume 20ºC

  int   sorig[] = {1, 2, 2,  2, 12, 3, 3, 6, 6,  9,  9, 12};
  int   torig[] = {2, 6, 9, 12,  3, 4, 5, 7, 8, 10, 11, 13};
  float dorig[] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.075, 0.075, 0.075, 0.1, 0.1, 0.075, 0.075};
  // float lorig[] = {1.3 , 19.3, 5.27, 131,  4.20, 0.69,  53.8,  0.420, 5.23,1.3, 0.39,  420};
  float lorig[] = {100 , 100, 100, 100,  100, 100,  100,  100, 100, 100, 100,  100};
  float rough[] = {0.00015, 0.0000015, 0.000015, 0.0000015, 0.0015, 0.0000015, 0.000015, 0.0000015, 0.00015, 0.0000015, 0.000015, 0.0015};
  // float outfw[] = {3.6, 2.7, 0, 23.1, 0.5, 1.8, 3.6};
  float outfw[] = {0.02, 0.03, 0.05, 0.01, 0.02, 0.01, 0.02};
  float height_reservoir = 70;
  int num_pipes = 12;

  //Create graph
  Graph *g = graph_new(NULL, num_pipes, sorig, torig);

  graph_set_fluid_viscosity(g, water_viscosity);
  graph_set_fluid_density(g, water_density);
  graph_set_friction_model(g, friction_model_churchill);

  //Store vectors for sporadic use
  Node **node_v = graph_get_nodes(g);
  // Pipe **pipe_v = graph_get_pipes(g);

  //Set pipe geometries
  graph_set_diameters(g, dorig);
  graph_set_roughness(g, rough);
  graph_set_lengths(g, lorig);

  //Compute mass conservation matrix
  // graph_compute_mass_conservation_matrix(g);

  //Set the outflow measured at each house
  int num = graph_get_n_output_nodes(g);
  printf("%d demand nodes\n", num);
  for (int i = 0; i < num; i++){
    Node *n = graph_get_nth_output_node(g, i);
    node_set_flowrate_measured(n, outfw[i]);
  }

  //Propagate to inflow
  graph_set_inflow_evenly(g);
  graph_inflow_calc_to_real(g);
  graph_outflow_real_to_calc(g);

  //Set inflow parametrers
  node_set_height(node_v[1], height_reservoir); //Node is 25 meters high
  float input_pressure = node_input_compute_pressure(node_v[1]);
  node_set_pressure_calculated(node_v[1], input_pressure);
  // node_set_fluid_velocity(node_v[1], velocity_reservoir);


  ////Compute friction
  graph_backpropagate_flowrate(g);
  graph_propagate_pressure(g);


  //View result
  graph_print(g);


  //Check if there are leaks
  if (! graph_detect_leaks(g)){
    printf("No leaks detected\n");
  } else {
    printf("Leaks detected!!!\n");
  }
}

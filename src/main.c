#include <graph.h>
#include <stdlib.h>

#include <stdio.h>
#include <stdbool.h>

#include <fluid_mechanics.h>

#include <lodepng.h>

#include <graphviz/cgraph.h>

int main(int argc, char *argv[]){
  //INITIAL DATA:
  float water_viscosity = 0.08903; float water_density = 997.08;   //Asume 20ºC
  // float water_viscosity = 1; float water_density = 1;   //Asume 20ºC

  int   sorig[] = {1, 2, 2,  2, 12, 3,
                   3, 6, 6,  9,  9, 12};
  int   torig[] = {2, 6, 9, 12,  3, 4,
                   5, 7, 8, 10, 11, 13};
  float dorig[] = {0.15,  0.15, 0.15, 0.15, 0.15, 0.075,
                   0.075, 0.075, 0.1, 0.1, 0.075, 0.075};
  // float lorig[] = {1.3 , 19.3, 5.27, 131,  4.20, 0.69,  53.8,  0.420, 5.23,1.3, 0.39,  420};
  float lorig[] = {100 , 100, 100, 100,  100, 100,
                   100,  100, 100, 100, 100,  100};
  float rough[] = {0.00015, 0.0000015, 0.000015, 0.0000015, 0.0015, 0.0000015,
                   0.000015, 0.0000015, 0.00015, 0.0000015, 0.000015, 0.0015};
  // float outfw[] = {3.6, 2.7, 0, 23.1, 0.5, 1.8, 3.6};
  float outfw[] = {0.02, 0.03, 0.05, 0.01, 0.02, 0.01, 0.02};
  float height_reservoir = 70;
  int num_pipes = 12;

  //Aux data:
  srand(420);

  //Create graph
  Graph *g = graph_new(NULL, num_pipes, sorig, torig);

  //Set physics variables
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
  graph_set_output_flowrates(g, outfw);

  //Propagate to inflow
  graph_set_inflow_evenly(g);
  graph_inflow_calc_to_real(g);
  graph_outflow_real_to_calc(g);

  //Set inflow parametrers
  node_set_height(node_v[1], height_reservoir); //Node is 70 meters high
  float input_pressure = node_input_compute_pressure(node_v[1]);
  node_set_pressure_calculated(node_v[1], input_pressure);


  //Compute friction and pressures
  graph_backpropagate_flowrate(g);
  graph_propagate_pressure(g);


  ////View result
  //printf("GRAPH WITHOUT LEAKS: \n");
  //graph_print(g);
  //printf("\n\n");


  int nleaks = 1;
  Leaks *leaks_gen = graph_generate_random_leaks(g, nleaks);
  if (leaks_gen == NULL){
    printf("Error generating random leaks...\n");
  }

  //Add leak outflow to input node
  graph_add_leaks_to_inflow(g);


  ////View result
  //printf("GRAPH WITH LEAKS: \n");
  //graph_print(g);

  //Check if there are leaks
  _Bool has_leaks = graph_has_leaks(g);
  if (! has_leaks){
    printf("No leaks detected\n");
    printf("Nothing to do... Exiting\n");
    return 0;
  }

  printf("Leaks detected!!!\n");
  printf("Will try to optimize measurements...\n");

  Leaks *leaks_calc = graph_find_leaks(g);
  if (leaks_calc == NULL){
    printf("Couldnt find leaks\n");
  }

  graph_plot(g, 1920, 1080);

  graph_destroy(g);
  leaks_destroy(leaks_calc);
}

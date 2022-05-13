#include <fluid_mechanics.h>
#include <math.h>

//FRICTION MODELS:----------------------------

float friction_model_churchill(float diam, float rough, float dens, float visc, float vel){
  float Re = compute_reynolds_number(vel, dens, diam, visc);
  float Ca = pow((-2.457*log(pow((7.0/Re), 0.9) + 0.27*(rough/diam))), 16.0);
  float Cb = pow((37530.0/Re), 16.0);

  float fd = 8*pow(pow(8.0/Re, 12) + pow(Ca + Cb, -1.5) , 1.0/12.0);
  return fd;
}

//Stokes model ignores friction
float friction_model_stokes(float diam, float rough, float dens, float visc, float vel){
  return 0;
}


//Misc
//u = velocity, d = density, L = characteristic linear dimension, v = viscosity
float compute_reynolds_number(float vel, float dens, float lineardim, float visc){
  return ((dens*vel*lineardim)/visc);
}

float calculate_pressure_drop(float vel, float diam, float length, float friction, float dens){
  float drop = friction * length/diam * dens/2 * pow(vel, 2);
  // float drop = (friction * length * dens * pow(vel, 2)) / (2 * 9.806 * diam);
  return drop;
  // return ((pressure_in + drop));
}

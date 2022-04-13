#ifndef FLUID_MECHANICS_H
#define FLUID_MECHANICS_H

//Friction Models:
typedef float (*FrictionModel)(float diam,
                               float rough,
                               float dens,
                               float visc,
                               float vel);

float friction_model_churchill(float diam, float rough, float dens, float visc, float vel);
float friction_model_stokes(float diam, float rough, float dens, float visc, float vel);


float compute_reynolds_number(float u, float d, float l, float v);

float calculate_pressure_drop(float vel, float diam, float length, float friction, float dens);

#endif //FLUID_MECHANICS_H

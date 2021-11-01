/**
 * @file main.cpp
 * @author Mit Bailey (mitbailey99@gmail.com)
 * @brief Satellite simulation program.
 * @version See Git tags for version information.
 * @date 2021.10.31
 * 
 * Program for the simulation of a CubeSat in low Earth orbit. Takes into account the effects of residual magnetorquer
 * dipole moments using the world magnetic model and satellite aerodynamics to calculate satellite disturbance torques.
 * The program then uses this information to produce the orientation of the satellite, represented by quaternions, once
 * active mangetorquer dipole moments are applied towards maintaining a nominal attitude. 
 * 
 * Adapted from a Simulink model using MatLab created by Sunip K. Mukherjee.
 * Adapts code from the World Magnetic Model GEOMAG DRIVER created by Dr. John Quinn, Stefan Maus, and Manoj Nair.
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <stdio.h>
#include "wmm_geomag.h"
#include "meb_debug.h"

int main ()
{
    // LLA in degrees and meters, time in decimal years (i.e. 2021.7).
    double latitude = 42.0;
    double longitude = 76.0;
    double altitude = 410.0;
    double time = 2021.0;

    dbprintlf("INPUTS");
    dbprintlf("===============");
    dbprintlf("LLA:  %.02f deg, %.02f deg, %.04f km", latitude, longitude, altitude);
    dbprintlf("Time: %f", time);

    magfield_t results = wmm_main(latitude, longitude, altitude, time);

    printf("\n");
    dbprintlf("MAGFIELD RESULTS");
    dbprintlf("===============");
    dbprintlf("F: %f nT", results.F);
    dbprintlf("H: %f nT", results.H);
    dbprintlf("X: %f nT", results.X);
    dbprintlf("Y: %f nT", results.Y);
    dbprintlf("Z: %f nT", results.Z);

    return 0;
}
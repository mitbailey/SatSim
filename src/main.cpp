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
#include <math.h>
#include "wmm_geomag.h"
#include "meb_debug.h"
#include "CoordGeodetic.h"
#include "CoordTopocentric.h"
#include "SGP4.h"
#include "satsim.hpp"

#define PERFORM_TESTS

#define APPROX_EQUALS(var, ideal, error) ({(var > ideal - error) && (var < ideal + error);})

int test_wmm();
int test_sgp4();

int main()
{
#ifdef PERFORM_TESTS
    if (test_sgp4() == 1 && test_wmm() == 1)
    {
        dbprintlf(GREEN_FG "Tests passed.");
    }
#endif // PERFORM_TESTS

    // /*** SGP4 Test ***/
    // SGP4 *satellite = new SGP4(Tle(ISS_TLE[0], ISS_TLE[1]));
    // // Observer *ground_station = new Observer(GS_LAT, GS_LON, ELEV);
    // DateTime time_now = DateTime::Now(true);
    // Eci sat_position = satellite->FindPosition(time_now);
    // CoordGeodetic sat_lla = sat_position.ToGeodetic();

    // dbprintlf("Current Position of Satellite");
    // dbprintlf("Latitude: %f", sat_lla.latitude TO_DEGREES);
    // dbprintlf("Longitude: %f", sat_lla.longitude TO_DEGREES);
    // dbprintlf("Altitude: %f", sat_lla.altitude);

    // /*** WMM Test ***/

    // // LLA in degrees and meters, time in decimal years (i.e. 2021.7).
    // double latitude = 42.0;
    // double longitude = 76.0;
    // double altitude = 410.0;
    // double time = 2021.0;

    // dbprintlf("INPUTS");
    // dbprintlf("===============");
    // dbprintlf("LLA:  %.02f deg, %.02f deg, %.04f km", latitude, longitude, altitude);
    // dbprintlf("Time: %f", time);

    // magfield_t results = wmm_main(latitude, longitude, altitude, time);

    // printf("\n");
    // dbprintlf("MAGFIELD RESULTS");
    // dbprintlf("===============");
    // dbprintlf("F: %f nT", results.F);
    // dbprintlf("H: %f nT", results.H);
    // dbprintlf("X: %f nT", results.X);
    // dbprintlf("Y: %f nT", results.Y);
    // dbprintlf("Z: %f nT", results.Z);

    return 0;
}

int test_sgp4()
{
    const char Test_TLE[2][70] = {"1 25544U 98067A   21305.49174769  .00006352  00000-0  12389-3 0  9990",
                                  "2 25544  51.6454  23.8015 0003616 161.6486  16.2104 15.48875649309872"};

    SGP4 *satellite = new SGP4(Tle(ISS_TLE[0], ISS_TLE[1]));
    DateTime time_now = DateTime::Now(true);
    Eci sat_position = satellite->FindPosition(time_now);
    CoordGeodetic sat_lla = sat_position.ToGeodetic();

    dbprintlf("* SGP4 UNIT TEST *");
    dbprintlf("Is this the current LLA of the ISS?");
    dbprintlf("http://www.satflare.com/track.asp#TOP");
    dbprintlf("LAT (deg): %.02f", sat_lla.latitude TO_DEGREES);
    dbprintlf("LON (deg): %.02f", sat_lla.longitude TO_DEGREES);
    dbprintlf("ALT (km):  %.04f", sat_lla.altitude);

    return 1;
}

int test_wmm()
{
    /*** WMM Test ***/

    // LLA in degrees and meters, time in decimal years (i.e. 2021.7).
    double latitude = 42.0;
    double longitude = 76.0;
    double altitude = 410.0;
    double time = 2021.0;

    magfield_t results = wmm_main(latitude, longitude, altitude, time);

    if (APPROX_EQUALS(results.X, 21009.499560, 0.00001)
        && APPROX_EQUALS(results.Y, 1429.307425, 0.00001)
        && APPROX_EQUALS(results.Z, 39479.132185, 0.00001))
    {
        return 1;
    }
    else
    {
        printf("\n");
        dbprintlf("EXPECTED RESULTS");
        dbprintlf("===============");
        dbprintlf("X: 21009.499560 nT");
        dbprintlf("Y: 1429.307425 nT");
        dbprintlf("Z: 39479.132185 nT");

        printf("\n");
        dbprintlf("ACTUAL RESULTS");
        dbprintlf("===============");
        dbprintlf("X: %f nT", results.X);
        dbprintlf("Y: %f nT", results.Y);
        dbprintlf("Z: %f nT", results.Z);
        return 0;
    }
}
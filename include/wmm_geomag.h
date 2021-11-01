/**
 * @file wmm_geomag.h
 * @author Mit Bailey (mitbailey99@gmail.com)
 * @brief 
 * @version See Git tags for version information.
 * @date 2021.10.31
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef WMM_GEOMAG_H
#define WMM_GEOMAG_H

#define SUPPRESS_RESULTS
#define NaN log(-1.0)

typedef struct
{
    double F;
    double H;
    double X;
    double Y;
    double Z;
    double D_deg;
    double D_min;
    double I_deg;
    double I_min;
    double dF;
    double dH;
    double dX;
    double dY;
    double dZ;
    double dD_min;
    double dI_min;
} magfield_t;

int my_isnan(double d);

magfield_t wmm_main(double latitude, double longitude, double altitude, double argtime);

static void E0000(int IENTRY, int *maxdeg, double alt,double glat,double glon, double time, double *dec, double *dip, double *ti, double *gv);

// static void E0000_GEOMAG(int *maxdeg, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);

// static void E0000_GEOMG1(int *maxdeg, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);

void geomag(int *maxdeg);

void geomg1(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);

char geomag_introduction(double epochlowlim);

#endif // WMM_GEOMAG_H
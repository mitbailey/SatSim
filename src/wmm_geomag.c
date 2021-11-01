/**
 * @file wmm_geomag.c
 * @author Mit Bailey (mitbailey99@gmail.com)
 * @brief
 * @version See Git tags for version information.
 * @date 2021.10.31
 *
 * Original by Dr. John Quinn with additions by Stefan Maus and Manoj Nair.
 *
 * @copyright Copyright (c) 2021
 *
 */

/* PROGRAM MAGPOINT (GEOMAG DRIVER) */
/************************************************************************
     Contact Information

     Software and Model Support
        National Geophysical Data Center
        NOAA EGC/2
        325 Broadway
        Boulder, CO 80303 USA
        Attn: Manoj Nair or Stefan Maus
        Phone:  (303) 497-4642 or -6522
        Email:  Manoj.C.Nair@Noaa.gov or Stefan.Maus@noaa.gov
        Web: http://www.ngdc.noaa.gov/geomag/WMM/

     Sponsoring Government Agency
       National Geospatial-Intelligence Agency
           PRG / CSAT, M.S. L-41
           3838 Vogel Road
           Arnold, MO 63010
           Attn: Craig Rollins
           Phone:  (314) 263-4186
           Email:  Craig.M.Rollins@Nga.Mil

      Original Program By:
        Dr. John Quinn
        FLEET PRODUCTS DIVISION, CODE N342
        NAVAL OCEANOGRAPHIC OFFICE (NAVOCEANO)
        STENNIS SPACE CENTER (SSC), MS 39522-5001

        3/25/05 Version 2.0 Stefan Maus corrected 2 bugs:
         - use %c instead of %s for character read
         - help text: positive inclination is downward
        1/29/2010 Version 3.0 Manoj Nair
        Converted floating variables from single precision to double
        Changed : height above AMSL (WGS84) to Height above WGS84 Ellipsoid
        Removed the NaN forcing at the geographic poles
        A new function "my_isnan" for improved portablility

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wmm_geomag.h"
#include "meb_debug.h"

int my_isnan(double d)
{
    return (d != d); /* IEEE: only NaN is not equal to itself */
}

// #define NaN log(-1.0)

// static void E0000(int IENTRY, int *maxdeg, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);
// void geomag(int *maxdeg);
// void geomg1(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);
// char geomag_introduction(double epochlowlim);

// Degrees, degrees, kilometers, decimal years.
magfield_t wmm_main(double latitude, double longitude, double altitude, double argtime)
{
    int warn_H, warn_H_strong, warn_P;
    static int maxdeg;
    static double altm, dlat, dlon;
    static double ati, adec, adip;
    static double alt, time, dec, dip, ti, gv;
    static double time1, dec1, dip1, ti1;
    static double time2, dec2, dip2, ti2;
    char answer, ans;
    char decd[7], dipd[7], d_str[81], modl[20];
    char goodbye[81];
    double x1, x2, y1, y2, z1, z2, h1, h2;
    double ax, ay, az, ah;
    double rTd = 0.017453292;
    double epochlowlim, epochuplim;
    double epochrange = 5.0;
    double warn_H_val, warn_H_strong_val;
    double dmin, imin, ddeg, ideg;
    FILE *wmmtemp;

    wmmtemp = fopen("WMM/WMM.COF", "r");
    if (wmmtemp == NULL)
    {
        fprintf(stderr, "Error opening model file WMM.COF\n");
        exit(1);
    }

    fgets(d_str, 80, wmmtemp);
    if (sscanf(d_str, "%lf%s", &epochlowlim, modl) < 2)
    {
        fprintf(stderr, "Invalid header in model file WMM.COF\n");
        exit(1);
    }

    fclose(wmmtemp);
    strcpy(goodbye, "\n -- End of WMM Point Calculation Program -- \n\n");
    // ans = geomag_introduction(epochlowlim);

    maxdeg = 12;
    warn_H = 0;
    warn_H_val = 99999.0;
    warn_H_strong = 0;
    warn_H_strong_val = 99999.0;
    warn_P = 0;

    geomag(&maxdeg);

    // printf("\n\n\nENTER LATITUDE IN DECIMAL DEGREES ");
    // printf("\n(North latitude positive, South latitude negative \n");
    // printf("i.e. 25.5 for 25 degrees 30 minutes north.) \n");
    // scanf("%lf%*[^\n]", &dlat);
    dlat = latitude;
    // getchar();

    // printf("ENTER LONGITUDE IN DECIMAL DEGREES");
    // printf("(East longitude positive, West negative \n");
    // printf("i.e.- 100.0 for 100.0 degrees west.)\n");
    // scanf("%lf%*[^\n]", &dlon);
    dlon = longitude;
    // getchar();

    // printf("ENTER ALTITUDE IN KILOMETERS ABOVE WGS84 ELLIPSOID\n");
    // scanf("%lf%*[^\n]", &altm);
    altm = altitude;
    // getchar();
    alt = altm;

    epochuplim = epochlowlim + epochrange;
    // printf("ENTER TIME IN DECIMAL YEAR (%-7.2lf - %-7.2lf)\n", epochlowlim, epochuplim);
    // scanf("%lf%*[^\n]", &time);
    time = argtime;
    // getchar();

    geomg1(alt, dlat, dlon, time, &dec, &dip, &ti, &gv);
    time1 = time;
    dec1 = dec;
    dip1 = dip;
    ti1 = ti;
    time = time1 + 1.0;

    /*COMPUTE X, Y, Z, AND H COMPONENTS OF THE MAGNETIC FIELD*/

    x1 = ti1 * (cos((dec1 * rTd)) * cos((dip1 * rTd)));
    y1 = ti1 * (cos((dip1 * rTd)) * sin((dec1 * rTd)));
    z1 = ti1 * (sin((dip1 * rTd)));

    magfield_t results = {0};
    results.X = x1;
    results.Y = y1;
    results.Z = z1;

    return results;
}
/*************************************************************************/

static void E0000(int IENTRY, int *maxdeg, double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv)
{
    static int maxord, i, icomp, n, m, j, D1, D2, D3, D4;
    static double c[13][13], cd[13][13], tc[13][13], dp[13][13], snorm[169],
        sp[13], cp[13], fn[13], fm[13], pp[13], k[13][13], pi, dtr, a, b, re,
        a2, b2, c2, a4, b4, c4, epoch, gnm, hnm, dgnm, dhnm, flnmj, otime, oalt,
        olat, olon, dt, rlon, rlat, srlon, srlat, crlon, crlat, srlat2,
        crlat2, q, q1, q2, ct, st, r2, r, d, ca, sa, aor, ar, br, bt, bp, bpp,
        par, temp1, temp2, parp, bx, by, bz, bh;
    static char model[20], c_str[81], c_new[5];
    static double *p = snorm;
    char answer;

    FILE *wmmdat;

    switch (IENTRY)
    {
    case 0:
        goto GEOMAG;
    case 1:
        goto GEOMG1;
    }

GEOMAG:
    wmmdat = fopen("WMM/WMM.COF", "r");
    if (wmmdat == NULL)
    {
        fprintf(stderr, "Error opening model file WMM.COF\n");
        exit(1);
    }

    /* INITIALIZE CONSTANTS */
    maxord = *maxdeg;
    sp[0] = 0.0;
    cp[0] = *p = pp[0] = 1.0;
    dp[0][0] = 0.0;
    a = 6378.137;
    b = 6356.7523142;
    re = 6371.2;
    a2 = a * a;
    b2 = b * b;
    c2 = a2 - b2;
    a4 = a2 * a2;
    b4 = b2 * b2;
    c4 = a4 - b4;

    /* READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS */
    c[0][0] = 0.0;
    cd[0][0] = 0.0;

    fgets(c_str, 80, wmmdat);
    if (sscanf(c_str, "%lf%s", &epoch, model) < 2)
    {
        fprintf(stderr, "Invalid header in model file WMM.COF\n");
        exit(1);
    }

S3:
    if (fgets(c_str, 80, wmmdat) == NULL)
        goto S4;

    /* CHECK FOR LAST LINE IN FILE */
    for (i = 0; i < 4 && (c_str[i] != '\0'); i++)
    {
        c_new[i] = c_str[i];
        c_new[i + 1] = '\0';
    }
    icomp = strcmp("9999", c_new);
    if (icomp == 0)
        goto S4;
    /* END OF FILE NOT ENCOUNTERED, GET VALUES */
    sscanf(c_str, "%d%d%lf%lf%lf%lf", &n, &m, &gnm, &hnm, &dgnm, &dhnm);

    if (n > maxord)
        goto S4;
    if (m > n || m < 0.0)
    {
        fprintf(stderr, "Corrupt record in model file WMM.COF\n");
        exit(1);
    }

    if (m <= n)
    {
        c[m][n] = gnm;
        cd[m][n] = dgnm;
        if (m != 0)
        {
            c[n][m - 1] = hnm;
            cd[n][m - 1] = dhnm;
        }
    }
    goto S3;

    /* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
S4:
    *snorm = 1.0;
    fm[0] = 0.0;
    for (n = 1; n <= maxord; n++)
    {
        *(snorm + n) = *(snorm + n - 1) * (double)(2 * n - 1) / (double)n;
        j = 2;
        for (m = 0, D1 = 1, D2 = (n - m + D1) / D1; D2 > 0; D2--, m += D1)
        {
            k[m][n] = (double)(((n - 1) * (n - 1)) - (m * m)) / (double)((2 * n - 1) * (2 * n - 3));
            if (m > 0)
            {
                flnmj = (double)((n - m + 1) * j) / (double)(n + m);
                *(snorm + n + m * 13) = *(snorm + n + (m - 1) * 13) * sqrt(flnmj);
                j = 1;
                c[n][m - 1] = *(snorm + n + m * 13) * c[n][m - 1];
                cd[n][m - 1] = *(snorm + n + m * 13) * cd[n][m - 1];
            }
            c[m][n] = *(snorm + n + m * 13) * c[m][n];
            cd[m][n] = *(snorm + n + m * 13) * cd[m][n];
        }
        fn[n] = (double)(n + 1);
        fm[n] = (double)n;
    }
    k[1][1] = 0.0;

    otime = oalt = olat = olon = -1000.0;
    fclose(wmmdat);
    return;

    /*************************************************************************/

GEOMG1:

    dt = time - epoch;
    if (otime < 0.0 && (dt < 0.0 || dt > 5.0))
    {
        printf("\n\n WARNING - TIME EXTENDS BEYOND MODEL 5-YEAR LIFE SPAN");
        printf("\n CONTACT NGDC FOR PRODUCT UPDATES:");
        printf("\n         National Geophysical Data Center");
        printf("\n         NOAA EGC/2");
        printf("\n         325 Broadway");
        printf("\n         Boulder, CO 80303 USA");
        printf("\n         Attn: Manoj Nair or Stefan Maus");
        printf("\n         Phone:  (303) 497-4642 or -6522");
        printf("\n         Email:  Manoj.C.Nair@Noaa.Gov");
        printf("\n         or");
        printf("\n         Stefan.Maus@noaa.gov");
        printf("\n         Web: http://www.ngdc.noaa.gov/geomag/WMM/");
        printf("\n\n EPOCH  = %.3lf", epoch);
        printf("\n TIME   = %.3lf", time);
        printf("\n Do you wish to continue? (y or n) ");
        scanf("%c%*[^\n]", &answer);
        getchar();
        if ((answer == 'n') || (answer == 'N'))
        {
            printf("\n Do you wish to enter more point data? (y or n) ");
            scanf("%c%*[^\n]", &answer);
            getchar();
            if ((answer == 'y') || (answer == 'Y'))
                goto GEOMG1;
            else
                exit(0);
        }
    }

    pi = 3.14159265359;
    dtr = pi / 180.0;
    rlon = glon * dtr;
    rlat = glat * dtr;
    srlon = sin(rlon);
    srlat = sin(rlat);
    crlon = cos(rlon);
    crlat = cos(rlat);
    srlat2 = srlat * srlat;
    crlat2 = crlat * crlat;
    sp[1] = srlon;
    cp[1] = crlon;

    /* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
    if (alt != oalt || glat != olat)
    {
        q = sqrt(a2 - c2 * srlat2);
        q1 = alt * q;
        q2 = ((q1 + a2) / (q1 + b2)) * ((q1 + a2) / (q1 + b2));
        ct = srlat / sqrt(q2 * crlat2 + srlat2);
        st = sqrt(1.0 - (ct * ct));
        r2 = (alt * alt) + 2.0 * q1 + (a4 - c4 * srlat2) / (q * q);
        r = sqrt(r2);
        d = sqrt(a2 * crlat2 + b2 * srlat2);
        ca = (alt + d) / r;
        sa = c2 * crlat * srlat / (r * d);
    }
    if (glon != olon)
    {
        for (m = 2; m <= maxord; m++)
        {
            sp[m] = sp[1] * cp[m - 1] + cp[1] * sp[m - 1];
            cp[m] = cp[1] * cp[m - 1] - sp[1] * sp[m - 1];
        }
    }
    aor = re / r;
    ar = aor * aor;
    br = bt = bp = bpp = 0.0;
    for (n = 1; n <= maxord; n++)
    {
        ar = ar * aor;
        for (m = 0, D3 = 1, D4 = (n + m + D3) / D3; D4 > 0; D4--, m += D3)
        {
            /*
               COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
               AND DERIVATIVES VIA RECURSION RELATIONS
            */
            if (alt != oalt || glat != olat)
            {
                if (n == m)
                {
                    *(p + n + m * 13) = st * *(p + n - 1 + (m - 1) * 13);
                    dp[m][n] = st * dp[m - 1][n - 1] + ct * *(p + n - 1 + (m - 1) * 13);
                    goto S50;
                }
                if (n == 1 && m == 0)
                {
                    *(p + n + m * 13) = ct * *(p + n - 1 + m * 13);
                    dp[m][n] = ct * dp[m][n - 1] - st * *(p + n - 1 + m * 13);
                    goto S50;
                }
                if (n > 1 && n != m)
                {
                    if (m > n - 2)
                        *(p + n - 2 + m * 13) = 0.0;
                    if (m > n - 2)
                        dp[m][n - 2] = 0.0;
                    *(p + n + m * 13) = ct * *(p + n - 1 + m * 13) - k[m][n] * *(p + n - 2 + m * 13);
                    dp[m][n] = ct * dp[m][n - 1] - st * *(p + n - 1 + m * 13) - k[m][n] * dp[m][n - 2];
                }
            }
        S50:
            /*
                TIME ADJUST THE GAUSS COEFFICIENTS
            */
            if (time != otime)
            {
                tc[m][n] = c[m][n] + dt * cd[m][n];
                if (m != 0)
                    tc[n][m - 1] = c[n][m - 1] + dt * cd[n][m - 1];
            }
            /*
                ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
            */
            par = ar * *(p + n + m * 13);
            if (m == 0)
            {
                temp1 = tc[m][n] * cp[m];
                temp2 = tc[m][n] * sp[m];
            }
            else
            {
                temp1 = tc[m][n] * cp[m] + tc[n][m - 1] * sp[m];
                temp2 = tc[m][n] * sp[m] - tc[n][m - 1] * cp[m];
            }
            bt = bt - ar * temp1 * dp[m][n];
            bp += (fm[m] * temp2 * par);
            br += (fn[n] * temp1 * par);
            /*
                SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
            */
            if (st == 0.0 && m == 1)
            {
                if (n == 1)
                    pp[n] = pp[n - 1];
                else
                    pp[n] = ct * pp[n - 1] - k[m][n] * pp[n - 2];
                parp = ar * pp[n];
                bpp += (fm[m] * temp2 * parp);
            }
        }
    }
    if (st == 0.0)
        bp = bpp;
    else
        bp /= st;
    /*
        ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
        GEODETIC COORDINATES
    */
    bx = -bt * ca - br * sa;
    by = bp;
    bz = bt * sa - br * ca;
    /*
        COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
        TOTAL INTENSITY (TI)
    */
    bh = sqrt((bx * bx) + (by * by));
    *ti = sqrt((bh * bh) + (bz * bz));
    *dec = atan2(by, bx) / dtr;
    *dip = atan2(bz, bh) / dtr;
    /*
        COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
        GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
        (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

        OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
    */
    *gv = -999.0;
    if (fabs(glat) >= 55.)
    {
        if (glat > 0.0 && glon >= 0.0)
            *gv = *dec - glon;
        if (glat > 0.0 && glon < 0.0)
            *gv = *dec + fabs(glon);
        if (glat < 0.0 && glon >= 0.0)
            *gv = *dec + glon;
        if (glat < 0.0 && glon < 0.0)
            *gv = *dec - fabs(glon);
        if (*gv > +180.0)
            *gv -= 360.0;
        if (*gv < -180.0)
            *gv += 360.0;
    }
    otime = time;
    oalt = alt;
    olat = glat;
    olon = glon;
    return;
}

/*************************************************************************/

void geomag(int *maxdeg)
{
    E0000(0, maxdeg, 0.0, 0.0, 0.0, 0.0, NULL, NULL, NULL, NULL);
}

/*************************************************************************/

void geomg1(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv)
{
    E0000(1, NULL, alt, glat, glon, time, dec, dip, ti, gv);
}

/*************************************************************************/

char geomag_introduction(double epochlowlim)
{
    char help;
    static char ans;

    printf("\n\n Welcome to the World Magnetic Model (WMM) %4.0lf C-Program\n\n", epochlowlim);
    printf("            --- Version 3.0, January 2010 ---\n\n");
    printf("\n This program estimates the strength and direction of ");
    printf("\n Earth's main magnetic field for a given point/area.");
    printf("\n Enter h for help and contact information or c to continue.");
    printf("\n >");
    scanf("%c%*[^\n]", &help);
    getchar();

    if ((help == 'h') || (help == 'H'))
    {
        printf("\n Help information ");

        printf("\n The World Magnetic Model (WMM) for %7.2lf", epochlowlim);
        printf("\n is a model of Earth's main magnetic field.  The WMM");
        printf("\n is recomputed every five (5) years, in years divisible by ");
        printf("\n five (i.e. 2010, 2015).  See the contact information below");
        printf("\n to obtain more information on the WMM and associated software.");
        printf("\n ");
        printf("\n Input required is the location in geodetic latitude and");
        printf("\n longitude (positive for northern latitudes and eastern ");
        printf("\n longitudes), geodetic altitude in meters, and the date of ");
        printf("\n interest in years.");

        printf("\n\n\n The program computes the estimated magnetic Declination");
        printf("\n (D) which is sometimes called MAGVAR, Inclination (I), Total");
        printf("\n Intensity (F or TI), Horizontal Intensity (H or HI), Vertical");
        printf("\n Intensity (Z), and Grid Variation (GV). Declination and Grid");
        printf("\n Variation are measured in units of degrees and are considered");
        printf("\n positive when east or north.  Inclination is measured in units");
        printf("\n of degrees and is considered positive when pointing down (into");
        printf("\n the Earth).  The WMM is reference to the WGS-84 ellipsoid and");
        printf("\n is valid for 5 years after the base epoch.");

        printf("\n\n\n It is very important to note that a  degree and  order 12 model,");
        printf("\n such as WMM, describes only the long  wavelength spatial magnetic ");
        printf("\n fluctuations due to  Earth's core.  Not included in the WMM series");
        printf("\n models are intermediate and short wavelength spatial fluctuations ");
        printf("\n that originate in Earth's mantle and crust. Consequently, isolated");
        printf("\n angular errors at various  positions on the surface (primarily over");
        printf("\n land, incontinental margins and  over oceanic seamounts, ridges and");
        printf("\n trenches) of several degrees may be expected.  Also not included in");
        printf("\n the model are temporal fluctuations of magnetospheric and ionospheric");
        printf("\n origin. On the days during and immediately following magnetic storms,");
        printf("\n temporal fluctuations can cause substantial deviations of the geomagnetic");
        printf("\n field  from model  values.  If the required  declination accuracy  is");
        printf("\n more stringent than the WMM  series of models provide, the user is");
        printf("\n advised to request special (regional or local) surveys be performed");
        printf("\n and models prepared. Please make requests of this nature to the");
        printf("\n National Geospatial-Intelligence Agency (NGA) at the address below.");

        printf("\n\n\n Contact Information");

        printf("\n  Software and Model Support");
        printf("\n	National Geophysical Data Center");
        printf("\n	NOAA EGC/2");
        printf("\n	325 Broadway");
        printf("\n	Boulder, CO 80303 USA");
        printf("\n	Attn: Susan McLean or Stefan Maus");
        printf("\n	Phone:  (303) 497-6478 or -6522");
        printf("\n	Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov ");

        printf("\n\n\n Continue with program? (y or n) ");
        scanf("%c%*[^\n]", &ans);
        getchar();
    }
    else
    {
        ans = 'y';
    }

    return (ans);
}

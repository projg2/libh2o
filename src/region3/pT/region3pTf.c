/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Specific Volume
 * as a Function of Pressure and Temperature v(p,T)
 * for Region 3 of the IAPWS Industrial Formulation 1997 for the
 * Thermodynamic Properties of Water and Steam */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	-0.251756547792325E-07, +0.601307193668763E-05,
	-0.100615977450049E-02, +0.999969140252192E+00,
	+0.214107759236486E+01, -0.165175571959086E+02,
	-0.141987303638727E-02, +0.269251915156554E+01,

	+0.349741815858722E+02, -0.300208695771783E+02,
	-0.131546288252539E+01, -0.839091277286169E+01,
	+0.181545608337015E-09, -0.591099206478909E-03,
	+0.152115067087106E+01, +0.252956470663225E-04,

	+0.100726265203786E-14, -0.149774533860650E+01,
	-0.793940970562969E-09, -0.150290891264717E-03,
	+0.151205531275133E+01, +0.470942606221652E-05,
	+0.195049710391712E-12, -0.911627886266077E-08,

	+0.604374640201265E-03, -0.225132933900136E-15,
	+0.610916973582981E-11, -0.303063908043404E-06,
	-0.137796070798409E-04, -0.919296736666106E-03,
	+0.639288223132545E-09, +0.753259479898699E-06,

	-0.400321478682929E-12, +0.756140294351614E-08,
	-0.912082054034891E-11, -0.237612381140539E-07,
	+0.269586010591874E-04, -0.732828135157839E-10,
	+0.241995578306660E-09, -0.405735532730322E-03,

	+0.189424143498011E-09, -0.486632965074563E-09
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 10, 12, 14, 16, 18, 20, 22, 24, 28, 32
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8, 9,
	9, 9, 10, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 15, 15,
	16, 17
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3
};

static const int J[] = {
	0,

	6, 7, 8, 9, 10, 11, 8, 10, 11, 12, 9, 10, 4, 7, 9, 6, 2, 10, 3, 5,
	10, 3, 1, 2, 5, 0, 1, 2, 3, 5, 1, 2, 0, 1, 0, 1, 3, 0, 0, 5, 0, 0
};

static const double vstar = 0.0064; /* [m³/kg] */
static const double pstar = 40; /* [MPa] */
static const double Tstar = 730; /* [K] */

double h2o_region3f_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = poly_value(sqrt(pi - 0.587), theta - 0.891,
			I, Ipows, 0, 18, 0,
			J, Jpows, 9, 13, 0,
			n, 42);

	return pow4(sum) * vstar;
}

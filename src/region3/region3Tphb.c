/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Revised Supplementary Release on Backward Equations for the Functions
 * T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial
 * Formulation 1997 for the Thermodynamic Properties of Water and Steam
 * s. 3.3: Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	+0.323254573644920E-4, -0.127575556587181E-3,
	-0.475851877356068E-3, +0.156183014181602E-2,
	+0.105724860113781E+0, -0.858514221132534E+2,
	+0.724140095480911E+3, +0.296475810273257E-2,

	-0.592721983365988E-2, -0.126305422818666E-1,
	-0.115716196364853E+0, +0.849000969739595E+2,
	-0.108602260086615E-1, +0.154304475328851E-1,
	+0.750455441524466E-1, +0.252520973612982E-1,

	-0.602507901232996E-1, -0.307622221350501E+1,
	-0.574011959864879E-1, +0.503471360939849E+1,
	-0.925081888584834E+0, +0.391733882917546E+1,
	-0.773146007130190E+2, +0.949308762098587E+4,

	-0.141043719679409E+7, +0.849166230819026E+7,
	+0.861095729446704E+0, +0.323346442811720E+0,
	+0.873281936020439E+0, -0.436653048526683E+0,
	+0.286596714529479E+0, -0.131778331276228E+0,

	+0.676682064330275E-2
};

static const double Ipows[] = {
	-12, -10, -8, -6, -4, -3, -2, -1, 0, 1, 3, 5, 6, 8
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
	3, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7, 7, 7, 7,
	8, 8, 9, 10, 11, 12, 13
};

static const double Jpows[] = {
	0, 1, 2, 4, 5, 6, 10, 12, 14, 16
};

static const int J[] = {
	0,

	0, 1, 0, 1, 4, 6, 7, 0, 1, 2, 3, 6,
	0, 1, 2, 0, 1, 4, 0, 3, 2, 3, 5, 6, 8, 9,
	0, 2, 1, 1, 1, 1, 1
};

static const double Tstar = 860; /* [K] */
static const double pstar = 100; /* [MPa] */
static const double hstar = 2800; /* [kJ/kg] */

double h2o_region3b_T_ph(double p, double h)
{
	double pi = p / pstar;
	double eta = h / hstar;

	return poly_value(pi + 0.298, eta - 0.720,
			I, Ipows, 8, 14, 0,
			J, Jpows, 0, 10, 0,
			n, 33) * Tstar;
}

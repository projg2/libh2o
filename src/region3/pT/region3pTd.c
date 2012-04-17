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

	-0.452484847171645E-09, +0.315210389538801E-04,
	-0.214991352047545E-02, +0.508058874808345E+03,
	-0.127123036845932E+08, +0.115371133120497E+13,
	-0.197805728776273E-15, +0.241554806033972E-10,

	-0.156481703640525E-05, +0.277211346836625E-02,
	-0.203578994462286E+02, +0.144369489909053E+07,
	-0.411254217946539E+11, +0.623449786243773E-05,
	-0.221774281146038E+02, -0.689315087933158E+05,

	-0.195419525060713E+08, +0.316373510564015E+04,
	+0.224040754426988E+07, -0.436701347922356E-05,
	-0.404213852833996E-03, -0.348153203414663E+03,
	-0.385294213555289E+06, +0.135203700099403E-06,

	+0.134648383271089E-03, +0.125031835351736E+06,
	+0.968123678455841E-01, +0.225660517512438E+03,
	-0.190102435341872E-03, -0.299628410819229E-01,
	+0.500833915372121E-02, +0.387842482998411E+00,

	-0.138535367777182E+04, +0.870745245971773E+00,
	+0.171946252068742E+01, -0.326650121426383E-01,
	+0.498044171727877E+04, +0.551478022765087E-02
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 3
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5,
	5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 11
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16
};

static const int J[] = {
	0,

	4, 6, 7, 9, 10, 12, 0, 2, 4, 6, 8, 9, 11, 3, 7, 8, 9, 6, 8, 1, 2, 5, 7, 0,
	1, 7, 2, 4, 0, 1, 0, 1, 5, 0, 2, 0, 6, 0
};

static const double vstar = 0.0029; /* [m³/kg] */
static const double pstar = 40; /* [MPa] */
static const double Tstar = 690; /* [K] */

double h2o_region3d_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.559, theta - 0.939,
			I, Ipows, 9, 12, 0,
			J, Jpows, 0, 13, 0,
			n, 38);

	return pow4(sum) * vstar;
}

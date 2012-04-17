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

	+0.106905684359136E+01, -0.148620857922333E+01,
	+0.259862256980408E+15, -0.446352055678749E-11,
	-0.566620757170032E-06, -0.235302885736849E-02,
	-0.269226321968839E+00, +0.922024992944392E+01,

	+0.357633505503772E-11, -0.173942565562222E+02,
	+0.700681785556229E-05, -0.267050351075768E-03,
	-0.231779669675624E+01, -0.753533046979752E-12,
	+0.481337131452891E+01, -0.223286270422356E+22,

	-0.118746004987383E-04, +0.646412934136496E-02,
	-0.410588536330937E-09, +0.422739537057241E+20,
	+0.313698180473812E-12, +0.164395334345040E-23,
	-0.339823323754373E-05, -0.135268639905021E-01,

	-0.723252514211625E-14, +0.184386437538366E-08,
	-0.463959533752385E-01, -0.992263100376750E+14,
	+0.688169154439335E-16, -0.222620998452197E-10,
	-0.540843018624083E-07, +0.345570606200257E-02,

	+0.422275800304086E+11, -0.126974478770487E-14,
	+0.927237985153679E-09, +0.612670812016489E-13,
	-0.722693924063497E-11, -0.383669502636822E-03,
	+0.374684572410204E-03, -0.931976897511086E+05,

	-0.247690616026922E-01, +0.658110546759474E+02
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 7, 8, 10, 12, 14, 18, 20, 22, 24, 32, 36
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 9, 9,
	9, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 13, 14, 14, 15, 15,
	16, 16
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 5, 10, 12
};

static const int J[] = {
	0,

	9, 10, 13, 5, 7, 8, 9, 9, 4, 9, 6, 7, 8, 3, 8, 14, 5, 6, 3, 13, 2,
	0, 3, 5, 1, 2, 5, 12, 0, 1, 2, 3, 11, 0, 1, 0, 0, 2, 1, 4, 1, 2
};

static const double vstar = 0.0041; /* [m³/kg] */
static const double pstar = 25; /* [MPa] */
static const double Tstar = 660; /* [K] */

double h2o_region3i_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(sqrt(pi - 0.910), theta - 0.984,
			I, Ipows, 0, 17, 0,
			J, Jpows, 9, 15, 0,
			n, 42);

	return pow4(sum) * vstar;
}

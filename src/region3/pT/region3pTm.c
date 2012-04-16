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
	+0.000000000000000E00,

	+0.811384363481847E00, -0.568199310990094E04,
	-0.178657198172556E11, +0.795537657613427E32,
	-0.814568209346872E05, -0.659774567602874E08,
	-0.152861148659302E11, -0.560165667510446E12,

	+0.458384828593949E06, -0.385754000383848E14,
	+0.453735800004273E08, +0.939454935735563E12,
	+0.266572856432938E28, -0.547578313899097E10,
	+0.200725701112386E15, +0.185007245563239E13,

	+0.185135446828337E09, -0.170451090076385E12,
	+0.157890366037614E15, -0.202530509748774E16,
	+0.368193926183570E60, +0.170215539458936E18,
	+0.639234909918741E42, -0.821698160721956E15,

	-0.795260241872306E24, +0.233415869478510E18,
	-0.600079934586803E23, +0.594584382273384E25,
	+0.189461279349492E40, -0.810093428842645E46,
	+0.188813911076809E22, +0.111052244098768E36,

	+0.291133958602503E46, -0.329421923951460E22,
	-0.137570282536696E26, +0.181508996303902E28,
	-0.346865122768353E30, -0.211961148774260E38,
	-0.128617899887675E49, +0.479817895699239E65
};

static const double Ipows[] = {
		0, 1, 2, 3, 4, 5, 6, 8, 12, 14, 16, 20, 24, 28
};

static const int I[] = {
	0,

	0, 3, 7, 11, 1, 3, 4, 5, 1, 6, 2, 4, 9, 2, 5, 3, 0, 1, 1, 1, 13, 2,
	10, 0, 5, 0, 3, 4, 8, 10, 1, 7, 9, 0, 2, 3, 4, 7, 9, 12
};

static const double Jpows[] = {
	0, 1, 2, 5, 6, 7, 8, 10, 12, 14, 18, 20, 22, 24, 28, 32, 36
};

static const int J[] = {
	0,

	0, 0, 0, 2, 3, 3, 3, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 11, 11,
	12, 12, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 16, 16, 16, 16, 16,
	16, 16
};

static const double vstar = 0.0028; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3m_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(pi - 1.000, sqrt(sqrt(theta - 0.997)),
			I, Ipows, 0, 14, 0,
			J, Jpows, 0, 17, 0,
			n, 40) * vstar;
}

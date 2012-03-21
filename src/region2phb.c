/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6.3.1: The Backward Equation T(p, h) ... */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	+0.14895041079516E+04, +0.74307798314034E+03,
	-0.97708318797837E+02, +0.24742464705674E+01,
	-0.63281320016026E+00, +0.11385952129658E+01,
	-0.47811863648625E+00, +0.85208123431544E-02,
	+0.93747147377932E+00, +0.33593118604916E+01, /* [10] */

	+0.33809355601454E+01, +0.16844539671904E+00,
	+0.73875745236695E+00, -0.47128737436186E+00,
	+0.15020273139707E+00, -0.21764114219750E-02,
	-0.21810755324761E-01, -0.10829784403677E+00,
	-0.46333324635812E-01, +0.71280351959551E-04, /* [20] */

	+0.11032831789999E-03, +0.18955248387902E-03,
	+0.30891541160537E-02, +0.13555504554949E-02,
	+0.28640237477456E-06, -0.10779857357512E-04,
	-0.76462712454814E-04, +0.14052392818316E-04,
	-0.31083814331434E-04, -0.10302738212103E-05, /* [30] */

	+0.28217281635040E-06, +0.12704902271945E-05,
	+0.73803353468292E-07, -0.11030139238909E-07,
	-0.81456365207833E-13, -0.25180545682962E-10,
	-0.17565233969407E-17, +0.86934156344163E-14
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 0, 0, /* [8] */
	1, 1, 1, 1, 1, 1, 1, 1, /* [16] */
	2, 2, 2, 2,
	3, 3, 3, 3, /* [24] */
	4, 4, 4, 4, 4, 4,
	5, 5, 5, /* [33] */
	6,
	7, 7, /* [36] */
	9, 9
};

static const double Jpows[] = {
	0, 1, 2, 6, 8, 12, 18, 24, 28, 40
};

static const int J[] = {
	0,

	0, 1, 2, 5, 6, 7, 8, 9, /* [8] */
	0, 2, 3, 5, 6, 7, 8, 9, /* [16] */
	2, 4, 6, 9,
	1, 2, 5, 7, /* [24] */
	2, 5, 6, 7, 8, 9,
	6, 7, 9, /* [33] */
	8,
	2, 8, /* [36] */
	1, 9
};

static const double hstar = 2000; /* [kJ/kg] */

double h2o_region2b_T_ph(double p, double h) /* [MPa, kJ/kg] -> [K] */
{
	double eta = h / hstar;
	double etaexpr = eta - 2.6;
	double piexpr = p - 2;

	double sum = 0;

	int i;

	double pipowers[10], etapowers[10];

	fill_powers_incr(pipowers, 10, piexpr, 0);
	fill_powers(etapowers, Jpows, 0, 10, etaexpr, 0);

	for (i = 1; i <= 38; ++i)
	{
		double pipow = pipowers[I[i]];
		double etapow = etapowers[J[i]];

		sum += n[i] * pipow * etapow;
	}

	return sum;
}

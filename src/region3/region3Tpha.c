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

	-0.133645667811215E-6, +0.455912656802978E-5,
	-0.146294640700979E-4, +0.639341312970080E-2,
	+0.372783927268847E+3, -0.718654377460447E+4,
	+0.573494752103400E+6, -0.267569329111439E+7,

	-0.334066283302614E-4, -0.245479214069597E-1,
	+0.478087847764996E+2, +0.764664131818904E-5,
	+0.128350627676972E-2, +0.171219081377331E-1,
	-0.851007304583213E+1, -0.136513461629781E-1,

	-0.384460997596657E-5, +0.337423807911655E-2,
	-0.551624873066791E+0, +0.729202277107470E+0,
	-0.992522757376041E-2, -0.119308831407288E+0,
	+0.793929190615421E+0, +0.454270731799386E+0,

	+0.209998591259910E+0, -0.642109823904738E-2,
	-0.235155868604540E-1, +0.252233108341612E-2,
	-0.764885133368119E-2, +0.136176427574291E-1,
	-0.133027883575669E-1
};

static const double Ipows[] = {
	-12, -10, -8, -5, -3, -2, -1, 0, 1, 3, 4, 10, 12
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2,
	3, 4, 5, 5, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 10, 12, 14, 16, 20, 22
};

static const int J[] = {
	0,

	0, 1, 2, 6, 9, 10, 11, 12, 1, 5, 8, 0, 2, 4, 7, 2,
	0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5
};

static const double Tstar = 760; /* [K] */
static const double pstar = 100; /* [MPa] */
static const double hstar = 2300; /* [kJ/kg] */

double h2o_region3a_T_ph(double p, double h)
{
	double pi = p / pstar;
	double eta = h / hstar;

	double piexpr = pi + 0.240;
	double etaexpr = eta - 0.615;

	double sum = 0;

	int i;

	double pipowers[13], etapowers[13];

	fill_powers(pipowers, Ipows, 7, 13, piexpr, 0);
	fill_powers(etapowers, Jpows, 0, 13, etaexpr, 0);

	for (i = 1; i <= 31; ++i)
	{
		double pipow = pipowers[I[i]];
		double etapow = etapowers[J[i]];

		sum += n[i] * pipow * etapow;
	}

	return sum * Tstar;
}

/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations p(h,s) for Region 3, Equations as
 * a Function of h and s for the Region Boundaries, and an Equation Tsat(h,s)
 * for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic
 * Properties of Water and Steam
 * s. 3.3: Backward Equations p(h,s) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E00,

	+0.770889828326934E01, -0.260835009128688E02,
	+0.267416218930389E03, +0.172221089496844E02,
	-0.293542332145970E03, +0.614135601882478E03,
	-0.610562757725674E05, -0.651272251118219E08,

	+0.735919313521937E05, -0.116646505914191E11,
	+0.355267086434461E02, -0.596144543825955E03,
	-0.475842430145708E03, +0.696781965359503E02,
	+0.335674250377312E03, +0.250526809130882E05,

	+0.146997380630766E06, +0.538069315091534E20,
	+0.143619827291346E22, +0.364985866165994E20,
	-0.254741561156775E04, +0.240120197096563E28,
	-0.393847464679496E30, +0.147073407024852E25,

	-0.426391250432059E32, +0.194509340621077E39,
	+0.666212132114896E24, +0.706777016552858E34,
	+0.175563621975576E42, +0.108408607429124E29,
	+0.730872705175151E44, +0.159145847398870E25,

	+0.377121605943324E41
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 14, 18, 20, 22, 24, 28, 32
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3,
	4, 4, 4, 4, 5, 6, 7, 8, 9, 9,
	10, 11, 12, 13, 13, 14, 15, 15, 16, 16
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 10, 14, 16, 22, 24, 28, 32, 36
};

static const int J[] = {
	0,

	0, 1, 5, 0, 3, 4, 7, 9, 6, 10, 0, 2, 3,
	0, 1, 4, 5, 13, 13, 12, 1, 14, 15,
	11, 13, 15, 10, 13, 15, 10, 15, 8, 13
};

static const double pstar = 99; /* [MPa] */
static const double hstar = 2300; /* [kJ/kg] */
static const double sstar = 4.4; /* [kJ/kgK] */

double h2o_region3a_p_hs(double h, double s)
{
	double eta = h / hstar;
	double etaexpr = eta - 1.01;

	double sigma = s / sstar;
	double sigmaexpr = sigma - 0.75;

	double sum = 0;

	int i;

	double etapowers[17], sigmapowers[16];

	fill_powers(etapowers, Ipows, 0, 17, etaexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 16, sigmaexpr, 0);

	for (i = 1; i <= 33; ++i)
	{
		double etapow = etapowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * etapow * sigmapow;
	}

	return sum * pstar;
}

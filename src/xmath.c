/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "xmath.h"

#include <assert.h>
#include <math.h>

static void fill_powers(double* powers, const double* exponents,
		int zeropos, int count, double expr, int deriv)
{
	int i;
	double tmp = expr;

	powers[zeropos + deriv] = 1;
	powers[zeropos + deriv + 1] = tmp;

	for (i = zeropos + deriv + 2; i < count; ++i)
	{
		if (exponents[i] - 1 == exponents[i - 1])
			tmp *= expr;
		else
			tmp = pow(expr, exponents[i] - deriv);
		powers[i] = tmp;
	}

	tmp = 1.0;
	for (i = zeropos + deriv - 1; i >= 0; --i)
	{
		if (exponents[i] + 1 == exponents[i + 1])
			tmp /= expr;
		else
			tmp = pow(expr, exponents[i] - deriv);
		powers[i] = tmp;
	}
}

double poly_value(double x,
		int min, int max, int deriv,
		const double n[])
{
	double sum;

	const int maxn = 1 - min + max;

	int i, j;

	sum = n[maxn];

	for (j = deriv; j > 0; --j)
		sum *= maxn + min - j;

	for (i = maxn - 1; i >= 1; --i)
	{
		double coeff = n[i];

		for (j = deriv; j > 0; --j)
			coeff *= i + min - j;

		sum *= x;
		sum += coeff;
	}

	if (min - deriv != 0)
		sum *= pow(x, min - deriv);

	return sum;
}

double twoarg_poly_value(double x1, double x2,
		const int I[], const double Ipows[], int Ipowzero,
		int Ipowlen, int x1der,
		const int J[], const double Jpows[], int Jpowzero,
		int Jpowlen, int x2der,
		const double n[], int nlen)
{
	double sum = 0;

	int i;

	double x1powers[20], x2powers[34];

	assert(x1der >= 0 && x1der <= 2);
	assert(x2der >= 0 && x2der <= 2);
	assert(Ipowlen <= 20);
	assert(Jpowlen <= 34);

	fill_powers(x1powers, Ipows, Ipowzero, Ipowlen, x1, x1der);
	fill_powers(x2powers, Jpows, Jpowzero, Jpowlen, x2, x2der);

	for (i = 1; i <= nlen; ++i)
	{
		double x1pow = x1powers[I[i]];
		double x2pow = x2powers[J[i]];

		double memb = n[i] * x1pow * x2pow;
		if (x1der >= 1)
			memb *= Ipows[I[i]];
		if (x1der == 2)
			memb *= Ipows[I[i]] - 1;
		if (x2der >= 1)
			memb *= Jpows[J[i]];
		if (x2der == 2)
			memb *= Jpows[J[i]] - 1;

		sum += memb;
	}

	return sum;
}

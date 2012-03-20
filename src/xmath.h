/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_XMATH_H
#define _H2O_XMATH_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

/* using pow(arg, N) is basically equivalent
 * but it will be optimized out only with -ffast-math;
 * this way, optimization always takes place. */

static inline double pow2(double arg)
{
	return arg * arg;
}

static inline double pow4(double arg)
{
	return pow2(pow2(arg));
}

static inline double quadr_value(double a, double b, double c, double x)
{
	return (a * x + b) * x + c;
}

static inline void fill_powers(double* powers, const double* exponents,
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

static inline void fill_powers_incr(double* powers, int count, double expr)
{
	int i;
	double tmp = expr;

	powers[0] = 1;
	powers[1] = tmp;

	for (i = 2; i < count; ++i)
	{
		tmp *= expr;
		powers[i] = tmp;
	}
}

static inline void fill_powers_decr(double* powers, int min, double expr)
{
	int i;
	double tmp = 1.0;

	for (i = -1; i >= min; --i)
	{
		tmp /= expr;
		powers[i] = tmp;
	}
}

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_XMATH_H*/

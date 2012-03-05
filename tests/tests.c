/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdio.h>
#include <math.h>

#include "tests.h"

int exit_status;

void _check(double result, double expected, double precision, const char* call)
{
	double difference = fabs(expected - result);

	if (difference >= precision)
	{
		fprintf(stderr, "[FAIL] %s = %.9e, while %.9e expected.\n",
				call, result, expected);
		exit_status++;
	}
	else
		fprintf(stderr, "[ OK ] %s = %.9e.\n",
				call, result);
}

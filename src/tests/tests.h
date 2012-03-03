/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

extern int exit_status;

void _check(double result, double expected, double precision, const char* call);

#define CHECK(call, expected, precision) _check(call, expected, precision, #call)

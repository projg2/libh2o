/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_COMMON_H
#define _H2O_COMMON_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

enum h2o_region
{
	H2O_REGION_OUT_OF_RANGE,

	H2O_REGION1,
	H2O_REGION2,
	H2O_REGION3,
	H2O_REGION4,
	H2O_REGION5,

	H2O_REGION_MAX
};

double h2o_b23_p_T(double T); /* [K] -> [MPa] */
double h2o_b23_T_p(double p); /* [MPa] -> [K] */

enum h2o_region h2o_region_pT(double p, double T); /* [MPa, K] */

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_COMMON_H*/

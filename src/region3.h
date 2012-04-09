/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION3_H
#define _H2O_REGION3_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region3_p_rhoT(double rho, double T);
double h2o_region3_u_rhoT(double rho, double T);
double h2o_region3_s_rhoT(double rho, double T);
double h2o_region3_h_rhoT(double rho, double T);

/* special use functions */

enum h2o_region3_subregion
{
	H2O_REGION3A,
	H2O_REGION3B,

	H2O_REGION3_MAX
};

enum h2o_region3_subregion h2o_region3_subregion_ph(double p, double h);
enum h2o_region3_subregion h2o_region3_subregion_ps(double p, double s);

double h2o_region3_b3ab_h_p(double p);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION3_H*/

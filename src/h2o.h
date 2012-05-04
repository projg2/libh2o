/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_H2O_H
#define _H2O_H2O_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#include "boundaries.h"

typedef struct
{
	enum h2o_region region;

/*private:*/
	double _arg1, _arg2;
} h2o_t;

h2o_t h2o_new_pT(double p, double T);
h2o_t h2o_new_Tx(double T, double x);
h2o_t h2o_new_px(double p, double x);
h2o_t h2o_new_ph(double p, double h);
h2o_t h2o_new_ps(double p, double s);
h2o_t h2o_new_hs(double h, double s);
h2o_t h2o_new_rhoT(double rho, double T);

double h2o_get_p(const h2o_t state);
double h2o_get_T(const h2o_t state);
double h2o_get_x(const h2o_t state);
double h2o_get_rho(const h2o_t state);

double h2o_get_v(const h2o_t state);
double h2o_get_u(const h2o_t state);
double h2o_get_h(const h2o_t state);
double h2o_get_s(const h2o_t state);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_H2O_H*/

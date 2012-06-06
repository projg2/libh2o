/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION4_H
#define _H2O_REGION4_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region4_p_T(double T);
double h2o_region4_T_p(double p);

double h2o_region4_v_Tx(double T, double x);
double h2o_region4_u_Tx(double T, double x);
double h2o_region4_s_Tx(double T, double x);
double h2o_region4_h_Tx(double T, double x);
double h2o_region4_cp_Tx(double T, double x);
double h2o_region4_cv_Tx(double T, double x);
double h2o_region4_w_Tx(double T, double x);

double h2o_region4_x_Ts(double T, double s);
double h2o_region4_x_Th(double T, double h);
double h2o_region4_T_hs(double h, double s);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION4_H*/

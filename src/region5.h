/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION5_H
#define _H2O_REGION5_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region5_v_pT(double p, double T);
double h2o_region5_u_pT(double p, double T);
double h2o_region5_s_pT(double p, double T);
double h2o_region5_h_pT(double p, double T);
double h2o_region5_cp_pT(double p, double T);
double h2o_region5_cv_pT(double p, double T);
double h2o_region5_w_pT(double p, double T);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION5_H*/

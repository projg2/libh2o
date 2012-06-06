/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#pragma once

#ifndef _H2O_REGION2_H
#define _H2O_REGION2_H 1

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

double h2o_region2_v_pT(double p, double T);
double h2o_region2_u_pT(double p, double T);
double h2o_region2_s_pT(double p, double T);
double h2o_region2_h_pT(double p, double T);
double h2o_region2_cp_pT(double p, double T);
double h2o_region2_cv_pT(double p, double T);
double h2o_region2_w_pT(double p, double T);

double h2o_region2_T_ph(double p, double h);
double h2o_region2_T_ps(double p, double s);
double h2o_region2_p_hs(double h, double s);

/* meta-stable vapor region functions */

double h2o_region2_meta_v_pT(double p, double T);
double h2o_region2_meta_u_pT(double p, double T);
double h2o_region2_meta_s_pT(double p, double T);
double h2o_region2_meta_h_pT(double p, double T);
double h2o_region2_meta_cp_pT(double p, double T);
double h2o_region2_meta_cv_pT(double p, double T);
double h2o_region2_meta_w_pT(double p, double T);

/* special use functions */

enum h2o_region2_subregion
{
	H2O_REGION2A,
	H2O_REGION2B,
	H2O_REGION2C,

	H2O_REGION2_MAX
};

enum h2o_region2_subregion
	h2o_region2_subregion_ph(double p, double h);
enum h2o_region2_subregion
	h2o_region2_subregion_ps(double p, double s);
enum h2o_region2_subregion
	h2o_region2_subregion_hs(double h, double s);

double h2o_region2_b2bc_p_h(double h);
double h2o_region2_b2bc_h_p(double p);
double h2o_region2_b2ab_h_s(double s);

double h2o_region2a_T_ph(double p, double h);
double h2o_region2b_T_ph(double p, double h);
double h2o_region2c_T_ph(double p, double h);

double h2o_region2a_T_ps(double p, double s);
double h2o_region2b_T_ps(double p, double s);
double h2o_region2c_T_ps(double p, double s);

double h2o_region2a_p_hs(double h, double s);
double h2o_region2b_p_hs(double h, double s);
double h2o_region2c_p_hs(double h, double s);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION2_H*/

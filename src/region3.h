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

double h2o_region3_T_ph(double p, double h);
double h2o_region3_v_ph(double p, double h);
double h2o_region3_T_ps(double p, double s);
double h2o_region3_v_ps(double p, double s);
double h2o_region3_p_hs(double h, double s);

/* special use functions */

enum h2o_region3_subregion
{
	H2O_REGION3A,
	H2O_REGION3B,

	H2O_REGION3_MAX
};

enum h2o_region3_subregion h2o_region3_subregion_ph(double p, double h);
enum h2o_region3_subregion h2o_region3_subregion_ps(double p, double s);
enum h2o_region3_subregion h2o_region3_subregion_hs(double h, double s);

double h2o_region3_b3ab_h_p(double p);

double h2o_region3a_T_ph(double p, double h);
double h2o_region3b_T_ph(double p, double h);
double h2o_region3a_v_ph(double p, double h);
double h2o_region3b_v_ph(double p, double h);

double h2o_region3a_T_ps(double p, double s);
double h2o_region3b_T_ps(double p, double s);
double h2o_region3a_v_ps(double p, double s);
double h2o_region3b_v_ps(double p, double s);
double h2o_region3a_p_hs(double h, double s);
double h2o_region3b_p_hs(double h, double s);

double h2o_region3_psat_h(double h);
double h2o_region3_psat_s(double s);

double h2o_region3a_v_pT(double p, double T);
double h2o_region3b_v_pT(double p, double T);
double h2o_region3c_v_pT(double p, double T);
double h2o_region3d_v_pT(double p, double T);
double h2o_region3e_v_pT(double p, double T);
double h2o_region3g_v_pT(double p, double T);
double h2o_region3h_v_pT(double p, double T);
double h2o_region3k_v_pT(double p, double T);
double h2o_region3l_v_pT(double p, double T);
double h2o_region3q_v_pT(double p, double T);
double h2o_region3r_v_pT(double p, double T);
double h2o_region3s_v_pT(double p, double T);
double h2o_region3t_v_pT(double p, double T);
double h2o_region3u_v_pT(double p, double T);
double h2o_region3v_v_pT(double p, double T);
double h2o_region3w_v_pT(double p, double T);
double h2o_region3x_v_pT(double p, double T);
double h2o_region3y_v_pT(double p, double T);
double h2o_region3z_v_pT(double p, double T);

#ifdef __cplusplus
};
#endif /*__cplusplus*/

#endif /*_H2O_REGION3_H*/

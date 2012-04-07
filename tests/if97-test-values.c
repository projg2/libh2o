/* libh2o -- water & steam properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <stdio.h>
#include <math.h>

#include "boundaries.h"
#include "region1.h"
#include "region2.h"
#include "region3.h"
#include "region4.h"
#include "region5.h"

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

#define CHECK(call, expected, precision) _check(call, expected, precision, #call)

int main(void)
{
	CHECK(h2o_region4_p_T(300), 0.353658941E-2, 1E-11);
	CHECK(h2o_region4_p_T(500), 0.263889776E+1, 1E-8);
	CHECK(h2o_region4_p_T(600), 0.123443146E+2, 1E-7);

	CHECK(h2o_region4_T_p(0.1), 0.372755919E3, 1E-6);
	CHECK(h2o_region4_T_p(1.0), 0.453035632E3, 1E-6);
	CHECK(h2o_region4_T_p(10.), 0.584149488E3, 1E-6);

	CHECK(h2o_b23_p_T(0.62315E3), 0.165291643E2, 1E-7);
	CHECK(h2o_b23_T_p(0.165291643E2), 0.62315E3, 1E-6);

	CHECK(h2o_region1_v_pT(3., 300), 0.100215168E-2, 1E-11);
	CHECK(h2o_region1_v_pT(80, 300), 0.971180894E-3, 1E-12);
	CHECK(h2o_region1_v_pT(3., 500), 0.120241800E-2, 1E-11);
	CHECK(h2o_region1_h_pT(3., 300), 0.115331273E+3, 1E-6);
	CHECK(h2o_region1_h_pT(80, 300), 0.184142828E+3, 1E-6);
	CHECK(h2o_region1_h_pT(3., 500), 0.975542239E+3, 1E-6);
	CHECK(h2o_region1_u_pT(3., 300), 0.112324818E+3, 1E-6);
	CHECK(h2o_region1_u_pT(80, 300), 0.106448356E+3, 1E-6);
	CHECK(h2o_region1_u_pT(3., 500), 0.971934985E+3, 1E-6);
	CHECK(h2o_region1_s_pT(3., 300), 0.392294792E+0, 1E-9);
	CHECK(h2o_region1_s_pT(80, 300), 0.368563852E+0, 1E-9);
	CHECK(h2o_region1_s_pT(3., 500), 0.258041912E+1, 1E-8);

	CHECK(h2o_region1_T_ph(3., 500.), 0.391798509E+3, 1E-6);
	CHECK(h2o_region1_T_ph(80, 500.), 0.378108626E+3, 1E-6);
	CHECK(h2o_region1_T_ph(80, 1500), 0.611041229E+3, 1E-6);

	CHECK(h2o_region1_T_ps(3., 0.5), 0.307842258E+3, 1E-6);
	CHECK(h2o_region1_T_ps(80, 0.5), 0.309979785E+3, 1E-6);
	CHECK(h2o_region1_T_ps(80, 3.0), 0.565899909E+3, 1E-6);

	CHECK(h2o_region1_p_hs(1E-3, 0.0), 0.9800980612E-3, 1E-12);
	CHECK(h2o_region1_p_hs(90E0, 0.0), 0.9192954727E+2, 1E-7);
	CHECK(h2o_region1_p_hs(15E2, 3.4), 0.5868294423E+2, 1E-7);

	CHECK(h2o_region2_v_pT(35E-4, 300), 0.394913866E+2, 1E-7);
	CHECK(h2o_region2_v_pT(35E-4, 700), 0.923015898E+2, 1E-7);
	CHECK(h2o_region2_v_pT(30E+0, 700), 0.542946619E-2, 1E-11);

	CHECK(h2o_region2_h_pT(35E-4, 300), 0.254991145E+4, 1E-5);
	CHECK(h2o_region2_h_pT(35E-4, 700), 0.333568375E+4, 1E-5);
	CHECK(h2o_region2_h_pT(30E+0, 700), 0.263149474E+4, 1E-5);

	CHECK(h2o_region2_u_pT(35E-4, 300), 0.241169160E+4, 1E-5);
	CHECK(h2o_region2_u_pT(35E-4, 700), 0.301262819E+4, 1E-5);
	CHECK(h2o_region2_u_pT(30E+0, 700), 0.246861076E+4, 1E-5);

	CHECK(h2o_region2_s_pT(35E-4, 300), 0.852238967E+1, 1E-8);
	CHECK(h2o_region2_s_pT(35E-4, 700), 0.101749996E+2, 1E-7);
	CHECK(h2o_region2_s_pT(30E+0, 700), 0.517540298E+1, 1E-8);

	CHECK(h2o_region2_meta_v_pT(1.0, 450), 0.192516540E+0, 1E-9);
	CHECK(h2o_region2_meta_v_pT(1.0, 440), 0.186212297E+0, 1E-9);
	CHECK(h2o_region2_meta_v_pT(1.5, 450), 0.121685206E+0, 1E-9);

	CHECK(h2o_region2_meta_h_pT(1.0, 450), 0.276881115E+4, 1E-5);
	CHECK(h2o_region2_meta_h_pT(1.0, 440), 0.274015123E+4, 1E-5);
	CHECK(h2o_region2_meta_h_pT(1.5, 450), 0.272134539E+4, 1E-5);

	CHECK(h2o_region2_meta_u_pT(1.0, 450), 0.257629461E+4, 1E-5);
	CHECK(h2o_region2_meta_u_pT(1.0, 440), 0.255393894E+4, 1E-5);
	CHECK(h2o_region2_meta_u_pT(1.5, 450), 0.253881758E+4, 1E-5);

	CHECK(h2o_region2_meta_s_pT(1.0, 450), 0.656660377E+1, 1E-8);
	CHECK(h2o_region2_meta_s_pT(1.0, 440), 0.650218759E+1, 1E-8);
	CHECK(h2o_region2_meta_s_pT(1.5, 450), 0.629170440E+1, 1E-8);

	CHECK(h2o_region2_b2bc_p_h(0.3516004323E4), 0.1E3, 1E-6);
	CHECK(h2o_region2_b2bc_h_p(0.1E3), 0.3516004323E4, 1E-5);
	CHECK(h2o_region2_b2ab_h_s(7), 0.3376437884E4, 1E-5);

	CHECK(h2o_region2a_T_ph(1E-3, 3000), 0.534433241E3, 1E-6);
	CHECK(h2o_region2a_T_ph(3.00, 3000), 0.575373370E3, 1E-6);
	CHECK(h2o_region2a_T_ph(3.00, 4000), 0.101077577E4, 1E-5);
	CHECK(h2o_region2b_T_ph(5.00, 3500), 0.801299102E3, 1E-6);
	CHECK(h2o_region2b_T_ph(5.00, 4000), 0.101531583E4, 1E-5);
	CHECK(h2o_region2b_T_ph(25.0, 3500), 0.875279054E3, 1E-6);
	CHECK(h2o_region2c_T_ph(40.0, 2700), 0.743056411E3, 1E-6);
	CHECK(h2o_region2c_T_ph(60.0, 2700), 0.791137067E3, 1E-6);
	CHECK(h2o_region2c_T_ph(60.0, 3200), 0.882756860E3, 1E-6);

	CHECK(h2o_region2a_T_ps(0.1, 7.50), 0.399517097E3, 1E-6);
	CHECK(h2o_region2a_T_ps(0.1, 8.00), 0.514127081E3, 1E-6);
	CHECK(h2o_region2a_T_ps(2.5, 8.00), 0.103984917E4, 1E-5);
	CHECK(h2o_region2b_T_ps(8.0, 6.00), 0.600484040E3, 1E-6);
	CHECK(h2o_region2b_T_ps(8.0, 7.50), 0.106495556E4, 1E-5);
	CHECK(h2o_region2b_T_ps(90., 6.00), 0.103801126E4, 1E-5);
	CHECK(h2o_region2c_T_ps(20., 5.75), 0.697992849E3, 1E-6);
	CHECK(h2o_region2c_T_ps(80., 5.25), 0.854011484E3, 1E-6);
	CHECK(h2o_region2c_T_ps(80., 5.75), 0.949017998E3, 1E-6);

	CHECK(h2o_region2a_p_hs(2800, 6.5), 0.1371012767E+1, 1E-9);
	CHECK(h2o_region2a_p_hs(2800, 9.5), 0.1879743844E-2, 1E-12);
	CHECK(h2o_region2a_p_hs(4100, 9.5), 0.1024788997E+0, 1E-10);
	CHECK(h2o_region2b_p_hs(2800, 6.0), 0.4793911442E+1, 1E-9);
	CHECK(h2o_region2b_p_hs(3600, 6.0), 0.8395519209E+2, 1E-8);
	CHECK(h2o_region2b_p_hs(3600, 7.0), 0.7527161441E+1, 1E-9);
	CHECK(h2o_region2c_p_hs(2800, 5.1), 0.9439202060E+2, 1E-8);
	CHECK(h2o_region2c_p_hs(2800, 5.8), 0.8414574124E+1, 1E-9);
	CHECK(h2o_region2c_p_hs(3400, 5.8), 0.8376903879E+2, 1E-8);

	CHECK(h2o_region3_p_rhoT(500, 650), 0.255837018E2, 1E-7);
	CHECK(h2o_region3_p_rhoT(200, 650), 0.222930643E2, 1E-7);
	CHECK(h2o_region3_p_rhoT(500, 750), 0.783095639E2, 1E-7);
	CHECK(h2o_region3_h_rhoT(500, 650), 0.186343019E4, 1E-5);
	CHECK(h2o_region3_h_rhoT(200, 650), 0.237512401E4, 1E-5);
	CHECK(h2o_region3_h_rhoT(500, 750), 0.225868845E4, 1E-5);
	CHECK(h2o_region3_u_rhoT(500, 650), 0.181226279E4, 1E-5);
	CHECK(h2o_region3_u_rhoT(200, 650), 0.226365868E4, 1E-5);
	CHECK(h2o_region3_u_rhoT(500, 750), 0.210206932E4, 1E-5);
	CHECK(h2o_region3_s_rhoT(500, 650), 0.405427273E1, 1E-8);
	CHECK(h2o_region3_s_rhoT(200, 650), 0.485438792E1, 1E-8);
	CHECK(h2o_region3_s_rhoT(500, 750), 0.446971906E1, 1E-8);

	CHECK(h2o_region4_T_hs(1800, 5.3), 0.3468475498E3, 1E-7);
	CHECK(h2o_region4_T_hs(2400, 6.0), 0.4251373305E3, 1E-7);
	CHECK(h2o_region4_T_hs(2500, 5.5), 0.5225579013E3, 1E-7);

	CHECK(h2o_region5_v_pT(.5, 1500), 0.138455090E+1, 1E-08);
	CHECK(h2o_region5_v_pT(30, 1500), 0.230761299E-1, 1E-10);
	CHECK(h2o_region5_v_pT(30, 2000), 0.311385219E-1, 1E-10);

	CHECK(h2o_region5_h_pT(.5, 1500), 0.521976855E+4, 1E-5);
	CHECK(h2o_region5_h_pT(30, 1500), 0.516723514E+4, 1E-5);
	CHECK(h2o_region5_h_pT(30, 2000), 0.657122604E+4, 1E-5);

	CHECK(h2o_region5_u_pT(.5, 1500), 0.452749310E+4, 1E-5);
	CHECK(h2o_region5_u_pT(30, 1500), 0.447495124E+4, 1E-5);
	CHECK(h2o_region5_u_pT(30, 2000), 0.563707038E+4, 1E-5);

	CHECK(h2o_region5_s_pT(.5, 1500), 0.965408875E+1, 1E-8);
	CHECK(h2o_region5_s_pT(30, 1500), 0.772970133E+1, 1E-8);
	CHECK(h2o_region5_s_pT(30, 2000), 0.853640523E+1, 1E-8);

	CHECK(h2o_b14_h_s(1.0), 0.3085509647E3, 1E-7);
	CHECK(h2o_b14_h_s(2.0), 0.7006304472E3, 1E-7);
	CHECK(h2o_b14_h_s(3.0), 0.1198359754E4, 1E-6);

	return exit_status;
}

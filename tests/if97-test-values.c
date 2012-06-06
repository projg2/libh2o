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

int tests_done = 0;
int tests_failed = 0;

void _check(double result, double expected, double precision, const char* call)
{
	double difference = fabs(expected - result);

	if (difference >= precision)
	{
		fprintf(stderr, "[FAIL] %s = %.9e, while %.9e expected.\n",
				call, result, expected);
		++tests_failed;
	}
	else
		fprintf(stderr, "[ OK ] %s = %.9e.\n",
				call, result);

	++tests_done;
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
	CHECK(h2o_region1_cp_pT(3., 300), 0.417301218E+1, 1E-8);
	CHECK(h2o_region1_cp_pT(80, 300), 0.401008987E+1, 1E-8);
	CHECK(h2o_region1_cp_pT(3., 500), 0.465580682E+1, 1E-8);
	CHECK(h2o_region1_w_pT(3., 300), 0.150773921E+4, 1E-5);
	CHECK(h2o_region1_w_pT(80, 300), 0.163469054E+4, 1E-5);
	CHECK(h2o_region1_w_pT(3., 500), 0.124071337E+4, 1E-5);

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
	CHECK(h2o_region2_cp_pT(35E-4, 300), 0.191300162E+1, 1E-8);
	CHECK(h2o_region2_cp_pT(35E-4, 700), 0.208141274E+1, 1E-8);
	CHECK(h2o_region2_cp_pT(30E+0, 700), 0.103505092E+2, 1E-7);
	CHECK(h2o_region2_w_pT(35E-4, 300), 0.427920172E+3, 1E-6);
	CHECK(h2o_region2_w_pT(35E-4, 700), 0.644289068E+3, 1E-6);
	CHECK(h2o_region2_w_pT(30E+0, 700), 0.480386523E+3, 1E-6);

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
	CHECK(h2o_region2_meta_cp_pT(1.0, 450), 0.276349265E+1, 1E-8);
	CHECK(h2o_region2_meta_cp_pT(1.0, 440), 0.298166443E+1, 1E-8);
	CHECK(h2o_region2_meta_cp_pT(1.5, 450), 0.362795578E+1, 1E-8);
	CHECK(h2o_region2_meta_w_pT(1.0, 450), 0.498408101E+3, 1E-6);
	CHECK(h2o_region2_meta_w_pT(1.0, 440), 0.489363295E+3, 1E-6);
	CHECK(h2o_region2_meta_w_pT(1.5, 450), 0.481941819E+3, 1E-6);

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
	CHECK(h2o_region3_cp_rhoT(500, 650), 0.138935717E2, 1E-7);
	CHECK(h2o_region3_cp_rhoT(200, 650), 0.446579342E2, 1E-7);
	CHECK(h2o_region3_cp_rhoT(500, 750), 0.634165359E1, 1E-8);
	CHECK(h2o_region3_w_rhoT(500, 650), 0.502005554E3, 1E-6);
	CHECK(h2o_region3_w_rhoT(200, 650), 0.383444594E3, 1E-6);
	CHECK(h2o_region3_w_rhoT(500, 750), 0.760696041E3, 1E-6);

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
	CHECK(h2o_region5_cp_pT(.5, 1500), 0.261609445E+1, 1E-8);
	CHECK(h2o_region5_cp_pT(30, 1500), 0.272724317E+1, 1E-8);
	CHECK(h2o_region5_cp_pT(30, 2000), 0.288569882E+1, 1E-8);
	CHECK(h2o_region5_w_pT(.5, 1500), 0.917068690E+3, 1E-6);
	CHECK(h2o_region5_w_pT(30, 1500), 0.928548002E+3, 1E-6);
	CHECK(h2o_region5_w_pT(30, 2000), 0.106736948E+4, 1E-5);

	CHECK(h2o_b13_h_s(3.7), 0.1632525047E4, 1E-6);
	CHECK(h2o_b13_h_s(3.6), 0.1593027214E4, 1E-6);
	CHECK(h2o_b13_h_s(3.5), 0.1566104611E4, 1E-6);
	CHECK(h2o_b14_h_s(1.0), 0.3085509647E3, 1E-7);
	CHECK(h2o_b14_h_s(2.0), 0.7006304472E3, 1E-7);
	CHECK(h2o_b14_h_s(3.0), 0.1198359754E4, 1E-6);
	CHECK(h2o_b3a4_h_s(3.8), 0.1685025565E4, 1E-6);
	CHECK(h2o_b3a4_h_s(4.0), 0.1816891476E4, 1E-6);
	CHECK(h2o_b3a4_h_s(4.2), 0.1949352563E4, 1E-6);
	CHECK(h2o_b2ab4_h_s(7.0), 0.2723729985E4, 1E-6);
	CHECK(h2o_b2ab4_h_s(8.0), 0.2599047210E4, 1E-6);
	CHECK(h2o_b2ab4_h_s(9.0), 0.2511861477E4, 1E-6);
	CHECK(h2o_b2c3b4_h_s(5.5), 0.2687693850E4, 1E-6);
	CHECK(h2o_b2c3b4_h_s(5.0), 0.2451623609E4, 1E-6);
	CHECK(h2o_b2c3b4_h_s(4.5), 0.2144360448E4, 1E-6);

	CHECK(h2o_b23_T_hs(2600, 5.10), 0.7135259364E+3, 1E-7);
	CHECK(h2o_b23_T_hs(2700, 5.15), 0.7685345532E+3, 1E-7);
	CHECK(h2o_b23_T_hs(2800, 5.20), 0.8176202120E+3, 1E-7);

	CHECK(h2o_region3_b3ab_h_p(25), 0.2095936454E+4, 1E-6);

	CHECK(h2o_region3a_T_ph(20., 1700), 0.6293083892E+3, 1E-7);
	CHECK(h2o_region3a_T_ph(50., 2000), 0.6905718338E+3, 1E-7);
	CHECK(h2o_region3a_T_ph(100, 2100), 0.7336163014E+3, 1E-7);
	CHECK(h2o_region3b_T_ph(20., 2500), 0.6418418053E+3, 1E-7);
	CHECK(h2o_region3b_T_ph(50., 2400), 0.7351848618E+3, 1E-7);
	CHECK(h2o_region3b_T_ph(100, 2700), 0.8420460876E+3, 1E-7);

	CHECK(h2o_region3a_v_ph(20., 1700), 0.1749903962E-2, 1E-12);
	CHECK(h2o_region3a_v_ph(50., 2000), 0.1908139035E-2, 1E-12);
	CHECK(h2o_region3a_v_ph(100, 2100), 0.1676229776E-2, 1E-12);
	CHECK(h2o_region3b_v_ph(20., 2500), 0.6670547043E-2, 1E-12);
	CHECK(h2o_region3b_v_ph(50., 2400), 0.2801244590E-2, 1E-12);
	CHECK(h2o_region3b_v_ph(100, 2700), 0.2404234998E-2, 1E-12);

	CHECK(h2o_region3a_T_ps(20., 3.8), 0.6282959869E+3, 1E-7);
	CHECK(h2o_region3a_T_ps(50., 3.6), 0.6297158726E+3, 1E-7);
	CHECK(h2o_region3a_T_ps(100, 4.0), 0.7056880237E+3, 1E-7);
	CHECK(h2o_region3b_T_ps(20., 5.0), 0.6401176443E+3, 1E-7);
	CHECK(h2o_region3b_T_ps(50., 4.5), 0.7163687517E+3, 1E-7);
	CHECK(h2o_region3b_T_ps(100, 5.0), 0.8474332825E+3, 1E-7);

	CHECK(h2o_region3a_v_ps(20., 3.8), 0.1733791463E-2, 1E-12);
	CHECK(h2o_region3a_v_ps(50., 3.6), 0.1469680170E-2, 1E-12);
	CHECK(h2o_region3a_v_ps(100, 4.0), 0.1555893131E-2, 1E-12);
	CHECK(h2o_region3b_v_ps(20., 5.0), 0.6262101987E-2, 1E-12);
	CHECK(h2o_region3b_v_ps(50., 4.5), 0.2332634294E-2, 1E-12);
	CHECK(h2o_region3b_v_ps(100, 5.0), 0.2449610757E-2, 1E-12);

	CHECK(h2o_region3a_p_hs(1700, 3.8), 0.2555703246E+2, 1E-8);
	CHECK(h2o_region3a_p_hs(2000, 4.2), 0.4540873468E+2, 1E-8);
	CHECK(h2o_region3a_p_hs(2100, 4.3), 0.6078123340E+2, 1E-8);
	CHECK(h2o_region3b_p_hs(2600, 5.1), 0.3434999263E+2, 1E-8);
	CHECK(h2o_region3b_p_hs(2400, 4.7), 0.6363924887E+2, 1E-8);
	CHECK(h2o_region3b_p_hs(2700, 5.0), 0.8839043281E+2, 1E-8);

	CHECK(h2o_region3_psat_h(1700), 0.1724175718E2, 1E-8);
	CHECK(h2o_region3_psat_h(2000), 0.2193442957E2, 1E-8);
	CHECK(h2o_region3_psat_h(2400), 0.2018090839E2, 1E-8);
	CHECK(h2o_region3_psat_s(3.8), 0.1687755057E2, 1E-8);
	CHECK(h2o_region3_psat_s(4.2), 0.2164451789E2, 1E-8);
	CHECK(h2o_region3_psat_s(5.2), 0.1668968482E2, 1E-8);

	CHECK(h2o_region3a_v_pT(50.000, 630.00), 0.1470853100E-2, 1E-12);
	CHECK(h2o_region3a_v_pT(80.000, 670.00), 0.1503831359E-2, 1E-12);
	CHECK(h2o_region3b_v_pT(50.000, 710.00), 0.2204728587E-2, 1E-12);
	CHECK(h2o_region3b_v_pT(80.000, 750.00), 0.1973692940E-2, 1E-12);
	CHECK(h2o_region3c_v_pT(20.000, 630.00), 0.1761696406E-2, 1E-12);
	CHECK(h2o_region3c_v_pT(30.000, 650.00), 0.1819560617E-2, 1E-12);
	CHECK(h2o_region3d_v_pT(26.000, 656.00), 0.2245587720E-2, 1E-12);
	CHECK(h2o_region3d_v_pT(30.000, 670.00), 0.2506897702E-2, 1E-12);
	CHECK(h2o_region3e_v_pT(26.000, 661.00), 0.2970225962E-2, 1E-12);
	CHECK(h2o_region3e_v_pT(30.000, 675.00), 0.3004627086E-2, 1E-12);
	CHECK(h2o_region3f_v_pT(26.000, 671.00), 0.5019029401E-2, 1E-12);
	CHECK(h2o_region3f_v_pT(30.000, 690.00), 0.4656470142E-2, 1E-12);
	CHECK(h2o_region3g_v_pT(23.600, 649.00), 0.2163198378E-2, 1E-12);
	CHECK(h2o_region3g_v_pT(24.000, 650.00), 0.2166044161E-2, 1E-12);
	CHECK(h2o_region3h_v_pT(23.600, 652.00), 0.2651081407E-2, 1E-12);
	CHECK(h2o_region3h_v_pT(24.000, 654.00), 0.2967802335E-2, 1E-12);
	CHECK(h2o_region3i_v_pT(23.600, 653.00), 0.3273916816E-2, 1E-12);
	CHECK(h2o_region3i_v_pT(24.000, 655.00), 0.3550329864E-2, 1E-12);
	CHECK(h2o_region3j_v_pT(23.500, 655.00), 0.4545001142E-2, 1E-12);
	CHECK(h2o_region3j_v_pT(24.000, 660.00), 0.5100267704E-2, 1E-12);
	CHECK(h2o_region3k_v_pT(23.000, 660.00), 0.6109525997E-2, 1E-12);
	CHECK(h2o_region3k_v_pT(24.000, 670.00), 0.6427325645E-2, 1E-12);
	CHECK(h2o_region3l_v_pT(22.600, 646.00), 0.2117860851E-2, 1E-12);
	CHECK(h2o_region3l_v_pT(23.000, 646.00), 0.2062374674E-2, 1E-12);
	CHECK(h2o_region3m_v_pT(22.600, 648.60), 0.2533063780E-2, 1E-12);
	CHECK(h2o_region3m_v_pT(22.800, 649.30), 0.2572971781E-2, 1E-12);
	CHECK(h2o_region3n_v_pT(22.600, 649.00), 0.2923432711E-2, 1E-12);
	CHECK(h2o_region3n_v_pT(22.800, 649.70), 0.2913311494E-2, 1E-12);
	CHECK(h2o_region3o_v_pT(22.600, 649.10), 0.3131208996E-2, 1E-12);
	CHECK(h2o_region3o_v_pT(22.800, 649.90), 0.3221160278E-2, 1E-12);
	CHECK(h2o_region3p_v_pT(22.600, 649.40), 0.3715596186E-2, 1E-12);
	CHECK(h2o_region3p_v_pT(22.800, 650.20), 0.3664754790E-2, 1E-12);
	CHECK(h2o_region3q_v_pT(21.100, 640.00), 0.1970999272E-2, 1E-12);
	CHECK(h2o_region3q_v_pT(21.800, 643.00), 0.2043919161E-2, 1E-12);
	CHECK(h2o_region3r_v_pT(21.100, 644.00), 0.5251009921E-2, 1E-12);
	CHECK(h2o_region3r_v_pT(21.800, 648.00), 0.5256844741E-2, 1E-12);
	CHECK(h2o_region3s_v_pT(19.100, 635.00), 0.1932829079E-2, 1E-12);
	CHECK(h2o_region3s_v_pT(20.000, 638.00), 0.1985387227E-2, 1E-12);
	CHECK(h2o_region3t_v_pT(17.000, 626.00), 0.8483262001E-2, 1E-12);
	CHECK(h2o_region3t_v_pT(20.000, 640.00), 0.6227528101E-2, 1E-12);
	CHECK(h2o_region3u_v_pT(21.500, 644.60), 0.2268366647E-2, 1E-12);
	CHECK(h2o_region3u_v_pT(22.000, 646.10), 0.2296350553E-2, 1E-12);
	CHECK(h2o_region3v_v_pT(22.500, 648.60), 0.2832373260E-2, 1E-12);
	CHECK(h2o_region3v_v_pT(22.300, 647.90), 0.2811424405E-2, 1E-12);
	CHECK(h2o_region3w_v_pT(22.150, 647.50), 0.3694032281E-2, 1E-12);
	CHECK(h2o_region3w_v_pT(22.300, 648.10), 0.3622226305E-2, 1E-12);
	CHECK(h2o_region3x_v_pT(22.110, 648.00), 0.4528072649E-2, 1E-12);
	CHECK(h2o_region3x_v_pT(22.300, 649.00), 0.4556905799E-2, 1E-12);
	CHECK(h2o_region3y_v_pT(22.000, 646.84), 0.2698354719E-2, 1E-12);
	CHECK(h2o_region3y_v_pT(22.064, 647.05), 0.2717655648E-2, 1E-12);
	CHECK(h2o_region3z_v_pT(22.000, 646.89), 0.3798732962E-2, 1E-12);
	CHECK(h2o_region3z_v_pT(22.064, 647.15), 0.3701940010E-2, 1E-12);

	CHECK(h2o_region3ab_T_p(40.0), 0.6930341408E3, 1E-7);
	CHECK(h2o_region3cd_T_p(25.0), 0.6493659208E3, 1E-7);
	CHECK(h2o_region3ef_T_p(40.0), 0.7139593992E3, 1E-7);
	CHECK(h2o_region3gh_T_p(23.0), 0.6498873759E3, 1E-7);
	CHECK(h2o_region3ij_T_p(23.0), 0.6515778091E3, 1E-7);
	CHECK(h2o_region3jk_T_p(23.0), 0.6558338344E3, 1E-7);
	CHECK(h2o_region3mn_T_p(22.8), 0.6496054133E3, 1E-7);
	CHECK(h2o_region3op_T_p(22.8), 0.6500106943E3, 1E-7);
	CHECK(h2o_region3qu_T_p(22.0), 0.6456355027E3, 1E-7);
	CHECK(h2o_region3rx_T_p(22.0), 0.6482622754E3, 1E-7);
	CHECK(h2o_region3uv_T_p(22.3), 0.6477996121E3, 1E-7);
	CHECK(h2o_region3wx_T_p(22.3), 0.6482049480E3, 1E-7);

	if (tests_failed == 0)
		fprintf(stderr, "%d tests done. All tests suceeded.\n", tests_done);
	else
		fprintf(stderr, "%d of %d tests failed.\n", tests_failed, tests_done);

	return tests_failed ? 1 : 0;
}

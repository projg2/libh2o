#include <stdio.h>

#include "src/region1.h"
#include "src/saturation.h"

int main(void)
{
	double T;

	for (T = 300; T < 500; T += 0.004)
	{
		printf("%.11e\n", h2o_region1_v_pT(3, T));
	}

	return 0;
}

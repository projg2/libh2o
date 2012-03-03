#include <stdio.h>

#include "lowlevel/saturation.h"

int main(void)
{
	printf("%.13e\n", h2o_saturation_p_T(300.0));
	printf("%.13e\n", h2o_saturation_p_T(500.0));
	printf("%.13e\n", h2o_saturation_p_T(600.0));

	return 0;
}

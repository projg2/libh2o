#include <stdio.h>

double steam_saturation_p_T(double T);

int main(void)
{
	printf("%.13e\n", steam_saturation_p_T(300.0));
	printf("%.13e\n", steam_saturation_p_T(500.0));
	printf("%.13e\n", steam_saturation_p_T(600.0));

	return 0;
}

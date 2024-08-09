#include "de.h"
#include <math.h>
#include <cstdio>
#include <ctime>

#define TAM 1



static real model(real *x)
{
	  real fit = 0;
	  fit = x[0]*x[0] - x[0] - 12;
		if(fit < 0){fit *= -1;}
		return fit;
}

int main()
{
	int n_gens = 1000000;
	int pop_size = 30;
	int num_par = TAM;

	real l_bounds[num_par];
	real h_bounds[num_par];

	for(int i = 0; i < num_par; i++)
	{
	 	l_bounds[i] = -10;
	 	h_bounds[i] =  10;
	}

	printf("n_gens: %d, pop_size: %d, num_par: %d\n", n_gens, pop_size, num_par);
	real *x = differential_evolution(false, model, l_bounds, h_bounds, n_gens, 0.7, 0.5, 1.0, pop_size, num_par, BEST1BIN, time(NULL), 0, 0.01);

	printf("%ld, ",  sizeof(x));
	for(int i = 0; i < num_par; i++)
	{
        	printf("%lf, ",  x[i]);
    }

	printf("\n%f", model(x));

	return 0;
}
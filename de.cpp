#include "de.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int stop(real *a, int n, real atol, real tol)
{
    real std = stddev(a, n);
    real aux = atol + tol * abs(mean(a, n));
    //printf("\n %lf <= %lf", std , aux);
    if(std <= aux)
    {
        return 0;
    }
    return 1;
}

real stddev(real *a, int n)
{
  int i;
  real fsum=0.0, m;
  m = mean(a,n);

  for(i=0;i< n;i++)
  {
    fsum = fsum + (a[i]-m)*(a[i]-m);
  }

  return(sqrt(fsum/n));
}


real mean(real *a, int n)
{
  int i;
  real sum=0.0;

  for(i=0;i< n;i++)
  {
    sum = sum + a[i];
  }

  return(sum/n);
}

real **load_pop(const char *filename, int pop_size, int num_par) {

	real **filedata = NULL;

	filedata = (real**)malloc(sizeof(real*)*pop_size);

	for (int i = 0; i < pop_size; i++) {
		filedata[i] = (real*)malloc(sizeof(real)*num_par);
	}


    size_t len = 0;

    FILE *fp;
    fp = fopen(filename, "r");

    if(fp == NULL) {
        fprintf(stderr, "Error reading file %s\n", filename);
        return NULL;
    }

    real data;
	char ch;
	for (int i = 0; i < pop_size; i++) {
		for(int j = 0; j < num_par; j++) {
			fscanf(fp, "%lf", &data);
			filedata[i][j] = data;
		}

		while((ch = fgetc(fp)) != '\n');
		ch=fgetc(fp);
    }

    fclose(fp);

    return filedata;
}

static real rnd_uni(long *idum) {

    long j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    real temp;

    if(*idum <= 0) {
        if(-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for(j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if(*idum < 0)
                *idum += IM1;
            if(j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if(*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if(idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1)
        iy += IMM1;
    if((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

static int progress(bool from_file, void *instance, const real *x, const real *g, const real fx,
                    const real xnorm, const real gnorm, const real step, int n, int k,
                    int ls) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

real *differential_evolution(bool from_file, evalfn *evaluate, real *inibound_l, real *inibound_h, int genmax, real CR,
                             real F_min, real F_max, int NP, int D, int strategy, long seed, real atol, real tol) {

    /*-----Initialize random number generator-----------------------------*/

    real F = F_min;
    long rnd_uni_init = -(long)seed; /* initialization of rnd_uni() */

    real *cost = (real*)malloc(sizeof(real)*NP);

    real *best = (real*)malloc(sizeof(real)*D);

    real *bestit = (real*)malloc(sizeof(real)*D);
    real *tmp = (real*)malloc(sizeof(real)*D);

    real ***pold, ***pnew, ***pswap;

  	real **c = (real **)malloc(sizeof(real*)*NP);

    real **d = (real**)malloc(sizeof(real*)*NP);
    for(int i = 0; i < NP; i++) {
        d[i] = (real*)malloc(sizeof(real)*D);
    }

    real *original_ind = (real*)malloc(sizeof(real)*D);

    /*------Initialization------------------------------------------------*/
    /*------Right now this part is kept fairly simple and just generates--*/
    /*------random numbers in the range [-initfac, +initfac]. You might---*/
    /*------want to extend the init part such that you can initialize-----*/
    /*------each parameter separately.------------------------------------*/

	if(from_file) {
		c = load_pop("population.txt", NP, D);
	}

	else {
		c = (real **)malloc(sizeof(real*)*NP);
		for(int i = 0; i < NP; i++) {
			c[i] = (real*)malloc(sizeof(real)*D);
		}
	}


    for(int i = 0; i < NP; i++) {
        for(int j = 0; j < D; j++) {
			if(from_file) {
				c[i][j] = inibound_l[j] + c[i][j] * (inibound_h[j] - inibound_l[j]);
			}
			else {
				c[i][j] = inibound_l[j] + rnd_uni(&rnd_uni_init) * (inibound_h[j] - inibound_l[j]);
			}
        }

        cost[i] = evaluate(c[i]); /* obj. funct. value */
    }


    real cmin = cost[0];
    real trial_cost;
    int imin = 0;

    for(int i = 1; i < NP; i++) {
        if(cost[i] < cmin) {
            cmin = cost[i];
            imin = i;
        }
    }

    memcpy(best, c[imin], D * sizeof(real));
    memcpy(bestit, c[imin], D * sizeof(real));

    pold = &c; /* old population (generation G)   */
    pnew = &d; /* new population (generation G+1) */

    int r1, r2, r3, r4, r5; /* placeholders for random indexes    */

    /*=======================================================================*/
    /*=========Iteration loop================================================*/
    /*=======================================================================*/

    int gen = 0;
    int n, L;
    real cmean, cvar;

    while((gen < genmax) && (stop(cost, NP, 0.0, 0.01))) {
        gen++;
        imin = 0;

        for(int i = 0; i < NP; i++) {

            do { /* Endless loop for NP < 2 !!!     */
                r1 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while(r1 == i);

            do { /* Endless loop for NP < 3 !!!     */
                r2 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r2 == i) || (r2 == r1));

            do { /* Endless loop for NP < 4 !!!     */
                r3 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r3 == i) || (r3 == r1) || (r3 == r2));

            do { /* Endless loop for NP < 5 !!!     */
                r4 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r4 == i) || (r4 == r1) || (r4 == r2) || (r4 == r3));

            do { /* Endless loop for NP < 6 !!!     */
                r5 = (int)(rnd_uni(&rnd_uni_init) * NP);
            } while((r5 == i) || (r5 == r1) || (r5 == r2) || (r5 == r3) || (r5 == r4));

            /*=======Choice of strategy===============================================================*/
            /*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
            /*=======DE :  stands for Differential Evolution==========================================*/
            /*=======x  :  a string which denotes the vector to be perturbed==========================*/
            /*=======y  :  number of difference vectors taken for perturbation of x===================*/
            /*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
            /*                                                                                        */
            /*=======There are some simple rules which are worth following:===========================*/
            /*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
            /*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
            /*=======3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=*/
            /*           happens.                                                                     */
            /*=======4)  If you increase NP, F usually has to be decreased============================*/
            /*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/

            /*=======EXPONENTIAL CROSSOVER============================================================*/

            /*-------DE/best/1/exp--------------------------------------------------------------------*/
            /*-------Our oldest strategy but still not bad. However, we have found several------------*/
            /*-------optimization problems where misconvergence occurs.-------------------------------*/
            if(strategy == BEST1EXP) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));

                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));

                memcpy(original_ind, (*pold)[r1], D * sizeof(real));

            }
            /*-------DE/rand/1/exp-------------------------------------------------------------------*/
            /*-------This is one of my favourite strategies. It works especially well when the-------*/
            /*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
            /*-------as a first guess.---------------------------------------------------------------*/
            else if(strategy == RAND1EXP) /* strategy DE1 in the techreport */
            {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                memcpy(original_ind, tmp, D * sizeof(real));

                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = (*pold)[r1][n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));


            }
            /*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
            /*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
            /*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
            /*-------should play around with all three control variables.----------------------------*/
            else if(strategy == RANDBEST1EXP) /* similiar to DE2 but generally better */
            {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                memcpy(original_ind, tmp, D * sizeof(real));

                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = tmp[n] + F * (bestit[n] - tmp[n]) + F * ((*pold)[r1][n] - (*pold)[r2][n]);
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }
            /*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
            else if(strategy == BEST2EXP) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                memcpy(original_ind, tmp, D * sizeof(real));

                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }
            /*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
            else if(strategy == RAND2EXP) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                L = 0;
                do {
                    tmp[n] = (*pold)[r5][n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    n = (n + 1) % D;
                    L++;
                } while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
            }

            /*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

            /*-------DE/best/1/bin--------------------------------------------------------------------*/
            else if(strategy == BEST1BIN) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                memcpy(original_ind, tmp, D * sizeof(real));

                n = (int)(rnd_uni(&rnd_uni_init) * D);
                real F_J = F_min + rnd_uni(&rnd_uni_init) * (F_max - F_min);
                //real F_J = F_min;
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = bestit[n] + F_J * ((*pold)[r2][n] - (*pold)[r3][n]);
                    }
                    n = (n + 1) % D;
                }

            }
            /*-------DE/rand/1/bin-------------------------------------------------------------------*/
            else if(strategy == RAND1BIN) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = (*pold)[r1][n] + F * ((*pold)[r2][n] - (*pold)[r3][n]);
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
            else if(strategy == RANDBEST1BIN) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = tmp[n] + F * (bestit[n] - tmp[n]) + F * ((*pold)[r1][n] - (*pold)[r2][n]);
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/best/2/bin--------------------------------------------------------------------*/
            else if(strategy == BEST2BIN) {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] = bestit[n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    }
                    n = (n + 1) % D;
                }
            }
            /*-------DE/rand/2/bin--------------------------------------------------------------------*/
            else {
                memcpy(tmp, (*pold)[i], D * sizeof(real));
                n = (int)(rnd_uni(&rnd_uni_init) * D);
                for(L = 0; L < D; L++) /* perform D binomial trials */
                {
                    if((rnd_uni(&rnd_uni_init) < CR) || L == (D - 1)) /* change at least one parameter */
                    {
                        tmp[n] =
                            (*pold)[r5][n] + ((*pold)[r1][n] + (*pold)[r2][n] - (*pold)[r3][n] - (*pold)[r4][n]) * F;
                    }
                    n = (n + 1) % D;
                }
            }

            for(int j = 0; j < D; j++) { //----and bounce back----------------------------------------
                if(tmp[j] < inibound_l[j]) {
                    tmp[j] = inibound_l[j] + rnd_uni(&rnd_uni_init) * (original_ind[j] - inibound_l[j]);
                }
                if(tmp[j] > inibound_h[j]) {
                    tmp[j] = inibound_h[j] + rnd_uni(&rnd_uni_init) * (original_ind[j] - inibound_h[j]);
                }
            }

            /*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/

            trial_cost = evaluate(tmp); /* Evaluate new vector in tmp[] */

            if(trial_cost <= cost[i]) /* improved objective function value ? */
            {
                cost[i] = trial_cost;
                memcpy((*pnew)[i], tmp, D * sizeof(real));
                if(trial_cost < cmin)  /* Was this a new minimum? */
                {                      /* if so...*/
                    cmin = trial_cost; /* reset cmin to new low...*/
                    imin = i;
                    memcpy(best, tmp, D * sizeof(real));
                }
            } else {
                memcpy((*pnew)[i], (*pold)[i], D * sizeof(real)); /* replace target with old value */
            }
        } /* End mutation loop through pop. */

        memcpy(bestit, best, D * sizeof(real)); /* Save best population member of current iteration */

        /* swap population arrays. New generation becomes old one */

        pswap = pold;
        pold = pnew;
        pnew = pswap;

        /*----Compute the energy variance (just for monitoring purposes)-----------*/

        cmean = 0.; /* compute the mean value first */
        for(int j = 0; j < NP; j++) {
            cmean += cost[j];
        }
        cmean = cmean / NP;

        cvar = 0.; /* now the variance              */
        for(int j = 0; j < NP; j++) {
            cvar += (cost[j] - cmean) * (cost[j] - cmean);
        }
        cvar = cvar / (NP - 1);
    }
    /*=======================================================================*/
    /*=========End of iteration loop=========================================*/
    /*=======================================================================*/
    printf("\tit: %d\n", gen);
    return best;
}

/*-----------End of main()------------------------------------------*/
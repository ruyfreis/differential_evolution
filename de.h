/*------Constants for rnd_uni()--------------------------------------------*/
#include <stdbool.h>
#include <math.h>


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

#define BEST1EXP 1
#define RAND1EXP 2
#define RANDBEST1EXP 3
#define BEST2EXP 4
#define RAND2EXP 5
#define BEST1BIN 6
#define RAND1BIN 7
#define RANDBEST1BIN 8
#define BEST2BIN 9

#define real float

typedef real evalfn(real *);

real *differential_evolution(bool from_file, evalfn *evaluate, real *inibound_l, real *inibound_h, int genmax, real CR,
                             real F_min, real F_max, int NP, int D, int strategy, long seed, real atol, real tol);

//statistics functions//
int stop(real *a, int n, real atol, real tol);
real stddev(real *a, int n);
real mean(real *a, int n);
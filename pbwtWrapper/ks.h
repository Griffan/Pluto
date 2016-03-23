//
// Created by Fan Zhang on 12/23/15.
//

#ifndef PLUTO_KS_H
#define PLUTO_KS_H



//#include <R.h>
//#include <Rinternals.h>
//#include <Rmath.h>		/* constants */
#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif
#include <vector>
double ks_test(std::vector<int> &x, std::vector<int> &y, int EXACT=-1);
double K(int n, double d);
void m_multiply(double *A, double *B, double *C, int m);
void m_power(double *A, int eA, double *V, int *eV, int m, int n);
double psmirnov2x(double *x, int m, int n);
#endif //PLUTO_KS_H

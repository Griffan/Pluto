
/* ks.c
   Compute the asymptotic distribution of the one- and two-sample
   two-sided Kolmogorov-Smirnov statistics, and the exact distributions
   in the two-sided one-sample and two-sample cases.
*/


#include <math.h>
#include <stdlib.h>
#include "ks.h"
#include "pbwt/utils.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <limits>
#include <algorithm>
#define NA_INTEGER INT_MIN
double R_PosInf = INFINITY;
double R_NegInf = -INFINITY;

typedef union
{
    double value;
    unsigned int word[2];
} ieee_double;
/* gcc had problems with static const on AIX and Solaris
   Solaris was for gcc 3.1 and 3.2 under -O2 32-bit on 64-bit kernel */
#ifdef _AIX
#define CONST
#elif defined(sparc) && defined (__GNUC__) && __GNUC__ == 3
#define CONST
#else
#define CONST const
#endif

#ifdef WORDS_BIGENDIAN
static CONST int hw = 0;
static CONST int lw = 1;
#else  /* !WORDS_BIGENDIAN */
static CONST int hw = 1;
static CONST int lw = 0;
#endif /* WORDS_BIGENDIAN */

static double R_ValueOfNA(void)
{
    /* The gcc shipping with Fedora 9 gets this wrong without
     * the volatile declaration. Thanks to Marc Schwartz. */
    volatile ieee_double x;
    x.word[hw] = 0x7ff00000;
    x.word[lw] = 1954;
    return x.value;
}

int R_finite(double x)
{
#ifdef HAVE_WORKING_ISFINITE
    return isfinite(x);
#else
    return (!isnan(x) & (x != R_PosInf ) & (x != R_NegInf));
#endif
}
static double R_pow_di(double x, int n)
{
    const double  NA_REAL = R_ValueOfNA();
    double xn = 1.0;

    if (isnan(x)) return x;
    if (n == NA_INTEGER) return NA_REAL;

    if (n != 0) {
        if (!R_finite(x)) return pow(x, (double)n);

        BOOL is_neg = (n < 0);
        if(is_neg) n = -n;
        for(;;) {
            if(n & 01) xn *= x;
            if(n >>= 1) x *= x; else break;
        }
        if(is_neg) xn = 1. / xn;
    }
    return xn;
}

/* Two-sample two-sided asymptotic distribution */
static void
pkstwo(int n, double *x, double tol)
{
/* x[1:n] is input and output
 *
 * Compute
 *   \sum_{k=-\infty}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = 1 + 2 \sum_{k=1}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^\infty \exp(-(2k-1)^2\pi^2/(8x^2))
 *
 * See e.g. J. Durbin (1973), Distribution Theory for Tests Based on the
 * Sample Distribution Function.  SIAM.
 *
 * The 'standard' series expansion obviously cannot be used close to 0;
 * we use the alternative series for x < 1, and a rather crude estimate
 * of the series remainder term in this case, in particular using that
 * ue^(-lu^2) \le e^(-lu^2 + u) \le e^(-(l-1)u^2 - u^2+u) \le e^(-(l-1))
 * provided that u and l are >= 1.
 *
 * (But note that for reasonable tolerances, one could simply take 0 as
 * the value for x < 0.2, and use the standard expansion otherwise.)
 *
 */
    double new_v, old_v, s, w, z;
    int i, k, k_max;

    k_max = (int) sqrt(2 - log(tol));

    for(i = 0; i < n; i++) {
	if(x[i] < 1) {
	    z = - (M_PI_2 * M_PI_4) / (x[i] * x[i]);
	    w = log(x[i]);
	    s = 0;
	    for(k = 1; k < k_max; k += 2) {
		s += exp(k * k * z - w);
	    }
	    x[i] = s / M_1_SQRT_2PI;
	}
	else {
	    z = -2 * x[i] * x[i];
	    s = -1;
	    k = 1;
        old_v = 0;
        new_v = 1;
	    while(fabs(old_v - new_v) > tol) {
            old_v = new_v;
            new_v += 2 * s * exp(z * k * k);
		s *= -1;
		k++;
	    }
	    x[i] = new_v;
	}
    }
}

/* Two-sided two-sample */
double psmirnov2x(double *x, int m, int n)
{
    double md, nd, q, *u, w;
    int i, j;

    if(m > n) {
	i = n; n = m; m = i;
    }
    md = (double) m;
    nd = (double) n;
    /*
       q has 0.5/mn added to ensure that rounding error doesn't
       turn an equality into an inequality, eg abs(1/2-4/5)>3/10 

    */
    q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
    u = (double *) calloc(n + 1, sizeof(double));

    for(j = 0; j <= n; j++) {
	u[j] = ((j / nd) > q) ? 0 : 1;
    }
    for(i = 1; i <= m; i++) {
	w = (double)(i) / ((double)(i + n));
	if((i / md) > q)
	    u[0] = 0;
	else
	    u[0] = w * u[0];
	for(j = 1; j <= n; j++) {
	    if(fabs(i / md - j / nd) > q) 
		u[j] = 0;
	    else
		u[j] = w * u[j] + u[j - 1];
	}
    }
    return u[n];
}

double
K(int n, double d)
{
    /* Compute Kolmogorov's distribution.
       Code published in
	 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
	 "Evaluating Kolmogorov's distribution".
	 Journal of Statistical Software, Volume 8, 2003, Issue 18.
	 URL: http://www.jstatsoft.org/v08/i18/.
    */

   int k, m, i, j, g, eH, eQ;
   double h, s, *H, *Q;

   /* 
      The faster right-tail approximation is omitted here.
      s = d*d*n; 
      if(s > 7.24 || (s > 3.76 && n > 99)) 
          return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   */
   k = (int) (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = (double*) calloc(m * m, sizeof(double));
   Q = (double*) calloc(m * m, sizeof(double));
   for(i = 0; i < m; i++)
       for(j = 0; j < m; j++)
	   if(i - j + 1 < 0)
	       H[i * m + j] = 0;
	   else
	       H[i * m + j] = 1;
   for(i = 0; i < m; i++) {
       H[i * m] -= R_pow_di(h, i + 1);
       H[(m - 1) * m + i] -= R_pow_di(h, (m - i));
   }
   H[(m - 1) * m] += ((2 * h - 1 > 0) ? R_pow_di(2 * h - 1, m) : 0);
   for(i = 0; i < m; i++)
       for(j = 0; j < m; j++)
	   if(i - j + 1 > 0)
	       for(g = 1; g <= i - j + 1; g++)
		   H[i * m + j] /= g;
   eH = 0;
   m_power(H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];
   for(i = 1; i <= n; i++) {
       s = s * i / n;
       if(s < 1e-140) {
	   s *= 1e140;
	   eQ -= 140;
       }
   }
   s *= R_pow_di(10.0, eQ);
   free(H);
   free(Q);
   return(s);
}

void
m_multiply(double *A, double *B, double *C, int m)
{
    /* Auxiliary routine used by K().
       Matrix multiplication.
    */
    int i, j, k;
    double s;
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++) {
	    s = 0.;
	    for(k = 0; k < m; k++)
		s+= A[i * m + k] * B[k * m + j];
	    C[i * m + j] = s;
	}
}

void
m_power(double *A, int eA, double *V, int *eV, int m, int n)
{
    /* Auxiliary routine used by K().
       Matrix power.
    */
    double *B;
    int eB , i;

    if(n == 1) {
	for(i = 0; i < m * m; i++)
	    V[i] = A[i];
	*eV = eA;
	return;
    }
    m_power(A, eA, V, eV, m, n / 2);
    B = (double*) calloc(m * m, sizeof(double));
    m_multiply(V, V, B, m);
    eB = 2 * (*eV);
    if((n % 2) == 0) {
	for(i = 0; i < m * m; i++)
	    V[i] = B[i];
	*eV = eB;
    }
    else {
	m_multiply(A, B, V, m);
	*eV = eA + eB;
    }
    if(V[(m / 2) * m + (m / 2)] > 1e140) {
	for(i = 0; i < m * m; i++)
	    V[i] = V[i] * 1e-140;
	*eV += 140;
    }
    free(B);
}

/* Two-sided two-sample */
//SEXP pSmirnov2x(SEXP statistic, SEXP snx, SEXP sny)
//{
//    int nx = asInteger(snx), ny = asInteger(sny);
//    double st = asReal(statistic);
//    return ScalarReal(psmirnov2x(&st, nx, ny));
//}
//
///* Two-sample two-sided asymptotic distribution */
//SEXP pKS2(SEXP statistic, SEXP stol)
//{
//    int n = LENGTH(statistic);
//    double tol = asReal(stol);
//    SEXP ans = duplicate(statistic);
//    pkstwo(n, REAL(ans), tol);
//    return ans;
//}
//
//
///* The two-sided one-sample 'exact' distribution */
//SEXP pKolmogorov2x(SEXP statistic, SEXP sn)
//{
//    int n = asInteger(sn);
//    double st = asReal(statistic), p;
//    p = K(n, st);
//    return ScalarReal(p);
//}

double pkstwo_wrapper(int n, double* statistic, long tol = 1e-06) {

    std::vector<double> p;
    std::vector<int> index;
    for (int i = 0; i <n; ++i) {
        if(isnan(statistic[i]))
            continue;
        else if(statistic[i]>0) {
            index.push_back(i);
            p.push_back(0);
        }
    }
    if(n!=1) { std::cerr<<"only designed to handle 1 p value at a time!"<<std::endl;exit(EXIT_FAILURE);}
    if (index.size()>0)
        pkstwo(n,statistic,tol);
    return *statistic;
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

double ks_test(std::vector<int> &x, std::vector<int> &y, int EXACT) {
    std::vector<int> tmpx;
    for (int i = 0; i < x.size(); ++i) {
        if (!isnan(x[i]))
            tmpx.push_back(x[i]);
    }
    x = tmpx;
    double n = x.size();
    if (n < 1) {
        fprintf(stderr, "not enough 'x' data");
        exit(EXIT_FAILURE);
    }
    double PVAL = std::numeric_limits<double>::max();
    double STATISTIC=-1.0;
    /*if (isdigit(y[0]))*/ {
        tmpx.clear();
        for (int i = 0; i < y.size(); ++i) {
            if (!isnan(y[i]))
                tmpx.push_back(y[i]);
        }
        y = tmpx;
        double nx = n;
        double ny = y.size();
        if (ny < 1) {
            fprintf(stderr, "not enough 'y' data");
            exit(EXIT_FAILURE);
        }
        if(EXACT==-1) {
            EXACT = (nx*ny < 10000);
        }

        std::string METHOD("Two-sample Kolmogorov-Smirnov test");
        bool TIES = FALSE;
        n = nx * ny / (nx + ny);

        std::vector<int> w(x);
        std::vector<size_t> wOrder;
        std::vector<double> z;
        std::move(y.begin(), y.end(), std::back_inserter(w));

        //std::sort(w.begin(),w.end());
        wOrder=sort_indexes(w);
//        for (auto i: wOrder) {
//            cout << wOrder[i] << endl;
//        }
        int decline_cycle=0;
        int previous_sign=0;
        int num_switch=0;
        for (int j = 0; j <wOrder.size() ; ++j) {
//            if(j==0)
//            {
//                if(wOrder[j]+1 <= nx) {
//                    z.push_back(1. / (nx+ny));
//                    previous_sign=1;
//                }else {
//                    z.push_back(-1. / (nx+ny));
//                    previous_sign=-1;
//                }
//            }
//            else {
//                double addOn=0;
//                if(wOrder[j] + 1 <= nx)
//                {
//                    if(previous_sign==1) {
//                        decline_cycle++;
//                    }
//                    else
//                    {
//                        num_switch++;
//                        decline_cycle=0;
//                    }
//                    addOn=1. / (nx+ny) * pow(.5,decline_cycle);
//                    previous_sign=1;
//                }
//                else
//                {
//                    if(previous_sign==-1) decline_cycle++;
//                    else
//                    {
//                        num_switch++;
//                        decline_cycle=0;
//                    }
//                    addOn=-1./(nx+ny) * pow(.5,decline_cycle);
//                    previous_sign=-1;
//                }
//                z.push_back(z.back() + addOn);
//            }
            if(j==0)
            {
                if(wOrder[j]+1 <= nx) {
                    z.push_back(1. / nx);
                }else {
                    z.push_back(-1. / ny);
                }
            }
            else {
                double addOn=0;
                if(wOrder[j] + 1 <= nx)
                {
                    addOn=1. /nx;
                }
                else
                {
                    addOn=-1./ny;
                }
                z.push_back(z.back() + addOn);
            }
        }
        std::sort(w.begin(),w.end());

        std::vector<int> tmpW=w;

        auto last=std::unique(tmpW.begin(),tmpW.end());
        tmpW.erase(last,tmpW.end());

        if (tmpW.size() < (nx + ny)) {
            if (EXACT) {
                fprintf(stderr,"[Waring]cannot compute exact p-value with ties\n");
                EXACT=0;
            }
            else fprintf(stderr,"[Waring]p-value will be approximate in the presence of ties\n");
            std::vector<double> tmp;
            for (int i = 0; i <w.size()-1; ++i) {
                if((w[i+1]-w[i])!=0)
                    tmp.push_back(z[i]);
            }
            tmp.push_back(z[nx+ny-1]);
            z=tmp;

            TIES =TRUE;
        }

        /*STATISTIC < -
        switch (alternative, two.sided = max(abs(z)),
                             greater = max(z), less = -min(z))
            nm_alternative < -
        switch (alternative, two.sided = "two-sided",
                             less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
            if (exact && (alternative == "two.sided") && !TIES)
                PVAL < -1 - .Call(C_pSmirnov2x, STATISTIC, n.x, n.y)*/
        //here we only consider two-sided(equal to)

        for (int k = 0; k < z.size(); ++k) {
            if(STATISTIC < fabs(z[k]))
                STATISTIC=fabs(z[k]);
        }
        /*correction for switch numbers*/
        //STATISTIC*=10./(num_switch+1e-6);
        if(EXACT && !TIES)
            PVAL = 1 - psmirnov2x(&STATISTIC,nx,ny);
    }

//    names(STATISTIC) < -
//    switch (alternative, two.sided = "D",
//                         greater = "D^+", less = "D^-")
        if (PVAL==std::numeric_limits<double>::max()) {
//            pkstwo < -function(x, tol = 1e-06)
//            {
//                if (is.numeric(x))
//                    x < -as.double(x)
//                else stop("argument 'x' must be numeric")
//                p < -rep(0, length(x))
//                p[is.na(x)] < -NA
//                IND < -which(!is.na(x) & (x > 0))
//                if (length(IND))
//                    p[IND] < - .Call(C_pKS2, p = x[IND], tol)
//                p
//            }
//            PVAL < -ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) *
//                                                                  STATISTIC), exp(-2 * n * STATISTIC ^ 2))
            STATISTIC=sqrt(double(n)) * STATISTIC;
            PVAL=1 - pkstwo_wrapper(1,&STATISTIC,1e-06);
        }
    PVAL = std::min(1., std::max(0., PVAL));

    return (PVAL);
}
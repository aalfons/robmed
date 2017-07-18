#include <S.h>
#include "R_ext/Rdynload.h"
#include "robust.h"

void R_init_mypkg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


#ifdef USING_R
void F77_SUB(fseedi)(void)
{
  long x = 100;
  seed_in(&x);
}


void F77_SUB(fseedo)(void)
{
  long x = 100;
  seed_out(&x);
}


void F77_SUB(marriage)(Sint* a, Sint* b, Sint* c, Sfloat* d, Sint* e, Sfloat* f,
                       Sint* g, Sint* h, Sint* i, Sint* j, Sint* k, Sint* l)
{
	error("Genetic sampling not implemented.  Please report this to the package maintainer.");
}


void F77_SUB(getrandind)(Sint* n, Sint* np, Sint* maxslen, Sint* ntind, Sint* ni)
{
	error("getrandind not implemented.  Please report this to the package maintainer.");
}


void F77_SUB(xerror)(const char* msg, Sint* n, Sint* p, Sint* i)
{
  error(msg);
}
#endif

void F77_SUB(roblibrunif)(Sfloat* px)
{
#ifdef USING_R
	*px = unif_rand();
#else
	*px = unif_rand(S_evaluator);
#endif
}


Sfloat F77_SUB(robliberfc)(Sfloat* px)
//Sfloat erfc1(Sint ind, Sfloat x)
{
    Sint ind = 0;
		Sfloat x = *px;

    /* Initialized data */

    static Sfloat c__ = .564189583547756;
    static Sfloat a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static Sfloat b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static Sfloat p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static Sfloat q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static Sfloat r__[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static Sfloat s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    Sfloat ret_val, d__1;

    /* Local variables */
    Sfloat e, t, w, ax, bot, top;

/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */

/*                     ABS(X) .LE. 0.5 */

    ax = fabs(x);
    if (ax > 0.5) {
	goto L10;
    }
    t = x * x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
    ret_val = 0.5 - x * (top / bot) + 0.5;
    if (ind != 0) {
	ret_val = exp(t) * ret_val;
    }
    return ret_val;

/*                  0.5 .LT. ABS(X) .LE. 4 */

L10:
    if (ax > 4.0) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = top / bot;
    goto L40;

/*                      ABS(X) .GT. 4 */

L20:
    if (x <= -5.6f) {
	goto L50;
    }
    if (ind != 0) {
	goto L30;
    }
    if (x > 100.0) {
	goto L60;
    }
/*    if (x * x > -exparg(1)) {
	goto L60;
    }*/

L30:
/* Computing 2nd power */
    d__1 = 1.0 / x;
    t = d__1 * d__1;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    ret_val = (c__ - t * top / bot) / ax;

/*                      FINAL ASSEMBLY */

L40:
    if (ind == 0) {
	goto L41;
    }
    if (x < 0.0) {
	ret_val = exp(x * x) * 2.0 - ret_val;
    }
    return ret_val;
L41:
    w = (Sfloat) (x) * (Sfloat) (x);
    t = w;
    e = w - t;
    ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
    if (x < 0.0) {
	ret_val = 2.0 - ret_val;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE NEGATIVE X */

L50:
    ret_val = 2.0;
    if (ind != 0) {
	ret_val = exp(x * x) * 2.0;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE POSITIVE X */
/*                       WHEN IND = 0 */

L60:
    ret_val = 0.0;
    return ret_val;
} /* erfc1 */


Sfloat F77_SUB(robliberf)(Sfloat* px)
//Sfloat rl_erf(Sfloat x)
{
    Sfloat x = *px;

    /* Initialized data */

    static Sfloat c__ = .564189583547756;
    static Sfloat a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static Sfloat b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static Sfloat p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static Sfloat q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static Sfloat r__[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static Sfloat s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    Sfloat ret_val;

    /* Local variables */
    Sfloat t, x2, ax, bot, top;

/* ----------------------------------------------------------------------- */
/*             EVALUATION OF THE REAL ERROR FUNCTION */
/* ----------------------------------------------------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
    ax = fabs(x);
    if (ax > 0.5) {
	goto L10;
    }
    t = x * x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
    ret_val = x * (top / bot);
    return ret_val;

L10:
    if (ax > 4.0) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = 0.5 - exp(-x * x) * top / bot + 0.5;
    if (x < 0.0) {
	ret_val = -ret_val;
    }
    return ret_val;

L20:
    if (ax >= 5.8) {
	goto L30;
    }
    x2 = x * x;
    t = 1.0 / x2;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    ret_val = (c__ - top / (x2 * bot)) / ax;
    ret_val = 0.5 - exp(-x2) * ret_val + 0.5;
    if (x < 0.0) {
	ret_val = -ret_val;
    }
    return ret_val;

L30:
    ret_val = x > 0 ? 1 : -1;
    return ret_val;
} /* erf__ */



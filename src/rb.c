/* Robust bootstrap for MM-regression estimates
 * 2001 - Matias Salibian-Barrera
 * School of Mathematics and Statistics
 * Carleton University
 * 1125 Colonel By Drive
 * Ottawa, ON - K1S 5B6 - Canada
 * matias@math.carleton.ca
 * Reference: Salibian-Barrera and Zamar (2000) "Bootstrapping 
 * robust estimates of regression". To appear in The Annals of
 * Statistics.  
 */

#include <S.h>
#include "robust.h"
#include <math.h>

#define MAX_ITER 500
#define NA 99999.99
/* #define max(a,b) ((a)>(b)?(a):(b)) */
#define TOL_INVERSE 1e-8

void rl_rb_rand(Sfloat *X, Sfloat *y, Sint *N, Sint *P, Sint *M, Sfloat *ours, 
                Sfloat *Beta_m, Sfloat *Beta_s, Sfloat *Scale, Sint *Seed, 
                Sfloat *C, Sfloat *Psi_c, Sint *chi_fn, Sint *psi_fn,
                Sfloat *Beta)
{
	/* chi_fn = 1 -> Bisquare
	 * chi_fn = 2 -> Optimal
	 * psi_fn = 1 -> Bisquare
	 * psi_fn = 2 -> Optimal
	 */

	Sfloat rl_Psi_reg(Sfloat,Sfloat,Sint);
	Sfloat rl_Psi_reg_prime(Sfloat,Sfloat,Sint);
	Sfloat rl_Chi_prime(Sfloat,Sfloat,Sint);
	Sfloat rl_Chi(Sfloat,Sfloat,Sint);
	void rl_sampler_i(Sint, Sint *);
	Sint rl_inverse(Sfloat **,Sfloat **, Sint);
	void rl_reset_mat(Sfloat **, Sint , Sint );
	void rl_reset_vec(Sfloat *, Sint);
	void rl_vec_vec_symmetric(Sfloat **, Sfloat *, Sint);
	void rl_scalar_mat(Sfloat **, Sfloat, Sfloat **, Sint, Sint);
	void rl_scalar_vec(Sfloat *, Sfloat, Sfloat *, Sint);
	void rl_sum_mat(Sfloat **,Sfloat **, Sfloat **, Sint, Sint);
	void rl_sum_vec(Sfloat *, Sfloat *, Sfloat *, Sint);
	void rl_dif_vec(Sfloat *, Sfloat *, Sfloat *, Sint);
	void rl_mat_vec(Sfloat **, Sfloat *, Sfloat *, Sint, Sint);
	void rl_mat_mat(Sfloat **, Sfloat **, Sfloat **, Sint, Sint, Sint);

	Sfloat **xb; 
	Sfloat **x, **x2, **x3, **x4; 
	Sfloat *Fi, *res, *res_s, *w, scale;
	Sfloat *v, *v2, *v_aux, *yb;
	Sfloat u, u2, s, c, Psi_constant, beta;
	Sint n, p, m, seed, *indices;
	register Sint i, j, k;

	beta = *Beta;
	c = *C;
	Psi_constant = *Psi_c;
	n = *N;
	p = *P;
	m = *M;
	seed = *Seed;

	indices = Salloc(n, Sint);
	v = Salloc(p, Sfloat);
	v2 = Salloc(p, Sfloat);
	v_aux = Salloc(p, Sfloat);
	yb = Salloc(n, Sfloat);
	x = Salloc(n, Sfloat*);
	xb = Salloc(n, Sfloat*);
	Fi = Salloc(n, Sfloat);
	res = Salloc(n, Sfloat);
	res_s = Salloc(n, Sfloat);
	w = Salloc(n, Sfloat);
	x2 = Salloc(p, Sfloat*);
	x3 = Salloc(p, Sfloat*);
	x4 = Salloc(p, Sfloat*);

	for(i = 0; i < n; i++) {
		x[i] = Salloc(p, Sfloat);
		xb[i] = Salloc(p + 1, Sfloat);
	};

	for(i = 0; i < p; i++) {
		x2[i] = Salloc(p, Sfloat);
		x3[i] = Salloc(p, Sfloat);
		x4[i] = Salloc(p, Sfloat);
	};

	/* copy X into x for easier handling */
	for(i = 0; i < n; i++)
		for(j=0;j<p;j++)
			x[i][j]=X[j*n+i];

	scale = *Scale;

	/* get M-fitted values in Fi */
	rl_mat_vec(x,Beta_m,Fi,n,p);

	/* get residuals of M-est in res */
	rl_dif_vec(y,Fi,res,n);

	/* get S-fitted values in res_s */
	rl_mat_vec(x,Beta_s,res_s,n,p);

	/* get residuals of S-est in res_s */
	rl_dif_vec(y,res_s,res_s,n);

	/* set auxiliary matrices to zero */
	rl_reset_mat(x3,p,p);
	rl_reset_mat(x4,p,p);
	rl_reset_vec(v,p);

	u2 = 0.0;

	/* calculate correction matrix */
	for(i=0;i<n;i++) {
		u = res[i] / scale;
		w[i] = rl_Psi_reg(u,Psi_constant,*psi_fn)/res[i];
		rl_vec_vec_symmetric(x2,x[i],p);
		rl_scalar_mat(x2, rl_Psi_reg_prime(u, Psi_constant, *psi_fn), x2, p, p);
		rl_sum_mat(x3,x2,x3,p,p);
		rl_vec_vec_symmetric(x2,x[i],p);
		rl_scalar_mat(x2,w[i],x2,p,p);
		rl_sum_mat(x4,x2,x4,p,p);
		rl_scalar_vec(x[i],rl_Psi_reg_prime(u,Psi_constant,*psi_fn)*u,v_aux,p);
		rl_sum_vec(v,v_aux,v,p);
		u2 += rl_Chi_prime(u, c, *chi_fn) * u;
	};

	rl_scalar_vec(v, beta * (Sfloat) (n-p) * scale / u2 , v, p);
	rl_inverse(x3, x2, p);
	rl_mat_mat(x2, x4, x3, p, p, p);
	rl_mat_vec(x2, v, v2, p, p);
	rl_scalar_mat(x3, scale, x3, p, p);

	/* the correction matrix is now in x3 */
	/* the correction vector is now in v2 */

	/* srand(seed); */

	/* start the bootstrap replications */
	for(i = 0; i < m; i++) {

		rl_sampler_i(n, indices);

		/* get pseudo observed y's */
		for(j = 0; j < n; j++) /* xb[j][p] = */ 
			yb[j] = y[indices[j]];

		for(j = 0; j < n ; j++) 
			for(k = 0; k < p; k++) {
				xb[j][k] = x[indices[j]][k]; 
			};
	
	/* calculate robust bootsrap */
		rl_reset_vec(v, p); /* v <- 0 */ 
		rl_reset_mat(x2, p, p); /* x2 <- 0 */ 
		s = 0.0;
		for(j=0;j<n;j++) {
			rl_scalar_vec(xb[j], yb[j]* w[indices[j]], v_aux, p);
			rl_sum_vec(v, v_aux, v, p);
			rl_vec_vec_symmetric(x4, xb[j], p);
			rl_scalar_mat(x4, w[indices[j]], x4, p, p);
			rl_sum_mat(x2, x4, x2, p, p); 
			s += rl_Chi(res_s[indices[j]] / scale , c, *chi_fn);
		};

	s = s * scale / beta / (Sfloat) (n - p); 
	rl_inverse(x2, x4, p);		/* x4 <- x2^-1 */
	rl_mat_vec(x4, v, v_aux, p, p);	/* v_aux <- x4 * v */ 
	rl_dif_vec(v_aux, Beta_m, v_aux, p); 	/* v_aux <- v_aux - beta_m */

	/* v has the robust bootstrapped vector, correct it */ 
	rl_mat_vec(x3, v_aux, v, p, p);	/* v <- x3 * v_aux */ 
	rl_scalar_vec(v2, s-scale, v_aux, p);
	rl_sum_vec(v_aux, v, v, p);

	/* store the betas (roblib-wise!) */
	for(j = 0; j < p; j++) 
		ours[j * m + i] = v[j];
	};

	return;
}
 

void rl_rb_fixed(Sfloat *X, Sfloat *y, Sint *N, Sint *P, Sint *M, Sfloat *ours, 
                 Sfloat *Beta_m, Sfloat *Beta_s, Sfloat *Scale, Sint *Seed,
                 Sfloat *C, Sfloat *Psi_c, Sint *chi_fn, Sint *psi_fn,
                 Sfloat *Beta)
{
	/* chi_fn = 1 -> Bisquare
	 * chi_fn = 2 -> Optimal
	 * psi_fn = 1 -> Bisquare
	 * psi_fn = 2 -> Optimal
	 */

	void rl_reset_mat(Sfloat **, Sint , Sint );
	void rl_reset_vec(Sfloat *, Sint);
	Sfloat rl_Psi_reg(Sfloat,Sfloat,Sint);
	Sfloat rl_Psi_reg_prime(Sfloat,Sfloat,Sint);
	Sfloat rl_Chi_prime(Sfloat,Sfloat,Sint);
	Sfloat rl_Chi(Sfloat,Sfloat,Sint);
	void rl_sampler_i(Sint, Sint *);
	Sint rl_inverse(Sfloat **,Sfloat **, Sint);
	void rl_vec_vec_symmetric(Sfloat **, Sfloat *, Sint);
	void rl_scalar_mat(Sfloat **, Sfloat, Sfloat **, Sint, Sint);
	void rl_scalar_vec(Sfloat *, Sfloat, Sfloat *, Sint);
	void rl_sum_mat(Sfloat **,Sfloat **, Sfloat **, Sint, Sint);
	void rl_sum_vec(Sfloat *, Sfloat *, Sfloat *, Sint);
	void rl_dif_vec(Sfloat *, Sfloat *, Sfloat *, Sint);
	void rl_mat_vec(Sfloat **, Sfloat *, Sfloat *, Sint, Sint);
	void rl_mat_mat(Sfloat **, Sfloat **, Sfloat **, Sint, Sint, Sint);

	Sfloat **x, **x2, **x3, **x4;
	Sfloat *Fi, *res, *res_s, *w, scale;
	Sfloat *v, *v2, *v_aux, *yb, beta;
	Sfloat u, u2, s, c, Psi_constant;
	Sint n, p, m, seed, *indices;
	register Sint i, j;

	beta = *Beta;
	scale = *Scale;
	c = *C;
	Psi_constant = *Psi_c;
	n = *N;
	p = *P;
	m = *M;
	seed = *Seed;

	indices = Salloc(n, Sint);
	v = Salloc(p, Sfloat);
	v2 = Salloc(p, Sfloat);
	v_aux = Salloc(p, Sfloat);
	yb = Salloc(n, Sfloat);
	x = Salloc(n, Sfloat*);
	Fi = Salloc(n, Sfloat);
	res = Salloc(n, Sfloat);
	res_s = Salloc(n, Sfloat);
	w = Salloc(n, Sfloat);
	x2 = Salloc(p, Sfloat*);
	x3 = Salloc(p, Sfloat*);
	x4 = Salloc(p,Sfloat*);

	for(i = 0; i < n; i++) 
		x[i] = Salloc(p + 1, Sfloat);

	for(i = 0; i < p; i++) {
		x2[i] = Salloc(p, Sfloat);
		x3[i] = Salloc(p, Sfloat);
		x4[i] = Salloc(p, Sfloat);
	};

	/* copy X into x for easier handling */
	for(i=0;i<n;i++) 
		for(j=0;j<p;j++)
			x[i][j]=X[j*n+i];

	/* get M-fitted values in Fi */
	rl_mat_vec(x,Beta_m,Fi,n,p);

	/* get residuals of M-est in res */
	rl_dif_vec(y,Fi,res,n);

	/* get S-fitted values in res_s */
	rl_mat_vec(x,Beta_s,res_s,n,p);

	/* get residuals of S-est in res_s */
	rl_dif_vec(y,res_s,res_s,n);

	/* set auxiliary matrices to zero */
	rl_reset_mat(x3,p,p);
	rl_reset_mat(x4,p,p);
	rl_reset_vec(v,p);
	u2 = 0.0;

	/* calculate correction matrix */
	for(i=0;i<n;i++) {
		u = res[i]/scale ;
		w[i] = rl_Psi_reg(u,Psi_constant,*psi_fn)/res[i];
		rl_vec_vec_symmetric(x2,x[i],p);
		rl_scalar_mat(x2,rl_Psi_reg_prime(u,Psi_constant,*psi_fn),
								 x2,p,p);
		rl_sum_mat(x3,x2,x3,p,p);
		rl_vec_vec_symmetric(x2,x[i],p);
		rl_scalar_mat(x2,w[i],x2,p,p);
		rl_sum_mat(x4,x2,x4,p,p);
		rl_scalar_vec(x[i],rl_Psi_reg_prime(u,Psi_constant,*psi_fn)*u,v_aux,p);
		rl_sum_vec(v,v_aux,v,p);
		u2 += rl_Chi_prime(u,c, *chi_fn) * u;
	};

	rl_scalar_vec(v, beta * (Sfloat) (n-p) * scale / u2 , v, p);
	rl_inverse(x3, x2, p);
	rl_mat_mat(x2, x4, x3, p, p, p);
	rl_mat_vec(x2, v, v2, p, p);
	rl_scalar_mat(x3, scale, x3, p, p);


	/* the correction matrix is now in x3 */
	/* the correction vector is now in v2 */

	/* srand(seed); */

	/* start the bootstrap replications */
	for(i=0;i<m;i++) {
		rl_sampler_i(n,indices);

		/* get pseudo observed y's */
		for(j=0;j<n;j++) yb[j] = x[j][p] = Fi[j] + res[indices[j]];

		/* calculate the robust bootstrap estimate */
		rl_reset_vec(v,p); /* v <- 0 */
		rl_reset_mat(x2,p,p);	/* x2 <- 0 */
		s = 0.0;

		for(j=0;j<n;j++) {
			rl_scalar_vec(x[j],yb[j]*w[indices[j]],v_aux,p);
			rl_sum_vec(v,v_aux,v,p);
			rl_vec_vec_symmetric(x4,x[j],p);
			rl_scalar_mat(x4,w[indices[j]],x4,p,p);
			rl_sum_mat(x2,x4,x2,p,p);
			s += rl_Chi(res_s[indices[j]] / scale , c, *chi_fn);
		};

		s = s * scale / beta / (Sfloat) (n-p);
		rl_inverse(x2,x4,p);		/* x4 <- x2^-1 */
		rl_mat_vec(x4,v,v_aux,p,p);	/* v_aux <- x4 * v */
		rl_dif_vec(v_aux,Beta_m,v_aux,p); 	/* v_aux <- v_aux - beta_m */

		/* v has the robust bootstrapped vector, correct it */
		rl_mat_vec(x3,v_aux,v,p,p);	/* v <- x3 * v_aux */
		rl_scalar_vec(v2,s-scale,v_aux,p);
		rl_sum_vec(v_aux,v,v,p);

		/* store the betas */
		for(j=0;j<p;j++) 
			ours[j*m+i]=v[j];
	};

	return;
}


Sfloat rl_Chi(Sfloat x, Sfloat c, Sint chi_fn)
{ 
	/* chi_fn = 1 -> Bisquare
	 * chi_fn = 2 -> Optimal
	 */

	if( chi_fn == 1) {
		Sfloat t;
		if( fabs(x) > c ) return(1.0);
		else {
			t = x / c;
			return( 3.0*t*t - 3.0*t*t*t*t + t*t*t*t*t*t );
		}
	}

	else {

		Sfloat g1, g2, g3, g4, h1, h2, h3, h4;
		Sfloat tmp,ans;
		Sfloat tmp2, tmp4, tmp6, tmp8;
		g1 = -1.944;
		g2 = 1.728;
		g3 = -.312;
		g4 = 0.016;
		h1 = g1 / 2.;
		h2 = g2 / 4.;
		h3 = g3 / 6.;
		h4 = g4 / 8.;
		tmp = x / c;
		tmp2 = tmp * tmp;
		tmp4 = tmp2 * tmp2;
		tmp6 = tmp2 * tmp4;
		tmp8 = tmp6 * tmp2;

		if( fabs(tmp) > 3.0 )  ans = 3.25 * c * c;

		else if(fabs(tmp) <= 2.0) ans = x * x / 2.;

		else
			ans = c * c * (1.792 + h1 * tmp2 + h2 * tmp4 +
										 h3 * tmp6 + h4 * tmp8);
		return(ans);
	};
}


Sfloat rl_Psi_reg(Sfloat x, Sfloat c, Sint psi_fn)
{
	/* psi_fn = 1 -> Bisquare
	 * psi_fn = 2 -> Optimal
	 */

	if(psi_fn==2) {
		Sfloat g1, g2, g3, g4;
		Sfloat tmp,ans;
		Sfloat tmp2, tmp3, tmp5, tmp7;
		g1 = -1.944;
		g2 = 1.728;
		g3 = -.312;
		g4 = 0.016;
		tmp = x / c;
		tmp2 = tmp * tmp;
		tmp3 = tmp2 * tmp;
		tmp5 = tmp3 * tmp2;
		tmp7 = tmp5 * tmp2;

		if( fabs(tmp) > 3.0 )  ans = 0.;

		else if(fabs(tmp) <= 2.0) ans = x;

		else
			ans = c * (g1 * tmp + g2 * tmp3 +
								 g3 * tmp5 + g4 * tmp7);
		return(ans);
	}

	else {
		if (fabs(x) > c) return(0.0);

		else return( x / c * (1.0-(x/c)*(x/c))*
									(1.0-(x/c)*(x/c))  );
	};
}


Sfloat rl_Chi_prime(Sfloat x, Sfloat c, Sint chi_fn)
{
	/* chi_fn = 1 -> Bisquare
	 * chi_fn = 2 -> Optimal
	 */

	Sfloat g1, g2, g3, g4;
	Sfloat tmp,ans;
	Sfloat tmp2, tmp3, tmp5, tmp7;
	Sfloat t;
	g1 = -1.944;
	g2 = 1.728;
	g3 = -.312;
	g4 = 0.016;
	tmp = x / c;
	tmp2 = tmp * tmp;
	tmp3 = tmp2 * tmp;
	tmp5 = tmp3 * tmp2;
	tmp7 = tmp5 * tmp2;

	if( chi_fn == 1 ) {
		if( fabs(x) > c ) return(0.0);
		else { t = x / c ;
		return( 6.0*t*(1 - t*t) * (1-t*t) / c );
		}
	}

	else {
		if( fabs(tmp) > 3.0 )  ans = 0.;

		else if(fabs(tmp) <= 2.0) ans = x;

		else
			ans = c * (g1 * tmp + g2 * tmp3 +
								 g3 * tmp5 + g4 * tmp7);
		return(ans);
	};
}


Sfloat rl_Psi_reg_prime(Sfloat x, Sfloat c, Sint psi_fn)
{
	/* psi_fn = 1 -> Bisquare
	 * psi_fn = 2 -> Optimal
	 */

	if( psi_fn == 1) {
		if (fabs(x) > c) return(0.0);

		else return( ( 1.0 - (x/c)*(x/c) ) * 
								 ( 1.0 - 5.0 * x * x / c / c ) / c );
	}

	else {
		Sfloat g1, g2, g3, g4;
		Sfloat tmp,ans;
		Sfloat tmp2, tmp4, tmp6;
		g1 = -1.944;
		g2 = 1.728;
		g3 = -.312;
		g4 = 0.016;
		tmp = x / c;
		tmp2 = tmp * tmp;
		tmp4 = tmp2 * tmp2;
		tmp6 = tmp4 * tmp2;

		if ( fabs(tmp) > 3) ans = 0.0;

		else if( fabs(tmp)<2) ans = 1.;

		else
			ans = g1 + 3. * g2 * tmp2 +
				5. * g3 * tmp4 + 7. * g4 * tmp6;

		return(ans);
	};
}


void rl_sampler_i(Sint n, Sint *x)
{
	/* function to get a random sample of
	 * indices (0 to n-1)
	 * *x receives the output
	 * rand() returns an integer between 0 and RAND_MAX
	 */

  /* Sint i;
   * for(i=0;i<n;i++) 
   * x[i] = (Sint) ( (Sfloat) rand() / RAND_MAX * (Sfloat) (n-1) );
   */

  long ignored = 0;
  Sint i = 0;
  Sfloat u = 0.0;

  seed_in(&ignored);

  for(i = 0; i < n; i++) {
    u = unif_rand();
    x[i] = (Sint) (u * ((Sfloat) (n - 1)));
  }

  seed_out(&ignored);
}


void rl_sum_mat(Sfloat **a, Sfloat **b, Sfloat **c, Sint n, Sint m)
{
	register Sint i,j;
	for(i = 0; i < n; i++)
		for(j = 0; j < m; j++) 
			c[i][j] = a[i][j] + b[i][j];
}


void rl_vec_vec_symmetric(Sfloat **a, Sfloat *v1, Sint n)
{
	register Sint i,j;

	/* take advantage of symmetry  */
	for(i = 0; i < n; i++)
		for(j = i; j < n; j++)
			a[i][j] = a[j][i] = v1[i] * v1[j];
}


void rl_scalar_mat(Sfloat **a, Sfloat b, Sfloat **c, Sint n, Sint m)
{
	register Sint i,j;
	for(i = 0; i < n; i++)
		for(j = 0; j < m; j++)
			c[i][j]  = b * a[i][j];
}


void rl_scalar_vec(Sfloat *a, Sfloat b, Sfloat *c, Sint n)
{
	register Sint i;
	for(i = 0; i < n; i++)
		c[i] = b * a[i];
}


void rl_sum_vec(Sfloat *a, Sfloat *b, Sfloat *c, Sint n)
{
	register Sint i;
	for(i = 0; i < n; i++)
		c[i] = a[i] + b[i];
}


void rl_dif_vec(Sfloat *a, Sfloat *b, Sfloat *c, Sint n)
{
	register Sint i;
	for(i = 0; i < n; i++)
		c[i] = a[i] - b[i];
}


void rl_mat_vec(Sfloat **a, Sfloat *b, Sfloat *c, Sint n, Sint m)
{
	register Sint i,j; 
	for(i = 0; i < n; i++) 
		for(c[i] = 0, j = 0; j < m; j++)
			c[i] += a[i][j] * b[j];
}


void rl_mat_mat(Sfloat **a, Sfloat **b, Sfloat **c, Sint n, 
		Sint m, Sint l)
{
	register Sint i, j, k;
	for(i = 0; i < n; i++) 
		for(j = 0; j < l; j++) {
			c[i][j] = 0.0; 
			for(k = 0; k < m; k++)
				c[i][j] += a[i][k] * b[k][j];
		};
}


Sint rl_inverse(Sfloat **a, Sfloat **b, Sint n)
{
	Sint rl_lu(Sfloat **, Sint *, Sfloat *);
	register Sint i, j, k;
	Sfloat **c, *e;

	c = Salloc(n, Sfloat*);
	e = Salloc(n, Sfloat);

	for(i = 0; i < n; i++) 
		c[i] = Salloc(n + 1, Sfloat);

	for(i = 0; i < n; i++) {   /* i-th column */

		for(j=0;j<n;j++)
			for(k = 0; k < n; k++)
				c[j][k] = a[j][k];

		for(j = 0; j < i; j++)
			c[j][n] = 0.0;

		c[i][n] = 1.0;

		for(j = i + 1; j < n; j++)
			c[j][n] = 0.0;

		if( rl_lu(c,&n,e) == 1)
			return(1);

		for(j=0;j<n;j++)
			b[j][i] = e[j];
	};

	return(0);
}


void rl_reset_mat(Sfloat **a, Sint n, Sint m)
{
	register Sint i, j;
	for(i = 0; i < n; i++)
		for(j = 0; j < m; j++)
			a[i][j] = 0.0;
}


void rl_reset_vec(Sfloat *a, Sint n)
{
	register Sint i;
	for(i = 0; i < n; i++)
		a[i] = 0.0;
}


Sint rl_lu(Sfloat **a,Sint *P, Sfloat *x)
{
	Sint *pp, p;
	register Sint i,j,k;
	Sfloat *kk,s;
	p = *P;

	pp = Salloc(p, Sint);

	/* pp vector storing the permutations */
	for(j=0;j<p;j++) { /* cols */
		pp[j]=j;
		for(i=j;i<p;i++)   /* rows */
			if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
				pp[j]=i;
		if ( pp[j] != j ) {  /* permute rows swpping pointers */
			kk=a[j];
			a[j]=a[pp[j]];
			a[pp[j]]=kk;
		};

  /* exit if singular matrix (det=0)
   * if (j,j) pivot ``is zero''  */

		if ( fabs(a[j][j]) < TOL_INVERSE ) {
			return(1);
		};

		for(k=(j+1);k<p;k++)
			a[k][j] = a[k][j] / a[j][j];

		for(k=(j+1);k<p;k++)
			for(i=(j+1);i<p;i++)
				a[k][i] = a[k][i] - a[k][j] * a[j][i];

	}; /* end for(j=0;... */

	for(i=0;i<p;i++) {
		s=0.0;
	  for(j=0;j<i;j++)
	    s += a[i][j] * x[j];
		x[i] = a[i][p] - s;
	};

	for(i=(p-1);i>=0;i--) {
		s=0;
	  for(j=(i+1);j<p;j++)
	    s += a[i][j] * x[j];
	  x[i] = (x[i] - s) / a[i][i];
	};
	return(0);
}




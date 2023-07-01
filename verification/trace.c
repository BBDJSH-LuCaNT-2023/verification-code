#define _GNU_SOURCE // for versionsort()
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <dirent.h>
#include <sys/wait.h>
#include <acb.h>
#include <zlib.h>
#include "kronecker.c"
extern void refine_regulator(arb_t,long,double,long);

const int32_t k=6,sqrtX=331662,nprocs=64;
const uint64_t prec=280;
const double theta = 7.0/64; // this is exact
static char temp_dir[32]="/tmp/trace_XXXXXX";

static int32_t n,nmax,proc,*H12;
static uint64_t X,**binom;
static uint32_t *factor;
static fmpz_t kfac;
static arb_t *trace,*lambda;
static arf_t *Lmin;
static struct dirent **mffile,**cgfile;
static int nmffiles,ncgfiles;

static void remove_temp_dir(void) {
	char buf[64];
	sprintf(buf,"/bin/rm -rf %s",temp_dir);
	system(buf);
}

// compute binomial coefficients and k!
static void binomial_table(void) {
	int32_t j,l,t1,t2;
	uint64_t *p;

	binom = (uint64_t **)malloc((2*k+1)*(sizeof(uint64_t *)
		+(2*k+1)*sizeof(uint64_t)));
	p = (uint64_t *)&binom[2*k+1];
	for (j=0;j<=2*k;j++)
		binom[j] = p, p += 2*k+1;

	binom[0][0] = 1;
	for (j=1;j<=2*k;j++)
		binom[0][j] = 0;
	for (l=1;l<=2*k;l++)
		for (j=0;j<=2*k;j++)
			if (j > l)
				binom[l][j] = 0;
			else {
				if (j == 0) t1 = 0; else t1 = binom[l-1][j-1];
				if (j == l) t2 = 0; else t2 = binom[l-1][j];
				binom[l][j] = t1+t2;
			}

	// compute k!
	fmpz_init(kfac);
	fmpz_set_ui(kfac,1);
	for (j=2;j<=k;j++)
		fmpz_mul_ui(kfac,kfac,j);
}

// compute X and table of first factors
static void factor_table(void) {
	uint32_t p,t;

	X = (uint64_t)sqrtX*sqrtX;
	factor = (uint32_t *)calloc(sqrtX,sizeof(factor[0]));
	factor--;
	for (p=sqrtX;p>=2;p--)
		for (t=p;t<=sqrtX;t+=p)
			factor[t] = p;
	factor[1] = 1;
}

static void hurwitz_table(void) {
	int32_t s,t,N;

	H12 = (int32_t *)malloc(2*(sqrtX+1)*sizeof(H12[0]));
	for (N=0;N<2*sqrtX;) {
		// N == 0 (mod 4)
		for (s=1;(t=N+1-s*s)>0;s+=2) {
			H12[N] -= H12[t];
			if (!(t % s)) H12[N] += (t/s-s)<<1;
		}
		if (!t) H12[N] -= (H12[0]+s);
		N += 3;

		// N == 3 (mod 4)
		for (s=0;(t=N-s*s)>0;s+=2)
			H12[N] -= H12[t];
		for (s=1;(t=N-s*s)>0;s+=2)
			if (!(t % s)) H12[N] += (t/s-s)<<1;
		H12[N] <<= 1;
		N++;
	}
}

// compute nmax and allocate trace, lambda and Lmin arrays
static void coeff_arrays(void) {
	int32_t m;

	nmax = (sqrtX-1)>>1;
	trace = (arb_t *)malloc((2*nmax+1)*sizeof(arb_t));
	trace += nmax;
	for (m=-nmax;m<=nmax;m++)
		if (m) arb_init(trace[m]);

	lambda = (arb_t *)malloc(nmax*sizeof(arb_t));
	lambda--;
	for (m=1;m<=nmax;m++) 
		arb_init(lambda[m]);
	arb_set_ui(lambda[1],1);

	Lmin = (arf_t *)malloc((2*nmax+1)*sizeof(arf_t));
	Lmin += nmax;
	for (m=-nmax;m<=nmax;m++) {
		if (!m) continue;
		arf_init(Lmin[m]);
		arf_one(Lmin[m]);
	}
}

static void h(arb_t res,const arb_t r,int rnonzero) {
	static int init;
	static arb_t A,Am2,Am2j,logA,tr;
	static acb_t p,s,tc;
	int32_t j,l;

	if (!init) {
		arb_init(A); arb_init(Am2); arb_init(Am2j);
		arb_init(logA); arb_init(tr);
		acb_init(p); acb_init(s); acb_init(tc);
		init = 1;
	}

	// compute A = (X+2n+sqrt(X^2+4nX))/2|n|
	arb_set_ui(A,X+4*n);
	arb_mul_ui(A,A,X,prec);
	arb_sqrt(A,A,prec);
	arb_add_ui(A,A,X+2*n,prec);
	arb_div_ui(A,A,2*abs(n),prec);

	// compute A^-2, log(A)
	arb_inv(Am2,A,prec);
	arb_sqr(Am2,Am2,prec);
	arb_log(logA,A,prec);

	arb_set_ui(Am2j,1);
	if (rnonzero) {
		arb_neg(Am2,Am2);
		acb_zero(s);
		for (j=0;j<=k;j++) {
			acb_set_fmpz(p,kfac);
			arb_set(acb_imagref(tc),r);
			for (l=-j;l<=k-j;l++) {
				// p /= (l+I*r);
				arb_set_si(acb_realref(tc),l);
				acb_div(p,p,tc,prec);
			}
			arb_mul_ui(tr,Am2j,binom[k][j],prec);
			acb_mul_arb(p,p,tr,prec);
			acb_add(s,s,p,prec);
			arb_mul(Am2j,Am2j,Am2,prec);
		}
		arb_zero(acb_realref(tc));
		arb_mul(acb_imagref(tc),r,logA,prec);
		acb_exp(tc,tc,prec);
		acb_mul(s,s,tc,prec);
		arb_set(res,acb_realref(s));
	} else {
		arb_zero(res);
		for (j=0;j<=k;j++) {
			arb_set(acb_realref(s),logA);
			for (l=j-k;l<=j;l++) {
				if (!l) continue;
				arb_set_ui(tr,1);
				arb_div_si(tr,tr,l,prec);
				arb_add(acb_realref(s),acb_realref(s),tr,prec);
			}
			arb_mul_ui(tr,acb_realref(s),binom[k][j],prec);
			arb_mul_ui(tr,tr,binom[k][j],prec);
			arb_mul(tr,tr,Am2j,prec);
			arb_add(res,res,tr,prec);
			arb_mul(Am2j,Am2j,Am2,prec);
		}
	}

	// multiply by 2*(A|n|/X)^k
	arb_mul_ui(tr,A,abs(n),prec);
	arb_div_ui(tr,tr,X,prec);
	arb_pow_ui(tr,tr,k,prec);
	arb_mul_2exp_si(tr,tr,1);
	arb_mul(res,res,tr,prec);
}

static void spectral_terms(const char *dir) {
	int32_t i,parity,m,p,t;
	arb_t r,e,x,y,temp;
	gzFile fp;
	char buf[256];

	arb_init(temp);
	arb_init(r); arb_init(e);
	arb_init(x); arb_init(y);

	for (i=proc;i<nmffiles;i+=nprocs) {
		sprintf(buf,"%s/%s",dir,mffile[i]->d_name);
		if ( !(fp=gzopen(buf,"r")) ) {
			fprintf(stderr,"failed to open %s for reading\n",buf);
			_exit(1);
		}

		// skip first line
		gzgets(fp,buf,sizeof(buf));

		// second line is target accuracy
		// we use this as error interval rather than reported accuracy,
		//  to err on the pessimistic side
		gzgets(fp,buf,sizeof(buf));
		t = atoi(buf);
		arb_set_ui(temp,1);
		arb_mul_2exp_si(temp,temp,-t);
		arb_neg(e,temp);
		arb_union(e,e,temp,prec);

		// parity
		gzgets(fp,buf,sizeof(buf));
		parity = atoi(buf);

		// number of coefficients
		gzgets(fp,buf,sizeof(buf));
		t = atoi(buf);
		if (t < nmax) {
			// ensure that we have enough coefficients
			// to process all discriminants up to X
			fprintf(stderr,"not enough coefficients in %s\n",mffile[i]->d_name);
			_exit(1);
		}

		// r
		gzgets(fp,buf,sizeof(buf));
		arb_set_str(x,buf,prec);
		gzgets(fp,buf,sizeof(buf));
		arb_set_str(y,buf,prec);
		arb_union(r,x,y,prec);

		// read prime power coeffs and compute composite ones
		for (m=2;m<=nmax;m++) {
			p = factor[m];
			for (t=m/p;factor[t]==p;t/=p);
			if (t == 1) {
				gzgets(fp,buf,sizeof(buf));
				arb_set_str(lambda[m],buf,prec);
				arb_add(lambda[m],lambda[m],e,prec);
			} else
				arb_mul(lambda[m],lambda[t],lambda[m/t],prec);
		}
		do gzgets(fp,buf,sizeof(buf)); while (!gzeof(fp));

		// check that reported error is within target
		arb_set_str(y,buf,prec);
		arb_neg(x,y);
		arb_union(x,x,y,prec);
		if (!arb_contains(e,x)) {
			fprintf(stderr,"reported error is greater than target in %s\n",
				mffile[i]->d_name);
			_exit(1);
		}
		gzclose(fp);

		for (n=-nmax;n<=nmax;n++) {
			if (!n) continue;
			// do this inside the loop since h depends implicitly on n
			h(temp,r,1);
			if (parity && n < 0) arb_neg(temp,temp);
			arb_mul(temp,temp,lambda[abs(n)],prec);
			arb_add(trace[n],trace[n],temp,prec);
		}
	}

	// add error terms if we have the last file
	if ((nmffiles-1-proc) % nprocs == 0) {
		arb_set_d(x,theta);
		// compute b(n) for all n
		for (m=2;m<=nmax;m++) {
			p = factor[m];
			for (t=m/p;factor[t]==p;t/=p);
			if (t == 1) {
				// compute b(m) for a prime power
				arb_log_ui(temp,(int64_t)m*p,prec);
				arb_mul(temp,temp,x,prec);
				arb_sinh(lambda[m],temp,prec);
				arb_log_ui(temp,p,prec);
				arb_mul(temp,temp,x,prec);
				arb_sinh(temp,temp,prec);
				arb_div(lambda[m],lambda[m],temp,prec);
			} else
				arb_mul(lambda[m],lambda[t],lambda[m/t],prec);
		}

		// add truncation error, with R = r from the last file
		arb_pow_ui(e,r,k-1,prec);
		arb_inv(e,e,prec);
		arb_mul_fmpz(e,e,kfac,prec);
		arb_div_ui(e,e,3*(k-1),prec);
		for (m=1;m<=nmax;m++) {
			arb_mul(r,e,lambda[m],prec);
			arb_neg(temp,r);
			arb_union(temp,temp,r,prec);
			arb_add(trace[m],trace[m],temp,prec);
			arb_add(trace[-m],trace[-m],temp,prec);
		}
	}
}

// integral of f(sqrt(y^2-min(4n,0)))
static void i3(arb_t res) {
	static int init;
	static arb_t temp;
	int32_t j;

	if (!init) {
		arb_init(temp);
		init = 1;
	}
	if (n < 0) {
		arb_set_ui(temp,X+4*n);
		arb_div_ui(temp,temp,X,prec);
		arb_pow_ui(res,temp,k,prec);
		arb_sqrt(temp,temp,prec);
		arb_mul(res,res,temp,prec);
	} else
		arb_set_ui(res,1);
	arb_mul_ui(res,res,2*sqrtX,prec);
	for (j=1;j<=k;j++) {
		arb_mul_ui(res,res,2*j,prec);
		arb_div_ui(res,res,2*j+1,prec);
	}
}

// compute integral of f(y)/(y+m) with input m and f = f(m)
// assumes X >= (a-n/a)^2
static void i4(arb_t res,int32_t m,arb_srcptr f) {
	static int init;
	static arb_t c,mc2,mc2j,temp,inner,cml;
	int32_t j,l;

	if (!init) {
		// allocate memory for temp variables
		arb_init(c); arb_init(mc2); arb_init(mc2j);
		arb_init(temp); arb_init(inner); arb_init(cml);
		init = 1;
	}

	// compute c = m/sqrt(X) and mc2 = -c^2
	arb_set_ui(c,m);
	arb_div_ui(c,c,sqrtX,prec);
	arb_sqr(mc2,c,prec);
	arb_neg(mc2,mc2);

	// compute f(m)*log((1+c^{-1})/2)
	arb_inv(temp,c,prec);
	arb_add_ui(temp,temp,1,prec);
	arb_mul_2exp_si(temp,temp,-1);
	arb_log(temp,temp,prec);
	arb_mul(res,f,temp,prec);

	// evaluate nested sum
	arb_set(mc2j,mc2);
	arb_inv(cml,c,prec);
	arb_zero(inner);

	for (j=l=1;j<=k;j++) {
		// add term to running inner sum
		arb_sub_ui(temp,cml,1,prec);
		arb_div_ui(temp,temp,l,prec);
		arb_add(inner,inner,temp,prec);
		arb_div(cml,cml,c,prec);
		l++;

		// subtract term from running inner sum
		arb_sub_ui(temp,cml,1,prec);
		arb_div_ui(temp,temp,l,prec);
		arb_sub(inner,inner,temp,prec);
		arb_div(cml,cml,c,prec);
		l++;

		// add one term to the j sum and compute next (-c^2)^j
		arb_mul_ui(temp,mc2j,binom[k][j],prec);
		arb_mul(temp,temp,inner,prec);
		arb_sub(res,res,temp,prec);
		arb_mul(mc2j,mc2j,mc2,prec);
	}
}

// integral of f(y)/(y^2+|D|)
static void i5(arb_t res,uint64_t absD) {
	static int init;
	static arb_t c,c2,c2j,temp,inner,c1m2l;
	int32_t j,l;

	if (!init) {
		// allocate memory for temp variables
		arb_init(c); arb_init(c2); arb_init(c2j);
		arb_init(temp); arb_init(inner); arb_init(c1m2l);
		init = 1;
	}

	// compute c = sqrt(|D|/X) and c2 = c^2
	arb_set_ui(c2,absD);
	arb_div_ui(c2,c2,X,prec);
	arb_sqrt(c,c2,prec);

	// compute (1+c^2)^k*arctan(1/c)
	arb_add_ui(temp,c2,1,prec);
	arb_pow_ui(res,temp,k,prec);
	arb_inv(c1m2l,c,prec);
	arb_atan(temp,c1m2l,prec);
	arb_mul(res,res,temp,prec);

	// evaluate nested sum
	arb_set(c2j,c2);
	arb_zero(inner);

	for (j=l=1;j<=k;j++) {
		// add term to running inner sum
		arb_div_si(temp,c1m2l,(2*l-1)*(2*(l&1)-1),prec);
		arb_add(inner,inner,temp,prec);
		arb_div(c1m2l,c1m2l,c2,prec);
		l++;

		// add one term to the j sum and compute next c^{2j}
		arb_mul_ui(temp,c2j,binom[k][j],prec);
		arb_mul(temp,temp,inner,prec);
		arb_sub(res,res,temp,prec);
		arb_mul(c2j,c2j,c2,prec);
	}
	arb_mul_2exp_si(res,res,1);
	arb_rsqrt_ui(temp,absD,prec);
	arb_mul(res,res,temp,prec);
}

// compute integral of 2a(f(v-n/v)-f(a-n/a))/(v^2-a^2)
//   with input a and f = f(a-n/a)
// assumes X >= (a-n/a)^2
void i6(arb_t res,int32_t a,arb_srcptr f) {
	static int init;
	static arb_t temp,b,b2,a2n,a2nm,b2lm1,lsum,jsum,nX;
	int32_t m,j,sign;

	if (!init) {
		arb_init(b); arb_init(b2); arb_init(b2lm1);
		arb_init(a2n); arb_init(a2nm); arb_init(nX);
		arb_init(lsum); arb_init(jsum); arb_init(temp);
		init = 1;
	}

	arb_zero(res);
	// handle the borderline case separately to avoid taking log of 0
	if (a >= abs(n/a) && a-n/a == sqrtX) return;

	// b = (sqrt(X)+sqrt(X+4*n))/(2*a);
	arb_set_ui(b,X+4*n);
	arb_sqrt(b,b,prec);
	arb_add_ui(b,b,sqrtX,prec);
	arb_div_ui(b,b,2*a,prec);

	// b2 = b^2
	arb_sqr(b2,b,prec);

	// a2n = -a^2/n
	arb_set_ui(a2n,(uint64_t)a*a);
	arb_div_si(a2n,a2n,-n,prec);

	// nX = n/X
	arb_set_si(nX,n);
	arb_div_ui(nX,nX,X,prec);

	for (sign=-1;sign<=1;sign+=2) {
		// flip the sign of m
		arb_inv(a2n,a2n,prec);
		arb_inv(b,b,prec);
		arb_inv(b2,b2,prec);

		arb_set(a2nm,a2n);
		arb_set(b2lm1,b);
		arb_zero(lsum);
		for (m=1;m<=k;m++) {
			// compute sum over j by Horner's method
			arb_zero(jsum);
			for (j=k;j>=0;j--) {
				arb_mul(jsum,jsum,nX,prec);
				arb_add_ui(jsum,jsum,binom[k][j]*binom[j+j][j+m],prec);
			}

			// compute next running sum over l
			arb_sub_ui(temp,b2lm1,1,prec);
			arb_div_ui(temp,temp,2*m-1,prec);
			arb_add(lsum,lsum,temp,prec);

			arb_mul(temp,jsum,lsum,prec);
			arb_mul(temp,temp,a2nm,prec);
			arb_add(res,res,temp,prec);

			// compute next (-a^2/n)^m and b^(2l-1)
			arb_mul(a2nm,a2nm,a2n,prec);
			arb_mul(b2lm1,b2lm1,b2,prec);
		}
	}

	// compute f(a-n/a)*log((b+1)/(b-1))
	arb_sub_ui(temp,b,1,prec);
	arb_add_ui(b,b,1,prec);
	arb_div(b,b,temp,prec);
	arb_log(b,b,prec);
	arb_mul(temp,f,b,prec);

	// final result
	arb_mul_2exp_si(res,res,1);
	arb_sub(res,res,temp,prec);
}

// 1/3 times integral of (f(0)-f(y))/y^2
static void i7(arb_t res) {
	int32_t j;
	arb_set_ui(res,2*k+1);
	arb_div_ui(res,res,3*sqrtX,prec);
	for (j=1;j<=k;j++) {
		arb_mul_ui(res,res,2*j,prec);
		arb_div_ui(res,res,2*j+1,prec);
	}
}

// 1/2 times integral of (f(y)+f(1/y)-f(0))/y
static void i8(arb_t res) {
	static int init;
	static arb_t temp;
	int32_t j;

	if (!init) {
		arb_init(temp);
		init = 1;
	}
	arb_log_ui(res,X,prec);
	for (j=1;j<=k;j++) {
		arb_set_ui(temp,j);
		arb_inv(temp,temp,prec);
		arb_sub(res,res,temp,prec);
	}
	arb_mul_2exp_si(res,res,-1);
}

static void S(arb_t res,int32_t m,arb_srcptr f) {
	static int init;
	static arb_t c1,c2,temp;
	int32_t r,p,pr,t;

	if (!init) {
		// c1 = Euler + log(4*Pi)
		// c2 = 2*Euler + log(2*Pi)
		arb_init(c1); arb_init(c2); arb_init(temp);
		arb_const_euler(temp,prec);
		arb_const_log_sqrt2pi(c2,prec);
		arb_add(c2,c2,temp,prec);
		arb_mul_2exp_si(c2,c2,1);
		arb_sub(c1,c2,temp,prec);
		arb_const_log2(temp,prec);
		arb_add(c1,c1,temp,prec);
		init = 1;
	}

	if (m) {
		arb_set(res,c1);
		for (t=m;t>1;t/=pr) {
			pr = p = factor[t], r = 1;
			while (factor[t/pr] == p)
				pr *= p, r++;
			arb_log_ui(temp,p,prec);
			arb_mul_ui(temp,temp,r*pr-(pr-1)/(p-1),prec);
			arb_div_ui(temp,temp,pr,prec);
			arb_add(res,res,temp,prec);
		}
		arb_mul(res,res,f,prec);
		i4(temp,m,f);
		arb_add(res,res,temp,prec);
	} else {
		arb_set(res,c2);
		i7(temp);
		arb_add(res,res,temp,prec);
		i8(temp);
		arb_add(res,res,temp,prec);
	}
}

static void psum(arb_t res,int32_t a) {
	static int init;
	static fmpz_t m2,z1,z2;
	static arb_t t,s;
	int32_t d=n/a,p,z;
	int64_t m;

	// z = (sqrtX+(int32_t)sqrtl((long double)(X+4*n)))/(2*a);
	z = sqrtX/a; // suffices since 4n < sqrt(X)
	if (!init) {
		arb_init(t);
		arb_init(s);
		init = 1;
	}

	arb_zero(res);
	for (p=2;p<=z;p++) {
		if (factor[p] != p) continue;
		arb_zero(s);
		for (m=p;m<=z;m*=p) {
			arb_set_si(t,d);
			arb_div_ui(t,t,m,prec);
			arb_sub_ui(t,t,a*m,prec);
			arb_sqr(t,t,prec);
			arb_div_ui(t,t,X,prec);
			arb_sub_ui(t,t,1,prec);
			if (arb_is_nonnegative(t)) continue;
			arb_neg(t,t);
			arb_pow_ui(t,t,k,prec);
			arb_div_ui(t,t,m,prec);
			arb_add(s,s,t,prec);
		}
		arb_log_ui(t,p,prec);
		arb_mul(t,s,t,prec);
		arb_add(res,res,t,prec);
	}
	arb_mul_2exp_si(res,res,1);
}

static void divisor_sum(arb_t res) {
	static int init;
	static arb_t temp,f;
	int32_t a,m,sigma0,sigma1;

	if (!init) {
		arb_init(temp);
		arb_init(f);
		init = 1;
	}

	arb_zero(res);
	sigma0 = sigma1 = 0;
	for (a=1;a<=abs(n);a++) {
		if (n % a) continue;

		// count of divisors and sum of divisors
		sigma0++;
		sigma1 += a;

		// compute f(a-n/a)
		m = abs(a-n/a);
		arb_set_ui(f,(int64_t)m*m);
		arb_div_ui(f,f,X,prec);
		arb_sub_ui(f,f,1,prec);
		arb_neg(f,f);
		arb_pow_ui(f,f,k,prec);

		// add up main terms
		psum(temp,a);
		arb_add(res,res,temp,prec);
		i6(temp,a,f);
		arb_add(res,res,temp,prec);
		S(temp,m,f);
		arb_add(res,res,temp,prec);
	}

	// res -= sigma_0(|n|)*h(0)/4
	h(temp,temp,0);
	arb_mul_ui(temp,temp,sigma0,prec);
	arb_mul_2exp_si(temp,temp,-2);
	arb_sub(res,res,temp,prec);

	// res -= sigma_{-1}(|n|)*i3
	i3(temp);
	arb_mul_ui(temp,temp,sigma1,prec);
	arb_div_ui(temp,temp,abs(n),prec);
	arb_sub(res,res,temp,prec);
}

// assumes n > 0
static void negative_disc_sum(arb_t res) {
	static int init;
	static arb_t temp;
	int32_t t,N;

	if (!init) {
		arb_init(temp);
		init = 1;
	}
	arb_zero(res);
	for (t=0;(N=4*n-t*t)>0;t++) {
		i5(temp,N);
		arb_mul_si(temp,temp,H12[N],prec);
		arb_div_ui(temp,temp,t?6:12,prec);
		arb_add(res,res,temp,prec);
	}
}

static void geometric_terms(void) {
	int32_t i;
	arb_t temp;
	arb_init(temp);
	srandom(42); // ensure all processes generate the same sequence
	for (i=0;;i++) {
		n = i*nprocs+((uint64_t)random()+proc)%nprocs+1;
		if (n > nmax) break;
		divisor_sum(temp); arb_sub(trace[n],trace[n],temp,prec);
		negative_disc_sum(temp); arb_sub(trace[n],trace[n],temp,prec);
		n = -n;
		divisor_sum(temp); arb_sub(trace[n],trace[n],temp,prec);
	}
	arb_clear(temp);
}

static void printarf(FILE *fp,const arf_t x) {
	static int init;
	static fmpz_t m,e;

	if (!init) {
		fmpz_init(m); fmpz_init(e);
		init = 1;
	}
	arf_get_fmpz_2exp(m,e,x);
	fmpz_fprint(fp,m); fprintf(fp,"\n");
	fmpz_fprint(fp,e); fprintf(fp,"\n");
}

static void readarf(FILE *fp,arf_t x) {
	static int init;
	static fmpz_t m,e;

	if (!init) {
		fmpz_init(m); fmpz_init(e);
		init = 1;
	}
	fmpz_fread(fp,m);
	fmpz_fread(fp,e);
	arf_set_fmpz_2exp(x,m,e);
}

static void print_traces(void) {
	int32_t m;
	arf_t a;
	FILE *fp;
	char fname[128];

	sprintf(fname,"%s/%d",temp_dir,(int)proc);
	if ( !(fp=fopen(fname,"w")) ) {
		fprintf(stderr,"failed to open %s for writing\n",fname);
		_exit(1);
	}

	arf_init(a);
	for (m=-nmax;m<=nmax;m++) {
		if (!m) continue;
		arb_get_ubound_arf(a,trace[m],prec);
		printarf(fp,a);
		printarf(fp,Lmin[m]);
	}
	arf_clear(a);
	fclose(fp);
}

static void collate(void) {
	int status;
	int32_t m;
	arf_t temp,emax;
	double d;
	FILE *fp;
	char fname[128];

	// wait for child processes to finish
	while (wait(&status) > 0)
		if (WIFEXITED(status) && WEXITSTATUS(status)) {
			fprintf(stderr,"unexpected exit status of child process\n");
			exit(1);
		}
	
	arf_init(temp); arf_init(emax);
	for (proc=0;proc<nprocs;proc++) {
		sprintf(fname,"%s/%d",temp_dir,(int)proc);
		if ( !(fp=fopen(fname,"r")) ) {
			fprintf(stderr,"failed to open %s for reading\n",fname);
			exit(1);
		}
		for (m=-nmax;m<=nmax;m++) {
			if (!m) continue;
			readarf(fp,temp);
			arf_add(arb_midref(trace[m]),arb_midref(trace[m]),temp,prec,ARF_RND_CEIL);
			readarf(fp,temp);
			arf_min(Lmin[m],Lmin[m],temp);
		}
		fclose(fp);
		//unlink(fname);
	}

	// compute worst case error
	for (m=-nmax;m<=nmax;m++) {
		if (!m) continue;
		arf_div(temp,arb_midref(trace[m]),Lmin[m],prec,ARF_RND_CEIL);
		arf_max(emax,emax,temp);
	}
	arf_root(temp,emax,k,prec,ARF_RND_CEIL);
	arf_neg(temp,temp);
	arf_add_ui(temp,temp,1,prec,ARF_RND_FLOOR);
	arf_mul_ui(temp,temp,X,prec,ARF_RND_FLOOR);
	printf("proved up to %ld\n",arf_get_si(temp,ARF_RND_FLOOR));

	arf_clear(temp); arf_clear(emax);
}

static void sub_from_trace(int64_t D,arb_srcptr a,arf_srcptr L,int32_t l) {
	int32_t t,m;

	t = 0;
	if (D > 4*nmax)
		t = (uint32_t)sqrtl((long double)(D-4*nmax-1))+1;
	t += ((D-t)&1);
	for (;(m=((int64_t)t*t-D)>>2)<=nmax;t+=2) {
		arb_submul_ui(trace[m],a,1+(t!=0),prec);
		if (l == 1) arf_min(Lmin[m],Lmin[m],L);
	}
}

static int64_t arithmetic_factor(int64_t d,int32_t l) {
	int64_t a;
	int32_t m,p,q;

	for (m=l,a=1;m>1;) {
		p = factor[m], q = 1;
		do m /= p, q *= p; while (factor[m] == p);
		a *= 1+(q-1)/(p-1)*(p-kronecker(d,p));
	}
	return a;
}

// process the next line of the file
static int do_one_line(char *buf) {
	static int init;
	static arb_t temp,L;
	static arf_t lbound;
	char *s;
	double r;
	int32_t l;
	int64_t h,d,D;

	if (!init) {
		arb_init(temp); arb_init(L);
		arf_init(lbound);
		init = 1;
	}
	d = atol(strtok(buf," "));  // discriminant
	if (d >= X) return 0;
	r = atof(strtok(NULL," ")); // regulator
	h = 1;
	for (s=strtok(NULL,"[, ]\n");s;s=strtok(NULL,", ]\n"))
		h *= atol(s);

	refine_regulator(L,d,r,prec);
	arb_mul_ui(L,L,2*h,prec);
	arb_rsqrt_ui(temp,d,prec);
	arb_mul(L,L,temp,prec);
	arb_get_lbound_arf(lbound,L,prec);

	for (l=1;(D=d*l*l)<X;l++) {
		arb_set_ui(temp,X-D);
		arb_div_ui(temp,temp,X,prec);
		arb_pow_ui(temp,temp,k,prec);
		arb_mul(temp,temp,L,prec);
		arb_mul_ui(temp,temp,arithmetic_factor(d,l),prec);
		arb_div_ui(temp,temp,l,prec);
		sub_from_trace(D,temp,lbound,l);
	}

	return 1;
}

static void class_number_sums(const char *dir) {
	gzFile fp;
	char buf[256];
	int32_t i;

	for (i=proc;i<ncgfiles;i+=nprocs) {
		sprintf(buf,"%s/%s",dir,cgfile[i]->d_name);
		if ( !(fp=gzopen(buf,"r")) ) {
			fprintf(stderr,"failed to open %s for reading\n",buf);
			_exit(1);
		}
		while (gzgets(fp,buf,sizeof(buf)) && do_one_line(buf));
		gzclose(fp);
	}
}

static int filter_hidden_files(const struct dirent *f) {
	return f->d_name[0] != '.';
}

int main(int argc,char *argv[]) {
	if (argc != 3
	|| (nmffiles=scandir(argv[1],&mffile,filter_hidden_files,alphasort)) < 0
	|| (ncgfiles=scandir(argv[2],&cgfile,filter_hidden_files,versionsort)) < 0) {
		printf("usage: %s <modular form directory> <class group directory>\n",argv[0]);
		return 0;
	}

	binomial_table();
	factor_table();
	hurwitz_table();
	coeff_arrays();

	mkdtemp(temp_dir);
	atexit(remove_temp_dir);
	signal(SIGTERM,exit);

	for (proc=0;proc<nprocs;proc++)
	if (!fork()) {
		spectral_terms(argv[1]);
		geometric_terms();
		class_number_sums(argv[2]);
		print_traces();
		_exit(0);
	}
	collate();
	return 0;
}

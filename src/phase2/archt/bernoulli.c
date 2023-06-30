#ifndef __BERNOULLI__
#ifdef  __BERNOULLI__
gcc -O2 bernoulli.c -lm -lgmp -lmpfr -lmpfi -DBERNOULLI_MAIN -o bernoulli
./bernoulli
exit

//
// Algorithms for computing Bernoulli numbers.
//

Copyright (C) 2011--2014 University of Bristol, UK
  and     (C) 2011--2014 University of Uppsala, Sweden
---------------------------------------------------------------
    This file is part of archt.

    archt is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    archt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with archt.  If not, see <http://www.gnu.org/licenses/>.


Authors
-------
Andrew R. Booker, Andreas Strombergsson, Holger Then


Version: 1.0.3 (15 May 2014)
----------------------------

)))))
#endif
#define __BERNOULLI__
#include"archt.c"
#include"primes.c"

//
// Akiyama-Tanigawa algorithm for Bernoulli numbers.
// Slow, but exact.
//
// References:
//   [1] S. Akiyama, Y. Tanigawa, The Ramanujan Journal 5 (2001) 327--351.
//   [2] M. Kaneko, Journal of Integer Sequences 3 (2000) 00.2.9.
//
// Algorithm 1: a[0,m] := 1/(1+m)  forall  0<=m<=N;
//              for(n=0; n<N; n++) {
//                a[n+1,m] := (m+1)*( a[n,m] - a[n,m+1])  forall  0<=m<N-n; }.
//
// Theorem 1: After running algorithm 1, we have a[n,0]==B_n forall 0<=n<=N,
//            where B_n is the n-th Bernoulli number, but B_1=1/2.
// Proof: See references [1] and [2].
//
// Lemma 1: Algorithm 1 yields a[1,m-1]==a[0,m] forall 0<m<=N.
// Proof: trivial.
//
// Algorithm 2: b[l] := 1/(1+l)  forall 0<=l<=N;
//              for(n=1; n<N; n++) {
//                for(l=N; l>n; l--) {
//                  b[l] := (l-n)*( b[l-1] - b[l]); } }.
//
// Theorem 2: Algorithm 2 results in b[n]==B_n forall 0<=n<=N.
// Poof: Write a[n,m]=:b[n+m] and m+n+1=:l in algorithm 1.
//       With the descending inner loop the evaluation order is kept straight.
//       Lemma 1 allows to start the outer loop at n=1.
//

void bern_aki(mpq_t b[],unsigned long N) { // Algorithm 2.
  unsigned long n,l;
  for(l=0; l<=N; l++) { mpq_init(b[l]); mpq_set_si(b[l],1l,l+1); }
  for(n=1; n<N; n++) {
    for(l=N; l>n; l--) {
      mpq_sub(b[l],b[l-1],b[l]);
      mpz_mul_ui(mpq_numref(b[l]),mpq_numref(b[l]),(l-n));
      mpq_canonicalize(b[l]); } } }

void clear_bern_aki(mpq_t b[],unsigned long N) {
  for(N++; N-->0; ) mpq_clear(b[N]); }

//
// Algorithm for computing Bernoulli numbers.
//
// References:
//   [1] H. Cohen and K. Belabas, the PARI-2.2.11.alpha source code.
//   [2] K. J. McGown, http://wstein.org/projects/168/kevin_mcgown/bernproj.pdf
//

int bernfrac_(mpz_t a,mpz_t d,const unsigned long m,mpfi_t K,
	      unsigned long *l,unsigned long **p,unsigned long *primelimit) {
  int e=0;
  long s;
  unsigned long i,j,N,Kprec;
  double lcKd;
  mpz_t a2;
  mpfi_t z,t;
  mpfr_t c,c2;

  /* check whether m is reasonable */
  if(m<2) {
    if(m==0) { mpz_set_si(a,1l); mpz_set_si(d,1l); /* B_0=1 */ }
    else { mpz_set_si(a,-1l); mpz_set_si(d,2l); /* B_1=-1/2 */ }
    return 0;
  }
  if(m%2==1) { mpz_set_si(a,0l); mpz_set_si(d,1l); /* B_m=0 */ return 0; }

  /* make sure there are enough primes tabulated */
  if(m+1>*primelimit) { e|=primes(l,p,m+1);
    if((e&3)==0) *primelimit=m+1; else { if((e&3)==3) *primelimit=0; return e; }
  }

  /* d = prod_{p-1|m} p */
  mpz_set_si(d,1l);
  for(i=0;i<*l && (j=**p-1)<=m; i++,(*p)++) {
    if(m%j==0) mpz_mul_ui(d,d,**p);
  } *p-=i;

  // bounds
  //   const = precision overhead >= 256
  //   K = 2m!/(2pi)^m < 8((m+1)/(2pi e))^(m+1/2)
  //   z = prod_{p<=N} (1-p^-m)^-1 < 2

  /* lcKd = log( const*K*d) */ mpz_get_d_2exp(&s,d);
  lcKd = (m+.5)*log((m+1)/(2*M_PI*M_E))+(s+11)*M_LN2;
  /* N = ceil( (const*K*d)^(1/(m-1))) */ N=ceil(exp(lcKd/(m-1)));
  /* Kprec = ceil( log( d*K*z)/log(2.)) */ Kprec=ceil(lcKd/M_LN2);

  /* 2^Kprec *= (1 + 2^-Kprec)^m < exp(m * 2^-Kprec) */
  Kprec+=ceil(m*exp(-M_LN2*Kprec)/M_LN2)+64; /* plus 64 bits extra */

  /* check whether Kprec is reasonable */
  if(Kprec<MPFR_PREC_MIN) Kprec=MPFR_PREC_MIN;
  if(Kprec>mpfi_get_prec(K)) return 4;

  /* make sure there are enough primes tabulated */
  if(N>*primelimit) { e|=primes(l,p,N);
    if((e&3)==0) *primelimit=N; else { if((e&3)==3) *primelimit=0; return e; }
  }

  /* z = prod_{p<=N} (1-p^-m)^-1 */
  mpfi_init2(z,(mp_prec_t)Kprec);
  mpfi_init2(t,(mp_prec_t)Kprec);
  mpfi_set_si(z,1l);
  for(i=0; i<*l && **p<=N; i++,(*p)++) {
    mpz_ui_pow_ui(a,**p,m);
    mpfi_set_z(t,a);
    mpfi_inv(t,t);
    mpfi_si_sub(t,1l,t);
    mpfi_mul(z,z,t);
  } *p-=i;
  mpfi_inv(z,z);

  /* add the error bound [0, 1/((m-1)*N^(m-1))] to z */
  mpfr_init2(c,(mpfr_prec_t)Kprec);
  mpfr_init2(c2,(mpfr_prec_t)Kprec);
  mpz_ui_pow_ui(a,N,(unsigned long)(m-1));
  mpz_mul_ui(a,a,(unsigned long)(m-1));
  mpfr_set_z(c2,a,MPFR_RNDZ);
  mpfr_si_div(c2,1l,c2,MPFR_RNDU);
  mpfr_set_si(c,0l,MPFR_RNDZ);
  mpfi_interv_fr(t,c,c2);
  mpfi_add(z,z,t);
  mpfi_clear(t);

  /* a = (-1)^(m/2+1) integer(d*K*z) */
  mpz_init(a2);
  mpfi_mul_z(z,z,d);
  mpfi_mul(z,K,z);
  mpfi_get_right(c2,z);
  mpfi_get_left(c,z);
  mpfr_get_z(a2,c2,MPFR_RNDD);
  mpfr_get_z(a,c,MPFR_RNDU);
  mpfr_clear(c2); mpfr_clear(c); mpfi_clear(z);
  if(mpz_cmp(a,a2)!=0) { mpz_clear(a2); return 8; } else mpz_clear(a2);
  if(m%4==0) mpz_neg(a,a);
  return 0;
}

//
// Algorithm for approximating Bernoulli numbers.
//
// Method:
//   Expressing the Bernoulli numbers in terms of the Riemann zeta function,
//   B_m=(-1)^(m/2+1)*K*zeta(m) with K=2m!/(2pi)^m, and computing zeta(m)
//   via its Euler product.
//
// Variables need to be initialized with a 24 bit precision overhead.
//

int bernreal_(mpfi_t b,const unsigned long m,mpfi_t K,
	      unsigned long *l,unsigned long **p,unsigned long *primelimit) {
  int e=0;
  unsigned long i,N,prec;
  double t;
  mpz_t a;
  mpfi_t z;
  mpfr_t c,c2;

  /* check whether m is reasonable */
  if(m<2) {
    mpfi_set_si(b,1l); /* B_0=1 */ if(m==1) mpfi_div_si(b,b,-2l); /* B_1=-1/2 */
    return 0;
  }
  if(m%2==1) { mpfi_set_si(b,0l); /* B_m=0 */ return 0; }

  /* N=ceil(2^(prec/(m-1)) */
  prec=mpfi_get_prec(K); if(prec<24+MPFR_PREC_MIN) return 64;
  t=prec*M_LN2/(m-1);
  if(t<10) N=(unsigned long)ceil(exp(t)); else return 16;
  N=(unsigned long)ceil(exp(t));

  /* make sure there are enough primes tabulated */
  if(N>*primelimit) { e|=primes(l,p,N);
    if((e&3)==0) *primelimit=N; else { if((e&3)==3) *primelimit=0; return e; }
  }

  /* z = prod_{p<=N} (1-p^-m)^-1 */
  mpfi_init2(z,(mp_prec_t)prec);
  mpz_init(a);
  mpfi_set_si(z,1l);
  for(i=0; i<*l && **p<=N; i++,(*p)++) {
    mpz_ui_pow_ui(a,**p,m);
    mpfi_set_z(b,a);
    mpfi_inv(b,b);
    mpfi_si_sub(b,1l,b);
    mpfi_mul(z,z,b);
  } *p-=i;
  mpfi_inv(z,z);

  /* add the error bound [0, 1/((m-1)*N^(m-1))] to z */
  mpfr_init2(c,(mpfr_prec_t)prec);
  mpfr_init2(c2,(mpfr_prec_t)prec);
  mpz_ui_pow_ui(a,N,(unsigned long)(m-1));
  mpz_mul_ui(a,a,(unsigned long)(m-1));
  mpfr_set_z(c2,a,MPFR_RNDZ);
  mpfr_si_div(c2,1l,c2,MPFR_RNDU);
  mpfr_set_si(c,0l,MPFR_RNDZ);
  mpfi_interv_fr(b,c,c2);
  mpfi_add(z,z,b);

  /* B_m=(-1)^(m/2+1)*K*z with K=2m!/(2pi)^m */
  if(m%4==0) mpfi_neg(z,z);
  mpfi_mul(b,K,z);

  /* check accuracy */
  mpfi_diam_rel(c,b);
  mpfr_log(c,c,MPFR_RNDU);
  t=mpfr_get_d(c,MPFR_RNDU)/M_LN2;
  mpfr_clear(c2); mpfr_clear(c); mpz_clear(a); mpfi_clear(z);
  if(t+prec>24) return 32;
  return 0;
}

//
// Tabulate Bernoulli numbers b[n]=B_{2n+2}/((2n+2)(2n+1)),
// where B_{2n} is the 2n'th Bernoulli number.
//

void clear_bern_tab(mpfi_ptr b,unsigned long N) { while(N-->0) mpfi_clear(b+N);}

int bern_tab(mpfi_ptr b,unsigned long N,unsigned long prec) {
  int e=0;
  unsigned long m=0,l=0,*p=(unsigned long*)0,primelimit=0;
  mpfi_t t,K,T;
  mpz_t a,d;
  N=2*N+1;
  if(prec<MPFR_PREC_MIN) prec=MPFR_PREC_MIN;
  if(prec>MPFR_PREC_MAX-24) prec=MPFR_PREC_MAX-24;
  mpz_init(a); mpz_init(d);
  mpfi_init2(t,(mp_prec_t)prec+24);
  mpfi_init2(K,(mp_prec_t)prec+24);
  mpfi_init2(T,(mp_prec_t)prec+24);
  mpfi_set_si(K,2l);
  mpfi_const_pi(T); mpfi_mul_2si(T,T,1l); mpfi_sqr(T,T); mpfi_inv(T,T);
  if(N>primelimit) { e|=primes(&l,&p,N);
    if((e&3)==0) primelimit=N; else if((e&3)==3) primelimit=0; }
  while(m+2<N) {
    mpfi_mul(K,K,T); m++; mpfi_mul_ui(K,K,m); m++; mpfi_mul_ui(K,K,m);
    mpfi_init2(b+(m/2-1),(mp_prec_t)prec);
    if((e&4)==0) {
      e|=bernfrac_(a,d,m,K,&l,&p,&primelimit);
      if(!e) { mpfi_set_z(t,a); mpfi_div_z(t,t,d);
#ifdef BERNOULLI_MAIN
	gmp_printf("bernfrac(%lu)= %Zd / %Zd\n",m,a,d);
#endif
      }
    }
    if(e==4) { e|=bernreal_(t,m,K,&l,&p,&primelimit); }
    mpfi_div_ui(b+(m/2-1),t,m*(m-1));
  }
  mpfi_clear(T); mpfi_clear(K); mpfi_clear(t); mpz_clear(d); mpz_clear(a);
  free(p);
  return(e&~4);
}

//
// Print info of the bernoulli routines.
//

void bern_mesg(int e) {
  primes_mesg(e);
  if(e&4) printf("bernfrac_: Insufficient precision in K.\n");
  if(e&8) printf("bernfrac_: Cast to integer failed.\n");
  if(e&16) printf("bernreal_: Try to increase precision and use bernfrac_.\n");
  if(e&32) printf("bernreal_: More than 24 bits accuracy lost.\n");
  if(e&64) printf("bernreal_: Precision too low.\n");
  if(e&128) printf("bern_tab & bern_aki: Results disagree.\n");
  if(e&256) printf("bernoulli: Out of memory.\n");
}

#ifdef BERNOULLI_MAIN

#define BernN 2001
mpq_t bern[BernN];

int main(int argc,char **argv) {
  int e=0;
  unsigned long i,m,M,N,d=53;
  mpfi_ptr b;
  if(argc==1)
    printf("Usage: %s -d%lu NNN, where NNN is a positive integer.\n",*argv,d);
  for(i=1; i<argc; i++) {
    if(**(argv+i)=='-' && *(*(argv+i)+1)=='d') d=atol(*(argv+i)+2),**(argv+i)=0;
    else if((**(argv+i)>='A' && **(argv+i)<='Z') ||
            (**(argv+i)>='a' && **(argv+i)<='z'))
      printf("%s should be a number.\n",*(argv+i)),**(argv+i)=0; }

  for(i=1; i<argc; i++) {
    if(**(argv+i)) {
      N=atol(*(argv+i));
      b=malloc(N*sizeof(mpfi_t));
      if(b!=(mpfi_ptr)0) {
	e|=bern_tab(b,N,d);
	if(e==0) {
	  M=N<<1;
	  printf("B_{0..%lu} computed at %lu bits of precision.\n",M,d);
	  if(M>=BernN) M=BernN-1;
	  printf("Compute B_{0..%lu} at infinite precision...\n",M);
	  bern_aki(bern,M);
	  for(m=2; m<=M; b++,m+=2) {
	    mpz_mul_ui(mpq_denref(bern[m]),mpq_denref(bern[m]),m*(m-1));
	    mpq_canonicalize(bern[m]);
	    mpfi_sub_q(b,b,bern[m]);
	    if(!mpfi_has_zero(b)) e|=128;
	  } b-=(m-2)>>1;
	  clear_bern_aki(bern,M); }
	clear_bern_tab(b,N); free(b);
      } else e|=256;
      if(!e) printf("Results agree.\n"); else bern_mesg(e);
    }
  }
  return(e!=0);
}

#endif
#endif

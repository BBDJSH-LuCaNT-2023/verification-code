#ifndef __KBESSEL__
#ifdef  __KBESSEL__
gcc -O2 kbessel.c -lm -lgmp -lmpfr -lmpfi -DKBESSEL_MAIN -o kbessel
./kbessel
exit

//
// Algorithm for computing kbessel.
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


Version: 1.0.4 (1 November 2014)
--------------------------------
** kbessel is a rigorous routine for computing exp(Pi*r/2)*K_{ir}(x) at an
   absolute precision of 2^-d.

** The routine is restricted to purely imaginary order and positive argument.

** It is possible to call init_kbessel() multiple times which is useful for
   computing kbessel at different precisions, e.g.,
        f1=init_kbessel(d1); f2=init_kbessel(d2); ...
        if(accuracy<=d1) kbessel(f1,r,x);
        else if(accuracy<=d2) kbessel(f2,r,x);
        else ...

)))))
#endif
#define __KBESSEL__
#include"archt.c"
#include"lngamma.c"

//
// clear_kbessel(f);
//
void clear_kbessel(mpfi_ptr f) {
  unsigned long *d,n; mpfi_ptr C,*g,tt;
  // Reference variables.
  d=(unsigned long*)(f+1); g=(mpfi_ptr*)(d+4); C=(mpfi_ptr)(g+1); tt=C+15;
  // Clear variables and free memory.
  for(n=12+*(d+2)*2; n--;) mpfi_clear(tt+n);
  for(n=15; n--;) mpfi_clear(C+n);
  mpfi_clear(f); clear_lngamma(*g); free(f);
}

//
// f=init_kbessel(d);  Sets the absolute precision of kbessel(r,x) to 1/2^d.
//
mpfi_ptr init_kbessel(unsigned long d) {
  unsigned long *dd,Nmax,n; mpfi_ptr f,C,*g,tt;
  // Allocate memory.
  Nmax=99999;
  f=malloc(4*sizeof(unsigned long)+(28+Nmax*2)*sizeof(mpfi_t)+sizeof(mpfi_ptr));
  if(f!=(mpfi_ptr)0) {
    // Make sure precision is reasonable.
    if(d<53) d=53; else if(d>MPFR_PREC_MAX-25) d=MPFR_PREC_MAX-25;
    // Reference variables.
    dd=(unsigned long*)(f+1); g=(mpfi_ptr*)(dd+4); C=(mpfi_ptr)(g+1); tt=C+15;
    // Set variables.
    *dd=d+1; *(dd+1)=*dd+24; *(dd+2)=0; *(dd+3)=Nmax;
    *g=init_lngamma(*(dd+1));
    if(*g==(mpfi_ptr)0) { free(f); return((mpfi_ptr)0); }
    // Init variables.
    mpfi_init2(f,*(dd+1));
    for(n=0; n<15; n++) mpfi_init2(C+n,*(dd+1));
    for(n=0; n<12; n++) mpfi_init2(tt+n,*(dd+1));
    // Set variable.
    mpfr_set_inf((mpfr_ptr)tt,1); mpfr_set_inf(1+(mpfr_ptr)tt,-1); // [inf,-inf]
    mpfi_const_pi(tt+2); // pi
    mpfr_set_ui((mpfr_ptr)(tt+1),1l,MPFR_RNDU);
    mpfr_div_2ui((mpfr_ptr)(tt+1),(mpfr_ptr)(tt+1),*dd+1,MPFR_RNDU);// (1/2)/2^d
    mpfr_div_ui(1+(mpfr_ptr)(tt+1),(mpfr_ptr)(tt+1),5l,MPFR_RNDU);
    mpfr_mul_ui(1+(mpfr_ptr)(tt+1),1+(mpfr_ptr)(tt+1),6l,MPFR_RNDU);// (3/5)/2^d
  }
  return f;
}

//
// e=kbessel(f,r,x);  Computes exp(r*pi/2)*K_{ir}(x) and stores the result in f.
//
int kbessel(mpfi_ptr f,mpfi_t r,mpfi_t x) {
  int e=0; unsigned long J,j,j1,N,n;
  double B,X,c0,c1,y0,y1,y2,y3,t,t2,t3,t4;
  unsigned long *d; mpfi_ptr C,*g,tt;
  // Reference variables.
  d=(unsigned long*)(f+1); g=(mpfi_ptr*)(d+4); C=(mpfi_ptr)(g+1); tt=C+15;
  // Make sure that reasonable arguments are given.
  mpfr_set_nan((mpfr_ptr)f); mpfr_set_nan(1+(mpfr_ptr)f);
  c0=mpfr_get_d((mpfr_ptr)r,MPFR_RNDD);
  c1=mpfr_get_d(1+(mpfr_ptr)r,MPFR_RNDU);
  y0=mpfr_get_d((mpfr_ptr)x,MPFR_RNDD);
  y1=mpfr_get_d(1+(mpfr_ptr)x,MPFR_RNDU);
  if(!(c1>=c0)) e|=1;
  if(!(y1>=y0)) e|=2;
  // Make sure that |r|>eps and |r|<inf.
  if(c1<0) { mpfi_neg(C,r); t=c0; c0=-c1; c1=-t; } else mpfi_set(C,r);
  if(!(c0>1e-11)) e|=4;
  if(!(c1<1e111)) e|=8;
  // Make sure that x>eps.
  if(!(y0>1e-111)) e|=16;
  if(e&31) return e;
  // If x>(*d+1)*log(2)+Pi/2*|r| then exp(r*pi/2)*K_{ir}(x)<(1/2)/2^d.
  mpfi_set(C+1,x); y3=(*d+1.)*.693147180560+1.57079632680*c1;
  if(y0>y3) { mpfr_set_ui((mpfr_ptr)f,0l,MPFR_RNDZ);
    mpfr_set(1+(mpfr_ptr)f,(mpfr_ptr)(tt+1),MPFR_RNDU); e|=32; return e; }
  if(y1>y3) { mpfr_set_d(1+(mpfr_ptr)(C+1),y3,MPFR_RNDU); y1=y3; e|=64; }
  // If x>2 use Poisson summation.
  y2=2; if(y0>y2) e|=128;

  // If cache is not up to date than update it.
  if(mpfr_cmp((mpfr_ptr)C,(mpfr_ptr)tt)!=0
     || mpfr_cmp(1+(mpfr_ptr)C,1+(mpfr_ptr)tt)!=0) {
    // Determine B up to a factor of 1/(2*pi), rounded to the safe side.
    B=1; t=log(.5*y3); if(t>B) B=t;
    t3=.500000000001/y2; t2=*d*1.38629436112+log(200*t3)+3.14159265359*c1;
    t=(B+t2)*t3; if(t>0) t=log(t);
    if(t>B) { do { B=log((t+(t-B)+t2)*t3); t=log((B+t2)*t3); } while(t>B); B=t;}
    // Determine N, rounded to the safe side.
    t=log(6.36619772368*t3)+*d*.693147180560;
    if(t>=0) t=c1+.636619772368*t; else t=c1+.636619772367*t;
    t=.159154943092*B*t; if(t>=*(d+3)) { N=*(d+3); e|=512; return e; }
    else if(t>0) N=ceil(t); else { N=0; e|=1024; }
    // Clear and init variables.
    n=12+*(d+2)*2; *(d+2)=N; j1=12+*(d+2)*2;
    while(n-->j1) mpfi_clear(tt+n);
    while(++n<j1) mpfi_init2(tt+n,*(d+1));
    // Set variables.
    mpfi_set(tt,C); // |r|
    mpfr_set_d((mpfr_ptr)(tt+3),2*M_PI/B,MPFR_RNDU);
    mpfr_set(1+(mpfr_ptr)(tt+3),(mpfr_ptr)(tt+3),MPFR_RNDU); // 1/B
    mpfi_set_ui(tt+4,1l); mpfi_div_2ui(tt+4,tt+4,1l); for(n=1; n<=N; n++) {
      mpfi_mul_ui(tt+5,tt+3,n); mpfi_add(tt+6,tt+5,tt); mpfi_sub(tt+5,tt+5,tt);
      mpfi_div_2ui(tt+5,tt+5,1l); mpfi_div_2ui(tt+6,tt+6,1l);
      lngamma(*g,tt+4,tt+5); mpfi_set(tt+7,*g+2); mpfi_set(tt+8,*g+3);
      lngamma(*g,tt+4,tt+6);
      mpfi_mul(tt+10+2*n,tt+7,*g+2); mpfi_mul(tt+11,tt+8,*g+3);
      mpfi_sub(tt+10+2*n,tt+10+2*n,tt+11); // Re gamma(..)*gamma(..)
      mpfi_mul(tt+11+2*n,tt+7,*g+3); mpfi_mul(tt+11,tt+8,*g+2);
      mpfi_add(tt+11+2*n,tt+11+2*n,tt+11); } // Im gamma(..)*gamma(..)
    mpfi_div_2ui(tt+5,tt,1l); lngamma(*g,tt+4,tt+5);
    mpfi_sqr(tt+11,*g+2); mpfi_sqr(tt+6,*g+3); mpfi_add(tt+11,tt+11,tt+6);
    mpfi_div_2ui(tt+11,tt+11,1l); // .5*gamma(1/2+ir/2)*gamma(1/2-ir/2)
    mpfi_mul(tt+4,tt+2,tt+5); mpfi_exp(tt+5,tt+4); mpfi_mul(tt+6,tt+2,tt+5);
    mpfi_mul(tt+9,tt+2,tt); mpfi_sinh(tt+10,tt+9);
    mpfi_div(tt+8,tt+6,tt+10); // pi*exp(pi*r/2)/sinh(pi*r)
    mpfi_div(tt+9,tt+2,tt+3); mpfi_mul_2ui(tt+9,tt+9,1l);
    mpfi_div(tt+7,tt+5,tt+9); mpfi_div_2ui(tt+7,tt+7,1l);// exp(pi*r/2)/(4*pi*B)
    mpfi_mul(tt+10,tt+9,tt); mpfi_exp(tt+4,tt+9);
    mpfi_cos(tt+5,tt+10); mpfi_sin(tt+6,tt+10);
    mpfi_mul(tt+5,tt+5,tt+4); mpfi_mul(tt+6,tt+6,tt+4); // exp((2*pi*B)*(1+ir))
    mpfi_sqr(tt+4,tt+4); // exp(4*pi*B)
    mpfi_set_ui(tt+9,1l); lngamma(*g,tt+9,tt);
    mpfi_set(tt+9,*g+2); mpfi_set(tt+10,*g+3); // gamma(1+ir)
    e|=256;
  } else { N=*(d+2); B=6.28318530717/mpfr_get_d(1+(mpfr_ptr)(tt+3),MPFR_RNDU); }

  if(e&128) {
    // Determine j1=J, rounded to the safe side.
    t2=*d*.693147180560-B+log(85/(.999999999999-exp(-6.28318530717*c0)));
    t3=1.35914091423*y1*exp(-B); t3=t3*t3+1e-300; t4=c0*c0;
    j1=2; do { j1<<=1; X=t3/((j1+1.)*sqrt((j1+1.)*(j1+1.)+t4));
    }while((X>.999999999999 || t2-log((1-X)*(j1+1.))+j1*log(X)>0) && j1<*(d+3));
    if(j1>4) j=j1>>1; else j=0; if(j1>=*(d+3)) { e|=2048; return e; }
    while(j1>j) { J=(j+j1)>>1; X=t3/((J+1.)*sqrt((J+1.)*(J+1.)+t4));
      if(X>.999999999999 || t2-log((1-X)*(J+1.))+J*log(X)>0) j=J+1; else j1=J; }
    // Compute kbessel.
    //   Re sum_{n};
    mpfi_div_2ui(C+1,C+1,1l); // x/2
    mpfi_log(C+3,C+1); // log(x/2)
    mpfi_mul(C+2,tt+3,C+3); // 1/B*log(x/2)
    mpfi_set(C+10,tt+11); // .5*gamma(1/2+ir/2)*gamma(1/2-ir/2)
    for(n=1; n<=N; n++) {
      mpfi_mul_ui(C+4,C+2,n); // n/B*log(x/2)
      mpfi_cos(C+8,C+4); mpfi_sin(C+9,C+4); // (x/2)^{in/B}
      mpfi_mul(C+12,C+8,tt+10+2*n); mpfi_mul(C+13,C+9,tt+11+2*n);
      mpfi_add(C+12,C+12,C+13);
      mpfi_add(C+10,C+10,C+12); } // Re sum_{n}
    mpfi_div(C+12,tt+7,C+1);
    mpfi_mul(f,C+12,C+10); // exp(pi*r/2)/(2*pi*B*x)*Re sum_{n}
    //   Im sum_{j};
    mpfi_sqr(C+2,C+1); // (x/2)^2
    mpfi_mul(C+12,C,C+3);
    mpfi_cos(C+4,C+12); mpfi_sin(C+5,C+12); // (x/2)^ir
    mpfi_set(C+6,tt+5); mpfi_set(C+7,tt+6); // exp((2*pi*B)*(1+ir))
    mpfi_set(C+8,tt+9); mpfi_set(C+9,tt+10); // gamma(1+ir)
    mpfi_set_ui(C+10,1l); // 1
    mpfi_sub_ui(C+3,C+6,1l); // Re exp((2*pi*B)*(1+ir))-1
    mpfi_mul(C+12,C+8,C+3); mpfi_mul(C+13,C+9,C+7); mpfi_sub(C+12,C+12,C+13);
    mpfi_mul(C+13,C+8,C+7); mpfi_mul(C+14,C+9,C+3); mpfi_add(C+13,C+13,C+14);
    mpfi_sqr(C+3,C+12); mpfi_sqr(C+14,C+13); mpfi_add(C+14,C+14,C+3);
    mpfi_mul(C+3,C+10,C+14);
    mpfi_mul(C+14,C+4,C+13); mpfi_mul(C+1,C+5,C+12); mpfi_sub(C+1,C+1,C+14);
    mpfi_div(C+11,C+1,C+3); // Im(zero-th summand)
    for(j=1; j<j1; j++) {
      mpfi_mul(C+4,C+4,C+2); mpfi_mul(C+5,C+5,C+2); // (x/2)^{ir+2j}
      mpfi_mul(C+6,C+6,tt+4); mpfi_mul(C+7,C+7,tt+4); // exp((2*pi*B)*(2j+1+ir))
      mpfi_mul_ui(C+12,C+8,j); mpfi_mul(C+13,C+9,C); mpfi_sub(C+12,C+12,C+13);
      mpfi_mul(C+13,C+8,C); mpfi_mul_ui(C+14,C+9,j); mpfi_add(C+13,C+13,C+14);
      mpfi_swap(C+8,C+12); mpfi_swap(C+9,C+13); // gamma(1+ir+j)
      mpfi_mul_ui(C+10,C+10,j); // j!
      mpfi_sub_ui(C+3,C+6,1l); // Re exp((2*pi*B)*(2j+1+ir))-1
      mpfi_mul(C+12,C+8,C+3); mpfi_mul(C+13,C+9,C+7); mpfi_sub(C+12,C+12,C+13);
      mpfi_mul(C+13,C+8,C+7); mpfi_mul(C+14,C+9,C+3); mpfi_add(C+13,C+13,C+14);
      mpfi_sqr(C+3,C+12); mpfi_sqr(C+14,C+13); mpfi_add(C+14,C+14,C+3);
      mpfi_mul(C+3,C+10,C+14);
      mpfi_mul(C+14,C+4,C+13); mpfi_mul(C+1,C+5,C+12); mpfi_sub(C+1,C+1,C+14);
      mpfi_div(C+14,C+1,C+3); // Im(j-th summand)
      mpfi_add(C+11,C+11,C+14); } // Im sum_{j}
    mpfi_mul(C+12,tt+8,C+11); // pi*exp(pi*r/2)/sinh(pi*r)*Im sum_{j}
    mpfi_add(f,f,C+12); // (..)*Re sum_{n} + (..)*Im sum_{j}
    //   +/-(3/5)/2^d;
    mpfi_increase(f,1+(mpfr_ptr)(tt+1));
  }
  else {
    // Determine j1=J, rounded to the safe side.
    t2=*d*.693147180560+log(22/(.999999999999-exp(-6.28318530717*c0)));
    t3=1.35914091423*y1; t3=t3*t3+1e-300; t4=c0*c0;
    j1=2; do { j1<<=1; X=t3/((j1+1.)*sqrt((j1+1.)*(j1+1.)+t4));
    }while((X>.999999999999 || t2-log((1-X)*(j1+1.))+j1*log(X)>0) && j1<*(d+3));
    if(j1>4) j=j1>>1; else j=0; if(j1>=*(d+3)) { e|=2048; return e; }
    while(j1>j) { J=(j+j1)>>1; X=t3/((J+1.)*sqrt((J+1.)*(J+1.)+t4));
      if(X>.999999999999 || t2-log((1-X)*(J+1.))+J*log(X)>0) j=J+1; else j1=J; }
    // Compute kbessel.
    //   Im sum_{j};
    mpfi_div_2ui(C+1,C+1,1l); // x/2
    mpfi_sqr(C+2,C+1); // (x/2)^2
    mpfi_log(C+3,C+1); // log(x/2)
    mpfi_mul(C+12,C,C+3);
    mpfi_cos(C+4,C+12); mpfi_sin(C+5,C+12); // (x/2)^ir
    mpfi_set(C+8,tt+9); mpfi_set(C+9,tt+10); // gamma(1+ir)
    mpfi_set_ui(C+10,1l); // 1
    mpfi_sqr(C+3,C+8); mpfi_sqr(C+14,C+9); mpfi_add(C+14,C+14,C+3);
    mpfi_mul(C+3,C+10,C+14);
    mpfi_mul(C+14,C+4,C+9); mpfi_mul(C+1,C+5,C+8); mpfi_sub(C+14,C+14,C+1);
    mpfi_div(C+11,C+14,C+3); // Im(zero-th summand)
    for(j=1; j<j1; j++) {
      mpfi_mul(C+4,C+4,C+2); mpfi_mul(C+5,C+5,C+2); // (x/2)^{ir+2j}
      mpfi_mul_ui(C+12,C+8,j); mpfi_mul(C+13,C+9,C); mpfi_sub(C+12,C+12,C+13);
      mpfi_mul(C+13,C+8,C); mpfi_mul_ui(C+14,C+9,j); mpfi_add(C+13,C+13,C+14);
      mpfi_swap(C+8,C+12); mpfi_swap(C+9,C+13); // gamma(1+ir+j)
      mpfi_mul_ui(C+10,C+10,j); // j!
      mpfi_sqr(C+3,C+8); mpfi_sqr(C+14,C+9); mpfi_add(C+14,C+14,C+3);
      mpfi_mul(C+3,C+10,C+14);
      mpfi_mul(C+14,C+4,C+9); mpfi_mul(C+1,C+5,C+8); mpfi_sub(C+14,C+14,C+1);
      mpfi_div(C+1,C+14,C+3); // Im(j-th summand)
      mpfi_add(C+11,C+11,C+1); } // Im sum_{j}
    mpfi_mul(f,tt+8,C+11); // pi*exp(pi*r/2)/sinh(pi*r)*Im sum_{j}
    //   +/-(1/2)/2^d;
    mpfi_increase(f,(mpfr_ptr)(tt+1));
  }

  return e;
}

//
// kbessel_mesg(e);  Print info of e=kbessel(f,r,x).
//
void kbessel_mesg(int e) {
  if(e&1) printf("kbessel: r=nan.\n");
  if(e&2) printf("kbessel: x=nan.\n");
  if(e&4) printf("kbessel: |r|<1e-11 not implemented.\n");
  if(e&8) printf("kbessel: |r|>1e111 not implemented.\n");
  if(e&16) printf("kbessel: x<1e-111 not implemented.\n");
  if(e&32) printf("kbessel: min x>>r.\n");
  if(e&64) printf("kbessel: max x>>r.\n");
  if(e&128) printf("kbessel: Poisson summation used.\n");
  if(e&256) printf("kbessel: Cache updated.\n");
  if(e&512) printf("kbessel: Cache overflow.\n");
  if(e&1024) printf("kbessel: No need to sum in N.\n");
  if(e&2048) printf("kbessel: J overflow.\n");
  if(e&~4095) printf("kbessel: Unsupported bits e=%x.\n",e&~4095);
}

#ifdef KBESSEL_MAIN

int main(int argc,char **argv) {
  int i,e,v=0;
  unsigned long d=53;
  mpfi_ptr f; mpfr_t R; mpfi_t r,x;
  if(argc==1) printf("Usage: %s -d%lu -v r x\n",*argv,d);
  for(i=1; i<argc; i++) {
    if(**(argv+i)=='-' && *(*(argv+i)+1)=='d') d=atol(*(argv+i)+2),**(argv+i)=0;
    else if(**(argv+i)=='-' && *(*(argv+i)+1)=='v') v=1,**(argv+i)=0;
    else if((**(argv+i)>='A' && **(argv+i)<='Z') ||
	    (**(argv+i)>='a' && **(argv+i)<='z'))
      printf("%s should be a number.\n",*(argv+i)),**(argv+i)=0; }

  f=init_kbessel(d); if(f!=(mpfi_ptr)0) { // Initialization.
    mpfi_init2(r,mpfi_get_prec(f));
    mpfi_init2(x,mpfi_get_prec(f));

    mpfr_init2(R,10l);
    for(i=1; i<argc; i++) {
      mpfi_set_ui(r,0l);
      if(**(argv+i)) {
        mpfi_set_str(r,*(argv+i),10); mpfi_revert_if_needed(r);
        printf("\nr= "); archt_out(stdout,r); printf("\n");
        mpfi_set_ui(x,0l);
        for(i++; i<argc; i++) {
          if(**(argv+i)) {
            mpfi_set_str(x,*(argv+i),10); mpfi_revert_if_needed(x);
            printf("x= "); archt_out(stdout,x); printf("\n"); break;
          }
        }

	e=kbessel(f,r,x); // Compute exp(r*pi/2)*K_{ir}(x)

	if(v) {
	  kbessel_mesg(e); // Print info (optional).
	  printf("accuracy +/- ");
	  mpfi_diam_abs(R,f); mpfr_out_str(stdout,10,2,R,MPFR_RNDU);
	  printf(" (");
	  mpfr_set_ui(R,1l,MPFR_RNDN); mpfr_div_2ui(R,R,d,MPFR_RNDN);
	  mpfr_out_str(stdout,10,2,R,MPFR_RNDN); printf(" was requested)\n"); }
        printf("exp(|r|*pi/2)*K_{ir}(x)= "); archt_out(stdout,f); printf("\n");
      }
    }
    mpfr_clear(R);
    mpfi_clear(x); mpfi_clear(r);
    clear_kbessel(f); // Clean up.

  } else printf("init_kbessel(%lu) failed.\n",d);
  return(0);
}

#endif
#endif

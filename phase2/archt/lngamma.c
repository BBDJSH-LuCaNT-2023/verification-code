#ifndef __LNGAMMA__
#ifdef  __LNGAMMA__
gcc -O2 lngamma.c -lm -lgmp -lmpfr -lmpfi -DLNGAMMA_MAIN -o lngamma
./lngamma
exit

//
// Algorithm for computing lngamma and gamma.
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
** lngamma is a rigorous routine for computing the principal branch of lngamma
   at an absolute precision of 2^-d, and gamma at a relative precision of 2^-d.

** The routine is as general as possible. It computes lngamma for any complex
   argument. The argument can be an interval on input, resulting in an
   appropriate interval for the result.

   The user might be aware that covering the branch cut results in wide
   intervals for the output. As an exercise you should think of whether you
   would want to call
        mpfi_set_str(Re_z,"-99.5"); mpfi_set_str(Im_z,"[-1e-99,0]");
        lngamma(f,Re_z,Im_z);
   or whether you would prefer
        mpfi_set_str(Re_z,"-99.5"); mpfi_set_str(Im_z,"[0,1e-99]");
        lngamma(f,Re_z,Im_z); mpfi_neg(f+1,f+1);

** It is possible to call init_lngamma() multiple times which is useful for
   computing lngamma at different precisions, e.g.,
        f1=init_lngamma(d1); f2=init_lngamma(d2); ...
        if(accuracy<=d1) lngamma(f1,Re_z,Im_z);
        else if(accuracy<=d2) lngamma(f2,Re_z,Im_z);
        else ...

)))))
#endif
#define __LNGAMMA__
#include"archt.c"
#include"bernoulli.c"

//
// clear_lngamma(f);
//
void clear_lngamma(mpfi_ptr f) {
  unsigned long *N,i;mpfi_ptr b; mpz_ptr B;
  // Reference variables.
  N=(unsigned long*)(f+7); N+=2; b=(mpfi_ptr)(N+1); B=(mpz_ptr)(b+*N+13);
  // Clear variables and free memory.
  mpz_clear(B+1); mpz_clear(B);
  for(i=*N+13; i--;) mpfi_clear(b+i); for(i=7; i--;) mpfi_clear(f+i); free(f);
}

//
// f=init_lngamma(d);
//
//     Sets the absolute precision of lngamma(z) to 1/2^d,
//     unless z is too close to a pole or too giant in magnitude.
mpfi_ptr init_lngamma(unsigned long d) {
  int e;
  unsigned long n;
  double c,j,lnj,k,r,x,y;
  unsigned long *N; mpfi_ptr f,b; mpz_ptr B;
  // Make sure precision is reasonable.
  if(d<53) d=53; else if(d>MPFR_PREC_MAX/1.177-64) d=MPFR_PREC_MAX/1.177-64;
  // Determine maximal n.
  r=.177*(d+1); c=M_SQRT2;
  j=M_PI*M_E*r/c; lnj=log(j); k=.5*(0.5723649429247001+log(r)+d*M_LN2);
  if(k>lnj) { x=j/M_E-1; do { y=x; x=(k-lnj)/(lnj-(1+.25/y)*log1p(y));
    } while(y>x+.1); if(!(y>=x)) return((mpfi_ptr)0);
    n=ceil(x+1); if(n<1) n=1; } else n=1;
  // Allocate memory.
  f=malloc(3*sizeof(unsigned long)+(20+n)*sizeof(mpfi_t)+2*sizeof(mpz_t));
  if(f!=(mpfi_ptr)0) {
    // Reference variables.
    N=(unsigned long*)(f+7); N+=2; b=(mpfi_ptr)(N+1); b+=12; B=(mpz_ptr)(b+n+1);
    // Set variables.
    *N=n; *(N-1)=d; *(N-2)=1.177*d+40; d+=ceil(log((double)n)/M_LN2);
    // Init variables.
    for(n=0; n<7; n++) mpfi_init2(f+n,*(N-2));
    mpfi_init2(b-12,53l); mpfi_init2(b-11,d); mpfi_init2(b-10,d);
    for(n=0; n<7; n++) mpfi_init2(b-9+n,*(N-1));
    mpfi_init2(b-2,*(N-2)); mpfi_init2(b-1,*(N-1)+2);
    mpz_init(B); mpz_init(B+1);
    // pi/2
    mpfi_const_pi(b-2); mpfi_div_2ui(b-2,b-2,1l);
    // 1/2*(log(2*pi)-1)
    mpfi_mul_2ui(b-1,b-2,2l); mpfi_log(b-1,b-1);
    mpfi_sub_ui(b-1,b-1,1l); mpfi_div_2ui(b-1,b-1,1l);
    // *(b+n)=B_{2n+2}/((2n+2)(2n+1)) for n=0..*N
    e=bern_tab(b,*N+1,*(N-1));
    if(e) { bern_mesg(e); clear_lngamma(f); f=(mpfi_ptr)0; }
  }
  return f;
}

//
// e=lngamma(f,Re_z,Im_z);  Computes lngamma(z) and stores the result in f.
//
int lngamma(mpfi_ptr f,mpfi_t Re_z,mpfi_t Im_z) {
  int e=0; long i; unsigned long m,n,nu; double c,j,lnj,k,r,x,y;
  unsigned long *d,*N; mpz_ptr B; mpfr_ptr R; mpfi_ptr sigma,sum,zz,iz,zp,tmp,b;
  // Reference variables.
  sigma=f+2; zz=sigma+1; d=(unsigned long*)(zz+4); d+=1; N=d+1;
  R=(mpfr_ptr)(N+1); sum=(mpfi_ptr)(R+2); iz=sum+2; zp=iz+2; tmp=zp+2; b=tmp+5;
  B=(mpz_ptr)(b+*N+1);

  // Apply reflection and recursion such that max(x)>=0, max(r)>=.177*(*d+1)
  x=mpfr_get_d(1+(mpfr_ptr)Re_z,MPFR_RNDU);
  c=mpfr_get_d((mpfr_ptr)Re_z,MPFR_RNDD);
  mpfi_mag(R,Im_z); y=mpfr_get_d(R,MPFR_RNDU);
  if(x<0) { j=1-c; c=1-x; x=j; e|=1; }
  j=.177*(*d+1); j*=j; if(x*x+y*y<j) {
    m=ceil(sqrt(j-y*y)-x); x+=m; c+=m; e|=2; } else m=0;
  r=sqrt(x*x+y*y);

  // Determine n.
  if(y<c) { c=1; e|=4; } else c=sqrt(2/(1+x/r));
  j=M_PI*M_E*r/c; lnj=log(j); k=.5*(0.5723649429247001+log(r)+*d*M_LN2);
  if(k>lnj) {
    x=j/M_E-1; do { y=x; x=(k-lnj)/(lnj-(1+.25/y)*log1p(y));
      if(y<x) { e|=32; continue; }
    } while(!(y<x+.1));
    n=ceil(x+1); if(n>*N) { n=*N; e|=16; }
  } else { n=1; e|=8; }

  // sigma=x;
  if((e&3)==0) sigma=Re_z;
  else if(e&1) mpfi_ui_sub(sigma,1+m,Re_z); else mpfi_add_ui(sigma,Re_z,m);
  // zz=|z|^2;
  mpfi_sqr(zz,sigma); mpfi_sqr(zz+1,Im_z); mpfi_add(zz,zz,zz+1);
  // zp=1/z;
  mpfi_inv(zp,zz); mpfi_neg(zp+1,zp);
  mpfi_mul(zp,zp,sigma); mpfi_mul(zp+1,zp+1,Im_z);
  // iz=1/(z*z);
  mpfi_sqr(iz,zp); mpfi_sqr(iz+1,zp+1); mpfi_sub(iz,iz,iz+1);
  mpfi_mul(iz+1,zp,zp+1); mpfi_mul_2ui(iz+1,iz+1,1l);

  // Stirling's asymptotic expansion.
  //   sum=sum_{nu=0}^{n-1} B_{2*nu+2}/((2*nu+2)*(2*nu+1)*z^(2*n+1));
  mpfi_set_ui(sum,0l); mpfi_set_ui(sum+1,0l);
  for(nu=0; nu<n; nu++) {
    //     sum+=b[nu]*zp;
    mpfi_mul(tmp,b+nu,zp); mpfi_add(sum,sum,tmp);
    mpfi_mul(tmp,b+nu,zp+1); mpfi_add(sum+1,sum+1,tmp);
    //     zp*=iz;
    mpfi_mul(tmp,zp,iz); mpfi_mul(tmp+1,zp+1,iz+1); mpfi_mul(tmp+2,zp,iz+1);
    mpfi_sub(zp,tmp,tmp+1); mpfi_mul(tmp,zp+1,iz); mpfi_add(zp+1,tmp,tmp+2);
  }
  //   R=|b[n]|*|zp|*c^(2n+2); sum+/-=R;
  //     R=|zp|;
  mpfi_get_right(R,zp); mpfr_sqr(R+1,R,MPFR_RNDU);
  mpfi_get_right(R,zp+1); mpfr_sqr(R,R,MPFR_RNDU);
  mpfr_add(R,R,R+1,MPFR_RNDU); mpfr_sqrt(R+1,R,MPFR_RNDU);
  //     R*=|b[n]|;
  mpfi_mag(R,b+n); mpfr_mul(R,R,R+1,MPFR_RNDU);
  //     R*=c^(2n+2);
  if((e&4)==0) {
    mpfi_sqrt(zz+1,zz);
    mpfr_div(R+1,(mpfr_ptr)sigma,1+(mpfr_ptr)(zz+1),MPFR_RNDD);
    mpfr_add_ui(R+1,R+1,1l,MPFR_RNDD); mpfr_ui_div(R+1,2l,R+1,MPFR_RNDU);
    mpfr_pow_ui(R+1,R+1,n+1,MPFR_RNDU); mpfr_mul(R,R,R+1,MPFR_RNDU);
  }
  //     sum+/-=R;
  mpfi_increase(sum,R); mpfi_increase(sum+1,R);

  //   f=sum+1/2*(log(2*pi)-1)+(z-1/2)*(log(z)-1);
  mpfi_add(f,sum,b-1);
  mpfi_log(zz+1,zz); mpfi_div_2ui(zz+1,zz+1,1l);
  if(mpfi_is_strictly_pos(sigma)) {
    mpfi_div(zz+3,Im_z,sigma); mpfi_atan(zz+2,zz+3); }
  else {
    mpfi_div(zz+3,sigma,Im_z); mpfi_atan(zz+2,zz+3); mpfi_neg(zz+2,zz+2);
    if(mpfi_is_strictly_pos(Im_z)) mpfi_add(zz+2,zz+2,b-2);
    else if(mpfi_is_strictly_neg(Im_z)) mpfi_sub(zz+2,zz+2,b-2);
    else e|=128; e|=64; }
  mpfi_sub_ui(zz+1,zz+1,1l); mpfi_sub_d(zz,sigma,.5);
  mpfi_mul(zz+3,zz,zz+1); mpfi_add(f,f,zz+3);
  mpfi_mul(zz+3,Im_z,zz+2); mpfi_sub(f,f,zz+3);
  mpfi_mul(zz+3,Im_z,zz+1); mpfi_add(f+1,sum+1,zz+3);
  mpfi_mul(zz+3,zz,zz+2); mpfi_add(f+1,f+1,zz+3);

  // Recursion formula.
  if(e&2) {
    //   i=0; zz=1;
    i=0; mpfi_set_ui(zz,1l); mpfi_set_ui(zz+1,0l);
    while(m-->0) {
      //   sigma-=1;
      mpfi_sub_ui(sigma,sigma,1l);
      //   zz*=sigma+I*Im_z;
      mpfi_mul(zz+2,zz,sigma); mpfi_mul(zz+3,zz+1,Im_z);
      mpfi_sub(zz+2,zz+2,zz+3); mpfi_mul(zz+3,zz,Im_z);
      mpfi_swap(zz,zz+2); mpfi_mul(zz+2,zz+1,sigma); mpfi_add(zz+1,zz+2,zz+3);
      //   |i-arg(zz)|<pi/2;
      if(i&1) { if((i&3)==1) {
	  if(mpfi_is_strictly_neg(zz)) i++;
	  else if(mpfi_is_strictly_pos(zz)) i--;
	  else if(!mpfi_is_strictly_pos(zz+1)) e|=256;
	} else {
	  if(mpfi_is_strictly_pos(zz)) i++;
	  else if(mpfi_is_strictly_neg(zz)) i--;
	  else if(!mpfi_is_strictly_neg(zz+1)) e|=256; }
      } else { if((i&3)==0) {
	  if(mpfi_is_strictly_pos(zz+1)) i++;
	  else if(mpfi_is_strictly_neg(zz+1)) i--;
	  else if(!mpfi_is_strictly_pos(zz)) e|=256;
	} else {
	  if(mpfi_is_strictly_neg(zz+1)) i++;
	  else if(mpfi_is_strictly_pos(zz+1)) i--;
	  else if(!mpfi_is_strictly_neg(zz)) e|=256; } } }
    //   f-=log(zz);
    mpfi_sqr(zz+2,zz); mpfi_sqr(zz+3,zz+1); mpfi_add(zz+2,zz+2,zz+3);
    mpfi_sqrt(zz+3,zz+2); mpfi_log(zz+2,zz+3); mpfi_sub(f,f,zz+2);
    if((e&256)==0) {
      if(i&1) {
	mpfi_div(zz+2,zz,zz+1); mpfi_atan(zz,zz+2); mpfi_add(f+1,f+1,zz); }
      else {
	mpfi_div(zz+2,zz+1,zz); mpfi_atan(zz,zz+2); mpfi_sub(f+1,f+1,zz); }
      mpfi_mul_si(zz,b-2,-i); mpfi_add(f+1,f+1,zz);
    } else { mpfr_set_inf(R,0); mpfi_increase(f+1,R); }
  }
  // Reflection formula.
  if(e&1) {
    // f=-conj(f);
    mpfi_neg(f,f);
    //   if(pi*|y|>(d+2)*log(2)) use large |y| asymptotics;
    if(M_PI*mpfr_get_d((mpfr_ptr)Im_z,MPFR_RNDD)>(*(d-1)+2)*M_LN2) {
      e|=512; //     f=-pi*y+log(2*pi)+I*pi*(x-1/2)-conj(f);
      mpfi_mul_2ui(zz,b-2,1l); mpfi_mul(zz+1,zz,Im_z); mpfi_sub(f,f,zz+1);
      mpfi_mul_2ui(zz+1,b-2,2l); mpfi_log(zz+2,zz+1); mpfi_add(f,f,zz+2);
      mpfi_sub_d(zz+2,Re_z,.5); mpfi_mul(zz+3,zz,zz+2); mpfi_add(f+1,f+1,zz+3);}
    else if(-M_PI*mpfr_get_d(1+(mpfr_ptr)Im_z,MPFR_RNDU)>(*(d-1)+2)*M_LN2) {
      e|=512; //     f=pi*y+log(2*pi)-I*pi*(x-1/2)-conj(f);
      mpfi_mul_2ui(zz,b-2,1l); mpfi_mul(zz+1,zz,Im_z); mpfi_add(f,f,zz+1);
      mpfi_mul_2ui(zz+1,b-2,2l); mpfi_log(zz+2,zz+1); mpfi_add(f,f,zz+2);
      mpfi_sub_d(zz+2,Re_z,.5); mpfi_mul(zz+3,zz,zz+2); mpfi_sub(f+1,f+1,zz+3);}
    else { //   else f-=log(sin(pi*z)/pi);
      mpfr_get_z(B,(mpfr_ptr)Re_z,MPFR_RNDD); // B=floor(x)
      mpfr_get_z(B+1,1+(mpfr_ptr)Re_z,MPFR_RNDD);
      mpfi_mul_2ui(zz-1,b-2,1l); // pi                                      [-1]
      mpfi_sub_z(zz+1,Re_z,B+1); mpfi_mul(zz,zz-1,zz+1); // pi*{x}          [ 0]
      mpfi_mul(zz+1,zz-1,Im_z); // pi*y                                     [ 1]
      mpfi_sin(zz+3,zz); mpfi_sqr(zz+2,zz+3); // sin^2(pi*{x})
      mpfi_sinh(zz+3,zz+1); mpfi_sqr(zz+3,zz+3); // sinh^2(pi*y)
      mpfi_add(zz+2,zz+2,zz+3); // |sin(pi*z)|^2
      mpfi_sqr(zz+3,zz-1); mpfi_div(zz+2,zz+2,zz+3); // (|sin(pi*z)|/pi)^2
      mpfi_log(zz+3,zz+2); mpfi_div_2ui(zz+3,zz+3,1l); // log|sin(pi*z)/pi|
      //     Re(f)-=log|sin(pi*z)/pi|;
      mpfi_sub(f,f,zz+3);
      //     arg(sin(pi*z))=sign(y)*(atan(|tanh(pi*y)|*cot(pi*{x}))-B*pi);
      mpfr_sub_z((mpfr_ptr)zz,(mpfr_ptr)Re_z,B,MPFR_RNDD); //     ...pi*{x} [ 0]
      mpfr_mul((mpfr_ptr)zz,(mpfr_ptr)(zz-1),(mpfr_ptr)zz,MPFR_RNDD);
      mpfi_tanh(zz+1,zz+1); //                              ...|tanh(pi*y)| [ 1]
      mpfr_set_ui(R,0l,MPFR_RNDZ);
      if(mpfr_greaterequal_p((mpfr_ptr)(zz+1),R)) e|=2048;
      else if(mpfr_less_p(1+(mpfr_ptr)(zz+1),R)){ mpfi_neg(zz+1,zz+1); e|=1024;}
      else { mpfr_neg((mpfr_ptr)(zz+1),(mpfr_ptr)(zz+1),MPFR_RNDN);
	mpfi_revert_if_needed(zz+1); e|=3072; } // sign(y)=+/-1
      mpfr_sin((mpfr_ptr)(zz+2),(mpfr_ptr)zz,MPFR_RNDD); //  ...sin(pi*{x}) [ 2]
      if(mpfr_lessequal_p((mpfr_ptr)(zz+2),R))
	mpfr_set_ui((mpfr_ptr)(zz+2),0l,MPFR_RNDZ);
      mpfr_sin(1+(mpfr_ptr)(zz+2),1+(mpfr_ptr)zz,MPFR_RNDU);
      if(mpfr_lessequal_p(1+(mpfr_ptr)(zz+2),R))
	mpfr_set_ui(1+(mpfr_ptr)(zz+2),0l,MPFR_RNDZ);
      mpfr_cos((mpfr_ptr)(zz+3),(mpfr_ptr)zz,MPFR_RNDU); //...cos(pi*{x})   [ 3]
      mpfr_cos(1+(mpfr_ptr)(zz+3),1+(mpfr_ptr)zz,MPFR_RNDD);//..cot(pi*{x}) [ 0]
      mpfr_div((mpfr_ptr)zz,(mpfr_ptr)(zz+3),(mpfr_ptr)(zz+2),MPFR_RNDU);
      mpfr_div(1+(mpfr_ptr)zz,1+(mpfr_ptr)(zz+3),1+(mpfr_ptr)(zz+2),MPFR_RNDD);
      if(mpfr_less_p((mpfr_ptr)(zz+3),R)) //    ...|tanh(pi*y)|*cot(pi*{x}) [ 2]
	mpfr_mul((mpfr_ptr)(zz+2),(mpfr_ptr)(zz+1),(mpfr_ptr)zz,MPFR_RNDU);
      else mpfr_mul((mpfr_ptr)(zz+2),1+(mpfr_ptr)(zz+1),(mpfr_ptr)zz,MPFR_RNDU);
      if(mpfr_less_p(1+(mpfr_ptr)(zz+3),R)) mpfr_mul(
	1+(mpfr_ptr)(zz+2),1+(mpfr_ptr)(zz+1),1+(mpfr_ptr)zz,MPFR_RNDD);
      else mpfr_mul(
        1+(mpfr_ptr)(zz+2),(mpfr_ptr)(zz+1),1+(mpfr_ptr)zz,MPFR_RNDD);
      mpfr_atan((mpfr_ptr)(zz+1),(mpfr_ptr)(zz+2),MPFR_RNDU); //  ...atan() [ 1]
      mpfr_atan(1+(mpfr_ptr)(zz+1),1+(mpfr_ptr)(zz+2),MPFR_RNDD);
      mpfr_mul_z((mpfr_ptr)(zz+3),1+(mpfr_ptr)(zz-1),B,MPFR_RNDD); //...B*pi [3]
      mpfr_mul_z(1+(mpfr_ptr)(zz+3),(mpfr_ptr)(zz-1),B+1,MPFR_RNDU);
      mpfr_sub(1+(mpfr_ptr)zz,(mpfr_ptr)(zz+1),(mpfr_ptr)(zz+3),MPFR_RNDU);
      mpfr_sub((mpfr_ptr)zz,1+(mpfr_ptr)(zz+1),1+(mpfr_ptr)(zz+3),MPFR_RNDU);
      //     Im(f)-=arg(sin(pi*z));
      if((e&3072)==2048) mpfi_sub(f+1,f+1,zz);
      else if((e&3072)==1024) mpfi_add(f+1,f+1,zz);
      else { mpfr_sub((mpfr_ptr)(f+1),(mpfr_ptr)(f+1),1+(mpfr_ptr)zz,MPFR_RNDD);
	mpfr_add(1+(mpfr_ptr)(f+1),1+(mpfr_ptr)(f+1),1+(mpfr_ptr)zz,MPFR_RNDU);}
    }
  }
  // gamma(z).
  mpfi_exp(f+4,f); mpfi_cos(f+5,f+1); mpfi_sin(f+6,f+1);
  mpfi_mul(f+2,f+4,f+5); mpfi_mul(f+3,f+4,f+6);
  return e;
}

//
// lngamma_mesg(e);  Print info of e=lngamma(f,Re_z,Im_z).
//
void lngamma_mesg(int e) {
  if(e&1) printf("lngamma: Reflection formula used.\n");
  if(e&2) printf("lngamma: Recursion formula used.\n");
  if(e&4) printf("lngamma: Lindelof bound used.\n");
  if(e&8) printf("lngamma: n=1 used.\n");
  if(e&16) printf("lngamma: Not enough Bernoulli numbers tabulated.\n");
  if(e&32) printf("lngamma: Approached n from the left.\n");
  if(e&64) printf("lngamma: %s includes zero.\n",(e&128)?"z":"Re(z)");
  if(e&256) printf("lngamma: Im(f)=[-inf,inf] after recursion.\n");
  if(e&512) printf("lngamma: Used large |y| asymptotics in reflection.\n");
  if(e&1024) printf("lngamma: Im(z)<0.\n");
  if(e&2048) printf("lngamma: Im(z)>=0.\n");
  if(e&~4095) printf("lngamma: Unsupported bits e=%x.\n",e&~4095);
}

#ifdef LNGAMMA_MAIN

int main(int argc,char **argv) {
  int i,e,v=0;
  unsigned long d=53;
  mpfi_ptr f; mpfr_t R; mpfi_t Re_z,Im_z;
  if(argc==1) printf("Usage: %s -d%lu -v Re_z Im_z\n",*argv,d);
  for(i=1; i<argc; i++) {
    if(**(argv+i)=='-' && *(*(argv+i)+1)=='d') d=atol(*(argv+i)+2),**(argv+i)=0;
    else if(**(argv+i)=='-' && *(*(argv+i)+1)=='v') v=1,**(argv+i)=0;
    else if((**(argv+i)>='A' && **(argv+i)<='Z') ||
	    (**(argv+i)>='a' && **(argv+i)<='z'))
      printf("%s should be a number.\n",*(argv+i)),**(argv+i)=0; }

  f=init_lngamma(d); if(f!=(mpfi_ptr)0) { // Initialization.
    mpfi_init2(Re_z,mpfi_get_prec(f));
    mpfi_init2(Im_z,mpfi_get_prec(f));

    mpfr_init2(R,10l);
    for(i=1; i<argc; i++) {
      mpfi_set_ui(Re_z,0l);
      if(**(argv+i)) {
	mpfi_set_str(Re_z,*(argv+i),10); mpfi_revert_if_needed(Re_z);
	printf("\nz= "); archt_out(stdout,Re_z);
	mpfi_set_ui(Im_z,0l);
	for(i++; i<argc; i++) {
	  if(**(argv+i)) {
	    mpfi_set_str(Im_z,*(argv+i),10); mpfi_revert_if_needed(Im_z);
	    printf(" + I*"); archt_out(stdout,Im_z); break;
	  }
	}
	printf("\n");

	e=lngamma(f,Re_z,Im_z); // Compute lngamma(z) and gamma(z).

	if(v) { lngamma_mesg(e); // Print info (optional).
	  printf("accuracy +/- ");
	  mpfi_diam_abs(R,f); mpfr_out_str(stdout,10,2,R,MPFR_RNDU);
	  printf(" +/- I*");
	  mpfi_diam_abs(R,f+1); mpfr_out_str(stdout,10,2,R,MPFR_RNDU);
	  printf(" (");
	  mpfr_set_ui(R,1l,MPFR_RNDN); mpfr_div_2ui(R,R,d,MPFR_RNDN);
	  mpfr_out_str(stdout,10,2,R,MPFR_RNDN); printf(" was requested)\n"); }
	printf("lngamma(z)= "); archt_out(stdout,f);
	printf(" + I*"); archt_out(stdout,f+1); printf("\n");
	printf("gamma(z)= "); archt_out(stdout,f+2);
	printf(" + I*"); archt_out(stdout,f+3); printf("\n");
      }
    }
    mpfr_clear(R);
    mpfi_clear(Im_z); mpfi_clear(Re_z);
    clear_lngamma(f); // Clean up.
  } else printf("init_lngamma(%lu) failed.\n",d);
  return(0);
}

#endif
#endif

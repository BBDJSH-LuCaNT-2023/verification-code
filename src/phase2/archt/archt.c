#ifndef archt_version
#ifdef  archt_version
gcc -O2 archt.c -lm -lgmp -lmpfr -lmpfi -DARCHT_MAIN -o archt
./archt
exit

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


)))))
#endif
#define archt_version "1.0.4"
#define archt_version_date "1 November 2014"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<mpfr.h>
#include<mpfi.h>
#include<mpfi_io.h>

// Avoid MPFI bugs.
//   Never use any mpfi_cmp routine,
//   Never use any mpfi_is_inside or mpfi_is_strictly_inside routine.
//   The other mpfi_is_ routines are Ok.
// Fix the mpfi_ui_sub routine.
#define mpfi_ui_sub(rop,op1,op2) (mpfi_sub_ui(rop,op2,op1),mpfi_neg(rop,rop),MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT)

#ifndef archt_nocheck
// Checking for memory leaks.
static long *archt_check_p=(long*)0;
void archt_check(void) {
  int i,k=0; if(archt_check_p!=(long*)0) {
    for(i=1; i<=6; i++) if(*(archt_check_p+6-i)-*(archt_check_p+6+i)!=0) k=1;
    if(k) { printf("memory leak= { "); for(i=1; i<=6; i++)
	printf("%ld ",*(archt_check_p+6-i)-*(archt_check_p+6+i)); printf("}\n");
    } } }
void archt_check_add(void *p,int k) {
  int i; if(archt_check_p==(long*)0) { archt_check_p=malloc(13*sizeof(long));
    if(archt_check_p!=(long*)0) { atexit(archt_check);
      for(i=0; i<13; i++) *(archt_check_p+i)=0; } }
  if(archt_check_p!=(long*)0) { if(p==(void*)0) k=0; *(archt_check_p+6-k)+=1; }
}
void *archt_calloc(size_t m,size_t n) {
  void *p; p=calloc(m,n); archt_check_add(p,1); return p; }
void *archt_malloc(size_t m) {
  void *p; p=malloc(m); archt_check_add(p,1); return p; }
void archt_free(void *p) {
  archt_check_add(p,-1); free(p); }
void *archt_realloc(void *p,size_t n) {
  archt_check_add(p,-1); p=realloc(p,n); archt_check_add(p,1); return p; }
void archt_mpz_init(mpz_ptr z) {
  mpz_init(z); archt_check_add((void*)z,2); }
void archt_mpz_init2(mpz_ptr z,unsigned long n) {
  mpz_init2(z,n); archt_check_add((void*)z,2); }
void archt_mpz_clear(mpz_ptr z) {
  archt_check_add((void*)z,-2); mpz_clear(z); }
void archt_mpq_init(mpq_ptr q) {
  mpq_init(q); archt_check_add((void*)q,3); }
void archt_mpq_clear(mpq_ptr q) {
  archt_check_add((void*)q,-3); mpq_clear(q); }
void archt_mpf_init(mpf_ptr f) {
  mpf_init(f); archt_check_add((void*)f,4); }
void archt_mpf_init2(mpf_ptr f,unsigned long n) {
  mpf_init2(f,n); archt_check_add((void*)f,4); }
void archt_mpf_clear(mpf_ptr f) {
  archt_check_add((void*)f,-4); mpf_clear(f); }
void archt_mpfr_init(mpfr_ptr f) {
  mpfr_init(f); archt_check_add((void*)f,5); }
void archt_mpfr_init2(mpfr_ptr f,unsigned long d) {
  mpfr_init2(f,d); archt_check_add((void*)f,5); }
void archt_mpfr_clear(mpfr_ptr f) {
  archt_check_add((void*)f,-5); mpfr_clear(f); }
void archt_mpfi_init(mpfi_ptr f) {
  mpfi_init(f); archt_check_add((void*)f,6); }
void archt_mpfi_init2(mpfi_ptr f,unsigned long d) {
  mpfi_init2(f,d); archt_check_add((void*)f,6); }
void archt_mpfi_clear(mpfi_ptr f) {
  archt_check_add((void*)f,-6); mpfi_clear(f); }
// Undefine some macros before redefining them.
#undef mpz_init
#undef mpz_init2
#undef mpz_clear
#undef mpq_init
#undef mpq_clear
#undef mpf_init
#undef mpf_init2
#undef mpf_clear
// The following are checked.
#define calloc archt_calloc
#define malloc archt_malloc
#define free archt_free
#define realloc archt_realloc
#define mpz_init archt_mpz_init
#define mpz_init2 archt_mpz_init2
#define mpz_clear archt_mpz_clear
#define mpq_init archt_mpq_init
#define mpq_clear archt_mpq_clear
#define mpf_init archt_mpf_init
#define mpf_init2 archt_mpf_init2
#define mpf_clear archt_mpf_clear
#define mpfr_init archt_mpfr_init
#define mpfr_init2 archt_mpfr_init2
#define mpfr_clear archt_mpfr_clear
#define mpfi_init archt_mpfi_init
#define mpfi_init2 archt_mpfi_init2
#define mpfi_clear archt_mpfi_clear
#endif

// An alternative to mpfi_out_str.
void archt_out(FILE *fptr,mpfi_ptr x) {
  long a2,b2,c2;
  double a,b,c;
  mpfr_t R;
  if(mpfr_number_p((mpfr_ptr)x) && mpfr_number_p(1+(mpfr_ptr)x) &&
     !mpfi_has_zero(x)) {
    mpfr_init2(R,53l);
    mpfi_diam_abs(R,x); b=mpfr_get_d_2exp(&b2,R,MPFR_RNDN); if(b<0) b=-b;
    mpfi_mid(R,x); a=mpfr_get_d_2exp(&a2,R,MPFR_RNDN); if(a<0) a=-a;
    if(b<.25) { b=a; b2=-mpfi_get_prec(x); } else b2-=a2;
    c=(log(b)+b2*M_LN2)/M_LN10; c2=floor(-c); if(c2<=0) c2=0;
    mpfr_clear(R);
  } else c2=0;
  mpfi_out_str(fptr,10,(size_t)c2+2,x);
}

#ifdef ARCHT_MAIN

int main() {
  printf("-----------------------------------------------------------------\n");
  printf("This is archt version %s (%s).\nA", archt_version,archt_version_date);
  printf("uthors: Andrew R. Booker, Andreas Strombergsson, and Holger Then.\n");
  printf("\nCopyright (C) 2011--2014 University of Bristol, UK,\n  and     (C");
  printf(") 2011--2014 University of Uppsala, Sweden.\n\nA");
  printf("rcHT is free software: you can redistribute it and/or modify it\nun");
  printf("der the terms of the GNU General Public License as published by\nth");
  printf("e Free Software Foundation, either version 3 of the License, or\n(a");
  printf("t your option) any later version.\n\narcht is distributed in the ho");
  printf("pe that it will be useful,\nbut WITHOUT ANY WARRANTY; without even ");
  printf("the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULA");
  printf("R PURPOSE.  See the\nGNU General Public License for more details.\n");
  printf("\nYou should have received a copy of the GNU General Public License");
  printf("\nalong with archt.  If not, see <http://www.gnu.org/licenses/>.\n-");
  printf("-----------------------------------------------------------------\n");
  printf("The current file was compiled and linked with\n");
  printf("  GMP version %s\n",gmp_version);
  printf("  MPFR version %s\n",mpfr_get_version());
  printf("  MPFI version %s\n",mpfi_get_version());
  printf("  archt version %s\nand uses ",archt_version);
  printf("%d bits per limb.\n-----------------------------\n",mp_bits_per_limb);
  return(0);
}

#endif
#endif

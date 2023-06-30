#ifndef __PRIMES__
#ifdef  __PRIMES__
gcc -O2 primes.c -lm -DPRIMES_MAIN -o primes
./primes
exit

//
// Euler sieve for computing primes.
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
Usage (example):
  #include"primes.c"
  ...
  int e; unsigned long l=0,*p=(unsigned long*)0,primelimit=3;
  e=primes(&l,&p,primelimit);
  ...
  free(p);

The returned value indicates:

  0: The routine completed successfully.
     l counts the number of primes up to primelimit,
     and p points to the list of these l primes.

  1: Primelimit is far too large causing overflow.
     l and p, including their content, are unaltered.

  2: Primelimit is this large that not enough memory could be allocated.
     l and p, including their content, are unaltered.

  3: System error or compiler bug.
     l and p are reset to l=0, p=(unsigned long*)0.


)))))
#endif
#define __PRIMES__

#include<limits.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int primes(unsigned long *l,unsigned long **p,const unsigned long primelimit) {
  unsigned long i,j,jp,jn,k,*m,*n;

  /* check whether primelimit is reasonable */
  if(primelimit<2) { *l=0; free(*p); *p=(unsigned long*)0; return 0; }
  if(primelimit-1>ULONG_MAX/sizeof(unsigned long)) return 1; /* overflow */

  /* allocate memory */
  n=(unsigned long*)malloc((primelimit-1)*sizeof(unsigned long));
  if(n==(unsigned long*)0) return 2; /* out of memory */
  m=(unsigned long*)realloc(*p,(primelimit-1)*sizeof(unsigned long));
  if(m==(unsigned long*)0) { free(n); return 2; } /* out of memory */

  /* create list of integers */
  *m++=0; for(i=2; i<primelimit; i++) *m++=*n++=1; *n++=0;
  m-=primelimit+1; n-=primelimit+1;
  /* if(2<=i && i<=primelimit) {
     actual_number=i; previous_number=i-*(m+i); next_number=i+*(n+i); } */

  /* Euler sieve */
  *l=(unsigned long)/*floor*/sqrt((double)primelimit);
  for(i=2; i<=*l; i+=*(n+i)) {
    for(k=i; k*i<=primelimit && *(n+k)>0; k+=*(n+k));
    for(; k*i>primelimit; k-=*(m+k));
    for(; k>=i; k-=*(m+k)) {
      j=k*i; /* remove the non-prime j */
      jp=j-*(m+j); jn=j+*(n+j);
      if(jn>j) *(n+jp)=*(m+jn)=jn-jp; else *(n+jp)=0; if(*(m+k)==0) break;
    }
    if(*(n+i)==0) break;
  }
  /* if(2<=i && i<=primelimit && isprime(i)) {
     actual_prime=i; previous_prime=i-*(m+i); next_prime=i+*(n+i); } */

  /* collapse the list down to the primes */
  for(i=2,m+=2,*l=0; ; i+=*(n+i)) {
    *m++=i; (*l)++; if(*(n+i)==0) break; } m-=*l;
  /* m points to the consecutive list of primes */

  /* free memory */
  free(n+2); *p=(unsigned long*)realloc(m,*l*sizeof(unsigned long));
  if(*p==(unsigned long*)0) { *l=0; free(m); return 3; } /* data lost */
  return 0;
}

void primes_mesg(int e) {
  if((e&3)==1) printf("primes: Overflow.\n");
  if((e&3)==2) printf("primes: Out of memory.\n");
  if((e&3)==3) printf("primes: Data lost.\n");
}

#ifdef PRIMES_MAIN

int main(int argc,char **argv) {
  unsigned long i,j;
  int e; unsigned long l=0,*p=(unsigned long*)0,primelimit;
  if(argc==1) printf("Usage: %s NNN, where NNN is a positive integer.\n",*argv);

  for(i=1; i<argc; i++) {
    primelimit=atol(*(argv+i));
    printf("primelimit=%lu ",primelimit);
    e=primes(&l,&p,primelimit);
    primes_mesg(e);
    printf("#{p}=%lu p=[",l);
    for(j=0;j<l;j++) printf("%lu%s",*(p+j),j<l-1?", ":""); printf("]\n");
  }
  free(p); return(0);
}

#endif
#endif

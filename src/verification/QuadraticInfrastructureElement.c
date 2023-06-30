// memore management

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "QuadraticInfrastructureElement.h"

//Function Defintions


// basic arithmetic, XGCD methods

inline void DivRem(long * q, long * r, long dividend, long divisor){
  *q = (long) floor((double) dividend / (double) divisor);
  *r = dividend - (*q) * divisor;  
}


inline void XGCD(long *d, long *s, long *t, long a, long b)
  {
    long  u, v, u0, v0, u1, v1, u2, v2, q, r;
    long aneg = 0, bneg = 0;

    if (a < 0) {
      a = -a;
      aneg = 1;
    }

    if (b < 0) {
      b = -b;
      bneg = 1;
    }

    u1=1; v1=0;
    u2=0; v2=1;
    u = a; v = b;

    while (v != 0) {
      DivRem(&q,&r,u,v);
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 = u1 - q*u2;
      v2 = v1 - q*v2;
      u1 = u0;
      v1 = v0;
    }

    if (aneg)
      u1 = -u1;

    if (bneg)
      v1 = -v1;

    *d = u;
    *s = u1;
    *t = v1;
  }



inline void XGCD_LEFT(long *d, long *s, long a, long b)
  {
    long  u, v, u0, u1, u2, q, r;
    long aneg = 0;

    u = a; v = b;

    if (u < 0) {
      u = -u;
      aneg = 1;
    }

    if (v < 0)
      v = -v;

    u1=1;
    u2=0;

    while (v) {
      DivRem(&q,&r,u,v);
      u = v;
      v = r;
      u0 = u2;
      u2 =  u1 - q*u2;
      u1 = u0;
    }

    if (aneg)
      u1 = -u1;

    *d = u;
    *s = u1;
  }



inline void XGCD_PARTIAL(long *R2, long *R1, long *C2, long *C1, long bound)
  {
    long  q, r;

    *C2 = 0;
    *C1 = -1;

#ifdef DEBUG    
    printf("\nPARTIAL:  R2=%ld, R1=%ld, bound=%ld\n",*R2,*R1,bound);
    printf("OC = %ld\n",*C2);
    printf("CC = %ld\n",*C1);
    fflush(stdout);
#endif

    while (*R1 > bound) {
      DivRem(&q,&r,*R2,*R1);
      *R2 = *R1;
      *R1 = r;

      r = *C2 - q * (*C1);
      *C2 = *C1;
      *C1 = r;

#ifdef DEBUG
      printf("*** q=%ld, C=%ld, R=%ld\n",q,*C1,*R1);
#endif
    }

#ifdef DEBUG
    printf("\n");
#endif
  }




// print
void print(QuadInfElement * qie) {
    printf("(%ld, %ld, %ld, ",qie->a,qie->b,qie->c);
    arb_printn(qie->distance,qie->precision/3,0);
    printf(")\n");
    fflush(stdout);
}



// test_ideal
// Checks whether the given ideal is in the order of the specified discriminant
void test_ideal(QuadInfElement * qie, char *msg) {
  long temp = qie->b * qie->b - (qie->a * qie->c << 2);
  if (temp != qie->Delta) {
    printf("Ideal error:  %s\n",msg);
    print(qie);
    fflush(stdout);
    exit(1);
  }
}



void normalize(QuadInfElement * qie) {
  static long rootD, a, b, c, a2, q, r, temp, nb;
  
  a = qie->a;
  b = qie->b;
  c = qie->c;
  rootD = qie->floor_sqrt_Delta;

  a2 = a << 1;
  if (a <= rootD) {
    temp = rootD - a2;
    if (!(temp < b && b <= rootD)) {
	temp = rootD - b;
	DivRem(&q,&r,temp,a2);

	nb = rootD - r;

	c += q * ((b + nb) >> 1);
	b = nb;
    }
  }
  else {
    if (!(b > -a && b <= a)) {
	DivRem(&q, &r, b, a2);

	if (r > a) {
	  r -= a2;
	  ++q;
	}

      c -= q* ((b + r) >> 1);
      b = r;
    }      
  }

  qie->b = b;
  qie->c = c;
}



void reduce (QuadInfElement * qie) {
  static long q, r, temp, a2, nb, na, s;
  static long a, b, c, rootD;

  a = qie->a;
  b = qie->b;
  c = qie->c;
  rootD = qie->floor_sqrt_Delta;

  qie->num_q = 0;

  // na = (D - b^2) / 4a
  c = -c;

  a2 = abs(a) << 1;
  temp = rootD - a2;
  if (temp < 0)
    ++temp;

#ifdef DEBUG
  printf("\nREDUCE:  num_q=%ld\n",qie->num_q);
  printf("a=%ld, b=%ld, c=%ld\n",a,b,c);
  printf("temp=%ld, rootD=%ld\n",temp,rootD);
#endif
  
  while ((abs(temp) >= b) || (b > rootD) || (a < 0)) {
    s = (a > 0) ? 0 : 1;

    // (rootD+b) = (2a)q + r
    a2 = a << 1;
    temp = rootD+b;
    if (s)  ++temp;
    DivRem(&q,&r,temp,a2);

    qie->qlist[qie->num_q] = q;
    ++qie->num_q;

    // nb = rootD + s - r;
    nb = rootD - r;
    if (s)  ++nb;

    // na = c - q * (nb - b)/2
    na = c - q*((nb - b) >> 1);
   
    b = nb;
    c = a;
    a = na;
  
    a2 = abs(a) << 1;
    temp = rootD - a2;
    if (temp < 0)
       ++temp;

#ifdef DEBUG
    printf("\nq=%ld\n",q);
    printf("a=%ld, b=%ld, c=%ld\n",a,b,c);
    printf("\ntemp=%ld, rootD=%ld\n",temp,rootD);
#endif
  }

  c = -c;
  //  assign_abc (qie, a, b, c);
  qie->a = a;
  qie->b = b;
  qie->c = c;
  
#ifdef DEBUG
  printf("\nDONE REDUCE:  num_q=%ld\n",qie->num_q);
  print(qie);
  test_ideal(qie,"reduce (after)");
  printf("\n");
#endif
}



void init(QuadInfElement * qie, long Delta, long precision) {
  qie->a = 0;
  qie->b = 0;
  qie->c = 0;

  qie->Delta = Delta;
  qie->floor_sqrt_Delta = (long)floor(sqrt((double)Delta));
  qie->NC_BOUND = (long) floor(sqrt(sqrt((double) Delta)));

  qie->num_q = 0;
  
  arb_init(qie->sqrt_Delta);
  arb_sqrt_ui(qie->sqrt_Delta,Delta,precision);

  arb_init(qie->distance);
  arb_zero(qie->distance);

  qie->precision = precision;
}


//
// call when done with a QuadInfElement or before reinitializing
//
void clear(QuadInfElement * qie) {
  arb_clear(qie->sqrt_Delta);
  arb_clear(qie->distance);
}


// assign_one
// Assigns the identity ideal to *qie
void assign_one(QuadInfElement * qie) {
  static long a, b, c;

  a = 1;
  
  if(qie->Delta & 3) {
    b = 1;
    c = (1 - qie->Delta) / (a << 2);
  }
  else {
    b = 0;
    c = (-qie->Delta) / (a << 2);
  }
  
  qie->a = a;
  qie->b = b;
  qie->c = c;
  normalize(qie);

  qie->num_q = 0;
  
  arb_zero(qie->distance);
}


 
void assign_abc(QuadInfElement * qie, long a, long b, long c) {
  qie->a = a;
  qie->b = b;
  qie->c = c;
  normalize(qie);
}



// rho:  Computation of the QuadraticIdealBase via a single step in the continued
// fraction expansion

void rho(QuadInfElement * qie){
  long a, b, c, Delta, rootD, precision;
  long q, r, a2, nb, na, temp, s;
  arb_t relative_generator;

  a = qie->a;
  b = qie->b;
  c = qie->c;
  Delta = qie->Delta;
  rootD = qie->floor_sqrt_Delta;
  precision = qie->precision;

  c = -c;

  s = (a > 0) ? 0 : 1;

  // (rootD-b) = (2a)*q + r
  a2 = a << 1;
  temp = rootD + b;
  if (s)
    ++temp;
  DivRem (&q, &r, temp, a2);

  // nb = rootD + s - r;
  nb = rootD - r;
  if (s)
    ++nb;

  // na = c -q*((nb - b) >> 1);
  na = c - q * ((nb - b) >> 1);

  // d += ln( (nb + rd) / 2a )
  arb_init(relative_generator);
  arb_add_si(relative_generator, qie->sqrt_Delta, nb, precision);
  arb_div_si(relative_generator, relative_generator, a2, precision);
  arb_log(relative_generator, relative_generator, precision);
  
  // Updating the distance
  arb_add(qie->distance, qie->distance, relative_generator, precision);
  
  qie->b = nb;
  qie->c = -a;
  qie->a = na;

  arb_clear(relative_generator);
  
#ifdef DEBUG
  test_ideal(qie,"rho");
#endif  
}



// inverse_rho:  Computation of the QuadraticIdealBase via a single backwards step in
// the continued fraction expansion

void inverse_rho(QuadInfElement * qie){
  long a, b, c, Delta, rootD, precision;
  long q, r, a2, nb, oa, temp;
  arb_t relative_generator;

  a = qie->a;
  b = qie->b;
  c = qie->c;
  Delta = qie->Delta;
  rootD = qie->floor_sqrt_Delta;
  precision = qie->precision;

  oa = a;
  a = -c;
  a2 = a << 1;

  // q = floor((rootD + b) / 2a)
  temp = rootD + b;
  DivRem (&q, &r, temp, a2);
  if (temp < 0 && r != 0)
     --q;

  // d -= ln( (b + rd) / 2 a )
  arb_init(relative_generator);
  arb_add_si(relative_generator, qie->sqrt_Delta, b, precision);
  arb_div_si(relative_generator, relative_generator, a2, precision);
  arb_log(relative_generator, relative_generator, precision);

  // Updating the distance
  arb_sub(qie->distance, qie->distance, relative_generator, precision);

  
  // b = 2aq - b
  nb = rootD - r;

  // c =  q*((nb - b)/2) - c
  c = q * ((nb - b) >> 1) - oa;

  qie->a = a;
  qie->b = nb;
  qie->c = c;
  if (a < 0)
    {
      qie->a = -a;
      qie->c = -c;
    }

  arb_clear(relative_generator);
  
#ifdef DEBUG
  test_ideal(qie,"inverse_rho");
#endif    
}


//
// adjust
// perform baby steps until distance is <= dist, as close as possible
//
void adjust(QuadInfElement * qie, double bound){
  arb_t arb_bound;
  arb_init(arb_bound);
  arb_set_d(arb_bound, bound);

#ifdef DEBUG
  printf("\n--> ADJUST:  bound=%lf\n",bound);
  printf("--> ");
  print(qie);
#endif
  while (arb_le(qie->distance, arb_bound)) {
    rho(qie);
#ifdef DEBUG
    printf("--> ");
    print(qie);
#endif
  }
  while (arb_gt(qie->distance, arb_bound)) {
    inverse_rho(qie);
#ifdef DEBUG
    printf("--> ");
    print(qie);
#endif
  }

#ifdef DEBUG
  test_ideal(qie,"adjust");
#endif  
}



//
// adjust_to_one
// perform baby steps until the ideal is (1)
//
void adjust_to_one (QuadInfElement * qie) {
#ifdef DEBUG
  printf("\nADJUST_TO_ONE\n");
  print(qie);
#endif  
  while (qie->a != 1 || arb_contains_zero(qie->distance)) {
    rho(qie);
#ifdef DEBUG
    print(qie);
#endif
  }

#ifdef DEBUG
  test_ideal(qie,"adjust_to_one");
  printf("\n");
#endif  
}

   
// print_cycle
// Checks whether the given ideal is in the order of the specified discriminant
void print_cycle(long Delta, long precision) {
  QuadInfElement qie;
  init(&qie,Delta,precision);
 
  // print forward cycle
  printf("\nPRINT_CYCLE:  Delta = %ld\n",Delta);
  printf("FORWARD:\n");
  assign_one(&qie);
  print(&qie);
  do {
    rho(&qie);
    print(&qie);
  } while (qie.a != 1);

  // print backward cycle, hoping to end at distance zero
  printf("\nBACKWARD:\n");
  do {
    inverse_rho(&qie);
    print(&qie);
  } while (qie.a != 1);

  clear(&qie);
}



//
// construct relative generator
// computes the coefficients a, b, d of the relative generator
// gamma = (a + b*sqrt(Delta))/d such that A^2 = (gamma) B
// Pass OB=1 and BB=0 
//
void
construct_relative_generator(long *rel_gen_a, long *rel_gen_b, long *rel_gen_d, QuadInfElement *qie, long OB, long BB, long S)
{
  static long NB;

#ifdef DEBUG
  printf("\n-->--> construct_relative_generator:\n");
  printf("OB=%ld, BB=%ld\n",OB,BB);
#endif
  
  for (long i=0; i<qie->num_q; ++i) {
    NB = qie->qlist[i]*BB + OB;
    OB = BB;
    BB = NB;

#ifdef DEBUG
    printf("&&& q=%ld, B=%ld\n",qie->qlist[i],BB);
#endif
  }

  *rel_gen_a = S*((OB * qie->a << 1) + BB * qie->b);
  *rel_gen_b = -S*BB;
  *rel_gen_d = (qie->a << 1);

#ifdef DEBUG
  printf("OB=%ld, BB=%ld, a=%ld, b=%ld, S=%ld\n",OB,BB,qie->a,qie->b,S);
  printf("relgen_a=%ld, relgen_b=%ld, relgen_d=%ld\n",*rel_gen_a,*rel_gen_b,*rel_gen_d);
  printf("-->--> done construct_relative_generator\n\n");
#endif
}


/*
//
// construct relative generator
// computes the coefficients a, b, d of the relative generator (a + b*sqrt(Delta))/d
//
void
construct_relative_generator(long *rel_gen_a, long *rel_gen_b, long *rel_gen_d, QuadInfElement *qie, long Ca, long S)
{
  static long OB, BB, NB;

  OB=1;
  BB=0;
  
#ifdef DEBUG
  printf("\n-->--> construct_relative_generator:\n");
  printf("num_q=%ld\n",qie->num_q);
  printf("OB=%ld\n",OB);
  printf("BB=%ld\n",BB);
#endif
  
  for (long i=0; i<qie->num_q; ++i) {
    NB = qie->qlist[i]*BB + OB;
    OB = BB;
    BB = NB;
    
#ifdef DEBUG  
    printf("&&& q=%ld, B=%ld\n",qie->qlist[i],BB);
#endif
  }

  *rel_gen_a = (OB * qie->a << 1) + BB * qie->b;
  *rel_gen_b = BB;
  *rel_gen_d = (S*Ca << 1);

#ifdef DEBUG
  printf("OB=%ld, BB=%ld, a=%ld, b=%ld, Ca=%ld, S=%ld\n",OB,BB,qie->a,qie->b,Ca,S);
  printf("relgen_a=%ld, relgen_b=%ld, relgen_d=%ld\n",*rel_gen_a,*rel_gen_b,*rel_gen_d);
  printf("-->--> done construct_relative_generator\n\n");
#endif
}
*/


void square(QuadInfElement * qie){
  static long a1, b1, c1, Ca, Cb, Cc;
  static long S = 1, v1, K, T;
  static long precision;
  
#ifdef DEBUG
  printf("\n--> SQUARE:  ");
  print(qie);
#endif
  
  a1 = qie->a;
  b1 = qie->b;
  c1 = qie->c;
  precision = qie->precision;

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT(&S, &v1, b1, a1);

  // K = -v1 c1 (mod L)
  K = -(v1 * c1);

  if (S != 1)
    {
      a1 /= S;
      c1 *= S;
    }

  K %= a1;
  if (K < 0)  K += a1;  

#ifdef DEBUG
  printf("S=%ld, v1=%ld, K=%ld\n",S,v1,K);
  printf("a1=%ld\n",a1);
#endif
 
  // N = L = a1;
  // T = NK
  T = a1 * K;

  // C.a = a1^2 / S^2 = N^2
  // C.b = b1 + 2 a1 K = b1 + 2 T
  // C.c = (S c1 + K (b1 + T)) / L;
  Ca = a1 * a1;
  Cb = b1 + (T << 1);
  Cc = (c1 + (K * (b1 + T))) / a1;

  // Update qie coefficients, distance = 2*dist
  qie->a = Ca;
  qie->b = Cb;
  qie->c = Cc;
  normalize(qie);
  arb_add(qie->distance, qie->distance, qie->distance, precision);

#ifdef DEBUG
  printf("--> square (before reduce): ");
  print(qie);
  test_ideal(qie,"square (before reduce)");
  printf("\n"); 
#endif

  // Reduce and get the coefficients of the relative generator
  reduce(qie);

  // rrelative_generator
  long rel_gen_a, rel_gen_b, rel_gen_d;
  arb_t relative_generator;

  //  construct_relative_generator(&rel_gen_a, &rel_gen_b, &rel_gen_d, qie, Ca, S);
  construct_relative_generator(&rel_gen_a, &rel_gen_b, &rel_gen_d, qie, 1, 0, S);

  // (rel_gen_a + rel_gen_b*sqrt(Delta)) / rel_gen_d
  arb_init(relative_generator);
  arb_mul_si(relative_generator, qie->sqrt_Delta, rel_gen_b, precision);
  arb_add_si(relative_generator, relative_generator, rel_gen_a, precision);
  arb_div_si(relative_generator, relative_generator, rel_gen_d, precision);
  arb_abs(relative_generator, relative_generator);
  arb_log(relative_generator, relative_generator, precision);

  // Update qie distance = 2*dist - log(rel_gen)
  arb_sub(qie->distance, qie->distance, relative_generator, precision);

  arb_clear(relative_generator);
  
#ifdef DEBUG
  printf("--> square (after reduce): ");
  print(qie);
  test_ideal(qie,"square (after)");
  printf("dist = ");
  arb_printn(relative_generator,qie->precision/3,0);
  printf("\n");
#endif
}


//
// nudupl()
//
void nudupl(QuadInfElement * qie){
  static long a1, b1, c1, Ca, Cb, Cc;
  static long S, v1, K, T;
  static long R1, R2, C1, C2, M2, temp;
  static long precision;

#ifdef DEBUG
  printf("\n--> NUDUPL:  ");
  print(qie);
#endif
  
  a1 = qie->a;
  b1 = qie->b;
  c1 = qie->c;
  precision = qie->precision;

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT (&S, &v1, b1, a1);

  // K = -v1 c1 (mod L)
  K = -(v1 * c1);

  if (S != 1)
    {
      a1 /= S;
      c1 *= S;
    }

  K %= a1;
  if (K < 0)  K += a1;  

#ifdef DEBUG
  printf("S=%ld, v1=%ld, K=%ld\n",S,v1,K);
  printf("a1=%ld, bound=%ld\n",a1,qie->NC_BOUND);
#endif
  
  // N = L = a1
  // check if NUCOMP steps are required
  if (a1 <= qie->NC_BOUND) {
#ifdef DEBUG
    printf("EASY!\n");
#endif
    
    // compute with regular squaring formula (result will be reduced)

    // T = NK
    T = a1 * K;

    // C.a = a1^2 / S^2 = N^2
    // C.b = b1 + 2 a1 K = b1 + 2 T
    // C.c = (S c1 + K (b1 + T)) / L;
    Ca = a1 * a1;
    Cb = b1 + (T << 1);
    Cc = (c1 + (K * (b1 + T))) / a1;

#ifdef DEBUG
    printf("Ca=%ld, Cb=%ld, Cc=%ld\n",Ca,Cb,Cc);
#endif
      
    qie->a = Ca;
    qie->b = Cb;
    qie->c = Cc;
    normalize(qie);
    arb_add(qie->distance, qie->distance, qie->distance, precision);

    C2 = 1;
    C1 = 0;
    
#ifdef DEBUG
    printf("--> nudupl (before reduce): ");
    print(qie);
    test_ideal(qie,"nudupl (before reduce)");
    printf("\n");  
#endif
  }
  else {
#ifdef DEBUG
    printf("HARD!\n");
#endif
    // use NUCOMP formulas

    // Execute partial reduction
    R2=a1; R1=K;
    XGCD_PARTIAL(&R2, &R1, &C2, &C1, qie->NC_BOUND);

#ifdef DEBUG
    printf("AFTER PARTIAL:  R2=%ld, R1=%ld, C2=%ld, C1=%ld, numq=%ld\n",R2,R1,C2,C1,qie->num_q);
#endif
    
    // M1 = R1

    // M2 = (R1 b1 - c1 S C1) / L
    M2 = (R1*b1 - c1*C1) / a1;

#ifdef DEBUG
    printf("M1=%ld, M2=%ld\n",R1,M2);
#endif
    
    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    Ca = R1*R1;
    temp = C1*M2;
    if (C1 > 0)
      Ca -= temp;
    else
      Ca = temp - Ca;
 
    // C.b = 2 (N R1 - C.a C2) / C1 - b1
    Cb = (a1*R1 + C2*Ca) << 1;
    Cb /= C1;
    Cb -= b1;

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = (Cb*Cb - qie->Delta) / Ca;
    Cc >>= 2;
      
#ifdef DEBUG
    printf("Ca=%ld, Cb=%ld, Cc=%ld\n",Ca,Cb,Cc);
#endif
  
    qie->a = Ca;
    qie->b = Cb;
    qie->c = Cc;
    arb_add(qie->distance, qie->distance, qie->distance, precision);

#ifdef DEBUG
    printf("--> nudupl (before reduce): ");
    print(qie);
    test_ideal(qie,"nudupl (before reduce)");
    printf("\n");  
#endif
  }
  
  // Reduce and get the coefficients of the relative generator
  reduce(qie);
#ifdef DEBUG
  printf("numq=%ld\n",qie->num_q);
#endif

  // rrelative_generator
  long rel_gen_a, rel_gen_b, rel_gen_d;
  arb_t relative_generator;
  
  arb_init(relative_generator);

  construct_relative_generator(&rel_gen_a, &rel_gen_b, &rel_gen_d, qie, labs(C2), labs(C1), S);
  
  // (rel_gen_a + rel_gen_b*sqrt(Delta)) / rel_gen_d
  arb_mul_si(relative_generator, qie->sqrt_Delta, rel_gen_b, precision);
  arb_add_si(relative_generator, relative_generator, rel_gen_a, precision);
  arb_div_si(relative_generator, relative_generator, rel_gen_d, precision);
  arb_abs(relative_generator, relative_generator);
  arb_log(relative_generator, relative_generator, precision);

  // Update qie distance = 2*dist - log(rel_gen)
  arb_sub(qie->distance, qie->distance, relative_generator, precision);
  
  arb_clear(relative_generator);
  
#ifdef DEBUG
  printf("--> nudupl (after reduce): ");
  print(qie);
  test_ideal(qie,"nudupl (after)");
  printf("dist = ");
  arb_printn(relative_generator,qie->precision/3,0);
  printf("\n");
#endif
}



//
// close
// compute an ideal with distance <= dist, as close as posible
//
void close_ideal(QuadInfElement * qie, double dist){
  long j = 0, k = 0, s = 1, i, ex;
  long bound = (long) floor(dist);

#ifdef DEBUG
  printf("\nCLOSE: Delta=%ld, dist=%lf, bound=%ld\n",qie->Delta,dist,bound);
  fflush(stdout);
#endif
  
  assign_one(qie);

  if (bound <= 1)
    return;

  // Compute binary expansion of ex (hi order to low order)
  ex = bound;
  while (!(ex == 1)) {
    j <<= 1;
    if (ex % 2 == 1) {
      ++j;
    }
    ex >>= 1;
    ++k;
  }

#ifdef DEBUG
  printf("\ns=%ld\n",s);
  fflush(stdout);
#endif
  
  adjust(qie, s);

  for (i = 1; i <= k; ++i) {
    s <<= 1;
    if (j % 2 == 1)
      ++s;

#ifdef DEBUG
    printf("\ns=%ld\n",s);
    fflush(stdout);
#endif

#ifdef USE_MULT
    square(qie);
#else
    nudupl(qie);
#endif
    
    adjust(qie, s);

    j >>= 1;
  }

#ifdef DEBUG
  printf("\ndist=%lf\n",dist);
  fflush(stdout);
#endif
  adjust(qie, dist);

#ifdef DEBUG
  printf("\nDONE CLOSE!  dist=%lf, ",dist);
  print(qie);
  test_ideal(qie,"nuclose");
#endif
}


// refine_regulator
// Checks whether the given ideal is in the order of the specified discriminant
void refine_regulator(arb_t new_reg, long Delta, double reg, long precision) {
    QuadInfElement qie;
    init(&qie,Delta,precision);
    
    close_ideal(&qie, floor(reg));
    adjust_to_one(&qie);

    arb_set(new_reg,qie.distance);

    clear(&qie);
}



bool test_refined_regulator(arb_t new_reg, double reg, long precision) {
  arb_t diff,max_diff;

  arb_init(max_diff);
  arb_set_d(max_diff,0.1);
  
  arb_init(diff);
  arb_set_d(diff,reg);
  arb_sub(diff,diff,new_reg,precision);
  arb_abs(diff,diff);

  bool good = (bool) arb_lt(diff,max_diff);
    
#ifdef DEBUG
  printf("\nTEST_REFINED_REGULATOR:\n");
  printf("double:   %lf\n",reg);
  printf("refined:  ");
  arb_printn(new_reg,precision/3,0);
  printf("\n");
  printf("abs(diff):  ");
  arb_printn(diff,precision/3,0);
  printf("\n");
  printf("max_diff:  ");
  arb_printn(max_diff,precision/3,0);
  printf("\n");
#endif  

  arb_clear(diff);
  arb_clear(max_diff);
  
  return good;
}

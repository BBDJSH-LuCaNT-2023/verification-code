#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define archt_nocheck 1
#include "kbessel.c"
#include "mpfi_io.h"

/* We'll have M<Mmax. */
#define Mmax 600
mpfi_t T[Mmax][Mmax],b[Mmax],w,v[Mmax],x[Mmax]; 
mpfi_t xx[Mmax]; 
mpfi_t TI[Mmax][Mmax];
mpfi_t S[Mmax][Mmax];
mpfr_t TIb[Mmax][Mmax];
mpfr_t xrdb[Mmax],x1b[Mmax],x1bb[Mmax],x2b[Mmax],vb[Mmax];
mpfr_t wdbA[2][Mmax],wdbB[2][Mmax];
mpfi_t wdiag[Mmax];
int M,N,oe;
mpfr_t Y;

/* wrapper variables and functions for calling Holger's archt library */
int kbes_pr;
static mpfi_ptr wrapper_f;
static mpfi_t wrapper_r;

static void my_kbessel_init(mpfr_t r,double ymax,double ymaxpr,double ymin) {
  double rr,z;

	mpfi_init_set_fr(wrapper_r,r);
  rr = mpfr_get_d(r,GMP_RNDN);
  z=sqrt(ymaxpr*ymaxpr/(rr*rr)-1);
  z=rr*(z-atan(z))*M_LOG2E;

  kbes_pr=(int)z+mpfr_get_default_prec()+15;
	wrapper_f = init_kbessel(kbes_pr);
}

static void my_kbessel_clear() {
	clear_kbessel(wrapper_f);
	mpfi_clear(wrapper_r);
}

static void mpfi_whittaker(mpfi_t s,mpfi_t y) {
	kbessel(wrapper_f,wrapper_r,y);
	mpfi_sqrt(s,y);
	mpfi_mul(s,s,wrapper_f);
}

void mpfr_myprint(mpfr_t x) {
  mpfr_out_str(stderr,10,0,x,GMP_RNDN);
}

void mpfr_myprint1(mpfr_t x,int n) {
  mpfr_out_str(stderr,10,n,x,GMP_RNDN);
}

void mpfi_myprint(mpfi_t x) {
  mpfr_out_str(stderr,10,0,&x->left,GMP_RNDD);
  fprintf(stderr,"\n");
  mpfr_out_str(stderr,10,0,&x->right,GMP_RNDU);
  fprintf(stderr,"\n");
}

/* Print midpoint and radius *** NOW SLOPPY! *** */
void mpfi_print1sloppy(mpfi_t x) {
  mpfr_t xmid;
  
  mpfr_init(xmid);
  mpfi_mid(xmid,x);
  mpfr_out_str(stderr,10,0,xmid,GMP_RNDN);
  fprintf(stderr,"  abs err<=");
  mpfr_sub(xmid,&x->right,xmid,GMP_RNDU);
  mpfr_out_str(stderr,10,4,xmid,GMP_RNDU);
  fprintf(stderr,"\n");
  mpfr_clear(xmid);
}

void print_endpoints(mpfi_t x) {
  mpfr_out_str(stdout,10,0,&x->left,GMP_RNDD);
  printf("\n");
  mpfr_out_str(stdout,10,0,&x->right,GMP_RNDU);
  printf("\n");
}

/**********************************************************************************************
 * For given exact 0<r<y, compute an upper bound on   
 *   (y^(1/2)/(y^2-r^2)^(1/4))*exp(-r*u(y/r))  where  u(q)=sqrt(q^2-1)-arctan(sqrt(q^2-1)).
 * 
 * We use the fact that both functions   (y^(1/2)/(y^2-r^2)^(1/4))   and   exp(-r*u(y/r))   
 * are increasing w.r.t. r (for any fixed y) and decreasing w.r.t. y (for any fixed r).
 **********************************************************************************************/
void explmajorantfunc(mpfr_t rop,mpfr_t y,mpfr_t r) {
  mpfr_t y2,r2,t1,t2,q,t,u;

  mpfr_init(y2); mpfr_init(r2); mpfr_init(t1); mpfr_init(t2); 
  mpfr_init(q); mpfr_init(t);  mpfr_init(u); 
  /* Set y2= l.b. (lower bound of) y^2 */
  mpfr_sqr(y2,y,GMP_RNDD);
  /* Set r2= u.b. (upper bound of) r^2 */
  mpfr_sqr(r2,r,GMP_RNDU);
  /* Set t2= u.b.  y^(1/2)/(y^2-r^2)^(1/4).  */
  mpfr_sqrt(t2,y,GMP_RNDU);
  mpfr_sub(t1,y2,r2,GMP_RNDD);
  mpfr_root(t1,t1,4,GMP_RNDD);
  mpfr_div(t2,t2,t1,GMP_RNDU);

  /* Set q= l.b.  y/r */
  mpfr_div(q,y,r,GMP_RNDD);
  /* Set u= l.b.  u(q)=sqrt(q^2-1)-arctan(sqrt(q^2-1)) */
  /* Here we use the fact that u(q) is an increasing function of q, and that
     t->t-arctan(t) is an increasing function of t. */
  mpfr_sqr(t,q,GMP_RNDD);
  mpfr_sub_ui(t,t,1,GMP_RNDD);
  mpfr_sqrt(t,t,GMP_RNDD);
  mpfr_atan(t1,t,GMP_RNDU);
  mpfr_sub(u,t,t1,GMP_RNDD);
  /* Set t1= u.b.  exp(-r*u(y/r)) */
  mpfr_mul(t1,u,r,GMP_RNDD);
  mpfr_neg(t1,t1,GMP_RNDU);
  mpfr_exp(t1,t1,GMP_RNDU);

  /* Compute result. */
  mpfr_mul(rop,t2,t1,GMP_RNDU);
  /* Clear variables. */
  mpfr_clear(y2); mpfr_clear(r2); mpfr_clear(t1); mpfr_clear(t2); 
  mpfr_clear(q); mpfr_clear(t); mpfr_clear(u); 
}


/**********************************************************************************************
 * Compute an upper bound on the majorant function \fM^\eta_{M,Y,R} defined in 
 * section 5 of the manual (for eta=7/64),  via Lemma 5.5.
 **********************************************************************************************/
void MMajorant_basic(mpfr_t rop,int M,mpfr_t Y,mpfr_t R) {
  mpfr_t pilb,t,fY,R2,fY2,s;

  mpfr_init(pilb); mpfr_init(t); mpfr_init(fY); 
  mpfr_init(R2); mpfr_init(fY2); mpfr_init(s); 
  /* Set fY= l.b. 2*Pi*(M+1)*Y. */
  mpfr_const_pi(pilb,GMP_RNDD);
  mpfr_mul_ui(t,pilb,2,GMP_RNDD);
  mpfr_mul_ui(t,t,M+1,GMP_RNDD);
  mpfr_mul(fY,t,Y,GMP_RNDD);
  /* Set t= u.b.  fY^(1/2)/(fY^2-R^2)^(1/4)*exp(-R*u(fY/R)) */
  explmajorantfunc(t,fY,R);
  /* Set s= u.b.  1+ (M+1)/sqrt(fY^2-R^2) */
  mpfr_sqr(R2,R,GMP_RNDU);
  mpfr_sqr(fY2,fY,GMP_RNDD);
  mpfr_sub(s,fY2,R2,GMP_RNDD);
  mpfr_sqrt(s,s,GMP_RNDD);
  mpfr_ui_div(s,1,s,GMP_RNDU);
  mpfr_mul_ui(s,s,M+1,GMP_RNDU);
  mpfr_add_ui(s,s,1,GMP_RNDU);
  /* Set rop= u.b.  3*t*s. */
  mpfr_mul_ui(rop,t,3,GMP_RNDU);
  mpfr_mul(rop,rop,s,GMP_RNDU);
  mpfr_clear(pilb); mpfr_clear(t); mpfr_clear(fY); 
  mpfr_clear(R2); mpfr_clear(fY2); mpfr_clear(s); 
}


/**********************************************************************************************
 * Compute an upper bound on the majorant function \fM^\eta_{M,Y,R} defined in 
 * section 5 of the manual (for eta=7/64), using (5.7).
 **********************************************************************************************/
void MMajorant(mpfr_t rop,int M,mpfr_t Y,mpfr_t R) {
  mpfr_t pilb,t,sqrtpihalf,t1,s,fY;
  int m;

  mpfr_init(pilb); mpfr_init(t); mpfr_init(sqrtpihalf);
  mpfr_init(t1); mpfr_init(s); mpfr_init(fY);
  MMajorant_basic(rop,M+5,Y,R);
  mpfr_const_pi(pilb,GMP_RNDD);
  /* Set sqrtpihalf = u.b.  sqrt(Pi/2) */
  mpfr_const_pi(t,GMP_RNDU);
  mpfr_mul_d(t,t,0.5,GMP_RNDU);
  mpfr_sqrt(sqrtpihalf,t,GMP_RNDU);
  for(m=M+1; m<=M+5; m++) {
    /* Set t= u.b.  sigma_{2*eta}(m)/m^{eta+1/2}. */
    KimSarnakbound(t,m);
    mpfr_sqrt_ui(t1,m,GMP_RNDD);
    mpfr_div(t,t,t1,GMP_RNDU);
    /* Set fY= l.b.  2*Pi*m*Y. */
    mpfr_mul_ui(t,pilb,2,GMP_RNDD);
    mpfr_mul_ui(t,t,m,GMP_RNDD);
    mpfr_mul(fY,t,Y,GMP_RNDD);
    /* Set s= u.b.  fY^(1/2)/(fY^2-R^2)^(1/4)*exp(-R*u(fY/R)) */
    explmajorantfunc(s,fY,R);
    /* Add u.b. of  sqrt(Pi/2)*t*s  to rop. */
    mpfr_mul(t,sqrtpihalf,t,GMP_RNDU);
    mpfr_mul(t,t,s,GMP_RNDU);
    mpfr_add(rop,rop,t,GMP_RNDU);
  }
  mpfr_clear(pilb); mpfr_clear(t); mpfr_clear(sqrtpihalf);
  mpfr_clear(t1); mpfr_clear(s); mpfr_clear(fY);
}

/**********************************************************************************************
 * Same as MMajorant_basic, but compute non-rigorously using doubles.
 * (This is only used when fixing appropriate choices of M,N.)
 **********************************************************************************************/
double MMajorant_basic_heuristic(int M,double Y,double R) {
  double a,a1,a2,a3,t,u;

  a=2*3.141592653589793238462643383*(M+1)*Y;
  a1=a*a-R*R;
  if (a1<=0) {fprintf(stderr,"MMajorant_basic_heuristic called with too small M.\n"); exit(1);}
  a2=sqrt(a1);
  a3=3*sqrt(a)/sqrt(a2)*(1+(M+1)/a2);
  t=a/R;
  u=sqrt(t*t-1);
  u=u-atan(u);
  return(a3*exp(-R*u));
}




/*
 * Compute the pullback  xp+I*yp  of  x+I*y,  working in precision kbes_pr.
 * Note: We don't guarantee that the box returned for  xp+I*yp  completely belongs to the 
 * (standard) fundamental domain for PSL(2,Z), but we do guarantee that there is some 
 * transformation  T  in PSL(2,Z) such that  T(x+I*y) = xp+I*yp.
 */
int PSL2Zpullback(mpfi_t xp,mpfi_t yp,mpfi_t x,mpfi_t y) {
  mpfi_t temp,temp2,Hi,a;
  mpfr_t H;
  int oldpr;

  oldpr=mpfr_get_default_prec();
  mpfr_set_default_prec(kbes_pr);
  mpfi_init(temp); mpfi_init(temp2); mpfi_init(Hi); mpfi_init(a); mpfr_init(H);
  mpfi_set(xp,x);
  mpfi_set(yp,y);
  while(1) {
    mpfi_set_d(temp,0.5);
    mpfi_sub(temp,temp,xp);
    mpfr_floor(H,&temp->left);
    mpfi_set_fr(Hi,H);
    /* Here we know that the interval Hi contains an integer, so the following certainly 
       applies a transformation from PSL(2,Z) to xp+I*yp!  */
    mpfi_add(xp,xp,Hi);
    /* Compute a=abs(xp+I*yp)^2. */
    mpfi_sqr(temp,xp);
    mpfi_sqr(temp2,yp);
    mpfi_add(a,temp,temp2);
    if (mpfi_cmp_d(a,0.99999)>0) {
      mpfi_clear(temp); mpfi_clear(temp2); mpfi_clear(Hi); mpfi_clear(a); mpfr_clear(H);
      mpfr_set_default_prec(oldpr);
      return 1;
    }
    /* Replace  xp+I*yp  by  -1/(xp+I*yp). */
    mpfi_div(xp,xp,a);
    mpfi_neg(xp,xp);
    mpfi_div(yp,yp,a);
  }
}

/*
 * Same as above (viz., compute the pullback  xp+I*yp  of  x+I*y), but heuristically,
 * working with doubles.
 */
int PSL2Zpullback_heur(double *xp,double *yp,double x,double y) {
  double a;
  int H;

  *xp=x; *yp=y;
  while(1) {
    H=floor(0.5-*xp);
    xp+=H;
    a=(*xp)*(*xp)+(*yp)*(*yp);
    if(a>0.99999) return 1;
    *xp=-*xp/a;
    *yp=*yp/a;
  }
}

/* Print the Frobenius norm of Andy's matrix "H" minus the identity matrix. 
   Here H in fact equals diag[W(2*Pi*2Y)/sqrt(2),W(2*Pi*3Y)/sqrt(3),...]*TI.
*/
void checkfrobeniusnorm(int MM) {
  mpfi_t t1,s;
  int j,k;

  mpfi_init(t1); mpfi_init(s); 
  mpfi_set_ui(s,0);
  for(j=0; j<MM; j++) {
    for(k=0; k<MM; k++) {
      mpfi_mul(t1,TI[j][k],wdiag[j+1]);
      /* Now t1 equals the j,k-entry of H^{-1}, NEGATED. */
      if(j==k) mpfi_add_ui(t1,t1,1);
      mpfi_sqr(t1,t1);
      mpfi_add(s,s,t1);
    }
  }
  mpfi_sqrt(s,s);
//  fprintf(stderr,"Frobenius norm of H minus identity matrix:\n");
//  mpfi_myprint(s);
  mpfi_clear(t1); mpfi_clear(s); 
}


/* Compute T, b, w, v for the given (exact!) r-value. */
void initsystem(mpfr_t r,mpfr_t r0,mpfr_t eps) {
  mpfi_t xj,xjp,yjp,tpxjp,tpyjp,mtpxjp,mtpyjp,Yi,twopi,ni2,tpxj,t1,t2,s,tpYi,mtpYi;
  mpfr_t t3,t4,r0eps;
  double Yd,xp,yp,maxyp;
  int j,k,m,pr;

  Yd=mpfr_get_d(Y,GMP_RNDD);
  /* Compute the largest y-value among the pullbacks. */
  maxyp=0.0;
  for(j=1; j<=N; j++) {
    PSL2Zpullback_heur(&xp,&yp,(j-0.5)/(2.0*N),Yd);
    if(yp>maxyp) maxyp=yp;    
  }
  /* We will only compute K_{ir}(y) for y satisfying  y <= 2*Pi*M*maxyp  and  y >= 2*Pi*Y. */

//  fprintf(stderr,"Calling kbessel_init.\n");
  my_kbessel_init(r,6.283185309*M*maxyp,6.283185309*M*Yd,6.2831*Yd);  
//  fprintf(stderr,"kbessel_init done; will use kbes_pr=%d.\n",kbes_pr);
  /* Set all T,b,v,w to 0. */
  for(k=0; k<M-1; k++) {
    for(m=0; m<M-1; m++) mpfi_set_ui(T[k][m],0);
    mpfi_set_ui(b[k],0);
    mpfi_set_ui(v[k],0);
  }
  mpfi_set_ui(w,0);
  /* Init our local mpfi_t's. */
  mpfi_init2(xj,kbes_pr); mpfi_init2(xjp,kbes_pr); mpfi_init2(yjp,kbes_pr); 
  mpfi_init(tpxjp); mpfi_init2(tpyjp,kbes_pr);
  mpfi_init(mtpxjp); mpfi_init2(mtpyjp,kbes_pr); 
  mpfi_init2(Yi,kbes_pr);  mpfi_init2(twopi,kbes_pr);
  mpfi_init(ni2);  mpfi_init2(tpxj,kbes_pr);  mpfi_init(t1);  mpfi_init(t2);  
  mpfi_init2(s,kbes_pr);
  mpfi_init2(tpYi,kbes_pr); mpfi_init2(mtpYi,kbes_pr);
  mpfr_init(t3); mpfr_init(t4); mpfr_init(r0eps);
  /* Set twopi=2*Pi */
  mpfi_const_pi(twopi);
  mpfi_mul_ui(twopi,twopi,2);

  mpfi_set_fr(Yi,Y);
  /* Set ni2=2/N. */
  mpfi_set_ui(ni2,2);
  mpfi_div_ui(ni2,ni2,N);
  for(j=1; j<=N; j++) {
    /*    fprintf(stderr,"j=%d of %d\n",j,N); */
    /* Set xj=(j-0.5)/(2N). */
    mpfi_set_si(xj,j+j-1);
    mpfi_div_ui(xj,xj,4*N);
    mpfi_mul(tpxj,xj,twopi);
    /* Set xjp,yjp so that   xjp+I*yjp  is the pullback of  xj+I*Y.  */
    PSL2Zpullback(xjp,yjp,xj,Yi);
    /* Note that we are not guaranteed that the pullback lies in the standard fundamental
       domain for PSL(2,Z), but the only thing we need for our equation system to be 
       valid is that yjp>=sqrt(3)/2=0.8660254..., so that the cutoff error bound is correct. */
    if (mpfi_cmp_d(yjp,0.867)<=0) {
      fprintf(stderr,"Error: Pullback is not safely above sqrt(3)/2.\n");
      exit(1);
    }
    mpfi_mul(tpyjp,yjp,twopi);
    mpfi_mul(tpxjp,xjp,twopi);
    for(m=1; m<=M; m++) {
      mpfi_mul_ui(mtpyjp,tpyjp,m);
      mpfi_mul_ui(mtpxjp,tpxjp,m);
      /* Set s=(2/N)*cos(2*Pi*m*xjp)*W_{ir}(2*Pi*m*yjp)/sqrt(m) */
      mpfi_whittaker(s,mtpyjp); 
      mpfi_set_ui(t1,m);
      mpfi_sqrt(t1,t1);
      mpfi_div(s,s,t1);
      if(oe) mpfi_sin(t2,mtpxjp); else mpfi_cos(t2,mtpxjp);
      mpfi_mul(s,s,t2);
      mpfi_mul(s,s,ni2);
      for(k=1; k<=M; k++) {
	/* Set t1=cos(2*Pi*k*xj). */
	mpfi_mul_ui(t1,tpxj,k);
	if(oe) mpfi_sin(t1,t1); else mpfi_cos(t1,t1);
	/* Set t2=s*t1=(2/N)*cos(2*Pi*m*xjp)*cos(2*Pi*k*xj)*W_{ir}(2*Pi*m*yjp)/sqrt(m) */
	mpfi_mul(t2,s,t1);
	if(k==1) {    /* Record t2 as contribution to either v or w. */
	  if(m>1) mpfi_sub(v[m-2],v[m-2],t2); else mpfi_sub(w,w,t2);
	} else {    /* Record t2 as contribution to either T or b. */
	  if(m>1) mpfi_add(T[k-2][m-2],T[k-2][m-2],t2); else mpfi_sub(b[k-2],b[k-2],t2);
	}
      }
    }
  }
  /* Add the "delta-term" contributions to T and w. */
  mpfi_mul(tpYi,Yi,twopi);
  for(m=1; m<=M; m++) {
    /* Set s=W_{ir}(2*Pi*m*Y)/sqrt(m) */
    mpfi_mul_ui(mtpYi,tpYi,m);
    mpfi_whittaker(s,mtpYi); 
    mpfi_set_ui(t1,m);
    mpfi_sqrt(t1,t1);
    mpfi_div(s,s,t1);
    mpfi_set(wdiag[m-1],s);   /* Used in checkfrobeniusnorm. */
    /* Record s as contribution to either T or w. */
    if(m>1)
      mpfi_sub(T[m-2][m-2],T[m-2][m-2],s);
    else
      mpfi_add(w,w,s);
  }
  my_kbessel_clear();
  /* At present our vector b gives b^\sharp(r), in the notation of the manual.
     Now add e^* to b, so as to get the true vector b(r). See Lemma 5.6.  */

  mpfr_sqrt_ui(t3,3,GMP_RNDD);
  mpfr_mul_d(t3,t3,0.5,GMP_RNDD);
  mpfr_add(r0eps,r0,eps,GMP_RNDU);
  MMajorant(t4,M,t3,r0eps);
  mpfr_mul_ui(t4,t4,2,GMP_RNDU);
  for(k=0; k<M-1; k++) mpfi_increase(b[k],t4);
  for(k=0; k<M-1; k++) {
    MMajorant(t4,2*N-k-3,Y,r0eps);
    mpfi_increase(b[k],t4);
  }
#if 0
  fprintf(stderr,"Now printing system, %d.\n",M-1);
  for(m=0; m<M-1; m++) for(k=0; k<M-1; k++) {
      mpfi_myprint(T[m][k]);
    }
  exit(1);
#endif
  mpfi_clear(xj); mpfi_clear(xjp); mpfi_clear(yjp); 
  mpfi_clear(tpxjp); mpfi_clear(tpyjp);
  mpfi_clear(mtpxjp); mpfi_clear(mtpyjp); 
  mpfi_clear(Yi);  mpfi_clear(twopi);
  mpfi_clear(ni2);  mpfi_clear(tpxj);  mpfi_clear(t1);  mpfi_clear(t2);  
  mpfi_clear(s);
  mpfi_clear(tpYi); mpfi_clear(mtpYi);
  mpfr_clear(t3); mpfr_clear(t4); mpfr_clear(r0eps);
}




/**********************************************************************************************
 * Invert a given real MM x MM-matrix T (global matrix of mpfi_t's), if possible.
 * If the function does not end with error, it returns a matrix TI of mpfi_t's with the 
 * property that for EVERY exact real matrix T' with entries inside the corresponding 
 * T-intervals, T' is invertible and its inverse has entries inside the corresponding 
 * TI-intervals.
 **********************************************************************************************/
int invertT(int MM) {
  int j,jj,k,maxp,e;
  mpfr_t maxv,t;
  mpfi_t s,q1,q2;
  /* Set TI = identity matrix. */
  for(k=0; k<MM; k++) 
    for(j=0; j<MM; j++)
      if(j==k) mpfi_set_ui(TI[j][k],1); else mpfi_set_ui(TI[j][k],0);
  /* Init local variables. */
  mpfr_init(t); mpfr_init(maxv); 
  mpfi_init(q1); mpfi_init(q2);  mpfi_init(s);

  /* For each column j=0,1,...,MM-1. */
  for(j=0; j<MM; j++) {
    mpfr_set_ui(maxv,0,GMP_RNDN);
    maxp=-1;
    for(k=j; k<MM; k++) {
      mpfi_mig(t,T[k][j]);
      if(mpfr_cmp(t,maxv)>0) {mpfr_set(maxv,t,GMP_RNDN); maxp=k;}
    }
    if (maxp==-1) {
      fprintf(stderr,"Linear dependence occurred when %d lines remained in the Gauss elimination.\n",MM-j);
      exit(1);
    }
    /*    maxp=j; */   /* Can be a loss, but makes life simple. */
    if(maxp>j) {
      /* Swap rows j and maxp. */
      for(k=0; k<MM; k++) {mpfi_swap(T[j][k],T[maxp][k]); mpfi_swap(TI[j][k],TI[maxp][k]);}
    }

    /* Set q1=1/T[j][j]. */
    mpfi_inv(q1,T[j][j]);
    /* Multiply row j with q1. */
    for(k=0; k<MM; k++) {
      if (k!=j) mpfi_mul(T[j][k],T[j][k],q1);
      mpfi_mul(TI[j][k],TI[j][k],q1);
    }
    /* Set T[j][j]=EXACTLY 1. */
    mpfi_set_ui(T[j][j],1);
    for(jj=0; jj<MM; jj++) 
      if(jj!=j) {
	/* Set q2=-T[jj][j]. */
	mpfi_neg(q2,T[jj][j]);
	/* Add q2 times row j to row jj. */
	for(k=0; k<MM; k++) {
	  if (k!=j) {
	    mpfi_mul(s,T[j][k],q2);
	    mpfi_add(T[jj][k],T[jj][k],s);
	  }
	  mpfi_mul(s,TI[j][k],q2);
	  mpfi_add(TI[jj][k],TI[jj][k],s);
	}
	/* Set T[jj][j]=EXACTLY 0. */
	mpfi_set_ui(T[jj][j],0);
      }
  }
  /* Clear variables. */
  mpfr_clear(t); mpfr_clear(maxv); 
  mpfi_clear(q1); mpfi_clear(q2);  mpfi_clear(s);
}


/**********************************************************************************************
 * Compute x = T^{-1} b  (store this in the global variable x)  and  c=w+vx.
 * The following should have been previously computed (global variables):
 *   TI (giving T^{-1}), b,w,v.
 **********************************************************************************************/
void computexc(mpfi_t c) {
  mpfi_t f;
  int k,m;
  mpfi_init(f);
  for(k=0; k<M-1; k++) {
    mpfi_set_ui(x[k],0);
    for(m=0; m<M-1; m++) {
      mpfi_mul(f,TI[k][m],b[m]);
      mpfi_add(x[k],x[k],f);
    }
  }
  mpfi_set(c,w);
  for(k=0; k<M-1; k++) {mpfi_mul(f,v[k],x[k]); mpfi_add(c,c,f);}
  mpfi_clear(f);
}


/**********************************************************************************************
 * "Compute" c(r^*), via bound on its absolute value; see Lemma 5.7.
 **********************************************************************************************/
void crstarbound(mpfi_t c,mpfr_t r0,mpfr_t eps) {
  mpfr_t t3,t4,r0eps;

  mpfr_init(t3); mpfr_init(t4); mpfr_init(r0eps);
  mpfr_add(r0eps,r0,eps,GMP_RNDU);
  mpfr_sqrt_ui(t3,3,GMP_RNDD);
  mpfr_mul_d(t3,t3,0.5,GMP_RNDD);
  MMajorant(t4,M,t3,r0eps);
  mpfr_mul_ui(t4,t4,2,GMP_RNDU);
  MMajorant(t3,2*N-2,Y,r0eps);
  mpfr_add(t3,t3,t4,GMP_RNDU);
  mpfi_set_ui(c,0);
  mpfi_increase(c,t3);
  mpfr_clear(t3); mpfr_clear(t4); mpfr_clear(r0eps);
}


/**********************************************************************************************
 * Compute an upper bound on \fW_r^{(j)}(y), the majorant function defined in (7.2), (7.3). 
 * Here we assume r>0, y>0 and j is 1 or 2.
 **********************************************************************************************/
void whittakerderbound(mpfr_t fW,mpfr_t r,mpfr_t y,int j) {
  mpfr_t y2,r2,r2q,q,t,t1,t2,t3,u;

  mpfr_init(y2); mpfr_init(r2); mpfr_init(r2q); mpfr_init(q);
  mpfr_init(t);  mpfr_init(t1); mpfr_init(t2); mpfr_init(t3); mpfr_init(u);

  /* Set y2= l.b. (lower bound of) y^2 */
  mpfr_sqr(y2,y,GMP_RNDD);
  /* Set r2= u.b. (upper bound of) r^2 */
  mpfr_sqr(r2,r,GMP_RNDU);
  /* Set r2q= u.b.  r^2+1/4 */
  mpfr_add_d(r2q,r2,0.25,GMP_RNDU);
  if(mpfr_cmp(y2,r2q)<0) {
    /* Set fW= u.b.  128r^(1/6) (if j=1) or  655r^(1/6) (if j=2). */
    mpfr_root(fW,r,6,GMP_RNDU);    
    if(j==1) mpfr_mul_ui(fW,fW,128,GMP_RNDU); else mpfr_mul_ui(fW,fW,655,GMP_RNDU); 
  } else {
    /* Set q= l.b.  y/r */
    mpfr_div(q,y,r,GMP_RNDD);
    /* Note here we know  y^2>=r^2+1/4; hence y>r and q>=1; but let's test this anyway. */
    if(mpfr_cmp_ui(q,1)<0) {fprintf(stderr,"q<1 although y^2>=r^2+1/4.\n"); exit(1);}
    /* Set u= l.b.  u(q)=sqrt(q^2-1)-arctan(sqrt(q^2-1)) */
    /* Here we use the fact that u(q) is an increasing function of q, and that
       t->t-arctan(t) is an increasing function of t. */
    mpfr_sqr(t,q,GMP_RNDD);
    mpfr_sub_ui(t,t,1,GMP_RNDD);
    mpfr_sqrt(t,t,GMP_RNDD);
    mpfr_atan(t1,t,GMP_RNDU);
    mpfr_sub(u,t,t1,GMP_RNDD);
    /* Set t= u.b.  sqrt(y)*exp(-r*u(y/r)) */
    mpfr_sqrt(t,y,GMP_RNDU);
    mpfr_mul(t1,r,u,GMP_RNDD);
    mpfr_neg(t1,t1,GMP_RNDU);
    mpfr_exp(t1,t1,GMP_RNDU);
    mpfr_mul(t,t,t1,GMP_RNDU);
    /* Set t2= u.b.  12r^(-1/3) (if j=1) or  55r^(-1/3) (if j=2). */
    mpfr_root(t2,r,3,GMP_RNDD);
    if(j==1) mpfr_ui_div(t2,12,t2,GMP_RNDU); else mpfr_ui_div(t2,55,t2,GMP_RNDU);
    /* Set t3= u.b.  7/(y^2-r^2)^(1/4) (if j=1) or  33/(y^2-r^2)^(1/4) (if j=2). */
    mpfr_sub(t1,y2,r2,GMP_RNDD);
    /* If mpfr behaves as supposed to, we should certainly have t1>=1/4 here! */
    if(mpfr_cmp_ui(t1,0)<=0) {fprintf(stderr,"Very strange rounding behaviour!\n"); exit(1);}
    mpfr_root(t3,t1,4,GMP_RNDD);
    if(j==1) mpfr_ui_div(t3,7,t3,GMP_RNDU); else mpfr_ui_div(t3,33,t3,GMP_RNDU);
    /* Set fW= u.b.  t*min(t2,t3) */
    if(mpfr_cmp(t2,t3)<=0) mpfr_mul(fW,t,t2,GMP_RNDU); else mpfr_mul(fW,t,t3,GMP_RNDU);
  }
  mpfr_clear(y2); mpfr_clear(r2); mpfr_clear(r2q); mpfr_clear(q);
  mpfr_clear(t);  mpfr_clear(t1); mpfr_clear(t2); mpfr_clear(t3); mpfr_clear(u);
}

/**********************************************************************************************
 * Compute upper bounds on \fW_{r0+eps}^{(j)}(m*Pi*sqrt(3))/sqrt(m) 
 * and    \fW_{r0+eps}^{(j)}(m*2*Pi*Y)/sqrt(m)   for m=1,2,...,M and j=1,2;
 * store these in the global variables  wdbA, wdbB.
 **********************************************************************************************/
void wdbinit(mpfr_t r0,mpfr_t eps) {
  mpfr_t r0eps,pilb,tA,tB,tAm,tBm,t,invsqrtm;
  int m,j;

  mpfr_init(r0eps); mpfr_init(pilb); 
  mpfr_init(tA); mpfr_init(tB); mpfr_init(tAm); mpfr_init(tBm);
  mpfr_init(t); mpfr_init(invsqrtm);
  /* Set r0eps= u.b.  r0+eps. */
  mpfr_add(r0eps,r0,eps,GMP_RNDU);
  /* Set  tA= l.b.  Pi*sqrt(3). */
  mpfr_const_pi(pilb,GMP_RNDD);
  mpfr_set_ui(tA,3,GMP_RNDD);
  mpfr_sqrt(tA,tA,GMP_RNDD);
  mpfr_mul(tA,tA,pilb,GMP_RNDD);
  /* Set  tB= l.b.  2*Pi*Y. */
  mpfr_mul(tB,pilb,Y,GMP_RNDD);
  mpfr_mul_ui(tB,tB,2,GMP_RNDD);

  for(m=1; m<=M; m++) {
    /* Set  tA= l.b.  m*Pi*sqrt(3)  and  tB= l.b.  m*2*Pi*Y. */
    mpfr_mul_ui(tAm,tA,m,GMP_RNDD);
    mpfr_mul_ui(tBm,tB,m,GMP_RNDD);
    /* Set invsqrtm= u.b.  1/sqrt(m). */
    mpfr_sqrt_ui(t,m,GMP_RNDD);
    mpfr_ui_div(invsqrtm,1,t,GMP_RNDU);
    for(j=1; j<=2; j++) {
      whittakerderbound(t,r0eps,tAm,j);
      mpfr_mul(wdbA[j-1][m-1],t,invsqrtm,GMP_RNDU);
      whittakerderbound(t,r0eps,tBm,j);
      mpfr_mul(wdbB[j-1][m-1],t,invsqrtm,GMP_RNDU);
    }
  }
  mpfr_clear(r0eps); mpfr_clear(pilb); 
  mpfr_clear(tA); mpfr_clear(tB); mpfr_clear(tAm); mpfr_clear(tBm);
  mpfr_clear(t); mpfr_clear(invsqrtm);
}

/**********************************************************************************************
 * For m a positive integer: Compute an upper bound on \sigma_{2*eta}(m)/m^eta,  
 * for eta=7/64.    (Our algorithm for finding the divisors d is extremely inefficient 
 * but this should not matter compared with other things.)
 **********************************************************************************************/
int KimSarnakbound(mpfr_t rop,int m) {
  int d;
  mpfr_t t,ld,f1;
  mpfr_init(t);  mpfr_init(ld);  mpfr_init(f1);
  mpfr_set_ui(rop,0,GMP_RNDU);
  /* Set f1= u.b.  7/32. */
  mpfr_set_ui(f1,7,GMP_RNDU);
  mpfr_div_ui(f1,f1,32,GMP_RNDU);
  /* For each divisor d of m. */
  for(d=1; d<=m; d++) 
    if(m%d==0) {
      /* Compute t= u.b.  d^(7/32). */
      mpfr_set_ui(ld,d,GMP_RNDU);
      mpfr_log(ld,ld,GMP_RNDU);
      mpfr_mul(t,ld,f1,GMP_RNDU);
      mpfr_exp(t,t,GMP_RNDU);
      mpfr_add(rop,rop,t,GMP_RNDU);
    }
  /* Compute t= u.b.  m^(-7/64). */
  mpfr_set_si(f1,-7,GMP_RNDU);
  mpfr_div_ui(f1,f1,64,GMP_RNDU);
  mpfr_mul(t,ld,f1,GMP_RNDU);
  mpfr_exp(t,t,GMP_RNDU);
  mpfr_mul(rop,rop,t,GMP_RNDU);
  mpfr_clear(t);  mpfr_clear(ld);  mpfr_clear(f1);
}


/**********************************************************************************************
 * Compute upper bounds on the absolute values of the entries of the vector x_1(r0,delta),
 * and store these in the global variable   x1b.
 * The following should have been previously computed (global variables):
 *   xrdb  -- bounds on the entries of x(r0+delta).
 *   TIb   -- bounds on the entries of T(r_0)^{-1}.
 *   wdbA, wdbB.
 **********************************************************************************************/
void x1bound() {
  int m,k;
  mpfr_t s,f,bmTx[Mmax];

  mpfr_init(s); mpfr_init(f);
  for(m=0; m<M-1; m++) mpfr_init(bmTx[m]);

  /* We use formula (9.1). First store in bmTx upper bounds on the absolute 
     values of the entries of   b1(r0,delta)-T1(r0,delta)x(r0+delta).       */
  mpfr_set(s,wdbA[0][0],GMP_RNDU); 
  mpfr_mul_ui(s,s,2,GMP_RNDU);
  for(m=0; m<M-1; m++) {
    mpfr_mul_ui(f,wdbA[0][m+1],2,GMP_RNDU);
    mpfr_mul(f,f,xrdb[m],GMP_RNDU);
    mpfr_add(s,s,f,GMP_RNDU);
  }
  for(k=0; k<M-1; k++) mpfr_set(bmTx[k],s,GMP_RNDU); 
  for(k=0; k<M-1; k++) {
    mpfr_mul(f,wdbB[0][k+1],xrdb[k],GMP_RNDU);
    mpfr_add(bmTx[k],bmTx[k],f,GMP_RNDU);
  }
  /* Next multiply with our bounds on T(r0)^(-1). */
  for(k=0; k<M-1; k++) {
    mpfr_set_ui(x1b[k],0,GMP_RNDU);
    for(m=0; m<M-1; m++) {
      mpfr_mul(f,TIb[k][m],bmTx[m],GMP_RNDU);
      mpfr_add(x1b[k],x1b[k],f,GMP_RNDU);
    }
  }
  /* Clear variables. */
  mpfr_clear(s); mpfr_clear(f);
  for(m=0; m<M-1; m++) mpfr_clear(bmTx[m]);
}


/**********************************************************************************************
 * Compute upper bounds on the absolute values of the entries of the vector x_2(r0,delta),
 * and store these in the global variable   x2b.
 * The following should have been previously computed (global variables):
 *   xrdb  -- bounds on the entries of x(r0+delta).
 *   x1b   -- bounds on the entries of x_1(r0,delta).
 *   TIb   -- bounds on the entries of T(r_0)^{-1}.
 *   wdbA, wdbB.
 **********************************************************************************************/
void x2bound() {
  int m,k;
  mpfr_t s,f,bmTx[Mmax];

  mpfr_init(s); mpfr_init(f);
  for(m=0; m<M-1; m++) mpfr_init(bmTx[m]);

  /* We use formula (9.2). First store in bmTx upper bounds on the absolute values of 
     the entries of   b2(r0,delta)-T'(r0)x_1(r0,delta)-T_2(r0,delta)x(r0+delta).       */
  mpfr_set(s,wdbA[1][0],GMP_RNDU); 
  for(m=0; m<M-1; m++) {
    mpfr_mul_ui(f,wdbA[0][m+1],2,GMP_RNDU);
    mpfr_mul(f,f,x1b[m],GMP_RNDU);
    mpfr_add(s,s,f,GMP_RNDU);
    mpfr_mul(f,wdbA[1][m+1],xrdb[m],GMP_RNDU);
    mpfr_add(s,s,f,GMP_RNDU);
  }
  for(k=0; k<M-1; k++) mpfr_set(bmTx[k],s,GMP_RNDU); 
  for(k=0; k<M-1; k++) {
    mpfr_mul(f,wdbB[0][k+1],x1b[k],GMP_RNDU);
    mpfr_add(bmTx[k],bmTx[k],f,GMP_RNDU);
    mpfr_mul_d(f,wdbB[1][k+1],0.5,GMP_RNDU);
    mpfr_mul(f,f,xrdb[k],GMP_RNDU);
    mpfr_add(bmTx[k],bmTx[k],f,GMP_RNDU);
  }
  /* Next multiply with our bounds on T(r0)^(-1). */
  for(k=0; k<M-1; k++) {
    mpfr_set_ui(x2b[k],0,GMP_RNDU);
    for(m=0; m<M-1; m++) {
      mpfr_mul(f,TIb[k][m],bmTx[m],GMP_RNDU);
      mpfr_add(x2b[k],x2b[k],f,GMP_RNDU);
    }
  }
  /* Clear variables. */
  mpfr_clear(s); mpfr_clear(f);
  for(m=0; m<M-1; m++) mpfr_clear(bmTx[m]);
}


/**********************************************************************************************
 * Compute an upper bound on the absolute values of c_2(r0,delta).
 * The following should have been previously computed (global variables):
 *   xrdb  -- bounds on the entries of x(r0+delta).
 *   x1b   -- bounds on the entries of x_1(r0,delta).
 *   x2b   -- bounds on the entries of x_2(r0,delta).
 *   TIb   -- bounds on the entries of T(r0)^{-1}.
 *   vb    -- bounds on the entries of v(r0).
 *   wdbA, wdbB.
 **********************************************************************************************/
void c2bound(mpfr_t c2b) {
  mpfr_t f;
  int m;
  mpfr_init(f);
  /* We use formula (9.5). First term to bound: w_2(r0,delta). */
  mpfr_mul_d(f,wdbB[1][0],0.5,GMP_RNDU);
  mpfr_add(c2b,f,wdbA[1][0],GMP_RNDU);
  /* Next add bound on v(r0)x_2(r0,delta). */
  for(m=0; m<M-1; m++) {
    mpfr_mul(f,vb[m],x2b[m],GMP_RNDU);
    mpfr_add(c2b,c2b,f,GMP_RNDU);
  }
  /* Next add bound on v'(r0)x_1(r0,delta). */
  for(m=0; m<M-1; m++) {
    mpfr_mul_ui(f,wdbA[0][m+1],2,GMP_RNDU);
    mpfr_mul(f,f,x1b[m],GMP_RNDU);
    mpfr_add(c2b,c2b,f,GMP_RNDU);
  }
  /* Finally add bound on v_2(r0,delta)x(r0+delta). */
  for(m=0; m<M-1; m++) {
    mpfr_mul(f,wdbA[0][m+1],xrdb[m],GMP_RNDU);
    mpfr_add(c2b,c2b,f,GMP_RNDU);
  }
  mpfr_clear(f);
}


void initglobalzoomvariables() {
  int j,k,m;
  
  /* Init T,v,b,w. */
  for(k=0; k<M-1; k++) {
    for(m=0; m<M-1; m++) mpfi_init(T[k][m]);
    mpfi_init(b[k]); 
    mpfi_init(v[k]); 
  }
  mpfi_init(w);
  /* Init TI, TIb */
  for(k=0; k<M-1; k++) for(j=0; j<M-1; j++) {
      mpfi_init(TI[j][k]); mpfr_init(TIb[j][k]);
    }
  /* Init x,vb,xrdb,x1b,x2b */
  for(k=0; k<M-1; k++) {
    mpfi_init(x[k]);
    mpfr_init(vb[k]);
    mpfr_init(xrdb[k]);
    mpfr_init(x1b[k]);
    mpfr_init(x1bb[k]);
    mpfr_init(x2b[k]);
  }
  /* Init wdbA,wdbB */
  for(j=0; j<2; j++) 
    for(m=0; m<M; m++) {
      mpfr_init(wdbA[j][m]); mpfr_init(wdbB[j][m]);
    }
  /* Init wdiag */
  for(m=0; m<M; m++) mpfi_init(wdiag[m]);
}

void clearglobalzoomvariables() {
  int j,k,m;
  
  /* Clear T,v,b,w. */
  for(k=0; k<M-1; k++) {
    for(m=0; m<M-1; m++) mpfi_clear(T[k][m]);
    mpfi_clear(b[k]);
    mpfi_clear(v[k]);
  }
  mpfi_clear(w);
  /* Clear TI, TIb */
  for(k=0; k<M-1; k++) for(j=0; j<M-1; j++) {
      mpfi_clear(TI[j][k]); mpfr_clear(TIb[j][k]);
    }
  /* Clear x,vb,xrdb,x1b,x2b */
  for(k=0; k<M-1; k++) {
    mpfi_clear(x[k]);
    mpfr_clear(vb[k]);
    mpfr_clear(xrdb[k]);
    mpfr_clear(x1b[k]);
    mpfr_clear(x1bb[k]);
    mpfr_clear(x2b[k]);
  }
  /* Clear wdbA,wdbB */
  for(j=0; j<2; j++) 
    for(m=0; m<M; m++) {
      mpfr_clear(wdbA[j][m]); mpfr_clear(wdbB[j][m]);
    }
  /* Clear wdiag */
  for(m=0; m<M; m++) mpfi_clear(wdiag[m]);
}



/**********************************************************************************************
 * Perform one step of the zooming. The input data r0, rtilde, eps are as in the manual.
 * We assume (without checking here) that  rtilde  is <= r0+eps.
 * The function returns an interval  rop  in which the eigenvalue must provably lie.
 * Note that the global variables M,N,Y and also the default precision must be chosen 
 * appropriately before calling this function.
 *********************************************************************************************/
void zoomstep(mpfi_t rop,mpfr_t r0,mpfr_t rtilde,mpfr_t eps) {
  mpfi_t cr0,crtilde,crstar,t1,t2,deltatilde,cpr0,deltastar;
  mpfr_t c2rstar,c2rtilde;
  int k,m;

  initglobalzoomvariables();
  mpfi_init(cr0); mpfi_init(crtilde); mpfi_init(crstar);
  mpfi_init(t1); mpfi_init(t2); mpfi_init(deltatilde); mpfi_init(cpr0);
  mpfi_init(deltastar); 
  mpfr_init(c2rstar);  mpfr_init(c2rtilde); 

  wdbinit(r0,eps);
  initsystem(r0,r0,eps);
//  fprintf(stderr,"initsystem done; now invert the matrix.\n");
  invertT(M-1);
//  fprintf(stderr,"invert done.\n");

  /* Some checks related to proving convergence. */
  checkfrobeniusnorm(M-1);

  computexc(cr0);
//  fprintf(stderr,"cr0:\n"); mpfi_myprint(cr0);
  /* Save x to xx; will use at the end. */
  for(k=0; k<M-1; k++) mpfi_set(xx[k],x[k]);
  /* Set TIb to give upper bounds on the absolute values of the entries of TI. */
  for(k=0; k<M-1; k++) for(m=0; m<M-1; m++) mpfi_mag(TIb[k][m],TI[k][m]);
  /* Set vb to give upper bounds on the absolute values of the entries of v. */
  for(k=0; k<M-1; k++) mpfi_mag(vb[k],v[k]);
  /* Set xrdb from the Kim-Sarnak bound. */
  for(m=0; m<M-1; m++) KimSarnakbound(xrdb[m],m+2);
  x1bound();
  /* Save x1b to x1bb; will use at the end. */
  for(k=0; k<M-1; k++) mpfr_set(x1bb[k],x1b[k],GMP_RNDU);
  x2bound();
  c2bound(c2rstar);
  initsystem(rtilde,r0,eps);
  invertT(M-1);

  computexc(crtilde);
  /* Set xrdb from our computed vector x. */
  for(m=0; m<M-1; m++) mpfi_mag(xrdb[m],x[m]);
  x1bound();
  x2bound();
  c2bound(c2rtilde);

//  fprintf(stderr,"cr0:\n"); mpfi_myprint(cr0);
//  fprintf(stderr,"crtilde:\n"); mpfi_myprint(crtilde);

  /* Set deltatilde=rtilde-r0 */
  mpfi_init_set_fr(deltatilde,rtilde);
  mpfi_sub_fr(deltatilde,deltatilde,r0);
  /* Compute cpr0=c'(r0). */
  mpfi_sub(t1,crtilde,cr0);
  mpfi_div(t1,t1,deltatilde);
//  fprintf(stderr,"Diff quot for c'(r): "); mpfr_myprint1(&t1->right,5);
//  fprintf(stderr," c2bounds: "); mpfr_myprint1(c2rtilde,5); fprintf(stderr," "); mpfr_myprint1(c2rstar,5);
  mpfi_div_fr(t2,t1,c2rtilde);
//  fprintf(stderr," ratio: "); mpfr_myprint1(&t2->left,5); fprintf(stderr,"\n");

//  mpfi_myprint(t1);

  mpfi_set_ui(t2,0);
  mpfi_increase(t2,c2rtilde);
  mpfi_mul(t2,deltatilde,t2);
  mpfi_sub(cpr0,t1,t2);

//  fprintf(stderr,"Interval for c'(r):\n");
//  mpfi_myprint(cpr0);

  /* Compute c(r^*) */
  crstarbound(crstar,r0,eps);
  /* Compute delta^* */
  mpfi_sub(t1,crstar,cr0);
  mpfi_set_ui(t2,0);
  mpfi_increase(t2,c2rstar);
  mpfi_mul_fr(t2,t2,eps);
  mpfi_add(t2,cpr0,t2);
//  fprintf(stderr,"Interval for c'(r) plus [eps times c_2(r_0,delta^*)-bound]:\n");
//  mpfi_myprint(t2);

  mpfi_div(deltastar,t1,t2);
  /* Compute rop */
  mpfi_add_fr(rop,deltastar,r0);
  /* Compute rigorous bounds on the solution vector at r=r^*. */
  for(m=0; m<M-1; m++) {
    mpfi_set_ui(t1,0);
    mpfi_increase(t1,x1bb[m]);
    mpfi_mul(t1,t1,deltastar);
    mpfi_add(xx[m],xx[m],t1);
  }

  mpfi_clear(cr0); mpfi_clear(crtilde); mpfi_clear(crstar);
  mpfi_clear(t1); mpfi_clear(t2);  mpfi_clear(deltatilde); mpfi_clear(cpr0);
  mpfi_clear(deltastar);
  mpfr_clear(c2rstar); mpfr_clear(c2rtilde); 

  clearglobalzoomvariables();
}




/********************************************************************************************** 
 * Set the global variables M and N adapted to computing in precision pr.
 * (using r0, Y).
 **********************************************************************************************/
void findMN(mpfr_t r0,int pr) {
  double R,YY;
  mpfr_t m,Y0;

  mpfr_init(m); mpfr_init(Y0);
  R=mpfr_get_d(r0,GMP_RNDD);  
  for(M=1; 2*3.141592653589793238462643383*(M+1)*0.5*sqrt(3.0)<=R; M++) ;
  mpfr_set_d(Y0,0.5*sqrt(3.0),GMP_RNDD);
  while(1) {
    MMajorant_basic(m,M,Y0,r0);   /* Previously we used here a floating point version of MMajorant;
				     however this fails when the majorant is e.g. < 10^(-700)!   */
    if(mpfr_get_exp(m)<-pr) break;
    M++;
  }
  N=M+1;
  YY=mpfr_get_d(Y,GMP_RNDD);
  while(2*3.141592653589793238462643383*(N+1)*YY<=R) N++;
  while(1) {
    MMajorant_basic(m,N,Y,r0);
    if(mpfr_get_exp(m)<-pr) break;
    N++;
  }
  mpfr_clear(m); mpfr_clear(Y0);
}


/* sm and sr are strings giving the midpoint and radius of an interval I;
   this routine initializes m and r in such a way that I is contained in
   [m-r,m+r].
*/

void mpfr_strstomidandradius(mpfr_t m,mpfr_t r,char *sm,char *sr) {
  char *a;  
  mpfr_t mu,x;
  mpfr_strtofr(m,sm,&a,10,GMP_RNDU);
  mpfr_init_set(mu,m,GMP_RNDU);
  mpfr_strtofr(m,sm,&a,10,GMP_RNDD);
  mpfr_init(x);
  mpfr_sub(x,mu,m,GMP_RNDU);
  mpfr_strtofr(r,sr,&a,10,GMP_RNDU);
  mpfr_add(r,r,x,GMP_RNDU);
}

int main(int argc,char *argv[]) {
  mpfr_t r0exact,r0,eps,eps1,rtilde,z,a1,d;
  mpfi_t newr;
  mpfi_t rr,a,b,y;
  char s[100],s1[100];
  int pr,pr1,endpr,extrapr,r0exactpr,j,k,m,sp,eig_ok;
  char *aa;

  if (argc!=5) {
    printf("Usage: hejhalzoom oe (0 or 1) r0 eps endpr (in base 2)\n");
    return 0;
  }
  sscanf(argv[1],"%d",&oe);
  sscanf(argv[4],"%d",&endpr);
  mpfr_set_default_prec(50);
  /* Read r0 to "full precision". */
  r0exactpr=(int)(strlen(argv[2])*3.33)+100;
  mpfr_init2(r0exact,r0exactpr);
  mpfr_strtofr(r0exact,argv[2],&aa,10,GMP_RNDN);
//  fprintf(stderr,"r0exactpr=%d\n",r0exactpr);
  /* The following is just to read eps crudely. Float doesn't work e.g. if eps=1e-900. */
  mpfr_init(eps); mpfr_init(d);
  mpfr_strtofr(eps,argv[3],&aa,10,GMP_RNDU);
  pr=-(int)mpfr_get_exp(eps);   /* This is roughly the number of known binary digits. */
  pr=2*pr+15;   /* Add some extra precision in this first iteration; this is found 
		   to substantially improve the worst allowable eps input. */
  if (pr>endpr+extrapr) pr=endpr+extrapr;
  mpfr_clear(eps); mpfr_clear(d);
  /* Now read again, in the appropriate precision. */
  mpfr_set_default_prec(pr);
  mpfr_init(r0); mpfr_init(eps);
  mpfr_strstomidandradius(r0,eps,argv[2],argv[3]);  /* *** CHANGE! ... */
  eig_ok=0;

  while(1) {
    mpfr_init(Y); mpfr_init(rtilde); mpfi_init(newr); mpfr_init(eps1);
    mpfr_set_d(Y,0.86,GMP_RNDN);
    findMN(r0,pr);
//    fprintf(stderr,"Prec=%d M=%d N=%d\n",pr,M,N);
    /* Set eps1=2^(-pr/2) and check whether eps<eps1. */
    mpfr_set_d(eps1,-0.5*pr,GMP_RNDN);
    mpfr_exp2(eps1,eps1,GMP_RNDN);
    if(mpfr_cmp(eps,eps1)<0) {
//      fprintf(stderr,"Special case...\n");
      mpfr_sub(rtilde,r0,eps1,GMP_RNDD);
    } else
      mpfr_sub(rtilde,r0,eps,GMP_RNDD);
    /* Since we SUBTRACT above we are sure that  rtilde is  < r0 and hence <r0+eps. */
    for(m=0; m<M-1; m++) mpfi_init(xx[m]);
    zoomstep(newr,r0,rtilde,eps);
    /* Now prepare for new iteration. New precision? */
    mpfr_init(d);
    mpfi_diam_abs(d,newr);
    pr1=-(int)mpfr_get_exp(d);  /* This is roughly the number of binary digits we now have. */
//    fprintf(stderr,"'Wasted precision': %d\n",pr-pr1);
    extrapr=pr-pr1+10;  /* Seems a good guess on extra precision needed. */
    pr=pr1;
    mpfr_clear(d);
    if (pr>endpr) {
      /* Repeat one final time, unless eig_ok==1 or r0 is already inside newr. */
      if(mpfi_is_inside_fr(r0,newr)>0 || eig_ok) {
printf("%d\n%d\n%d\n",endpr,oe,M);
	print_endpoints(newr);
	for(m=0;m<M-1;m++)
	  print_endpoints(xx[m]);
	return 0;
      } else {
	eig_ok=1;
//	fprintf(stderr,"Almost done; proved eigenvalue to be:\n");
//	mpfi_print1sloppy(newr);
//	fprintf(stderr,"which is within desired range. We iterate one last time with r0 in\nthis interval to get good estimates on the Fourier coefficients.\n");
      }
    }
    mpfr_clear(Y); mpfr_clear(rtilde); mpfr_clear(r0); mpfr_clear(eps); mpfr_clear(eps1);
    for(m=0; m<M-1; m++) mpfi_clear(xx[m]);
    pr=2*pr+5;
    if (pr>endpr+extrapr) pr=endpr+extrapr;
    mpfr_set_default_prec(pr);
    mpfr_init(r0); mpfr_init(eps); mpfr_init(eps1);
    mpfi_mid(r0,newr);
    mpfr_sub(eps,&newr->right,r0,GMP_RNDU);
    mpfr_sub(eps1,r0,&newr->left,GMP_RNDU);
    mpfr_max(eps,eps,eps1,GMP_RNDU);
    /* If r0exact is well within the endpoints of newr then use it instead! */
    mpfr_sub(eps1,r0,r0exact,GMP_RNDN);
    mpfr_abs(eps1,eps1,GMP_RNDU);
    mpfr_mul_ui(eps1,eps1,2,GMP_RNDU);
    if(mpfr_cmp(eps1,eps)<0) {
//      fprintf(stderr,"We use the r0 given by user instead of midpoint.\n");
      mpfr_set(r0,r0exact,GMP_RNDN);
    }
    mpfr_clear(eps1);
//    fprintf(stderr,"New r0:"); mpfr_myprint(r0);
//    fprintf(stderr,"\n");
//    fprintf(stderr," eps:"); mpfr_out_str(stderr,10,5,eps,GMP_RNDU); fprintf(stderr,"\n");
    mpfi_clear(newr);
  }
  return(0);
}

// This fft is for \sum_{j=0}^{N-1}f(j)cos(2*Pi*j*m/2N) or \sum_{j=0}^{N-1}f(j)sin(2*Pi*j*m/2N)
// Require N>2

  /* Prepare for fft */
init_fft(mpfi_t *s,mpfi_t *c,int n){
  int i;
  mpfi_t t,p;

  mpfi_init(t); mpfi_init(p);
  mpfi_const_pi(p); mpfi_div_ui(p,p,n);
  for (i=0; i<n; i++){
    mpfi_init(s[i]); mpfi_init(c[i]);
    mpfi_mul_si(t,p,i);
    mpfi_sin(s[i],t); mpfi_cos(c[i],t);
  }
  mpfi_clear(t);
  mpfi_clear(p);
}


fft(mpfi_t *in,mpfi_t *os,mpfi_t *oc,mpfi_t *s,mpfi_t *c, int n){
  int i,j,l,k1,k2,b;
  mpfi_t p,t,h;

  mpfi_init(p); mpfi_init(t); mpfi_init(h);
  for ( i=0; i<(2*n); i+=2 ){
    // k=bit_reverse(i).
    for ( k1=0,j=1; j<=n; j<<=1 ){
      k1<<=1; if ( i&j ) k1|=1; }
//printf(" k=%ld ",k1);
    // set up blocks. 
    mpfi_init(oc[i]); mpfi_init(oc[i+1]);
    mpfi_set(oc[i],in[k1]); mpfi_set(oc[i+1],in[k1]);
    mpfi_init(os[i]); mpfi_init(os[i+1]);
    mpfi_set_ui(os[i+1],0); mpfi_set_ui(os[i],0);
  }

  for ( i=2, b=n/2; i<=n; i<<=1, b>>=1 ){
    for ( j=0; j<b; j++ ){
      k1=i*j*2; k2=k1+i;
      for ( l=0; l<n; l+=b ){
        mpfi_mul(t,oc[k2],c[l]); mpfi_mul(p,os[k2],s[l]);
        mpfi_sub(t,t,p);
        mpfi_mul(p,oc[k2],s[l]); mpfi_mul(h,os[k2],c[l]);
        mpfi_add(p,p,h);
        mpfi_sub(oc[k2],oc[k1],t); mpfi_sub(os[k2],os[k1],p);
        mpfi_add(oc[k1],oc[k1],t); mpfi_add(os[k1],os[k1],p);
        k2++; k1++;
      }
    }
  }
  mpfi_clear(p);  mpfi_clear(h);  mpfi_clear(t);

}



complete_fft(mpfi_t *h_r,int n,int omega){
  int i,j;
  mpfi_t *fftc,*ffts,*ts,*tc,t,p,u;

  ffts=(mpfi_t *)malloc(2*n*sizeof(mpfi_t));
  fftc=(mpfi_t *)malloc(2*n*sizeof(mpfi_t));
  ts=(mpfi_t *)malloc(n*sizeof(mpfi_t));
  tc=(mpfi_t *)malloc(n*sizeof(mpfi_t));

  init_fft(ts,tc,n); 

  fft(h_r,ffts,fftc,ts,tc,n);
//for(i=0;i<2*n;i++){printf("TS[%ld]=",i);mpfi_out_str(stdout,10,15,ffts[i]);printf("\n");}

  mpfi_init(t); mpfi_init(p); mpfi_init(u);
  mpfi_const_pi(t); mpfi_div_ui(t,t,2*n);  //t=Pi/(2N).

  for ( i=1; i<=n; i++ ){
    mpfi_mul_si(u,t,i);
    mpfi_sin(p,u); mpfi_cos(u,u); //u=cos(Pi*i/2N) p=sin(Pi*i/2N);
    if ( omega ){
      // odd case
      mpfi_mul(p,fftc[i],p); mpfi_mul(u,ffts[i],u);
      mpfi_add((h_r)[i-1],u,p);
    }
    else {
      // even case
      mpfi_mul(u,fftc[i],u); mpfi_mul(p,ffts[i],p);
      mpfi_sub((h_r)[i-1],u,p);
    } 
  }
//for(i=0;i<n;i++){printf("h_r[%d]=",i);mpfi_out_str(stdout,10,15,h_r[i]);printf("\n");}
  for ( i=0; i<n; i++ ){
    mpfi_clear(fftc[i]); mpfi_clear(fftc[i+n]);
    mpfi_clear(ffts[i]); mpfi_clear(ffts[i+n]);
    mpfi_clear(tc[i]); mpfi_clear(ts[i]);
  }
  mpfi_clear(p);  mpfi_clear(u);  mpfi_clear(t);
  free(fftc); free(ffts); free(tc); free(ts);
}

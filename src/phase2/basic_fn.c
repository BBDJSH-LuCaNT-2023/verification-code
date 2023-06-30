/*************************************************
 F01: Truncation Error boundary M_{ir}(y,d)
*************************************************/ 
int error_M(mpfi_t *f,mpfi_t r,mpfi_t y,mpfi_t d){
  mpfi_t t1,t2;
  mpfi_init2(t1,mpfi_get_prec(r)); mpfi_init2(t2,mpfi_get_prec(r));

  mpfi_sqr(t1,y); mpfi_sqr(t2,r); mpfi_sub(t1,t1,t2); mpfi_sqrt(t1,t1);	// t1=sqrt(y^2-r^2);
  mpfi_div(t2,t1,r); mpfi_atan(t2,t2); mpfi_mul(t2,t2,r);	// t2=exp(r*(arctg(t1/r)-t1));
  mpfi_sub(t2,t2,t1); mpfi_exp(t2,t2);
  mpfi_div(t1,y,t1); mpfi_div(*f,t1,d); mpfi_add_ui(*f,*f,1); mpfi_mul(t2,t2,*f);	// t1=y/sqrt(y^2-r^2); t2*=(1+t1/d);
  mpfi_sqrt(t1,t1); mpfi_mul_d(t1,t1,2.203125); mpfi_mul(*f,t1,t2);	// t1=C*sqrt(t1); f*=t2*t1;

  mpfi_clear(t1);  mpfi_clear(t2);
  return 0;
}


/***************************************************************
 F02: if e=-1, 2^(-d)<x_left, return d;
      if e=0, 2^(-d)>x_right-x_left, return d;
      if e=1, 2^(-d)>x_right, return d;
      if e>1, x is an array, for 0<=i<abs(e), return the lasgest mpfi_size2(x[i],0);
      if e<-1, x is an array,for 0<=i<abs(e), return the smallest mpfi_size2(x[i],0);
***************************************************************/
int mpfi_size2(mpfi_t *x, int e){
  int i,j,k,p=0,q=99999;
  mpfi_t y; mpfi_init2(y,mpfi_get_prec(*x));
  mpfr_t t; mpfr_init2(t,mpfi_get_prec(*x));

  if ( abs(e)>1 )
    for ( i=0; i<abs(e); i++ ){
      mpfi_abs(y,x[i]);
      j=mpfi_size2(&y,0);
      if ( e>0 ){ if ( j>p ){ k=i; p=j; }}
      else { if ( j<q && j>0 ) { k=i; q=j; }}
    }
  else{
    mpfi_abs(y,*x);
    if ( e==0 ){ 
      mpfi_diam_abs(t,y); mpfr_log2(t,t,MPFR_RNDZ); 
//      j=mpfi_size2(&y,1); //if ( j>0 ) j=0;
      k=-(int)mpfr_get_si(t,MPFR_RNDZ);
    }
    else if ( e==1 ){ mpfi_get_right(t,y); mpfr_log2(t,t,MPFR_RNDU); k=-(int)mpfr_get_si(t,MPFR_RNDU); }
    else{ mpfi_get_left(t,y); mpfr_log2(t,t,MPFR_RNDD); k=-(int)mpfr_get_si(t,MPFR_RNDD); }
  }
  mpfi_clear(y);  mpfr_clear(t);
  return(k);
}


/***************************************************************
 F02A: e=1, right_endpoint = left_endpoint
       e=0, left_endpoint = right_endpoint = (left_endpoint+right_endpoint)/2
       e=-1,left_endpoiont = right_endpoint 
***************************************************************/
int mpfi_endpoint(mpfi_t x,int e){
  int i=0;  mpfr_t t; mpfr_init2(t,mpfi_get_prec(x));
  switch ( e ){
    case -1: mpfi_get_left(t,x); mpfi_set_fr(x,t); break;
    case 0: mpfi_mid(t,x); mpfi_set_fr(x,t); break;
    case 1: mpfi_get_right(t,x); mpfi_set_fr(x,t); break;
    default: i=-1;
  }
  mpfr_clear(t);
  return(i);
}


/*****************************************************************************
 F02B: compute log(v)
       e=2: return smallest j, v<2^j;
       e=1: return smallest j, v<=2^j;
       e=-1: return largest j, 2^j<=v;
       e=-2: return largest j, 2^j<v;
*****************************************************************************/
int int_size2(int v, int e){
  int i=1,j=0;
  if ( v<1 ) return(-1);
  else  while( i<v ){ i*=2; j++; }
  switch ( e ){
    case 2: if ( i==v ) j++; break;
    case -1: if ( i>v ) j--; break;
    case -2: j--; break;
    default: break;
  }
  return(j);
}


/*****************************************************************************
 F03: compute g(y,j)=a[j]*sqrt{y*j}*K_{ir}(j*y)*cos/sin(jx)
      If oe=-1, record = sqrt(y*n)*K_{ir}(y*n), return = mpfi_size2(record,-1$
      If oe=-2, record = min{ sqrt(y*j)*K_{ir}(y*j) }, for 0<j<=n. return = i$
      If oe=-3, record[i] = sqrt(y*j)*K_{ir}(y*j), for 0<j<=n. return = index$
      If oe=-4, record = min{ sqrt(y*j)*K_{ir}(y*j) }, for prime j and 0<j<=n$
      If oe=0, record = sum of g(y,j) for 1<=j<=n, using cosine; return max-i$
      If oe=1, record = sum of g(y,j) for 1<=j<=n, using sine; return max-int$
                return index j of max interval g(y,j) occured;
*****************************************************************************/
int Cn_sqrt_bessel(mpfi_t *record,mpfi_t *a,mpfi_t y,mpfi_t R,int n,mpfi_t x,int oe){
#if 1
  int i=0,j=0,k=99999,m,*p, l;
  mpfi_t t,u,v; l=mpfi_get_prec(y); mpfi_init2(t,l); mpfi_init2(u,l); mpfi_init2(v,l);
  mpfi_ptr f; f=init_kbessel(mpfr_get_default_prec());

  p=(int *)malloc(n*sizeof(int));
  if ( oe==-3 || oe==-2 || oe==-4 ){
    mpfi_set_si(u,1024);
    prime_decomposition(p,n);
    for ( i=0; i<n; i++ )
      if ( oe+4 || p[i]<=0 ){
        mpfi_mul_si(t,y,i+1); kbessel(f,R,t); mpfi_sqrt(v,t); mpfi_mul(t,v,f); mpfi_abs(v,t);
        if ( oe==-3 ) mpfi_set(record[i],t);
        else if ( mpfi_cmp(v,u)<0 ){ j=i; mpfi_set(u,v); }
      }
    if ( oe+3 ) mpfi_set(*record,u);
    free(p);
  }
  else if ( oe==-1 ){
    mpfi_mul_si(t,y,n); kbessel(f,R,t); mpfi_sqrt(v,t);
    mpfi_mul(*record,v,f); j=mpfi_size2(record,-1);
  }
  else if ( oe==0||oe==1 ){
    mpfi_set_si(*record,0);
    for ( i=0; i<n; i++ ){
      mpfi_mul_si(t,y,i+1); kbessel(f,R,t); mpfi_mul(f,a[i],f);
      mpfi_sqrt(v,t); mpfi_mul(v,v,f); mpfi_mul_si(t,x,i+1);
      if ( oe ) mpfi_sin(t,t);
      else  mpfi_cos(t,t);
      mpfi_mul(v,t,v); mpfi_add(*record,*record,v);
      if ( (m=mpfi_size2(&v,0))<k ){ k=m; j=i; }
    }
  }
  else printf("Incorrect value of oe.\n");
  clear_kbessel(f); mpfi_clear(t);mpfi_clear(v); mpfi_clear(u);
  return j;
}
#else
  int i=0,j=0,k=0,m;
  unsigned long *p=0, l;
  mpfi_t t,u,v; mpfi_init(t); mpfi_init(u); mpfi_init(v);
  mpfi_ptr f; f=init_kbessel(mpfr_get_default_prec());

  if ( mpfi_cmp_si(y,0)<0 ) return(-1);
  if ( oe==-3 || oe==-2 ){
    mpfi_set_si(u,1024);
    mpfi_sqrt(v,y);
    for ( i=0; i<n; i++ ){
      mpfi_mul_si(t,y,i+1); kbessel(f,R,t); mpfi_mul(t,v,f);
      if ( oe==-3 ) mpfi_set(record[i],t);
      mpfi_abs(t,t);
      if ( mpfi_cmp(t,u)<0 ){ j=i; mpfi_set(u,t); }
    }
    if ( oe==-2 ) mpfi_set(*record,u);
  }
  else if ( oe==-1 ){
    mpfi_mul_si(t,y,n); kbessel(f,R,t); mpfi_sqrt(v,y);
    mpfi_mul(*record,v,f); j=mpfi_size2(record,-1);
  }
  else if ( oe==0||oe==1 ){
    mpfi_set_si(*record,0);
    mpfi_sqrt(v,y);
    for ( i=0; i<n; i++ ){
      mpfi_mul_si(t,y,i+1); kbessel(f,R,t); mpfi_mul(f,a[i],f);
      mpfi_mul(f,v,f); mpfi_mul_si(t,x,i+1);
      if ( oe ) mpfi_sin(t,t);
      else  mpfi_cos(t,t);
      mpfi_mul(t,t,f); mpfi_add(*record,*record,t);
      if ( (m=mpfi_size2(&t,0))>k ){ k=m; j=i; }
    }
  }
  else printf("Incorrect value of oe.\n");
  clear_kbessel(f); mpfi_clear(t);mpfi_clear(v); mpfi_clear(u);
  return j;
}
#endif



/********************************************************************
 F04: The interval of Box[i]=interval of y_s, which need i terms in inner sum;
      Return the max number of terms needed for this R;
*********************************************************************/
int compute_Box(mpfi_t **box,mpfi_t R,mpfi_t end,int CN){
  int i,k,j,N=abs(CN); i=(int)mpfi_get_prec(R);
  mpfi_t con,f,y0;  mpfi_init2(f,i); mpfi_init2(y0,i); mpfi_init2(con,i);

  mpfi_set_si(f,3); mpfi_sqrt(f,f); mpfi_const_pi(con); mpfi_mul(con,con,f); //con=sqrt(3)*Pi

  k=N+1;
  do {	 /***** How many coeff needed for y=\sqrt(3)Pi *****/
    k--;
    mpfi_mul_si(y0,con,k+1);  error_M(&f,R,y0,con);
  } while ( mpfi_cmp(f,end)<0 );
  i=k+1;
  if ( k==N ){ mpfi_clear(y0); mpfi_clear(con); mpfi_clear(f); return(-1); }
  else if ( CN<0 ){ mpfi_clear(y0); mpfi_clear(con); mpfi_clear(f); return(i); }
  (*box)=(mpfi_t *)malloc(i*sizeof(mpfi_t));

  mpfr_t z; mpfr_init(z);
  mpfi_t s,t,e;  mpfi_init(s);mpfi_init(t); mpfi_init(e);

  /********** Find out every Box interval **************/
  mpfi_set_si(e,-mpfr_get_default_prec()); mpfi_exp2(e,e);	//e=2^(-default_prec())
  mpfi_get_left(z,con);
  while ( k>0 ){
    j=1; mpfi_set_fr(t,z);
    do {	/********* compute upper bound *********/
      mpfi_add_si(t,t,j);
      mpfi_mul_si(y0,t,k);  error_M(&f,R,y0,t);
      j*=2;
    } while ( mpfi_cmp(f,end)>0 );
    mpfi_set_si(s,j-1); 
    while( mpfi_cmp(s,e)>0 ){	/***** compute the right endpoint *****/
      mpfi_div_si(s,s,2);
      if ( mpfi_cmp(f,end)<0 ) mpfi_sub(t,t,s);
      else  mpfi_add(t,t,s);
      mpfi_mul_si(y0,t,k); error_M(&f,R,y0,t);
    }

    mpfi_put_fr(t,z); mpfi_get_right(z,t);
    mpfi_init((*box)[k]);  mpfi_set((*box)[k],t);
    k--;
  }
  mpfi_set_fr(con,z);
  mpfi_interv_si(t,0,1); mpfi_inv(t,t); mpfi_get_right(z,t);    //z=\infty
  mpfi_put_fr(con,z);
  mpfi_init((*box)[k]); mpfi_set((*box)[k],con);
  mpfr_clear(z); mpfi_clear(s);mpfi_clear(t); mpfi_clear(con); mpfi_clear(e); mpfi_clear(f); mpfi_clear(y0);
  return(i);
}


/*******************************************************************
 F05: Compute which Box does the Y value need.
*******************************************************************/
int which_box(mpfi_t *box,mpfi_t Y,int n){
  int i,j,k;
  i=(n+1)/2; j=(n+3)/4;
  do {
    k=mpfi_cmp(Y,box[i]);
    if ( k<0 ) i+=j; else if ( k>0 ) i-=j;
    if ( i>=n ) i=n-1;
    if ( i<0 ) i=0;
    j=(j+1)/2;
  } while ( k ); 
  return i;
}



/*******************************************************************
 F06: Compute Y s.t.error_M(R,Y*(N-M),Y*N)/(sqrt(Y)*Kir(Y*M))<end
*******************************************************************/
int compute_Y(mpfi_t *re,mpfi_t *Y,mpfi_t R,int N,int M,int en,int *Rnd,int CN,int c){
  int i,e,ex,j,rnd,con=c+int_size2(N,1);
  mpfi_t t,u,v,y,s,w,end;
  i=(int)mpfi_get_prec(*Y);
  mpfi_init2(u,i); mpfi_init2(v,i); mpfi_init2(t,i); mpfi_init2(y,i); mpfi_init2(s,i); mpfi_init2(w,i); mpfi_init2(end,i);

  mpfi_set_si(end,-en); mpfi_exp2(end,end); mpfi_div_d(end,end,10.03125);
  if ( N<M*2 ) return(-2);

  mpfi_div_si(w,R,16*M);  mpfi_mul_si(y,w,15);
  i=*Rnd;
/*
printf("\tmax=%d,end=",i);
mpfi_out_str(stdout,10,5,end);printf("\n");
*/
  do{	/******* compute the upper bound ******/
    mpfi_add(y,y,w);
    mpfi_mul_si(u,y,N-M); mpfi_mul_si(v,y,N); error_M(re,R,u,v);
    e=Cn_sqrt_bessel(&t,&v,y,R,M,t,-1);
    mpfi_mul(u,t,end);
    rnd=int_size2(compute_Box(&Y,R,u,-CN),1)+con;
    if ( e<0 ) e=0;
/*
printf("\tcompute_Y0 e=%d rnd=%d Y=",e,rnd);
mpfi_out_str(stdout,10,5,y);printf(" re=");
mpfi_out_str(stdout,10,5,*re);printf(" t=");
mpfi_out_str(stdout,10,5,t);printf("\n");
*/
  } while ( mpfi_cmp(*re,u)>0 && e+en+rnd<=i);

  mpfi_mul_d(t,w,16.5);
  mpfi_set(*Y,y);
  if ( e+en+rnd > i ) e=-1;
  else{
    if ( mpfi_cmp(y,t)>0 ){	/***** compute the Y value *****/
      mpfi_set_si(s,-i); mpfi_exp2(s,s);
      while ( mpfi_cmp(w,s)>0 ){
        mpfi_div_si(w,w,2);
        if ( mpfi_cmp(*re,u)>0 ) mpfi_add(y,y,w);
        else mpfi_sub(y,y,w);
        mpfi_mul_si(u,y,N-M); mpfi_mul_si(v,y,N); error_M(re,R,u,v);
        e=Cn_sqrt_bessel(&t,&v,y,R,M,t,-1);
        mpfi_mul(u,end,t);      
      }
    }
/*
printf("\tcompute_Y1 e=%d rnd=%d Y=",e,rnd);
mpfi_out_str(stdout,10,5,y);printf(" re=");
mpfi_out_str(stdout,10,5,*re);printf(" t=");
mpfi_out_str(stdout,10,5,t);printf("\n");
*/
    mpfi_div_si(w,y,128); mpfi_sub(y,y,w);
    do{
      mpfi_add(y,y,w);
      mpfi_mul_si(u,y,N-M); mpfi_mul_si(v,y,N); error_M(re,R,u,v);
      j=Cn_sqrt_bessel(&t,&v,y,R,M,t,-4); e=mpfi_size2(&t,-1); if ( e<0 ) e=0;
      mpfi_mul(u,t,end); *Rnd=int_size2(compute_Box(&Y,R,u,-CN),1)+con;
      if ( j+1==M ){ ex=e; rnd=*Rnd; }     
      else{
        ex=Cn_sqrt_bessel(&s,&v,y,R,M,t,-1);
        mpfi_mul(v,s,end);  rnd=int_size2(compute_Box(&Y,R,v,-CN),1)+con;
      }
/*
printf("\tcompute_Y2 e=%d at %d Rnd=%d Y=",e,j,*Rnd);
mpfi_out_str(stdout,10,5,y);printf(" re=");
mpfi_out_str(stdout,10,5,*re);printf(" t=");
mpfi_out_str(stdout,10,5,t);printf("\n");

printf("\tcompute_Y2 ex=%d rnd=%d Y=",ex,rnd);
mpfi_out_str(stdout,10,5,y);printf(" re=");
mpfi_out_str(stdout,10,5,*re);printf(" t=");
mpfi_out_str(stdout,10,5,s);printf("\n");
*/
    } while ( (mpfi_cmp(*re,u)>0 || e+*Rnd+en>i)&&(ex+rnd+en<=i) );
    mpfi_set(*Y,y); mpfi_set(*re,u);
    if ( (ex+rnd+en)>i ) e=-1;
  }

  mpfi_clear(end); mpfi_clear(u); mpfi_clear(v); 
  mpfi_clear(s); mpfi_clear(t); mpfi_clear(y); mpfi_clear(w);
  return(e);
}


/**********************************************************
 F07: Compute the smallest N, for a given Y. 
      error_M(R,Y*(N-M),Y*N)/sqrt(Y)Kir(YM)<2^-en
      Return Pr, which is -log2(sqrt(Y)*Kir(YM))+EndPr
***********************************************************/
int compute_N(mpfi_t *re,mpfi_t Y,mpfi_t R,int *N,int M,int en){
  int e;
  mpfi_t u,v,w,end;
  mpfi_init(u); mpfi_init(v); mpfi_init(w);mpfi_init(end);

  mpfi_set_si(end,-en); mpfi_exp2(end,end); mpfi_div_d(end,end,10.03125);

  *N=4; while ( *N<M ) *N*=2;
  mpfi_set_si(w,1);
  e=Cn_sqrt_bessel(re,&w,Y,R,M,u,-1);
  e=mpfi_size2(re,-1);
  do{
    *N*=2;
    mpfi_mul_si(u,Y,*N-M); mpfi_mul_si(v,Y,*N); error_M(&w,R,u,v);
    mpfi_div(v,w,*re);
  } while ( mpfi_cmp(v,end)>0 );

  mpfi_clear(u); mpfi_clear(v); mpfi_clear(w); mpfi_clear(end);
  return(e);
}

//printf(" K="); mpfi_out_str(stdout,10,8,*re); printf("\n");

/*****************************************************
 F08: Transfer z=>z* by -1/z and abs(Re(z*))<0.5
      The x,x_s,y_s are scaled with 2Pi
*****************************************************/ 
compute_Points(mpfi_t *x_s, mpfi_t *y_s, mpfi_t Y, int N, int j){
  mpfi_t u,v,c,t,w;
  int i;  i=mpfi_get_prec(Y);
  mpfi_init2(t,i); mpfi_init2(v,i); mpfi_init2(c,i); mpfi_init2(u,i); mpfi_init2(w,i); 

  mpfi_const_pi(c); mpfi_mul_si(c,c,2);	/*** c=2*Pi ***/
  mpfi_set_d(u,0.5+(double)j); mpfi_div_si(u,u,N);
  mpfi_div(v,Y,c);  

  mpfi_sqr(t,u); mpfi_sqr(w,v); mpfi_add(t,t,w);

  while( mpfi_cmp_si(t,1)<0 ){
    mpfi_neg(u,u); mpfi_div(u,u,t);
    mpfi_add_d(w,u,0.5); mpfi_sub_si(u,u,floor(mpfi_get_d(w)));
    mpfi_div(v,v,t);
    mpfi_sqr(t,u); mpfi_sqr(w,v); mpfi_add(t,t,w);
  }
  mpfi_mul(*x_s,u,c); mpfi_mul(*y_s,v,c);
  mpfi_clear(u); mpfi_clear(c); mpfi_clear(v); mpfi_clear(w); mpfi_clear(t);
}


/********************************************************
 F09: prime decomposite integer upto x
      p_d[i]: 0--prime; -p--prime power; p^e--composite;
********************************************************/
int prime_decomposition(int *p_d,int x){
  int i=1,j,k,l,m,n,counter=0;

  l=floor((double)sqrt(x));
  for ( j=0; j<x; j++ ) p_d[j]=0;
  while ( i<x ){ i++;
    if ( !p_d[i-1] ){
      j=m=1; n=i;
      while ( m&&(n<=x)&&(i<l) ){ m=n*i; if((m<x)&&(m>0)){n=m;j++;} else m=0; }
      counter+=j;
      for ( ; j>0; j--,n/=i ){  
        p_d[n-1]=-1*i;   /* prime power case */
        if ( n<l ){
          m=floor((double)x/n);  m*=n;  if ((x-m)>=n) m+=n;
          for ( k=(i+1)*n; k<m+1; k+=n ) 
            if ( !p_d[k-1] ) p_d[k-1]=n; }  /* composite case */
      }
      p_d[0]++;  p_d[i-1]=0;   /* prime case */
    }  /****** end if ******/
  }
  return(counter);
}


/*********************************************
 F10: make the coeff. multiplicative 
*********************************************/
void coeff_mult(mpfi_t *C,int x){
  int i,m,n,*p_d;

  p_d=(int *)malloc(x*sizeof(int));
  prime_decomposition(p_d,x);  
  mpfi_set_si(*C,1);
  for ( i=3; i<x ; i++ ){
    m=p_d[i];
    if ( m>0 ){        //composite numbers.
      n=(i+1)/m-1; 
      mpfi_mul(C[i],C[m-1],C[n]);
    }
    else if ( m ){     //prime powers.
      m=-m; n=(i+1)/m;
      mpfi_mul(C[i],C[m-1],C[n-1]);
      n/=m;
      mpfi_sub(C[i],C[i],C[n-1]);
    }
  }
  free(p_d);
}



/*********************************************
 F11: Compositee & Deposite coefficients. return the lest accurate index
      cd>0, A_n --> A_n*sqrt(n);
      cd=0, find the least accurate A_n
      cd<0, A_n --> A_n/sqrt(n).
*********************************************/
int coeff_dm(mpfi_t *C,int cd){
  int i,j=99999,k,m;
  mpfi_t t; mpfi_init(t);
  for ( i=0; i<abs(cd); i++ ){
    mpfi_set_si(t,i+1); mpfi_sqrt(t,t);
    if ( cd<0 ) mpfi_div(C[i],C[i],t);
    else if( cd>0 )mpfi_mul(C[i],C[i],t);
    if ( (k=mpfi_size2(&C[i],0))<j ){ j=k;m=i; }
  }
  mpfi_clear(t);
  return(m);
}




/*******************************************************
 F12: Compute the best N2 and Y2Pi       
*******************************************************/
int N2_and_Y2Pi(mpfi_t R,mpfi_t *Y2Pi,int M,int *Pr, int EndPr, int *ExPr,mpfi_t *bot,int *RndPr,int CN){
  int i,N2,max=*Pr,con=*ExPr;

  i=int_size2(M,2); N2=pow(2,i);
  do{
    N2*=2;
    if ( N2<0 ) break;
    *RndPr=max;
    *ExPr=compute_Y(bot,Y2Pi,R,N2,M,EndPr,RndPr,CN,con);
    *Pr=*ExPr+EndPr+*RndPr;
/*
printf("best...N2=%d Pr=%d Rnd=%d Ex=%d ",N2,*Pr,*RndPr,*ExPr);
mpfi_out_str(stdout,10,5,*Y2Pi); printf("\n");
*/
  } while( *Pr>max || *ExPr<0 );

  if ( N2>0 ) mpfi_endpoint(*Y2Pi,0);
  return(N2);
}




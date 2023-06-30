/**********************************************************************************
  * This page contains functions:
  *   readC_full(): read in pre-computed coefficients c[i], for all 1<=i<=M;
  *   writeC_full(): write c[i], for all 1<=i<=M, to file;
  *   readC_short(): read in pre-computed coefficients c[i], for prime and 
  *                  prime power i<=M;
  *   writeC_short(): write c[i] to file, for prime and prime power i<=M;
  * readIn() returns stage -- the number of H[i] computed
  *          if Pr==0, initialize R, C[i], with 0<=i<=CN
  *          if Pr!=0, initialize C[i], but do not initialize R
  Note: Y is not initialized in any condition. You need read in
        with Pr=0 to get computing precision, then initialize
        Y (also need to initialize R if you changed Pr). After
        that read again (with Pr!=0) to obtain the info of H[i]   
  Note: readC and writeC are agree on the file format.
**********************************************************************************/

/*********************** File Format (readC & writeC) ***************************
int grogram prec -- default mpfi precision in program
int precision ----- desired precision of coefficients
int styleclass ---- 0 or 1
int CN ------------ number of coefficients
mpfi R
mpfi C[1] --------- coefficients
...
mpfi C[CN]
***************************************************************/

int readC1(char *str,int *oe,int *CN,int *Pr,int *EndPr, int *N2,mpfi_t *size,mpfi_t *R,mpfi_t *Y,mpfi_t **C,mpfi_t **H){
  int i,stage,inpr;
  FILE *fp;
  mpfr_t le,ri;

  if ((fp=fopen(str,"r"))==NULL) return(-1);

  fscanf(fp,"%d",&inpr);
  fscanf(fp,"%d",EndPr);
  if ( *EndPr==0 || *EndPr==1 ) *oe=*EndPr;
  else  fscanf(fp,"%d",oe);
  fscanf(fp,"%d",CN);

  if ( !(*Pr) ){
    mpfr_set_default_prec((long int)inpr);
    mpfi_init(*R);
  }  
  *C=(mpfi_t *)malloc((*CN)*sizeof(mpfi_t));
  mpfr_init(le); mpfr_init(ri);  
  for(i=0; i<(*CN); i++){ 
    mpfi_init((*C)[i]);
    mpfr_inp_str(le,fp,10,GMP_RNDD);
    mpfr_inp_str(ri,fp,10,GMP_RNDU);
    mpfi_interv_fr((*C)[i],le,ri);
  }
  mpfi_set(*R,(*C)[0]);
  mpfi_set_si((*C)[0],1);
  if ( fscanf(fp,"%d",&stage)==EOF ) stage=0;
  else if ( stage>0 && !(*Pr) ){
    fscanf(fp,"%d",N2);
    mpfi_init(*size);
    mpfr_inp_str(le,fp,10,GMP_RNDD);
    mpfr_inp_str(ri,fp,10,GMP_RNDU);
    mpfi_interv_fr(*size,le,ri);
    mpfi_init(*Y);
    mpfr_inp_str(le,fp,10,GMP_RNDD);
    mpfr_inp_str(ri,fp,10,GMP_RNDU);
    mpfi_interv_fr(*Y,le,ri);
    *H=(mpfi_t *)malloc((*N2)*sizeof(mpfi_t));
    for(i=0; i<stage; i++){ 
        mpfr_inp_str(le,fp,10,GMP_RNDD);
        mpfr_inp_str(ri,fp,10,GMP_RNDU);
        mpfi_init((*H)[i]);
        mpfi_interv_fr((*H)[i],le,ri);
    }
  }
  fclose(fp);
  if ( !(*Pr) )*Pr=inpr;
  mpfr_clear(le); mpfr_clear(ri);
  return(stage);
}


int writeC1(char *str,int oe,int CN,int Pr,mpfi_t size,int N2,mpfi_t R,mpfi_t Y,mpfi_t *C,mpfi_t *H,int Stage){
  int d,i,k;
  FILE *fp;
  mpfr_t t; mpfr_init(t);

  k=mpfr_get_default_prec();
  fp=fopen(str,"w");
  fprintf(fp,"%d\n",k);         
  fprintf(fp,"%d\n",Pr);         
  fprintf(fp,"%d\n",oe);                 
  fprintf(fp,"%d\n",CN);
  
  d=ceil(log(2)/log(10)*k)+3;
  mpfi_get_left(t,R);  mpfr_out_str(fp,10,d,t,GMP_RNDD); fprintf(fp,"\n");
  mpfi_get_right(t,R); mpfr_out_str(fp,10,d,t,GMP_RNDU); fprintf(fp,"\n");
  for ( i=0; i<CN; i++ ){
    mpfi_get_left(t,C[i]);  mpfr_out_str(fp,10,d,t,GMP_RNDD);fprintf(fp,"\n");
    mpfi_get_right(t,C[i]); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,"\n");
  }

  if ( Stage>0 ){
    fprintf(fp,"%d\n",Stage);
    fprintf(fp,"%d\n",N2);
    mpfi_get_left(t,size);  mpfr_out_str(fp,10,d,t,GMP_RNDD);fprintf(fp,"\n");
    mpfi_get_right(t,size); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,"\n");
    mpfi_get_left(t,Y);  mpfr_out_str(fp,10,d,t,GMP_RNDD);fprintf(fp,"\n");
    mpfi_get_right(t,Y); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,"\n");
  }
  for ( i=0; i<Stage; i++ ){
    mpfi_get_left(t,H[i]); mpfr_out_str(fp,10,d,t,GMP_RNDD);fprintf(fp,"\n");
    mpfi_get_right(t,H[i]); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,"\n");
  }
  fclose(fp);
  mpfr_clear(t);
  return 0;
}

int readC0(char *str,int *oe,int *CN,int *Pr,int *EndPr, mpfi_t *R, mpfi_t **C){
  int i,inpr;
  FILE *fp;
  mpfr_t le,ri;

  if ((fp=fopen(str,"r"))==NULL) return(-1);

  fscanf(fp,"%d",&inpr);
  fscanf(fp,"%d",EndPr);
  if ( *EndPr==0 || *EndPr==1 ) *oe=*EndPr;
  else  fscanf(fp,"%d",oe);
  fscanf(fp,"%d",CN);

  if ( !(*Pr) ){ mpfr_set_default_prec((long int)inpr); mpfi_init(*R); *Pr=inpr; }  

  *C=(mpfi_t *)malloc((*CN)*sizeof(mpfi_t));
  mpfr_init(le); mpfr_init(ri);  
  for(i=0; i<(*CN); i++){ 
    mpfi_init((*C)[i]);
    mpfr_inp_str(le,fp,10,GMP_RNDD);
    mpfr_inp_str(ri,fp,10,GMP_RNDU);
    mpfi_interv_fr((*C)[i],le,ri);
  }
  mpfi_set(*R,(*C)[0]);
  mpfi_set_si((*C)[0],1);
  fclose(fp);
  mpfr_clear(le); mpfr_clear(ri);
  return(0);
}


int writeC0(char *str,int oe,int CN,int Pr,int EndPr, mpfi_t R, mpfi_t *C){
  int d,i,k;
  FILE *fp;
  mpfr_t t; mpfr_init2(t,(unsigned long)Pr);

  fp=fopen(str,"w");
  fprintf(fp,"%d\n",Pr);         
  fprintf(fp,"%d\n",EndPr);
  fprintf(fp,"%d\n",oe);                 
  fprintf(fp,"%d\n",CN);
  
//  i=ceil(log(2.0)/log(10.0)*(double)(EndPr+20));
  d=ceil(log(2.0)/log(10.0)*(double)Pr);
//  if ( i<k ) d=i; else d=k;
  mpfi_get_left(t,R);  mpfr_out_str(fp,10,d,t,GMP_RNDD); fprintf(fp,"\n");
  mpfi_get_right(t,R); mpfr_out_str(fp,10,d,t,GMP_RNDU); fprintf(fp,"\n");
  for ( i=0; i<CN; i++ ){
    mpfi_get_left(t,C[i]);  mpfr_out_str(fp,10,d,t,GMP_RNDD);fprintf(fp,"\n");
    mpfi_get_right(t,C[i]); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,"\n");
  }

  fclose(fp);
  mpfr_clear(t);
  return 0;
}



int writeGP(char *str,mpfi_t R,mpfi_t *C,int CN,int oe,int EndPr){ 
  int i,d=ceil(log(2)/log(10)*(double)EndPr);
  FILE *fp;
  mpfr_t t; mpfr_init(t);

  fp=fopen(str,"w");
  fprintf(fp,"prec=%d;\n",d); d+=3;
  fprintf(fp,"classtype=%d;\n",oe);        
  fprintf(fp,"R=vector(2);\n");
  mpfi_get_left(t,R); fprintf(fp,"R[1]="); mpfr_out_str(fp,10,d,t,GMP_RNDD); fprintf(fp,";\n");
  mpfi_get_right(t,R); fprintf(fp,"R[2]="); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,";\n");
  fprintf(fp,"C=vector(%d);\n",CN);
  for ( i=0; i<CN; i++ ){
    fprintf(fp,"C[%d]=vector(2);\n",i+1);
    mpfi_get_left(t,C[i]); fprintf(fp,"C[%d][1]=",i+1); mpfr_out_str(fp,10,d,t,GMP_RNDD); fprintf(fp,";\n");
    mpfi_get_right(t,C[i]); fprintf(fp,"C[%d][2]=",i+1); mpfr_out_str(fp,10,d,t,GMP_RNDU);fprintf(fp,";\n");
  }
  fclose(fp);
  mpfr_clear(t);
  return 0;
}



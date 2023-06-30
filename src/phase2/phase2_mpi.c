#include "time.h"
#include <string.h>
#include <mpi.h>
#include "kbessel.c"
#include "read_write.c"
#include "basic_fn.c"
#include "mpfi_fft.c"

#define INITTAG 1
#define DIETAG 2
#define WORKTAG 3
#define MaxDigits 128	/*** the max digits of mpfi type going to be used ***/
#define RoundErrorConstant 20
typedef struct{ int index; char left[MaxDigits]; char right[MaxDigits]; } my_type;

Pack_In(my_type *var, mpfi_t z,int i){
  mpfr_t x,y;  mpfr_init(x); mpfr_init(y);

  mpfi_get_left(x,z); mpfi_get_right(y,z);  
  mpfr_sprintf((*var).left,"%RDe",x);
  mpfr_sprintf((*var).right,"%RUe",y);
  (*var).index=i;
  mpfr_clear(x); mpfr_clear(y);
}

Pack_Out(my_type var,mpfi_t *z,int *i){
  mpfr_t x,y;  mpfr_init(x); mpfr_init(y);

  mpfr_set_str(x,var.left,10,MPFR_RNDU);
  mpfr_set_str(y,var.right,10,MPFR_RNDD);
  mpfi_interv_fr(*z,x,y);
  *i=var.index;
  mpfr_clear(x); mpfr_clear(y);
}



int main(int argc,char *argv[]){
  if (argc !=5 ){ printf("Usage: InFile OutFile Coef endpr\n"); return(1); }

  int rank;
  char In[128];  sscanf(argv[1],"%s",In);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if ( rank!=0 ) slave(In);
  else {	//begin master code

  showTime();
/*************************************************************************
 * N2=2N, Number of FFT; M=# of C[i] to compute;CN=# of C[i] known from 
   phase 1; Omega=0/1:even/odd form;Stage:current H[i];
*************************************************************************/
  int N2=4,Omega,M,CN,Stage,Pr=0,EndPr,ExPr,RndPr, i,j,k,m,nproc;
  char Out[128];
  mpfi_t t,R, Y2Pi, *C, *H;	// z^*=x^*=i*y^*; Y2Pi=2*Pi*Y,Y below fundamental Domain. H[i]=\sum_1^{L_j} g(m)*cos(2\pi m x_s);
  MPI_Status status;
  my_type var;

  sscanf(argv[2],"%s",Out);
  sscanf(argv[3],"%d",&M);
  sscanf(argv[4],"%d",&EndPr);  EndPr+=(int_size2(M,1)+1)/2;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  readC0(In,&Omega,&CN,&Pr,&i,&R,&C);
  coeff_mult(C,CN);
printf("CN=%d M=%d Pr=%d R=",CN,M,Pr);mpfi_out_str(stdout,10,8,R); printf("\n");
  i=coeff_dm(C,-CN);	// change C_n to C_n/sqrt(n)

  /****************** Preliminary computation *******************/

  if ( i>=EndPr && CN>=M ){ printf("No computation needed.\n"); Stop_Slave(nproc); return(1); }

  Pr=mpfi_size2(&R,0)-mpfi_size2(&R,1);	// Check the max possible Precision

  if ( Pr<EndPr ){ printf("Previous Precision is not enough.\n"); Stop_Slave(nproc); return 0; }

printf("Max_Pr=%d\n",Pr);

  mpfi_init(Y2Pi); mpfi_init(t);

  ExPr=RoundErrorConstant;
  N2=N2_and_Y2Pi(R, &Y2Pi, M, &Pr, EndPr,&ExPr,&t,&RndPr,CN);	// compute the best N2 and Y2Pi
  if ( N2<1 ){ Stop_Slave(nproc); return(0); }

  /************* Main Computation ************/

printf("Final N=%d ExPr=%d RndPr=%d Pr=%d \nY=",N2,ExPr,RndPr,Pr);
mpfi_out_str(stdout,10,15,Y2Pi);
printf(" size=%d\nR=",mpfi_size2(&Y2Pi,0));
mpfi_out_str(stdout,10,15,R);
printf(" size=%d\n",mpfi_size2(&R,0));

  /************* Send Initial Info. to slaves ************/

  Pack_In(&var,Y2Pi,Pr);
  for( i=1; i<nproc; i++ )
    MPI_Send(&var,sizeof(my_type),MPI_BYTE,i,INITTAG,MPI_COMM_WORLD);
  Pack_In(&var,t,N2);
  for( i=1; i<nproc; i++ )
    MPI_Send(&var,sizeof(my_type),MPI_BYTE,i,INITTAG,MPI_COMM_WORLD);

  /************* Main Loop, send jobs to slaves ************/

  mpfr_set_default_prec(Pr); 
  H=(mpfi_t*)malloc((N2/2)*sizeof(mpfi_t));
showTime();
  i=Stage=0; j=0; k=N2/2-i; if ( k>nproc-1 ) k=nproc-1;
  while ( i<N2/2+nproc-1 ){
    if ( i==N2/2+k ) j=k;
    if ( i>Stage+k-1 && i<N2/2+k ){
      MPI_Recv(&var, sizeof(my_type), MPI_BYTE, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      m=var.index;    mpfi_init(H[m]);
      Pack_Out(var,&(H[m]),&j);
      j=status.MPI_SOURCE;
    }
    else j++;
    if ( i>=N2/2 )
      MPI_Send(0, 0, MPI_INT, j, DIETAG, MPI_COMM_WORLD);
    else
      MPI_Send(&i,1,MPI_INT,j,WORKTAG,MPI_COMM_WORLD);
    i++;
  }	//end while
showTime();
  
i=mpfi_size2(H,-N2/2);
j=mpfi_size2(H,N2/2);
printf("The min H[%d]=",i);mpfi_out_str(stdout,10,8,H[i]);printf(" size=%d\n",mpfi_size2(&H[i],0));
printf("The max H[%d]=",j);mpfi_out_str(stdout,10,8,H[j]);printf(" size=%d\n",mpfi_size2(&H[j],0));

  /************* FFT & Finishing ************/
printf("Finishing...\n");
  After_Loop(H, Y2Pi, R, N2, Omega, M);

  i=coeff_dm(H,M);	// change C_n to C_n*sqrt(n)

printf("Multiply sqtr[M]\n");
i=mpfi_size2(H,-M); 
j=mpfi_size2(H,M);
printf("  The min C[%d]=",i);mpfi_out_str(stdout,10,8,H[i]);
printf(" size=%d\n",mpfi_size2(&H[i],0));
printf("  The max C[%d]=",j);mpfi_out_str(stdout,10,8,H[j]);
printf(" size=%d\n",mpfi_size2(&H[j],0));

  j=mpfi_size2(&H[i],0);
  EndPr-=(int_size2(M,1)+1)/2;

  if ( EndPr>j ) EndPr=j;	// the exact min prec

printf("Writing out... The min pre %d at %d to %s\n",EndPr,i,Out);
  writeC0(Out,Omega,M,Pr,EndPr,R,H);
  sprintf(Out,"%s.gp",Out);
//  writeGP(Out,R,H,M,Omega,EndPr);
  showTime();
}	//end of master code
  MPI_Finalize();
  return 0;
}



/*******************************************************
 F01: Stop all slave processors
*******************************************************/
Stop_Slave(n){ int i; for ( i=1; i<n; i++ ) MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD); }


/*******************************************************
 F03: check the max precision can get from phase1 data
*******************************************************/
int Max_Precision(mpfi_t *C,mpfi_t R,int CN){
  int i,j=0;
  mpfi_t t,u,v;  mpfi_init(t); mpfi_init(u); mpfi_init(v);

//  j=mpfi_size2(&R,0)-mpfi_size2(&R,1);  //Check 1
printf("  MaxPr="); mpfi_out_str(stdout,10,20,R);printf("\tInPr=%d\n",j);

  mpfi_set_si(v,3); mpfi_sqrt(v,v); mpfi_const_pi(t); mpfi_mul(v,t,v);  // v=\sqrt{3}*Pi
  mpfi_set_si(u,0); 
  i=Cn_sqrt_bessel(&t,C,v,R,CN,u,-2);  //Check 2
  if ( (i=mpfi_size2(&t,0))<j ) j=i;
printf("  MaxPr="); mpfi_out_str(stdout,10,20,t);printf("\tInPr=%d\n",i);
  mpfi_mul_si(u,v,CN+1); error_M(&t,R,u,v);
  i=mpfi_size2(&t,1);  //Check 3
  if ( i<j ) j=i;
printf("  MaxPr="); mpfi_out_str(stdout,10,20,t);printf("\tInPr=%d\n",i);

  mpfi_clear(t); mpfi_clear(u); mpfi_clear(v);
  return j;
}
  


/*******************************************************
 F05: All computation after the main loop
*******************************************************/
After_Loop(mpfi_t *H,mpfi_t Y2Pi,mpfi_t R,int N2,int Omega,int M){
  int i,j;
  mpfi_t *tm,v;  mpfi_init(v);

  tm=(mpfi_t *)malloc(M*sizeof(mpfi_t));
  for ( i=0; i<M; i++ ) mpfi_init(tm[i]);

  complete_fft(H,N2/2,Omega);  // FFt.

  Cn_sqrt_bessel(tm,&v,Y2Pi,R,M,v,-3);

printf("B Left-hand side precision\n");
i=mpfi_size2(tm,-M);
j=mpfi_size2(tm,M);
printf("  The min M[%d]=",i);mpfi_out_str(stdout,10,8,tm[i]);
printf(" size=%d\n",mpfi_size2(&tm[i],0));
printf("  The max M[%d]=",j);mpfi_out_str(stdout,10,8,tm[j]);
printf(" size=%d\n",mpfi_size2(&tm[j],0));

printf("After FFT, precision of H\n");
i=mpfi_size2(H,-M); 
j=mpfi_size2(H,M);
printf("  The min H[%d]=",i);mpfi_out_str(stdout,10,8,H[i]);
printf(" size=%d\n",mpfi_size2(&H[i],0));
printf("  The max H[%d]=",j);mpfi_out_str(stdout,10,8,H[j]);
printf(" size=%d\n",mpfi_size2(&H[j],0));

  for ( i=0; i<M; i++ ){
    mpfi_div_si(H[i],H[i],N2/4);
    mpfi_div(H[i],H[i],tm[i]);
    mpfi_clear(tm[i]);
  }

printf("Finel result...C=A/(N2/4)/B\n");
i=mpfi_size2(H,-M); 
j=mpfi_size2(H,M);
printf("  The min H[%d]=",i);mpfi_out_str(stdout,10,8,H[i]);
printf(" size=%d\n",mpfi_size2(&H[i],0));
printf("  The max H[%d]=",j);mpfi_out_str(stdout,10,8,H[j]);
printf(" size=%d\n",mpfi_size2(&H[j],0));

  mpfi_clear(v);
}


/****************************************
 F06: Slave code 
****************************************/
slave(char *str){
  mpfi_t *C, Y2Pi, R, *Box, x_s, y_s, t,v;
  int k=0,BoxSize,Omega,CN,N2,i;
  my_type var;
  MPI_Status status;
  
/**************** revieve precision and initialize variables ***************/
  MPI_Recv(&var, sizeof(my_type), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  if (status.MPI_TAG == DIETAG) {return;}
  readC0(str,&Omega,&CN,&k,&i,&R,&C);	// read in
  i=coeff_dm(C,-CN);	// C_n -> C_n/sqrt(n)
  mpfi_init(Y2Pi); mpfi_init(x_s); mpfi_init(y_s);

  Pack_Out( var, &Y2Pi,&k );
  mpfr_set_default_prec(var.index);	// set default prec
  mpfi_init(t);

/***************** revieve precisioin to compute Box ***********************/
  MPI_Recv(&var, sizeof(my_type), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  if (status.MPI_TAG == DIETAG){return;}
  Pack_Out( var, &t,&N2 );
  BoxSize=compute_Box( &Box, R, t, CN);

/*********************** main task *************************************/
  while (1) {
    MPI_Recv(&k, 1, MPI_INT, 0, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
    if (status.MPI_TAG == DIETAG) {break;}
    compute_Points(&x_s,&y_s,Y2Pi,N2,k);
    i=which_box(Box,y_s,BoxSize);
    Cn_sqrt_bessel(&t,C,y_s,R,i,x_s,Omega);

    Pack_In(&var,t,k);

    MPI_Send(&var, sizeof(my_type), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  }
  for ( i=1; i<BoxSize; i++ ) mpfi_clear(Box[i]);  
  mpfi_clear(x_s); mpfi_clear(y_s); mpfi_clear(R); mpfi_clear(Y2Pi);mpfi_clear(t);
// mpfi_clear(v);
}



showTime(){
  char *date; time_t timer;
  timer=time(NULL);
  date = asctime(localtime(&timer));
  printf("\nCurrent Date: %s", date);
//  putchar('\n');
}









#ifndef QUADRATIC_INFRASTRUCTURE_ELEMENT_H_
#define QUADRATIC_INFRASTRUCTURE_ELEMENT_H_

#include <stdbool.h>
#include <arb.h>

//QuadraticInfrastructureElement
typedef struct {
  //Data
  long Delta;
  long floor_sqrt_Delta;
  arb_t sqrt_Delta;
  long NC_BOUND;
  
  long a;
  long b;
  long c;

  long qlist[100];
  long num_q;
  
  arb_t distance;
  long precision;
} QuadInfElement;



//QuadInfElement functions
void print         (QuadInfElement * qie);

// for testing
void test_ideal    (QuadInfElement * qie, char *msg);

//
// normalization, reduction, assignment
//
void normalize     (QuadInfElement * qie);
void reduce        (QuadInfElement * qie);

void init          (QuadInfElement * qie, long Delta, long precision);
void clear         (QuadInfElement * qie);

void assign_one    (QuadInfElement * qie);
void assign_abc    (QuadInfElement * qie, long a, long b, long c);


//
// baby step and adjustment routines
//
void rho           (QuadInfElement * qie);
void rho_inverse   (QuadInfElement * qie);
void adjust        (QuadInfElement * qie, double bound);
void adjust_to_one (QuadInfElement * qie);
// for testing
void print_cycle(const long delta, long precision);

//
// arithmetic
//

void construct_relative_generator(long *rel_gen_a, long *rel_gen_b, long *rel_gen_d, QuadInfElement *qie, long OB, long BB, long S);

//void construct_relative_generator(long *rel_gen_a, long *rel_gen_b, long *rel_gen_d, QuadInfElement *qie, long Ca, long S);

void square        (QuadInfElement * qie);
void nudupl        (QuadInfElement * qie);

void close_ideal         (QuadInfElement * qie, double dist);

//
// main function for regulator refinement
//
void refine_regulator (arb_t new_reg, long Delta, double reg, long precision);
bool test_refined_regulator(arb_t new_reg, double reg, long precision);

#endif

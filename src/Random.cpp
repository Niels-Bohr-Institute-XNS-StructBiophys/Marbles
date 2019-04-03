#include "Random.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>


#define DEBUG
using namespace std;

/** class constructor */
RandomNumbers::RandomNumbers() {

#ifdef DEBUG /** if the define DEBUG is present fix the seed to 0. Generates a reproducible string of rands */
  srand(0);
#else /** otherwise, use the given time as seed */
  srand(time(NULL));
#endif

  /** allocate gsl arrays for rng calculations */
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, rand());

  srand((unsigned) 1);
}
//------------------------------------------------------------------------------

RandomNumbers::~RandomNumbers() {
  gsl_rng_free(r); /** frees gsl_array */
}
//------------------------------------------------------------------------------

double RandomNumbers::in_range( double min, double max ) {
  return min + ( max - min ) * uniform();
}

double RandomNumbers::uniform() {
  return gsl_rng_uniform_pos(r); /** calls gsl function for uniform random numbers in [0,1] */
}
//------------------------------------------------------------------------------

double RandomNumbers::in_range2( double min, double max ) {
  double R; /*Random number from uniform distribution in the open interval: (min, max)*/
  do{
  R=(double) rand()/RAND_MAX ;
  } while (R==0 || R==1);
  R=R*(max-min)+min;
  return R;
}

double RandomNumbers::gaussian( double tau ) {

  /** compute two uniforms */
  double r = uniform();
  double theta = uniform();

  return sqrt( - 2.*log(r) ) * sin(2.*M_PI*theta) * sqrt(tau);

}
//------------------------------------------------------------------------------

/*******************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

#include "Random.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>


//#define DEBUG
using namespace std;

/** class constructor */
RandomNumbers::RandomNumbers() {

//#ifdef DEBUG /** if the define DEBUG is present fix the seed to 0. Generates a reproducible string of rands */
//srand(0);
//#else /** otherwise, use the given time as seed */
//int tmp = 1559051515;//time(NULL);
//cout << tmp << endl;
srand(time(NULL));
//#endif

  /** allocate gsl arrays for rng calculations */
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, rand());
  //gsl_rng_set(r, -3.);

  //srand((unsigned) 1);
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

//VERY STRANGE WAY TO GENERATE A UNIT VECTOR
std::vector<double> RandomNumbers::vector2( double min_modulus, double max_modulus ) {

  std::vector<double> vec(3);
  double v1, v2, s, alpha, modulus;

  modulus = in_range2( min_modulus, max_modulus );

  do{
    v1 = 2. * rand() / RAND_MAX - 1.;
    v2 = 2. * rand() / RAND_MAX - 1.;
    s = v1 * v1 + v2 * v2;
  } while( s >= 1.0 );

  alpha = 2. * sqrt( 1. - s );
  vec[0] = modulus * ( 2. * s - 1. );
  vec[1] = modulus * alpha * v2;
  vec[2] = modulus * alpha * v1;

  //cout << modulus << " " << sqrt( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] ) << endl;

  return vec;

}

std::vector<double> RandomNumbers::vector3( double modulus ) {

  std::vector<double> vec(3);
  double v1, v2, s, alpha;

  do{
    v1 = 2. * rand() / RAND_MAX - 1.;
    v2 = 2. * rand() / RAND_MAX - 1.;
    s = v1 * v1 + v2 * v2;
  } while( s >= 1.0 );

  alpha = 2. * sqrt( 1. - s );
  vec[0] = modulus * ( 2. * s - 1. );
  vec[1] = modulus * alpha * v2;
  vec[2] = modulus * alpha * v1;

  //cout << modulus << " " << sqrt( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] ) << endl;

  return vec;

}

std::vector<double> RandomNumbers::vector( double modulus ) {

  std::vector<double> vec(3);
  double x, y, z, norm;

  x = in_range( -1, 1 );
  y = in_range( -1, 1 );
  z = in_range( -1, 1 );

  norm = sqrt( x*x + y*y + z*z );

  vec[0] = modulus * x / norm;
  vec[1] = modulus * y / norm;
  vec[2] = modulus * z / norm;

  return vec;
}

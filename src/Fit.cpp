#include <math.h>
#include <nlopt.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "Fit.h"

//using namespace std;

Fit::Fit() {
  setup = false;
}

Fit::~Fit() {
}

double Fit::chi2_background( const std::vector<double> &x, std::vector<double> &grad, void *data ) {

  double sum_squares = 0;
  double tmp;

  bck_data *d = (bck_data *) data;
  std::vector<double> exper = d->exper;
  std::vector<double> err = d->err;
  unsigned int len = d->len;

  for( unsigned int i = 0; i < len; i++ ) {
    tmp  = ( x[0] - exper[i] ) / err[i];
    sum_squares += tmp * tmp;
  }

  return sum_squares / ( len - 1 );
}

double Fit::chi2_intensity( const std::vector<double> &x, std::vector<double> &grad, void *data ) {

  double r, q, tmp, tmp2, exponent;
  double sum_squares = 0;
  int_data *d = (int_data *) data;

  std::vector<double> intensity( d->nq );
  fill(intensity.begin(),intensity.end(),0);

  for( int i = 0; i < d->nq; i++ ) {
    q = d->exp_q[i];
    exponent = x[0] * q * x[0] * q;
    r = exp( - exponent / 2. );

    for(int l = 0; l <= d->harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp = abs( r * d->alpha.at( i, l, m ) + d->beta.at( i, l, m ) );
        tmp *= tmp;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
      }
    }

    intensity[i] = intensity[i] * d->e_scattlen * d->e_scattlen * d->scale_factor + d->background;
    tmp2 = intensity[i] - d->ref_intensity[i];
    sum_squares += tmp2 * tmp2;

    //std::cout << sum_squares << std::endl;

  }

  return sum_squares;
}

// void Fit::setup_intensity_fit( std::vector<std::complex<double> > a, std::vector<std::complex<double> > b, std::vector<double> intens, std::vector<double> q, double escatt, double scale, double bck, unsigned int num_q, unsigned int ho  ) {
//
//   unsigned int harmonics_order;
//   std::vector<double> exp_q;
//   std::vector<double> ref_intensity;
//
//   datai.nq = num_q;
//   datai.harmonics_order = ho;
//   datai.e_scattlen = escatt;
//   datai.background = bck;
//   datai.scale_factor = scale;
//   datai.exp_q = q;
//   datai.ref_intensity = intens;
//   datai.alpha.resize_width( datai.nq );
//   datai.beta.resize_width( datai.nq );
//
//   datai.alpha.copy_from_vector( a, datai.nq, datai.harmonics_order + 1, datai.harmonics_order + 1 );
//   datai.beta.copy_from_vector( b, datai.nq, datai.harmonics_order + 1, datai.harmonics_order + 1 );

  // if( !setup ) {
  //   harmonics_order = ho;
  //   nq              = num_q;
  //   e_scattlen      = escatt;
  //   scale_factor    = scale;
  //   background      = bck;
  //
  //   alpha.resize_width( nq );
  //   beta.resize_width( nq );
  //   exp_q.resize( nq );
  //   intensity.resize( nq );
  //   ref_intensity.resize(nq);
  //
  //   ref_intensity = intens;
  //   exp_q = q;
  //
  //   alpha.copy_from_vector( a, nq, harmonics_order + 1, harmonics_order + 1 );
  //   beta.copy_from_vector( b, nq, harmonics_order + 1, harmonics_order + 1 );
  // }
  //
  // setup = true;
//}

// void Fit::parametrized_intensity( double xr ) {
//
//   double r, q, tmp, exponent;
//   std::vector<double> intensity( nq );
//
//   fill(intensity.begin(),intensity.end(),0);
//
//   for( int i = 0; i < nq; i++ ) {
//     q = exp_q[i];
//     exponent = xr * q * xr * q;
//     r = exp( - exponent / 2. );
//
//     for(int l = 0; l <= harmonics_order; l++ ) {
//       for(int m = 0; m <= l; m++ ) {
//         tmp = abs( r * alpha.at( i, l, m ) + beta.at( i, l, m ) );
//         tmp *= tmp;
//         intensity[i] += ( (m > 0) + 1. ) * tmp;
//       }
//     }
//
//     intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor + background;
//
//     std::cout << intensity[i] - ref_intensity[i] << std::endl;
//
//   }
// }

void Fit::fit_background( std::vector<std::vector<double> > rad, unsigned int bck_npoints_to_fit ) {

  unsigned int opt_dimension = 1; //number of parameters passed to the optimizer
  unsigned int j;
  unsigned int nq = rad.size();
  std::vector<double> exper(bck_npoints_to_fit ), err(bck_npoints_to_fit );

  //load scattering data from file
  //in the main code this will be not necessary and data will be passed by BeadModeling
  for( unsigned int i = 0; i < bck_npoints_to_fit; i++ ) {
    j = nq - i - 1;
    exper[i] = rad[j][1];
    err[i]   = rad[j][2];
  }

  std::vector<double> lb(opt_dimension);
  y_bck.resize(opt_dimension);

  datab.len   = bck_npoints_to_fit; //number of high-q points to fit the background
  datab.exper = exper; //experimental intensity
  datab.err   = err; //error on the experimental intensity

  //COBYLA is a gradient-less local optimizer. Works very well with such a simple problem
  nlopt::opt cobyla(nlopt::LN_COBYLA, opt_dimension); //initialize the optimization class

  lb[0] = 0.; //lower bound for the value of the background is 0
  y_bck[0] = 0.; //initial value of the background fed to the optimizer
  cobyla.set_lower_bounds(lb); //assign the lower bounds to the optimizer

  cobyla.set_min_objective( chi2_background, &datab ); //assign to the optimizer the function to optimize
  cobyla.set_xtol_rel(1e-10); //assign to optimizer the relative tollerance

  try {
    nlopt::result result = cobyla.optimize(y_bck, chi2_bck); //run the otimization
  } catch(std::exception &e) {
    //something went wrong! call an exception and exit
    std::cout << "nlopt failed: " << e.what() << std::endl;
    exit(-1);
  }

  std::cout << "found minimum at X^2(" << y_bck[0] << ") = " << chi2_bck << std::endl;
}

void Fit::fit_intensity( std::vector<std::complex<double> > a, std::vector<std::complex<double> > b, std::vector<double> intens, std::vector<double> q, double escatt, double scale, double bck, unsigned int num_q, unsigned int ho ) {

  unsigned int opt_dimension = 1; //number of parameters passed to the optimizer

  std::vector<double> lb(opt_dimension);
  y_int.resize(opt_dimension);

  datai.nq = num_q;
  datai.harmonics_order = ho;
  datai.e_scattlen = escatt;
  datai.background = bck;
  datai.scale_factor = scale;
  datai.exp_q = q;
  datai.ref_intensity = intens;
  datai.alpha.resize_width( datai.nq );
  datai.beta.resize_width( datai.nq );
  datai.alpha.copy_from_vector( a, datai.nq, datai.harmonics_order + 1, datai.harmonics_order + 1 );
  datai.beta.copy_from_vector( b, datai.nq, datai.harmonics_order + 1, datai.harmonics_order + 1 );

  //COBYLA is a gradient-less local optimizer. Works very well with such a simple problem
  nlopt::opt cobyla(nlopt::LN_COBYLA, opt_dimension); //initialize the optimization class

  lb[0] = 0.; //lower bound for the value of the background is 0
  y_int[0] = 0.; //initial value of the background fed to the optimizer
  cobyla.set_lower_bounds(lb); //assign the lower bounds to the optimizer

  cobyla.set_min_objective( chi2_intensity, &datai ); //assign to the optimizer the function to optimize
  cobyla.set_xtol_rel(1e-10); //assign to optimizer the relative tollerance

  try {
    nlopt::result result = cobyla.optimize(y_int, chi2_int); //run the otimization
  } catch(std::exception &e) {
    //something went wrong! call an exception and exit
    std::cout << "nlopt failed: " << e.what() << std::endl;
    exit(-1);
  }

  std::cout << std::setprecision(5) << "found minimum at X^2(" << y_int[0] << ") = " << chi2_int << std::endl;
}

#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include <cmath>

#define GSL_NEW

using namespace std;

int main() {

  int order = 17;
  int ntheta = 36;
  double theta_step = M_PI / ntheta;
  double c, theta;
  double legendre[ntheta][order+1];

  for( int m = 0; m <= order; m += 2) {
    for(int l = m; l <= ntheta; l += 2 ) {
      for( int t = 0; t < ntheta; t++ ) {
        theta = ( ntheta + 0.5 ) * theta_step;
        c = cos(theta);

        if( l == m ) {
#ifdef GSL_NEW
          gsl_sf_legendre_array( GSL_SF_LEGENDRE_SPHARM, order, c, &legendre[t][l] );
          cout << legendre[t][l] << endl;
#else
          gsl_sf_legendre_sphPlm_array(order, m, c, &legendre[t][l] );
          cout << legendre[t][l] << endl;
#endif
        }
      }
    }
  }
}

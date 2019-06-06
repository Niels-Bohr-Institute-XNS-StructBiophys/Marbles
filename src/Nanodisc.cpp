#include "Nanodisc.h"
#include <cmath>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>
#define pi 3.141593

using namespace std;

Nanodisc::Nanodisc() {
  //info about the file from which to load are not available at initialization time
  //therefore load_input() has to be called seperately
}

Nanodisc::~Nanodisc() {
  //nothing to be done here
}

void Nanodisc::load_input( const string& best_fit ) {

  /**
   * Reads from the WillItFit output file.
   * The number of lines to skip is hardcoded and needs to be re-checked everytime the WillItFit version changes.
   * Excess scattering lengths are in units of the electron scattering length, and are thus adimensional.
   * Volumes are in A^3.
   *
   * Inputs
   * ------
   * const strin& best_fit: path to the best fit file output from WillItFit.
   */

  volume_tests = true; //might want to ask the user about this!

  ifstream file( best_fit );
  string line;
  string d = "=", d_ = " ", d__ = ","; //possible delimiters
  double tmp;

  if( file.is_open() ) {

    skip_lines( file, 18 );
    hbelt                 = stod( parse_line( file, d_ ) );
    nlipids               = stod( parse_line( file, d_ ) );

    hbelt = 24;
    nlipids = 126.65951;

    skip_lines( file, 1 );
    wathead               = stod( parse_line( file, d_ ) );

    wathead = 0.;

    skip_lines( file, 2 );
    xrough                = stod( parse_line( file, d_ ) );
    cvbelt                = stod( parse_line( file, d_ ) );
    cvlipid               = stod( parse_line( file, d_ ) );
    cvprotein             = stod( parse_line( file, d_ ) );

    xrough = 4.609096;
    cvbelt = 0.906671;
    cvlipid = 1.034406;

    printf("%.5f\n", nlipids);

    cout << nlipids << " " << nlipids << " " << wathead << endl;
    cout << xrough << " " << cvbelt << " " << cvlipid << endl;

    //HARDCODED HERE!!!1
    //cvprotein = 1.;
    //CHANGE THIS!!

    skip_lines( file, 7 );
    cvwater               = stod( parse_line( file, d_ ) );

    cvwater = 1.;

    skip_lines( file, 1 );
    vertical_axis_endcaps = stod( parse_line( file, d_ ) );
    scale_endcaps         = stod( parse_line( file, d_ ) );

    skip_lines( file, 31 );
    hlipid                = stod( parse_line( file, d ) );
    hcore                 = stod( parse_line( file, d ) );
    hmethyl               = stod( parse_line( file, d ) );

    hlipid = 40.1603;
    hcore = 29.8818;
    hmethyl = 3.4992;

    //cout << hlipid << " " << hcore << " " << hmethyl << endl;

    skip_lines( file, 7 );
    radius_major          = stod( parse_line( file, d ) );
    radius_minor          = stod( parse_line( file, d ) );

    radius_major = 43.3604;
    radius_minor = 29.8501;

    skip_lines( file, 2 );
    width_belt            = stod( parse_line( file, d ) );

    width_belt = 8.1721;

    skip_lines( file, 3 );
    rho_h2o               = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;
    rho_d2o               = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;
    rho_head              = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;
    rho_alkyl             = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;
    rho_methyl            = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;
    rho_belt              = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;

    rho_head = 4.62e-11 / e_scatt_len;
    rho_belt = 3.31e-9 / e_scatt_len;

    skip_lines( file, 1 );
    rho_protein           = stod( parse_double_delimiter( file, d_, d__ ) ) / e_scatt_len;

    skip_lines( file, 2 );
    vh2o                  = stod( parse_double_delimiter( file, d_, d__ ) );
    vd2o                  = stod( parse_double_delimiter( file, d_, d__ ) );
    vhead                 = stod( parse_double_delimiter( file, d_, d__ ) );
    valkyl                = stod( parse_double_delimiter( file, d_, d__ ) );
    vmethyl               = stod( parse_double_delimiter( file, d_, d__ ) );
    vbelt                 = stod( parse_double_delimiter( file, d_, d__ ) );

    skip_lines( file, 1 );
    vprotein              = stod( parse_double_delimiter( file, d_, d__ ) );

    tmp = sqrt( scale_endcaps * scale_endcaps - 1. ) / scale_endcaps;
    vertical_axis_ellipsoid = vertical_axis_endcaps / ( 1. - tmp );

    //compute excess scattering length densities from the fit-corrected volumes
    rho_head    += wathead * rho_h2o;
    rho_head    /= ( vhead * cvlipid + wathead * vh2o );
    rho_alkyl   /= ( valkyl * cvlipid );
    rho_methyl  /= ( vmethyl * cvlipid );
    rho_belt    /= ( vbelt * cvbelt );
    rho_protein /= ( vprotein * cvprotein );
    rho_h2o     /= ( vh2o * cvwater );

    /** Three checks on volumes are ran to verify the consistency of the quantities loaded from the WillItFit output. */
    if( volume_tests ) {

      double vh1, vh2, va1, va2, vm1, vm2;

      vh1 = ( vhead * cvlipid + wathead * vh2o ) * nlipids;
      vh2 = ( hlipid - hcore ) * radius_minor * radius_major * M_PI;

      if( fabs( vh1 - vh2 ) > 0.5 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check 0 failed: inconsistent volume of lipid heads." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }

      va1 = valkyl * cvlipid * nlipids;
      va2 = ( hcore - hmethyl ) * radius_minor * radius_major * M_PI;

      if( fabs( va1 - va2 ) > 60 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check 1 failed: inconsistent volume of alkyl heads." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }

      vm1 = vmethyl * cvlipid * nlipids;
      vm2 = hmethyl * radius_minor * radius_major * M_PI;

      if( fabs( vm1 - vm2 ) > 0.1 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check 2 failed: inconsistent volume of methyl groups." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }
    }

    // cout << hbelt << " " << nlipids << endl;
    // cout << watheads << endl;
    // cout << xrough << " " << cvbelt << " " << cvlipids << " " << cvmp << endl;
    // cout << cvwater << endl;
    // cout << vertical_axis_endcaps << " " << scale_endcaps << endl;
    // cout << radius_major << " " << radius_minor << " " << width_belt << endl;
    // cout << rho_h2o << " " << rho_d2o << " " << rho_head << " " << rho_alkyl << endl;
    // cout << rho_methyl << " " << rho_belt << " " << rho_protein << endl;
    // cout << v_h2o << " " << v_d2o << " " << v_head << " " << v_alkyl << endl;
    // cout << v_methyl << " " << v_belt << " " << v_protein << endl;

  } else {
    cout << "Cannot open '" << best_fit << "'" << endl;
  }

  file.close();
}

void Nanodisc::flat_disc_form_factor( double a, double b, double L, double rho, double q, int index ) {

  /** Computes the orientationally averaged form factor of a disc for different values of the momentum:
    * psi = 2 * J1( q1 ) / q1 * sinc( q * L * cos(alpha/2) )
    * where
    * q1 = q * r * sin(alpha)
    * r = sqrt( a^2 * sin^2(phi) + b^2 * cos^2(phi) )
    * Stores the value in F.
    *
    * Inputs
    * ------
    * double a:   major semiaxis of the disc
    * double b:   minor semiaxis of the disc
    * double L:   height of the disc
    * double rho: excess scattering length density of the disc
    * double q:   value of the momentum
    * int index:  index on the momentum mesh
    */

  double sinc, sin_t, theta, phi, a_phi, b_phi, qr;
  vector<double> r( nphi );

  double theta_step = M_PI / ntheta;
  double phi_step = 2. * M_PI / nphi;
  double v_rho = rho * a * b * M_PI * L;

  //precompute values for phi to avoid over-computation
  for( unsigned int p = 0; p < nphi; p++ ) {
    phi    = ( p + 0.5 ) * phi_step;
    a_phi  = a * cos(phi);
    b_phi  = b * sin(phi);
    r[p]   = sqrt( a_phi * a_phi + b_phi * b_phi );
  }

  for( unsigned int t = 0; t < ntheta; t++ ) {

    theta = ( t + 0.5 ) * theta_step;
    sinc  = boost::math::sinc_pi( L/2. * q * cos(theta) );
    sin_t = sin(theta);

    for( unsigned int p = 0; p < nphi; p++ ) {
      qr = q * r[p] * sin_t;

      if( qr == 0 ) {
        F.add( index, t, p, v_rho * sinc );
      } else {
        F.add( index, t, p, 2. * v_rho * boost::math::cyl_bessel_j( 1, qr ) / qr * sinc );
      }
    }
  }
}

double Sinc(double x){
    double result;
    if(x==0)
        result=1.0;
    else
        result=sin(x)/x;
    return result;
}

void Nanodisc::flat_disc_form_factor2( double a, double b, double L, double rho, double q, int index ) {

    //cout << "# FORM FACTOR SHIT DIOCANE" << endl;
    //Søren Kynde 2012
    //This function calculates the form factor F(Q,theta,phi) of a cylinder
    //with elliptical cross-section with half axes a and b and height L. The
    //cylinder may be dispaced from the origin by a (real space) vector
    //(r0, theta0, phi0).

    double r;
    double cosr0q=0;
    //int ntheta,nphi;
    double sinc,theta,phi,tmp;
    double thetastep = M_PI/ntheta, phistep = 2. * M_PI / nphi;
    double volume=a*b*L*M_PI;

    // printf("%lf %lf %lf %lf %lf\n", thetastep, phistep, volume, rho, L);
    // printf("%lf %lf\n", a, b);
    // exit(-1);
    //
    // a = 44.521200;
    // b = 30.286500;

    //cout << thetastep << " " << phistep << " " << volume << " " << rho << " " << L << endl;
    //exit(-1);

    // for(ntheta=0; ntheta < Ntheta; ntheta++){
    //   for(nphi=0; nphi<Nphi; nphi++){
    //     F[ntheta][nphi] = 0.;
    //   }
    // }

    //printf("%d %d", Ntheta, Nphi);

    for(int t =0; t < ntheta; t++){
        theta=(t+.5)*thetastep;
    //printf("P: %g, %g, %g, %g, %g, %g, %g\n",Q,a,b,rho,L,theta,sinc);
        //sinc=gsl_sf_bessel_j0( L/2.*Q*cos(theta) );
        sinc=Sinc( L/2.*q*cos(theta) );
        //printf("%d %g %g\n", ntheta, sinc, sin(theta) );
        for(int p=0; p<nphi; p++){
            phi=(p+.5)*phistep;

            r=sin(theta)*sqrt( pow(a*cos(phi),2)+pow(b*sin(phi),2) );
            //if( ntheta == 0 ) printf("%g\n", sqrt( pow(a*cos(phi),2)+pow(b*sin(phi),2) ));
            //printf("%d %d %g\n", ntheta, nphi, Q*r );

            //printf("%g\n", Q*r);
            if(q*r==0) {
                tmp = volume*rho*sinc;
                F.add(index, t, p, tmp);
                //F[ntheta][nphi] += volume*pol(1.0*sinc,0.0);
            } else {
                tmp = volume*rho*2*j1(q*r)/(q*r)*sinc;
                F.add( index, t, p, tmp);
              }
                //F[ntheta][nphi]+=volume*pol(2*j1(Q*r)/(Q*r)*sinc,-r0*Q*cosr0q);

            //printf("%g\n", F[ntheta][nphi] );
        }
    }
}

double Nanodisc::PsiEllipticCylinderWithEndcaps(double q, double Alpha, double Beta, double MajorRadius,
  double MinorRadius, double Height, double ScaleFactorOfEndcaps, double VerticalAxisOfEndcaps)
{
    /// Declarations
    // Dummies
    //int i;
    //int j;
    int k;

    double Dummy1;
    double Dummy2;
    double Dummy3;

    double EquivalentRadiusInCylinder;
    double EquivalentRadiusInEndcaps;

    const int TSteps = 50;
    double TMin;
    double TMax;
    double TStepSize;
    double T;

    double ReturnValue;
    double SumOverT = 0.0f;

    // Assign values to arguments
    //ScaleFactorOfEndcaps = 1.5f;
    double MajorRadiusOfCurvatureOfEndcaps = MajorRadius * ScaleFactorOfEndcaps;
    double MinorRadiusOfCurvatureOfEndcaps = MinorRadius * ScaleFactorOfEndcaps;
    double ShiftOfEndcaps;

    //printf("%g %g\n", VerticalAxisOfEndcaps, ScaleFactorOfEndcaps);


    // Caps center inside cylidner
    ShiftOfEndcaps = - VerticalAxisOfEndcaps / MinorRadiusOfCurvatureOfEndcaps * sqrt(pow(MinorRadiusOfCurvatureOfEndcaps, 2) - pow(MinorRadius, 2));

    EquivalentRadiusInCylinder = sqrt(pow(MajorRadius * cos(Beta), 2) + pow(MinorRadius * sin(Beta), 2));
    EquivalentRadiusInEndcaps  = sqrt(pow(MajorRadiusOfCurvatureOfEndcaps * cos(Beta), 2) + pow(MinorRadiusOfCurvatureOfEndcaps * sin(Beta), 2));

    TMin = - ShiftOfEndcaps / VerticalAxisOfEndcaps;
    TMax = 1.0f;
    TStepSize = (TMax - TMin) / TSteps;

    Dummy1 = pi * MajorRadius * MinorRadius * Height *
        2.0f * sin(q * Height / 2.0f * cos(Alpha)) / (q * Height / 2.0f * cos(Alpha)) *
        j1(q * EquivalentRadiusInCylinder * sin(Alpha)) / (q * EquivalentRadiusInCylinder * sin(Alpha));

    SumOverT = 0.0f;

    for (k = 0; k < TSteps; ++k) {
        T = k * TStepSize + TMin;

            Dummy2 = cos(q * cos(Alpha) * (VerticalAxisOfEndcaps * T + ShiftOfEndcaps + Height / 2.0f)) *
            (1 - pow(T, 2)) * j1(q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0f - pow(T, 2))) /
            (q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0f - pow(T, 2)));

            SumOverT += Dummy2 * TStepSize;
    }

    Dummy3 = 4.0f * pi * MajorRadiusOfCurvatureOfEndcaps * MinorRadiusOfCurvatureOfEndcaps * VerticalAxisOfEndcaps * SumOverT;
    ReturnValue=Dummy1+Dummy3;
    if(isnan(SumOverT)){
        printf("NH2\n");
    }

    return ReturnValue;
}

void Nanodisc::disc_w_endcaps_form_factor( double a, double b, double L, double rho, double q, int index ) {

    //This function calculates the form factor F(Q,theta,phi) of a cylinder
    //with elliptical endcaps. The
    //cylinder may be dispaced from the origin by a (real space) vector
    //(r0, theta0, phi0).

    double theta, phi, tmp;
    double theta_step = M_PI / ntheta;
    double phi_step = 2 * M_PI / nphi;

    for( unsigned int t = 0; t < ntheta; t++ ) {
        theta = ( t + 0.5 ) * theta_step;

        for( unsigned int p = 0; p < nphi; p++ ) {
            phi = ( p + 0.5 ) * phi_step;

            tmp = rho * PsiEllipticCylinderWithEndcaps(q, theta, phi, a, b, L, scale_endcaps, vertical_axis_ellipsoid);
            F.add( index, t, p, tmp );
        }
    }
}

// double Nanodisc::expand_sh( int index ) {
//
//   double theta, phi;
//   double theta_step = M_PI/ ntheta;
//   double phi_step = 2. * M_PI / nphi;
//   vector<complex<double> > fm(ntheta);
//   vector<vector<complex<double> > > phase( harmonics_order + 1, vector<complex<double> >(nphi));// [harmonics_order+1][nphi];
//   int l,m,i,j;
//   double intensity = 0;
//   vector<double> sinth(ntheta), w(ntheta);
//   vector<vector<double> > legendre( ntheta, vector<double>(harmonics_order+1));
//
//   for( int t = 0; t < ntheta; t++ ) {
//       w[t]  = 0.;
//       theta = ( t + 0.5 ) * theta_step;
//
//       for( int l = 0; l < ntheta/2; l++ ) {
//         w[t] += 2./(ntheta/2)*1./(2*l+1)*sin((2*l+1)*theta);
//       }
//
//       //cout << w[t] << endl;
//   }
//   //fine here!!
//
//   for( int m = 0; m < harmonics_order + 1; m += 2 ) {
//
//     for( int t = 0; t < ntheta; t++ ) {
//       fm[t] = 0;
//
//       for( int p = 0; p < nphi; p++ ) {
//         phi = ( p + 0.5 ) * phi_step;
//
//         if( t == 0 ) {
//           phase[m][p] = pol( phi_step, - m * phi );
//         }
//
//         fm[t] += F.at( index, t, p ) * phase[m][p];
//       }
//     }
//
//     for( int l = m; l <= harmonics_order; l += 2 ) {
//       for( int t = 0; t < ntheta; t++ ) {
//         theta = ( t + 0.5 ) * theta_step;
//
//         if( l == m ) {
//           gsl_sf_legendre_sphPlm_array( harmonics_order, m, cos(theta), &legendre[t][l] );
//           sinth[t] = sin(theta);
//         }
//
//         complex<double> tmp = 1/sqrt(4*M_PI) * legendre[t][l] * w[t] * sinth[t] * fm[t];
//         //cout << legendre[t][l] << " " << sinth[t] << endl;
//         alpha.add( index, l, m, tmp );
//         //cout << real(alpha.at( index, l, m )) << " " << imag(alpha.at( index, l, m )) << endl;
//       }
//
//       intensity += ((m>0)+1)*pow( abs( alpha.at(index,l,m) ), 2 );
//     }
//   }
//   //
//   // for( int m = 0; m < harmonics_order + 1; m += 2 ) {
//   //   for( int l = m; l <= harmonics_order; l += 2 ) {
//   //     cout << real(alpha.at( index, l, m )) << " " << imag(alpha.at( index, l, m )) << endl;
//   //   }
//   // }
//
//
//   return intensity;
// }

double Nanodisc::expand_sh( int index ) {

  double theta, phi;
  double theta_step = M_PI/ ntheta;
  double phi_step = 2. * M_PI / nphi;
  double intensity = 0.;
  double sqrt_4pi_1 = 1. / sqrt( 4. * M_PI );
  complex<double> tmp;
  vector<complex<double> > fm(ntheta);
  vector<vector<complex<double> > > phase( harmonics_order + 1, vector<complex<double> >(nphi));
  vector<double> sin_t(ntheta), cos_t(ntheta), w(ntheta);
  vector<vector<double> > legendre( ntheta, vector<double>(harmonics_order+1));

  for( int i = 0; i < ntheta; i++ ) {
    for( int j = 0; j < harmonics_order+1; j++ ) {
      legendre[i][j] = 0.;
    }
  }


  for( unsigned int t = 0; t < ntheta; t++ ) {
      w[t]  = 0.;
      theta = ( t + 0.5 ) * theta_step;
      sin_t[t] = sin(theta);
      cos_t[t] = cos(theta);

      for( int l = 0; l < ntheta/2; l++ ) {
        w[t] += 2. / (ntheta / 2) * 1. / (2 * l + 1) * sin( (2 * l + 1) * theta );
      }
  }

  for( int m = 0; m < harmonics_order + 1; m += 2 ) {

    for( int t = 0; t < ntheta; t++ ) {
      fm[t] = {0., 0.};

      for( int p = 0; p < nphi; p++ ) {
        phi = ( p + 0.5 ) * phi_step;

        if( t == 0 ) {
          phase[m][p] = pol( phi_step, - m * phi );
        }

        fm[t] += F.at( index, t, p ) * phase[m][p];
        //cout << real(F.at( index, t, p )) << " " << imag(F.at( index, t, p )) << endl;
      }
    }

    for( int l = m; l <= harmonics_order; l += 2 ) {
      for( int t = 0; t < ntheta; t++ ) {

        //if( l == m ) gsl_sf_legendre_array( norm, harmonics_order, m, cos_t[t], &legendre[t][l] );

        if( l == m ) gsl_sf_legendre_sphPlm_array( harmonics_order, m, cos_t[t], &legendre[t][l] );

         tmp = sqrt_4pi_1 * legendre[t][l] * w[t] * sin_t[t] * fm[t];
         //cout << sqrt_4pi_1 << " " << cos_t[t] << " " << legendre[t][l] << " " << w[t] << " " << sin_t[t]  << " " << real(fm[t]) << " " << imag(fm[t]) << endl;
        //cout << legendre[t][l] << " " << sinth[t] << endl;
        alpha.add( index, l, m, tmp );
        //cout << real(alpha.at( index, l, m )) << " " << imag(alpha.at( index, l, m )) << endl;
      }

      intensity += ((m>0)+1)*pow( abs( alpha.at(index,l,m) ), 2 );
    }
  }

  return intensity;
}

complex<double> pol2(double r, double phi)
{
  if( phi == 0 )
    return { r, 0. };
  else
    return { r * cos(phi), r * sin(phi) };
}

double Nanodisc::expand_sh2( int index ) {
//Søren Kynde 2011
//This function calculates the coefficients Alm of the spherical harmonics
//expansion of an analytical form factor F(theta,phi) (phi is the azimuthal angle)
    double theta,phi;
    double thetastep=M_PI/(ntheta), phistep=2*M_PI/(nphi);
    //int ntheta,nphi;
    complex<double> fm[ntheta];//={0};
    complex<double> phase[harmonics_order+1][nphi];
    int l,m,i,j,p, t;
    double Int=0;
    double sinth[ntheta];
    double w[ntheta];//={0};
    double legendre[ntheta][harmonics_order+1];

    for(i=0;i<ntheta;i++){
        fm[i]={0.,0.};
    }

    for(i=0;i<ntheta;i++){
        w[i]=0;
    }
    for(j=0;j<ntheta;j++){ //calculate weights dtheta for theta integral
        theta=(j+.5)*thetastep;
        for(l=0;l<ntheta/2;l++){
            w[j]+=(double) 2./(ntheta/2)*1./(2*l+1)*sin((2*l+1)*theta);
        }

        //printf("%g\n", w[j]);
    }

    for(m=0;m<harmonics_order+1;m+=2){ //Because symetry of disc only the even harmonics contribute
        for(t=0;t<ntheta;t++){
            fm[t]=0.;
            for(p=0;p<nphi;p++){  //integration over phi
                phi=phistep*(p+.5);
                if(t==0) {
                    phase[m][p]=pol2(phistep,-m*phi);

                  }
                fm[t]+= F.at(index,t,p)*phase[m][p];  //fm(theta)= int_0^2pi [ F(theta,phi) exp(-m*phi) dphi ]
                //printf("%g %g\n", creal(F[ntheta][nphi]), cimag(F[ntheta][nphi]) );
                //printf("%g %g\n", creal(fm[ntheta]), cimag(fm[ntheta]));
                //printf("%g\n", F[ntheta][nphi]);
            }
            //printf("%g %g\n", creal(fm[ntheta]), cimag(fm[ntheta]));

        }
        for(l=m;l<=harmonics_order;l+=2){   //For disc symmetry only even harmonics contribute
            for(t=0;t<ntheta;t++){    //integration over theta
                theta=(t+.5)*thetastep;
                if(l==m){
                    gsl_sf_legendre_sphPlm_array(harmonics_order,m,cos(theta),&legendre[t][l]);//Calculate all Plm(cos(th)) for l=m ... to l=N;
                    sinth[t]=sin(theta);
                }
                alpha.add( index, l, m, 1/sqrt(4*M_PI)*legendre[t][l]*w[t]*sinth[t]*fm[t]);
                //alpha[l][m]+=1/sqrt(4*M_PI)*legendre[t][l]*w[t]*sinth[t]*fm[t];//alpha_lm=int_0^pi [ P_lm*sin(theta)*fm(theta)*dtheta ]
                //printf("%lf %lf %lf %lf %lf %lf\n", 1/sqrt(4*M_PI), cos(theta), legendre[t][l], w[t], sinth[t], real(fm[t]), imag(fm[t]) );
                //printf("%lf %lf\n", real(alpha.at(index,l,m)), imag(alpha.at(index,l,m)) );

            }
            Int+=((m>0)+1)*pow( abs( alpha.at(index,l,m) ), 2 ); //sum over l and m
            //// Intensity can be calculated here or outside the function in
            //seperate loop
        }
    }

    //cout << Int << endl;

    // for(m=0;m<Nh+1;m+=skip){
    //   for(l=m;l<=Nh;l+=skip){   //For disc symmetry only even harmonics contribute
    //       printf("%g %g\n", creal(alpha[l][m]), cimag(alpha[l][m]));
    //     }
    //   }

    return Int;
}

//compatible with previous version within 5e-4 relative error.
void Nanodisc::nanodisc_form_factor( vector<double> exp_q ) {

  double q;
  int dim = exp_q.size();

  F.resize_width( dim );
  //F.initialize(0);

  alpha.resize_width( dim );
  alpha.initialize(0);

    //cout << radius_major << " " << radius_minor << " " << fabs(hmethyl) << " " << rho_methyl - rho_alkyl << endl;

  //clock_t begin = clock();
  for( int i = 0; i < dim; i++ ) {

    F.initialize(0);

    q = exp_q[i];
    flat_disc_form_factor( radius_major, radius_minor, hlipid, rho_head - rho_h2o, q, i);
    flat_disc_form_factor( radius_major, radius_minor, hcore, rho_alkyl - rho_head, q, i);
    //disc_w_endcaps_form_factor( radius_major, radius_minor, hlipid, rho_head - rho_h2o, q, i);
    //disc_w_endcaps_form_factor( radius_major, radius_minor, hcore, rho_alkyl - rho_head, q, i);
    flat_disc_form_factor( radius_major, radius_minor, fabs(hmethyl), rho_methyl - rho_alkyl, q, i);
    flat_disc_form_factor( radius_major, radius_minor, hbelt, rho_h2o - rho_belt, q, i);
    flat_disc_form_factor( radius_major + width_belt, radius_minor + width_belt, hbelt, rho_belt - rho_h2o, q, i);

    // for( unsigned int t = 0; t < ntheta; t++ ) {
    //   for( unsigned int p = 0; p < nphi; p++ ) {
    //     cout << F.at(i, t, p) << endl;
    //  }
    // }

    double intensity = expand_sh2( i ); //uncomment
    //cout << exp_q[i] << " " << intensity << endl;
  }
  //exit(-1);
  //clock_t end = clock();
  //double elapsed_secs = (double)(end - begin) / CLOCKS_PER_SEC;
  //cout << "Average time per execution: " << (1.*elapsed_secs)/dim << endl;

}

double Nanodisc::get_radius_major() {
  return radius_major;
}

double Nanodisc::get_radius_minor() {
  return radius_minor;
}

double Nanodisc::get_scale_endcaps() {
  return scale_endcaps;
}

double Nanodisc::get_vertical_axis_ellipsoid() {
  return vertical_axis_ellipsoid;
}

double Nanodisc::get_rho_solvent() {
  return rho_h2o;
}

double Nanodisc::get_hlipid() {
  return hlipid;
}

double Nanodisc::get_hmethyl() {
  return hmethyl;
}

double Nanodisc::get_hcore() {
  return hcore;
}

double Nanodisc::get_rho_methyl() {
  return rho_methyl;
}

double Nanodisc::get_rho_alkyl() {
  return rho_alkyl;
}

double Nanodisc::get_rho_head() {
  return rho_head;
}

double Nanodisc::get_cvprotein() {
  return cvprotein;
}

double Nanodisc::get_xrough() {
  return xrough;
}

double Nanodisc::get_e_scatt_len() {
  return e_scatt_len;
}

complex<double> Nanodisc::get_alpha( int i, int l, int m ) {
  return alpha.at( i, l, m );
}

std::vector<complex<double> > Nanodisc::get_alpha_buffer() {
  return alpha.get_buffer();
}

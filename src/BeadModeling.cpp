#include "BeadModeling.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ctime>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <iomanip>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#define AV_RESIDUE_MASS 0.1 //average mass of a residue in kDa

using namespace std;

/*
 * TODO:
 * find better way to do update_statistics(). In particular, I don't like the fact that you have to pre-compute the distance matrix.
 */

BeadModeling::BeadModeling( const string& filename ) {

  input_file = filename;

  ifstream file( input_file );
  string line;
  string d = " ";

  sphere_generated = false;
  compute_scale = true;
  init = true;

  clash_distance = 1.8; //hardcoded because experimented. Might want to leave the choice open for users though.
  sequence = "";
  shift = 50.; //same here: might want to leave this free for the user to choose
  insertion = 14;
  T_strength = 5;
  H_strength = 20;

  if( file.is_open() ) {

    rad_file       = parse_line( file, d );           //path to the experimental .rad file
    best_fit       = parse_line( file, d );           //path to the WillItFit results
    sequence_file  = parse_line( file, d );           //path to the file with the protein sequence
    dmax           = stod( parse_line( file, d ) );   //dmax from fit
    npasses        = stoi( parse_line( file, d ) );   //total number of passes
    loops_per_pass = stoi( parse_line( file, d ) );   //number of performed loops per pass
    outdir         = parse_line( file, d );           //directory where results are saved
    lambda         = stoi( parse_line( file, d ) );
    connect        = stoi( parse_line( file, d ) );

    string mkdir = "mkdir " + outdir;
    system( mkdir.c_str() );

    mkdir = "mkdir " + outdir + "configurations/";
    system( mkdir.c_str() );

    mkdir = "mkdir " + outdir + "intensities/";
    system( mkdir.c_str() );

  } else {
    cerr << "Cannot open " << input_file << endl;
    exit(-1);
  }

  cout << "# File '" << input_file << "' loaded" << endl;

  load_FASTA(); //load sequence file
  cout << "# File '" << sequence_file << "' loaded" << endl;
  nresidues = sequence.length();

  load_rad(); //experimental file with SAXS or SANS data
  cout << "# File '" << rad_file << "' loaded " << endl;

  load_statistics(); //statistics needed for penalty function
  cout << "# Statistics loaded" << endl;

  nnnum = nnum1_ref.size();

  cout << endl;
  cout << "# SUMMARY OF PARAMETERS" << endl;
  cout << "# ---------------------" << endl;
  cout << "# Number of beads: " << nresidues << endl;
  cout << "# Radius of initial sphere: " << 2. * dmax / 3. << endl;
  cout << "# Number of passes: " << npasses << endl;
  cout << "# Loops per pass: " << loops_per_pass << endl;
  cout << "# Storing results in: '" << outdir << "'" << endl;

  nq = rad.size();
  beads.resize( nresidues );
  exp_q.resize( nq );
  beta.resize_width( nq );
  beta.initialize(0);
  beta_old.resize_width( nq );
  distances.resize( nresidues, nresidues );
  distances_old.resize( nresidues, nresidues );
  distances.initialize(0);
  intensity.resize( nq );
  intensity_old.resize( nq );
  roughness_chi2.resize( 5 );

  for( unsigned int i = 0; i < nq; i++ ) {
    exp_q[i] = rad[i][0];
  }

  nd.load_input_flat( best_fit );

  fit.fit_background( rad, 20 );
  fit.set_default_roughness( 4.6 );
}
//------------------------------------------------------------------------------

void BeadModeling::load_rad() {
  rad = load_matrix( rad_file, 3 );
}
//------------------------------------------------------------------------------

void BeadModeling::load_statistics() {

  string ndist_file = "include/statistics/ndist.dat";
  string nnum1_file = "include/statistics/nnum_5.3.dat";
  string nnum2_file = "include/statistics/nnum_6.8.dat";
  string nnum3_file = "include/statistics/nnum_8.3.dat";

  ndist_ref = load_vector_from_matrix( ndist_file, 1, 2 );
  nnum1_ref = load_vector_from_matrix( nnum1_file, 1, 3 );
  nnum2_ref = load_vector_from_matrix( nnum2_file, 1, 3 );
  nnum3_ref = load_vector_from_matrix( nnum3_file, 1, 3 );

  ndist.resize( ndist_ref.size() );
  nnum1.resize( nnum1_ref.size() );
  nnum2.resize( nnum2_ref.size() );
  nnum3.resize( nnum3_ref.size() );

  ndist_old.resize( ndist_ref.size() );
  nnum1_old.resize( nnum1_ref.size() );
  nnum2_old.resize( nnum2_ref.size() );
  nnum3_old.resize( nnum3_ref.size() );
}
//------------------------------------------------------------------------------

void BeadModeling::load_FASTA() {

    ifstream file( sequence_file );

    if( file.is_open() ) {
      skip_lines( file, 1 );
      string tmp;

      while( !file.eof() ) {
        getline( file, tmp );
        sequence.append( tmp );
      }
    } else {
      cerr << "Cannot open " << sequence_file << endl;
    }

}
//------------------------------------------------------------------------------

double BeadModeling::distance( unsigned const int i, unsigned const int j ) {

  double x = beads[i].x - beads[j].x;
  double y = beads[i].y - beads[j].y;
  double z = beads[i].z - beads[j].z;

  return sqrt( x * x + y * y + z * z );
}
//------------------------------------------------------------------------------

bool BeadModeling::bead_clash( unsigned const int i ) {

  bool clash = false;
  for( unsigned int j = 0; j < nresidues; j++ ) {
    if( beads[i].position_assigned == true && beads[j].position_assigned == true ) {
      if( distance(i,j) < clash_distance && j != i ) {
        clash = true;
        break;
      }
    }
  }

  return clash;
}
//------------------------------------------------------------------------------

//####################################################
//I REMOVED THE SHIFT!!!!!!! REMEMBER TO PUT IT BACK!!!
//####################################################

void BeadModeling::initial_configuration() {

  //if( !sphere_generated ) {

    double x, y, z, r, r2;
    bool clash;

    r = 2. * dmax / 3.; /** radius of the sphere **/
    r2 = r * r;

    for( unsigned int i = 0; i < nresidues; i++ ) {
      beads[i].assign_volume_and_scattlen( string(1, sequence[i]) );

      do {

          clash = false;

          do {

            x = rng.in_range2( -r, r );
            y = rng.in_range2( -r, r );
            z = rng.in_range2( -r, r );
            beads[i].assign_position( x, y, z /*+ shift*/ );

          } while( x*x +  y*y + z*z > r2 ); // condition that defines a sphere

          clash = bead_clash( i );

      } while( clash == true );
    }

  //} else {
  //  cout << "# NOTE! Skipping initial configuration because the the system is already set up.";
  //}

}
//------------------------------------------------------------------------------

void BeadModeling::load_initial_configuration( const string& filename ) {

  ifstream file( filename );
  vector<string> parsed;
  vector<double> x, y, z;
  vector<string> residue;

  if( file.is_open() ) {
    string tmp;

    while( !file.eof() ) {
      getline( file, tmp );

      boost::split( parsed, tmp, boost::is_any_of(" ") ); //equivalent of python split
      parsed.erase( remove( parsed.begin(), parsed.end(), "" ), parsed.end() ); //removes all occurrences of an empty string

      // for( int i = 0; i < parsed.size(); i++ ) {
      //   cout << parsed[i] << ", ";
      //   //cout << beads[i].x << " " << beads[i].y << " " << beads[i].z << endl;
      // }
      // cout << endl;

      if( parsed[0] == "ATOM" ) {
        x.push_back( stod( parsed[6] ) );
        y.push_back( stod( parsed[7] ) );
        z.push_back( stod( parsed[8] ) );
        residue.push_back( parsed[3] );
      }
    }

    nresidues = x.size();
    beads.resize( nresidues );

    for( int i = 0; i < nresidues; i++ ) {
      beads[i].x = x[i];
      beads[i].y = y[i];
      beads[i].z = z[i];
      beads[i].assign_volume_and_scattlen( residue[i] );
      //cout << beads[i].x << " " << beads[i].y << " " << beads[i].z << " " << residue[i] << endl;
    }

    cout << "# NRESIDUES: " << nresidues << endl;

  } else {
    cerr << "Cannot open " << sequence_file << endl;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::write_pdb( const string& filename ) {

    FILE *fil;
    char *amin = "rca";

    fil = fopen( filename.c_str(), "w" );

    for( unsigned int i = 0; i < nresidues; i++ ) {
        fprintf(fil,"ATOM   %4d  CA  %s   %4d   %8.3lf%8.3lf%8.3lf  1.00 18.20           N\n",i+1,amin,i+1,beads[i].x, beads[i].y,beads[i].z);
    }

    fclose(fil);
}
//------------------------------------------------------------------------------

void BeadModeling::write_xyz( const string& filename ) {

  ofstream pdb;

  pdb.open(filename);
  if( pdb.is_open() ) {

    pdb << nresidues << endl << endl;

    for( unsigned int i = 0; i < nresidues; i++ ) {
      pdb << "CA " << beads[i].x << " " << beads[i].y << " " << beads[i].z << endl;
    }

  }
}
//------------------------------------------------------------------------------

void BeadModeling::write_toy_intensity( const string& filename, vector<double> toy_err ) {

  ofstream int_file;

  int_file.open(filename);
  if( int_file.is_open() ) {
    for( unsigned int i = 0; i < nq; i++ ) {
      int_file << exp_q[i] << "\t" << intensity[i] << "\t" << toy_err[i] << endl;
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::write_intensity( const string& filename ) {

  ofstream int_file;

  int_file.open(filename);
  if( int_file.is_open() ) {
    for( unsigned int i = 0; i < nq; i++ ) {
      int_file << exp_q[i] << "\t" << intensity[i] << endl;
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::update_rho( int i ) {

  string model = "endcaps";

  double radius_major            = nd.get_radius_major();
  double radius_minor            = nd.get_radius_minor();
  double scale_endcaps           = nd.get_scale_endcaps();
  double vertical_axis_ellipsoid = nd.get_vertical_axis_ellipsoid();
  double rho_solvent             = nd.get_rho_solvent();
  double hlipid                  = nd.get_hlipid();
  double hmethyl                 = nd.get_hmethyl();
  double hcore                   = nd.get_hcore();
  double rho_alkyl               = nd.get_rho_alkyl();
  double rho_methyl              = nd.get_rho_methyl();
  double rho_head                = nd.get_rho_head();
  double cvprotein               = nd.get_cvprotein();

  double x, y, z, fz;

  x = beads[i].x;
  y = beads[i].y;
  z = beads[i].z;
  fz = fabs(z);

  if( model == "endcaps" ) {
    double a_endcaps               = radius_major * scale_endcaps;
    double a_endcaps_1             = 1. / a_endcaps;
    double b_endcaps               = radius_minor * scale_endcaps;
    double b_endcaps_1             = 1. / b_endcaps;
    double shift_endcaps           = - vertical_axis_ellipsoid / a_endcaps * sqrt( a_endcaps * a_endcaps - radius_major * radius_major );
    double c_endcaps_1             = 1. / vertical_axis_ellipsoid;
    double shift_z_core            = ( hcore / 2. + shift_endcaps ) * c_endcaps_1;
    double shift_z_lipid           = ( hlipid / 2. + shift_endcaps ) * c_endcaps_1;

    //nmethyl = 0;
    //nalkyl  = 0;
    //nhead   = 0;

    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp12;
    bool cnd1, cnd2, cnd3, cnd4;

    tmp1 = x * a_endcaps_1 * x * a_endcaps_1;
    tmp2 = y * b_endcaps_1 * y * b_endcaps_1;
    tmp3 = ( z * c_endcaps_1 - shift_z_core ) * ( z * c_endcaps_1 - shift_z_core );
    tmp4 = ( z * c_endcaps_1 + shift_z_core ) * ( z * c_endcaps_1 + shift_z_core );
    tmp5 = ( z * c_endcaps_1 - shift_z_lipid ) * ( z * c_endcaps_1 - shift_z_lipid );
    tmp6 = ( z * c_endcaps_1 + shift_z_lipid ) * ( z * c_endcaps_1 + shift_z_lipid );
    tmp12 = tmp1 + tmp2;

    cnd1 = ( z > 0 && (tmp12 + tmp3 < 1) );
    cnd2 = ( z < 0 && (tmp12 + tmp4 < 1) );
    cnd3 = ( z > 0 && (tmp12 + tmp5 < 1) );
    cnd4 = ( z < 0 && (tmp12 + tmp6 < 1) );

    //FOR DEBUGGING ONLY!!!
    //#####################
    if( i == 0 ) {
      beads[i].v = 165.18;
      beads[i].rho = 71;
    } else if( i == nresidues - 1 ) {
      beads[i].v = 147.14;
      beads[i].rho = 65;
    }
    //#####################
    // FOR DEBUGGING ONLY!!

    if( fz < hmethyl * .5 ) {
      beads[i].type = 3;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_methyl;
      nmethyl++;
    } else if( cnd1 || cnd2 || fz < hcore * .5 ) {
      beads[i].type = 2;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_alkyl;
      nalkyl++;
    } else if( cnd3 || cnd4 || fz < hlipid * .5 ) {
      beads[i].type = 1;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_head;
      nhead++;
    } else {
      beads[i].type = 0;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_solvent;
    }

  } else if( model == "flat" ) {

    //the nanodisc is symmetric with respect the horizontal axis of the methyl layer

      if( fz < hmethyl * .5 ) {
        beads[i].type = 3;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_methyl;
        nmethyl++;
      } else if( fz < hcore * .5 ) {
        beads[i].type = 2;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_alkyl;
        nalkyl++;
      } else if( fz < hlipid * .5 ) {
        beads[i].type = 1;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_head;
        nhead++;
      } else {
        beads[i].type = 0;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_solvent;
      }

  }

  if( beads[i].type_old == 3 ) {
    nmethyl--;
  } else if( beads[i].type_old == 2  ) {
    nalkyl--;
  } else if( beads[i].type_old == 1 ) {
    nhead--;
  }

}
//------------------------------------------------------------------------------

void BeadModeling::expand_sh( double q, int index, int i, int sign, int indice ) {

  double x, y, z, r, theta, phi;
  double sqrt_4pi = sqrt( 4. * M_PI );
  int status;
  double bessel[ harmonics_order + 1 ];
  vector<double> legendre( harmonics_order + 1 );
  complex<double> j(0,1), tmp, p;

  if( sign < 0 ) {
    x     = beads[i].x_old;
    y     = beads[i].y_old;
    z     = beads[i].z_old;
  } else {
    x     = beads[i].x;
    y     = beads[i].y;
    z     = beads[i].z;
  }

  r     = sqrt( x*x + y*y + z*z );
  theta = acos( z / r );
  phi   = acos( x / ( r * sin(theta) ) ) * sgn( y );

  status = gsl_sf_bessel_jl_array( harmonics_order, q * r, bessel ); // Calculate spherical bessel functions for l=0,..,Nh

  for( int m = 0; m <= harmonics_order; m++ ) {
    status = gsl_sf_legendre_sphPlm_array( harmonics_order, m, cos(theta), &legendre[m] ); //Calculate legendre polynomials P_l(cos(theta)) of degree l=m,..., up to Nh
    //Store the values in legendre[m],legendre[m+1],...,legendre[Nh]

    p = pol( 1., -m * phi);

    for( unsigned int l = m; l <= harmonics_order; l++ ) {
      tmp = sqrt_4pi * pow(j, l) * beads[i].rho_modified * bessel[l] * legendre[l] * p;
      //cout << tmp << endl;

      if( sign >= 0 ) {
        beta.add( index, l, m, tmp );
        //cout << beta.at( index, l, m ) << endl;
      } else {
        beta.add( index, l, m, -tmp );
      }
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::only_prot_intensity( bool toy_model ) {

  double r, q, tmp, exponent, I0;
  double e_scattlen = nd.get_e_scatt_len();
  double background = fit.get_background();//0.00045;

  fill(intensity.begin(),intensity.end(),0);

  for( int i = 0; i < nq; i++ ) {

    for(int l = 0; l <= harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp = abs( beta.at( i, l, m ) );
        tmp *= tmp;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
      }
    }

    //insert a check for the value of the correction factor: if it is similar to the value of the numerical density it's fine. Otherwise, suggest the user to proceed at its own risk.
    if( compute_scale ) {
      I0 = intensity[0] * e_scattlen * e_scattlen;

      if( toy_model ) scale_factor = 1000. / I0; //default normalization
      else scale_factor = rad[0][1]/I0;

    }
    compute_scale = false;

    intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor + background;
  }

}
//------------------------------------------------------------------------------

void BeadModeling::calc_intensity( vector<double> exp_q ) {

  double xr = 4.609096;//fit.get_rough();//4.6;//5.014;//nd.get_xrough();
  //double xr = 5.014;
  double r, q, tmp, exponent, I0;
  double e_scattlen = nd.get_e_scatt_len();
  double background = 0.000346;//fit.get_background();//7.8e-5; //TODO! Load this from WillItFit
  //double background = 7.8e-5;
  //double correction_factor = 2.409e15; //TODO! Understand how to compute this factor

  //intensity.resize( nq );
  fill(intensity.begin(),intensity.end(),0);

  for( int i = 1; i < nq; i++ ) {
    q = exp_q[i];
    exponent = xr * q * xr * q;
    r = exp( - exponent / 2. );

    for(int l = 0; l <= harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp = abs( r * nd.get_alpha( i, l, m ) ); //+ beta.at( i, l, m ) );
        tmp *= tmp;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
      }
    }

    //insert a check for the value of the correction factor: if it is similar to the value of the numerical density it's fine. Otherwise, suggest the user to proceed at its own risk.
    if( compute_scale ) {
      I0 = intensity[1] * e_scattlen * e_scattlen;
      //scale_factor = rad[0][1]/I0; //rescale the computed intensity to the experimental I[0]
      scale_factor = 0.0358436 / I0;
    }
    compute_scale = false;

    //scale_factor = 4e+16;

    //cout << "SCALE " << scale_factor << endl;

    intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor + background;
    //intensity[i] = intensity[i] * e_scattlen * e_scattlen * correction_factor + background;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::distance_matrix() {

  double tmp;

  for( unsigned int i = 0; i < nresidues; i ++ ) {
    for( unsigned int j = i+1; j < nresidues; j++ ) {
      tmp = distance( i, j );
      distances.set( i, j, tmp );
      distances.set( j, i, tmp );
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::update_statistics() {

  int count1, count2, count3;
  double d;

  fill( ndist.begin(), ndist.end(), 0.);
  fill( nnum1.begin(), nnum1.end(), 0.);
  fill( nnum2.begin(), nnum2.end(), 0.);
  fill( nnum3.begin(), nnum3.end(), 0.);

  for( unsigned int i = 0; i < nresidues; i++ ) {

    count1 = count2 = count3 = 0;

    for( unsigned int j = 0; j < nresidues; j++ ) {

      d = distances.at(i,j);
      if( d < 12. ) ndist[ (int)(d) ] += 1. / nresidues;
      if( d < 5.3 ) count1++;
      if( d < 6.8 ) count2++;
      if( d < 8.3 ) count3++;
    }

    if( count1 < nnnum ) nnum1[count1-1] += 1. / nresidues;
    if( count2 < nnnum ) nnum2[count2-1] += 1. / nresidues;
    if( count3 < nnnum ) nnum3[count3-1] += 1. / nresidues;
  }

}
//------------------------------------------------------------------------------

void BeadModeling::chi_squared() {

  double tmp, err;
  //double len = rad.size();
  X = 0.;

  for( unsigned int i = 0; i < nq; i++ ) {
    tmp = intensity[i] - rad[i][1];
    err = rad[i][2];
    X += tmp * tmp / (err * err);
  }
}
//------------------------------------------------------------------------------

void BeadModeling::type_penalty() {

  int tmp;
  tmp = nalkyl + nmethyl + nhead - insertion;
  //cout << "TEMP " << tmp << endl;

  // if( init_type_penalty ) {
  //   T = 2. * T_strength * tmp * tmp;
  // } else {
  if( tmp > 0 ) {
    T = 0;
  } else {
    T = T_strength * tmp * tmp;
  }
  //}

  //init_type_penalty = false;
}
//------------------------------------------------------------------------------

void BeadModeling::histogram_penalty() {

  double tmp1, tmp2, tmp3, tmp4;
  double nnum_len = nnum1_ref.size();
  double ndist_len = ndist_ref.size();
  H = 0;

  for( unsigned int i = 0; i < nnum_len; i++ ) {

    tmp1 = nnum1_ref[i] - nnum1[i];
    tmp2 = nnum2_ref[i] - nnum2[i];
    tmp3 = nnum3_ref[i] - nnum3[i];

    if( i < ndist_len ) {
      tmp4 = ndist_ref[i] - ndist[i];
    } else {
      tmp4 = 0.;
    }

    H += H_strength * ( tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3 ) + tmp4 * tmp4;
  }

  H *= lambda;
}
//------------------------------------------------------------------------------

void BeadModeling::recursive_connect( int i, int s, int *pop ) {

  for( unsigned int j = 0; j < nresidues; j++ ) {
    if( distance(i,j) < 5.81 && i != j ) {
      if( beads[j].burn == 0 ) {
          beads[j].burn = 1;
          pop[s]++;
          recursive_connect(j,s,pop);
      }
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::connect_penalty() {

  int pop[nresidues];// = {0};
  int i, s = 0, max;

  for( unsigned int i = 0; i < nresidues; i++ ) {
      beads[i].burn = 0;
      pop[i] = 0.;
  }

  for( unsigned i = 0; i < nresidues; i++ ) {
      if( beads[i].burn == 0 ) {
          beads[i].burn = 1;
          pop[s]++;
          recursive_connect(i,s,pop);
          s++;
      }
  }

  max = pop[0];
  for(unsigned int i = 1; i < nresidues; i++ ) {
      if( pop[i] >= max ) {
          max = pop[i];
      }
  }

  C = fabs( connect * log( (1. * nresidues) / max ) );
}
//------------------------------------------------------------------------------

void BeadModeling::penalty() {

  P = 0;
  histogram_penalty();
  chi_squared();
  //type_penalty();
  connect_penalty();

  if( init ) {

    /* initially the penalty function is computed in a strange way: needs to be further tested!
     * compute the histogram penalty
     * then compute a modified type penalty
     * set the connect penalty to 0
     * add 100000 to everything
     */

    // double tmp = nalkyl + nmethyl + nhead - insertion;
    // T = 2. * T_strength * tmp * tmp;
    // C = 0;
    // X = 0;
    P = 100000;

  }

  init = false;
  //P += X + H + T + C;
  P += X + H + C;
}
//------------------------------------------------------------------------------

void BeadModeling::save_old_config() {

  for( unsigned int i = 0; i < nresidues; i++ ) {
    beads[i].save_old_config();
  }

  ndist_old     = ndist;
  nnum1_old     = nnum1;
  nnum2_old     = nnum2;
  nnum3_old     = nnum3;
  intensity_old = intensity;
  nmethyl_old   = nmethyl;
  nalkyl_old    = nalkyl;
  nhead_old     = nhead;
  P_old         = P;

  distances_old.copy_from( distances );
  beta_old.copy_from( beta );
}

void BeadModeling::reject_move() {

  for( unsigned int i = 0; i < nresidues; i++ ) {
    beads[i].recover_old_config();
  }

  ndist     = ndist_old;
  nnum1     = nnum1_old;
  nnum2     = nnum2_old;
  nnum3     = nnum3_old;
  intensity = intensity_old;
  nmethyl   = nmethyl_old;
  nalkyl    = nalkyl_old;
  nhead     = nhead_old;
  P         = P_old;

  distances.copy_from( distances_old );
  beta.copy_from( beta_old );

}
//------------------------------------------------------------------------------

bool BeadModeling::inside_ellipse( int i, double a, double b ) {

  double x = beads[i].x;
  double y = beads[i].y;
  double tmp = x * x / (a * a) + y * y / ( b * b );

  return ( tmp < 1 );
}
//------------------------------------------------------------------------------

void BeadModeling::move_only_protein() {

  int s = 0, i, j;
  bool legal;
  save_old_config();

  double rmax, rmin, d2, z_ref;
  vector<double> vec(3);
  d2 = rng.in_range2( 1.8, 5.1 ); //5.1 seems quite random

  do {

    s++;
    legal = true;
    if( s == 1 || s == 1001 ) { /*Try the same set of beads 1000 times*/
      do {
        s = 1;
        i = (int)( rng.in_range2(0, nresidues) ); /*Pick a bead to be moved*/
        j = (int)( rng.in_range2(0, nresidues) ); /*Pick another bead. n is to be placed in contact with m*/
      } while( i == j );
    }

    vec = rng.vector3( d2 );
    beads[i].assign_position( beads[j].x + vec[0], beads[j].y + vec[1], beads[j].z + vec[2] );

    //cout << "YOLO1" << endl;

    if( legal ) {
      legal = ! bead_clash( i );
    }

  } while( legal == false );

  for( unsigned int k = 0; k < nq; k++ ) {
    expand_sh( exp_q[k], k, i, -1, 0 ); //subtract the contribution of the previous position
    expand_sh( exp_q[k], k, i, 1, 0 );  //update beta with the new position
  }

  //cout << "YOLO2" << endl;

  //update_rho( i );

  // for( unsigned int k = 0; k < nq; k++ ) {
  //   expand_sh( exp_q[k], k, i, 1, 0 );  //update beta with the new position
  // }

  only_prot_intensity( false );
  distance_matrix();
  update_statistics();

  penalty();

}
//------------------------------------------------------------------------------

void BeadModeling::move( int l ) {

  int s = 0, i, j;
  bool legal;
  save_old_config();

  //JUST FOR DEBUGGING PURPOSE
  double rmax, rmin, d2, z_ref;
  vector<double> vec(3);
  // I would expect these to be radius_minor and radius_major
  rmax = 42.6;
  rmin = 29.0;
  d2 = rng.in_range2( 1.8, 5.1 ); //5.1 seems quite random
  z_ref = 14; //seems wuite random too
  //END DEBUGGING PURPOSE

  do {

    s++;
    legal = true;
    if( s == 1 || s == 1001 ) { /*Try the same set of beads 1000 times*/
      do {
        s = 1;
        i = (int)( rng.in_range2(0, nresidues) ); /*Pick a bead to be moved*/
        j = (int)( rng.in_range2(0, nresidues) ); /*Pick another bead. n is to be placed in contact with m*/
      } while( i == j );
    }

    vec = rng.vector3( d2 );
    beads[i].assign_position( beads[j].x + vec[0], beads[j].y + vec[1], beads[j].z + vec[2] );

    if( legal ) {
      legal = ( fabs( beads[i].z ) > z_ref || inside_ellipse( i, rmax, rmin ) );
    }

    if( legal ) {
      legal = ! bead_clash( i );
    }

  } while( legal == false );

  for( unsigned int k = 0; k < nq; k++ ) {
    expand_sh( exp_q[k], k, i, -1, 0 ); //subtract the contribution of the previous position
  }

  update_rho( i );

  for( unsigned int k = 0; k < nq; k++ ) {
    expand_sh( exp_q[k], k, i, 1, 0 );  //update beta with the new position
  }

  calc_intensity( exp_q );
  distance_matrix();
  update_statistics();

  penalty();
}
//------------------------------------------------------------------------------

void BeadModeling::generate_toy_model( const string& pdb, double eta ) {

  string toy_intensity = "toy_intensity.dat";
  double rho_solvent = nd.get_rho_solvent();
  double av, std, gauss_err;
  vector<double> toy_err(nq);
  bool toy_model = true;

  cout << endl;
  cout << "# GENERATING TOY MODEL FROM PDB" << endl;
  cout << "# -----------------------------" << endl;

  load_initial_configuration( pdb );

  nmethyl = 0;
  nalkyl = 0;
  nhead = 0;

  for( unsigned int i = 0; i < nresidues; i++ ) {
    beads[i].type = 0;
    beads[i].rho_modified = beads[i].rho - beads[i].v * rho_solvent;
  }

  for( unsigned int j = 0; j < nq; j++ ) {
    for( unsigned int i = 0; i < nresidues; i++ ) {
      expand_sh( exp_q[j], j, i, 1, 0 );
    }
  }

  cout << "# Compute form factor: done!" << endl;

  only_prot_intensity( toy_model );
  cout << "# Compute intensity: done!" << endl;

  for( int i = 0; i < nq; i++ ) {
    av = intensity[i];
    std = sqrt( av ) * eta;
    intensity[i] = av + rng.gaussian( std );
    toy_err[i] = std;

    cout << intensity[i] << endl;

    //if( i > int(90./100 * nq) ) {
    //  std = av * 2;
    //}
    //toy_err[i] = 0.02;//fabs(std * rng.gaussian( std ));
  }

  // double scale = 0.1/intensity[0];
  //
  // for( int i = 0; i < nq; i++ ) {
  //
  //   intensity[i] *= scale;
  //   toy_err[i] *= scale;
  //   //cout << av << " " << sqrt(av) << " " << std << endl;
  // }

  write_toy_intensity( toy_intensity, toy_err );

  cout << "# Done!" << endl;
}
//------------------------------------------------------------------------------

void BeadModeling::SA_only_protein() {

  bool decreasing_p, metropolis, accept;
  int iterations = 1, rough_counter;
  double mean, sq_sum, stdev, scale_tmp;
  bool toy_model = false;

  double rho_solvent = nd.get_rho_solvent();

  cout << endl;
  cout << "# PRELIMINARIES" << endl;
  cout << "# -------------" << endl;

  initial_configuration();
  cout << "# Initial configuration set." << endl;
  //string xyz =
  //write_xyz( xyz );

  nmethyl = 0;
  nalkyl = 0;
  nhead = 0;

  for( unsigned int i = 0; i < nresidues; i++ ) {
    beads[i].type = 0;
    beads[i].rho_modified = beads[i].rho - beads[i].v * rho_solvent;
  }

  for( unsigned int j = 0; j < nq; j++ ) {
    for( unsigned int i = 0; i < nresidues; i++ ) {
      expand_sh( exp_q[j], j, i, 1, 0 );
    }
  }

  cout << "# Compute form factor: done!" << endl;
  only_prot_intensity( toy_model );

  distance_matrix();
  update_statistics();

  cout << "# Update statistics: done!" << endl;
  cout << "# Background: " << std::setprecision(2) << fit.get_background() << " (X^2_R = " << fit.get_bck_chi2() << ")" << endl;
  //cout << "# Starting roughness: " << fit.get_rough() << endl;
  cout << "# Scale factor (/1e15): " << scale_factor/1.e15 << ", " << intensity[1] << endl;

  cout << endl;
  cout << "# SIMULATED ANNEALING" << endl;
  cout << "# -------------------" << endl;
  cout << endl;

  penalty();

  //B = X/10; //effective temperature
  B = X / 10;
  // if( B > 10000 ) {
  //   B = 10000;
  // }
  ofstream penalty_file;

  penalty_file.open( outdir + "penalty.dat" );
  penalty_file << "#Iterations\tTemperature\tChi2\tType\tHistogram\tConnect\tTotal" << endl;

  for( unsigned int p = 0; p < npasses; p++ ) {

    int c = 0;
    int attempts = 0;

    cout << "# PASS " << p << endl;
    for( unsigned int l = 0; l < loops_per_pass; l++ ) {

      do {

        attempts++;

        //cout << "CIAO" << endl;

        move_only_protein();

        //cout << "HEY" << endl;

        decreasing_p = ( P < P_old );
        double tmp = rng.in_range2(0.,1.);
        metropolis   = ( exp( - (P - P_old)/B ) > tmp );
        accept = ( decreasing_p || metropolis );

        if( !accept ) {
          reject_move();
        }

      } while( accept == false );

      penalty_file << iterations << "\t" << B << "\t" << X << "\t" << T << "\t" << H << "\t" << C << "\t" << P << endl;
      iterations++;
    }

    scale_tmp = rad[0][1]/intensity[0];
    scale_factor *= scale_tmp;


    //cout << "# Statistics                    " << endl;
    cout << fixed << setprecision(2) << setfill('0');
    cout << setw(5) << "# Acceptance ratio:  " << (1.*loops_per_pass)/attempts << endl;
    cout << setw(5) << "# Temperature:       " << B << endl;
    cout << setw(5) << "# Chi squared:       " << X << endl;
    //cout << setw(5) << "# Type penalty:      " << T << endl;
    cout << setw(5) << "# Histogram penalty: " << H << endl;
    cout << setw(5) << "# Connect penalty:   " << C << endl;
    cout << setw(5) << "# Total penalty:     " << P << endl;
    //cout << setw(5) << "# Inserted beads:    " << nhead << " " << nalkyl << " " << nmethyl << endl;
    cout << setw(5) << "# I_exp[0]/I[0]:     " << setprecision(3) << scale_tmp << endl;

    // if( fit_rough && p%skip_passes == 0 && p != 0 ) {
    //   cout << setw(5) << "# Fitted roughness:  " << fit.get_rough() << " (X^2_R = " << fit.get_rough_chi2() << ")" << endl;
    // } else {
    //   cout << setw(5) << "# Fitted roughness:  " << fit.get_rough() << endl;
    // }
    cout << setw(5) << scientific << "# Scale factor:      " << scale_factor << endl;
    cout << endl;

    string xyz = outdir + "configurations/" + to_string(p) + ".xyz";
    string pdb = outdir + "configurations/" + to_string(p) + ".pdb";
    string calc_intensity = outdir + "intensities/" + to_string(p) + ".dat";
    write_xyz( xyz );
    write_pdb( pdb );
    write_intensity( calc_intensity );

    if( B > 0.0001 ) {
      B *= 0.9;
    }

    if( (1.*loops_per_pass)/attempts < 0.011 ) {
      cout << "\n # Acceptance has reached the lower limit. Execution will be stopped." << endl;
      exit(0);
    }
  }

  penalty_file.close();


}

void BeadModeling::test_flat() {

  bool decreasing_p, metropolis, accept;
  int iterations = 1, rough_counter;
  double mean, sq_sum, stdev, scale_tmp;

  cout << endl;
  cout << "# PRELIMINARIES" << endl;
  cout << "# -------------" << endl;

  initial_configuration();
  cout << "# Initial configuration set." << endl;

  nmethyl = 0;
  nalkyl = 0;
  nhead = 0;

  //############################
  //REMEMBER TO UNCOMMENT!!!!!!
  //############################
  // for( unsigned int i = 0; i < nresidues; i++ ) {
  //   update_rho( i );
  // }
  cout << "# Update scattering lengths: done!" << endl;

  nd.nanodisc_form_factor( exp_q );

  for( unsigned int j = 0; j < nq; j++ ) {
    for( unsigned int i = 0; i < nresidues; i++ ) {
      expand_sh( exp_q[j], j, i, 1, 0 );
    //exit(-1);
    }
  }

  cout << "# Compute form factor: done!" << endl;


  calc_intensity( exp_q );
  cout << "# Compute intensity: done!" << endl;

  for( int i = 0; i < nq; i++ ) {
    cout << exp_q[i] << " " << intensity[i] << endl;
  }
  exit(-1);


  distance_matrix();
  update_statistics();

  cout << "# Update statistics: done!" << endl;
  cout << "# Background: " << std::setprecision(2) << fit.get_background() << " (X^2_R = " << fit.get_bck_chi2() << ")" << endl;
  cout << "# Starting roughness: " << fit.get_rough() << endl;
  cout << "# Scale factor (/1e15): " << scale_factor/1.e15 << endl;

  cout << endl;
  cout << "# SIMULATED ANNEALING" << endl;
  cout << "# -------------------" << endl;
  cout << endl;

  bool fit_rough = true;
  int skip_passes = 10;
  //fit.set_default_roughness( 6.335 );

  penalty();

  B = X/10; //effective temperature

  ofstream penalty_file;

  penalty_file.open( outdir + "penalty.dat" );
  penalty_file << "#Iterations\tTemperature\tChi2\tType\tHistogram\tConnect\tTotal" << endl;

  for( unsigned int p = 0; p < npasses; p++ ) {

    int c = 0;
    int attempts = 0;

    cout << "# PASS " << p << endl;
    for( unsigned int l = 0; l < loops_per_pass; l++ ) {
      //cout << "#LOOP " << l << endl;
      //cout << "# Loops " << l << "/" << loops_per_pass << flush;
      do {

        attempts++;

        move( l );

        decreasing_p = ( P < P_old );
        double tmp = rng.in_range2(0.,1.);
        metropolis   = ( exp( - (P - P_old)/B ) > tmp );
        accept = ( decreasing_p || metropolis );

        if( !accept ) {
          reject_move();
        }

      } while( accept == false );
      //cout << "Loop done" << endl;

      int diff = nalkyl + nmethyl + nhead - insertion;

      //penalty_file << nalkyl << " " << nmethyl << " " << nhead << " " << diff << " " << T << endl;

      penalty_file << iterations << "\t" << B << "\t" << X << "\t" << T << "\t" << H << "\t" << C << "\t" << P << endl;
      iterations++;

      //fit.fit_intensity( nd.get_alpha_buffer(), beta.get_buffer(), rad, scale_factor, harmonics_order );
    }

    scale_tmp = rad[0][1]/intensity[0];
    scale_factor *= scale_tmp;

    //cout << fit_rough << " " << (fit_rough == true && (p > skip_passes)) << endl;
    //
    if( fit_rough && p%skip_passes == 0 && p != 0 ) {
      fit.fit_intensity( nd.get_alpha_buffer(), beta.get_buffer(), rad, scale_factor, harmonics_order );
    }
    //
    // if( fit.get_rough_chi2() <= 1. && fit.get_rough_chi2() != -1 ) {
    //   fit_rough = false;
    // }

    //cout << "# Statistics                    " << endl;
    cout << fixed << setprecision(2) << setfill('0');
    cout << setw(5) << "# Acceptance ratio:  " << (1.*loops_per_pass)/attempts << endl;
    cout << setw(5) << "# Temperature:       " << B << endl;
    cout << setw(5) << "# Chi squared:       " << X << endl;
    cout << setw(5) << "# Type penalty:      " << T << endl;
    cout << setw(5) << "# Histogram penalty: " << H << endl;
    cout << setw(5) << "# Connect penalty:   " << C << endl;
    cout << setw(5) << "# Total penalty:     " << P << endl;
    cout << setw(5) << "# Inserted beads:    " << nhead << " " << nalkyl << " " << nmethyl << endl;
    cout << setw(5) << "# I_exp[0]/I[0]:     " << setprecision(3) << scale_tmp << endl;

    if( fit_rough && p%skip_passes == 0 && p != 0 ) {
      cout << setw(5) << "# Fitted roughness:  " << fit.get_rough() << " (X^2_R = " << fit.get_rough_chi2() << ")" << endl;
    } else {
      cout << setw(5) << "# Fitted roughness:  " << fit.get_rough() << endl;
    }
    cout << setw(5) << scientific << "# Scale factor:      " << scale_factor << endl;
    cout << endl;

    string xyz = outdir + "configurations/" + to_string(p) + ".xyz";
    string pdb = outdir + "configurations/" + to_string(p) + ".pdb";
    string calc_intensity = outdir + "intensities/" + to_string(p) + ".dat";
    write_xyz( xyz );
    write_pdb( pdb );
    write_intensity( calc_intensity );

    if( B > 0.0001 ) {
      B *= 0.9;
    }

    if( (1.*loops_per_pass)/attempts < 0.011 ) {
      cout << "\n # Acceptance has reached the lower limit. Execution will be stopped." << endl;
      exit(0);
    }
  }

  penalty_file.close();

}
//------------------------------------------------------------------------------

BeadModeling::~BeadModeling() {
}
//------------------------------------------------------------------------------

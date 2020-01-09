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
#include <algorithm>

using namespace std;

BeadModeling::BeadModeling( const string& seq, const string& data, const string& out, int passes, int loops, double dm, double conn,
                            double lm, double sched, double clash, double maxd, double connd, double tr, int n0 ) {

  rad_file          = data;           //path to the experimental .rad file
  sequence_file     = seq;           //path to the file with the protein sequence
  dmax              = dm;   //dmax from fit
  npasses           = passes;   //total number of passes
  loops_per_pass    = loops;   //number of performed loops per pass
  outdir            = out;           //directory where results are saved
  lambda            = lm;
  connect           = conn;
  clash_distance    = clash; //hardcoded because experimented. Might want to leave the choice open for users though.
  max_distance      = maxd;
  conn_distance     = connd;
  t_ratio           = tr;
  schedule          = sched;
  ni0               = n0;
  convergence_temp  = 0.01;
  sequence          = "";
  shift             = 0.;
  sphere_generated  = false;
  compute_scale     = true;
  init              = true;
  with_nanodisc     = false;

  string mkdir = "mkdir " + outdir;
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "configurations/";
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "intensities/";
  system( mkdir.c_str() );

  load_FASTA(); //load sequence file
  cout << "# File '" << sequence_file << "' loaded" << endl;
  nresidues = sequence.length();

  load_rad(); //experimental file with SAXS or SANS data
  cout << "# File '" << rad_file << "' loaded " << endl;

  load_statistics(); //statistics needed for penalty function

  cout << endl;
  cout << "# SUMMARY OF PARAMETERS" << endl;
  cout << "# ---------------------" << endl;
  cout << "# Number of beads:          " << nresidues << endl;
  cout << "# Radius of initial sphere: " << dmax / 2. << endl;
  cout << "# Max number of passes:     " << npasses << endl;
  cout << "# Loops per pass:           " << loops_per_pass << endl;
  cout << "# q points I(0) estimation  " << ni0 << endl;
  cout << "# Storing results in:      '" << outdir << "'" << endl;

  nq = rad.size();
  nnnum = nnum1_ref.size();

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
  com.resize( 3 );

  for( unsigned int i = 0; i < nq; i++ ) {
    exp_q[i] = rad[i][0];
  }

  logfile();
}
//------------------------------------------------------------------------------

BeadModeling::BeadModeling( const string& seq, const string& data, const string& ft, const string& out, int passes, int loops, double dm, double conn,
                            double lm, double tm, int ins, double sched, double clash, double maxd, double connd, double tr, int n0, int ndt, double zs,
                            int qs_B, double convt, double convar ) {

  rad_file          = data;           //path to the experimental .rad file
  sequence_file     = seq;           //path to the file with the protein sequence
  best_fit          = ft;
  dmax              = dm;   //dmax from fit
  npasses           = passes;   //total number of passes
  loops_per_pass    = loops;   //number of performed loops per pass
  outdir            = out;           //directory where results are saved
  lambda            = lm;
  connect           = conn;
  insertion         = ins;
  T_strength        = tm;
  clash_distance    = clash;
  max_distance      = maxd;
  conn_distance     = connd;
  t_ratio           = tr;
  schedule          = sched;
  ni0               = n0;
  n_dtail           = ndt;
  shift             = zs;
  convergence_temp  = convt;
  convergence_accr  = convar;
  qs_to_fit         = qs_B;

  sequence          = "";
  sphere_generated  = false;
  compute_scale     = true;
  init              = true;
  with_nanodisc     = true;
  nano_model        = "flat";

  if( n_dtail > 0 ) {
    Rg = 1.93 * pow( n_dtail, 0.588 ); //R0 * N^nu, where nu is the Flory constant and R0 is obtained by fitting the behavior of many IDPs
  }

  string mkdir = "mkdir " + outdir + " >/dev/null 2>&1";
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "configurations/ >/dev/null 2>&1";
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "intensities/ >/dev/null 2>&1";
  system( mkdir.c_str() );

  load_FASTA(); //load sequence file
  //cout << "# File '" << sequence_file << "' loaded" << endl;
  nresidues = sequence.length();

  nd.load_input_flat( best_fit );
  //cout << "# File '" << best_fit << "' loaded" << endl;

  load_rad(); //experimental file with SAXS or SANS data
  //cout << "# File '" << rad_file << "' loaded " << endl;

  load_statistics(); //statistics needed for penalty function

  // cout << endl;
  // cout << "# SUMMARY OF PARAMETERS" << endl;
  // cout << "# ---------------------" << endl;
  // cout << "# Number of beads:          " << nresidues << endl;
  // cout << "# Radius of initial sphere: " << dmax / 2. << endl;
  // cout << "# Tail Rg:                  " << Rg << endl;
  // cout << "# Max number of passes:     " << npasses << endl;
  // cout << "# q points I(0) estimation  " << ni0 << endl;
  // cout << "# Loops per pass:           " << loops_per_pass << endl;

  nq = rad.size();
  nnnum = nnum1_ref.size();
  ref_mat_sum = 0.;

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

  ref_cmap.resize( insertion );
  for( auto &it : ref_cmap ) it.resize( insertion );

  for( unsigned int i = 0; i < nq; i++ ) {
    exp_q[i] = rad[i][0];
  }

  fit.fit_background( rad, qs_to_fit ); //ask the user for the number of points!
  fit.set_default_roughness( nd.get_xrough() ); //read it from file!!

  cout << "# PRELIMINARIES" << endl;
  cout << "# -------------" << endl;
  //cout << "# Background:               " << fit.get_background() << endl;
  cout << "# Results folder:          " << outdir << endl;
  //cout << "# Parameters summary:      " << outdir << "/parameters.log" << endl;

  logfile();
}
//------------------------------------------------------------------------------

void BeadModeling::load_rad() {
  rad = load_matrix( rad_file, 3 );
}
//------------------------------------------------------------------------------

void BeadModeling::logfile() {

  ofstream log;
  log.open(outdir + "parameters.log");

  log << "##########################" << endl;
  log << "# MARBLES RUN PARAMETERS #" << endl;
  log << "##########################" << endl;
  log << endl;

  log << "# INPUT/OUTPUT" << endl;
  log << "#-------------" << endl;
  log << "# Protein sequence:            " << sequence_file << endl;
  log << "# Experimental SAXS intensity: " << rad_file << endl;

  if( with_nanodisc ) {
    log << "# Nanodisc best fit file:      " << best_fit << endl;
  }

  log << "# Storing results in:          " << outdir << endl;

  log << endl;
  log << "# PROTEIN OPTIONS" << endl;
  log << "#----------------" << endl;
  log << "# Number of beads:             " << nresidues << endl;
  log << "# Radius of initial sphere:    " << dmax / 2. << endl;
  log << "# Bead clash distance:         " << clash_distance << endl;
  if( with_nanodisc ) {
    log << "# z-Shift:                     " << shift << endl;
    log << "# Initial optimization:        " << endl;
  }


  log << endl;
  if( with_nanodisc ) {
    log << "# NANODISC OPTIONS" << endl;
    log << "#----------------" << endl;
    log << "# Nanodisc model:              " << nano_model << endl;
    log << "# Inserted residues:           " << insertion << endl;
    log << "# Background estimated with:   " << qs_to_fit << " high-q points" << endl;
  }

  log << endl;
  log << "# PENALTY FUNCTION" << endl;
  log << "#-----------------" << endl;
  log << "# Connect distance:            " << conn_distance << endl;
  log << "# Connectivity strength:       " << connect << endl;
  log << "# Neighbour strength:          " << lambda << endl;

  if( with_nanodisc ) {
    log << "# Insertion strength:          " << T_strength << endl;
  }

  log << endl;
  log << "# MONTE CARLO OPTIONS" << endl;
  log << "#--------------------" << endl;
  log << "# Number of passes:            " << npasses << endl;
  log << "# Loops per pass:              " << loops_per_pass << endl;
  log << "# Maximum move length:         " << max_distance << endl;
  log << "# Initial temperature:         " << "X2/" << t_ratio << endl;
  log << "# Scheduling:                  " << schedule << endl;
  log << "# Convergence temperature:     " << convergence_temp << endl;
  log << "# q points for I(0) estimation " << ni0 << endl;

  log.close();
}
//------------------------------------------------------------------------------

void BeadModeling::load_statistics() {

  string ndist_file = "include/statistics/ndist.dat";
  //string ndist_file = "include/statistics/neigh.dat";
  string nnum1_file = "include/statistics/nnum_5.3.dat";
  string nnum2_file = "include/statistics/nnum_6.8.dat";
  string nnum3_file = "include/statistics/nnum_8.3.dat";

  ndist_ref = load_matrix( ndist_file, 3 );
  nnum1_ref = load_matrix( nnum1_file, 3 );
  nnum2_ref = load_matrix( nnum2_file, 3 );
  nnum3_ref = load_matrix( nnum3_file, 3 );

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

  double x, y, z;

  if( j == -1 ) {
    x = beads[i].x;
    y = beads[i].y;
    z = beads[i].z;
  } else {
    x = beads[i].x - beads[j].x;
    y = beads[i].y - beads[j].y;
    z = beads[i].z - beads[j].z;
  }

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

void BeadModeling::initial_configuration( double x0, double y0 ) {

    double x, y, z, r, r2;
    bool clash;

    r = dmax / 2.; /** radius of the sphere **/
    r2 = r * r;

    for( unsigned int i = 0; i < nresidues; i++ ) {
      beads[i].assign_volume_and_scattlen( string(1, sequence[i]) );

      do {

          clash = false;

          do {

            x = rng.in_range( -r, r );
            y = rng.in_range( -r, r );
            z = rng.in_range( -r, r );
            beads[i].assign_position( x + x0, y + y0, z + shift );

          } while( x*x +  y*y + z*z > r2 ); // condition that defines a sphere

          clash = bead_clash( i );

      } while( clash == true );
    }

}
//------------------------------------------------------------------------------

void BeadModeling::write_pdb( const string& filename ) {

    FILE *fil;

    fil = fopen( filename.c_str(), "w" );

    for( unsigned int i = 0; i < nresidues; i++ ) {
        const char *c = beads[i].res.c_str();
        fprintf(fil,"ATOM   %4d  CA  %s   %4d   %8.3lf%8.3lf%8.3lf  1.00 18.20           N\n",i+1,c,i+1,beads[i].x, beads[i].y,beads[i].z);
    }

    fclose(fil);
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

void BeadModeling::write_statistics( vector<double> stat, const string& filename ) {

  ofstream stat_file;

  stat_file.open(filename);
  if( stat_file.is_open() ) {
    for( unsigned int i = 0; i < stat.size(); i++ ) {
      stat_file << stat[i] << endl;
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::update_rho( int i ) {

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

  if( nano_model == "endcaps" ) {
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

  } else if( nano_model == "flat" ) {

    //the nanodisc is symmetric with respect the horizontal axis of the methyl layer

      if( fz < hmethyl * .5 && inside_ellipse( i, radius_major-3., radius_minor-3. ) ) {
        beads[i].type = 3;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_methyl;
        nmethyl++;
      } else if( fz < hcore * .5 && inside_ellipse( i, radius_major-3., radius_minor-3. ) ) {
        beads[i].type = 2;
        beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_alkyl;
        nalkyl++;
      } else if( fz < hlipid * .5 && inside_ellipse( i, radius_major-3., radius_minor-3. ) ) {
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

      if( sign >= 0 ) {
        beta.add( index, l, m, tmp );
      } else {
        beta.add( index, l, m, -tmp );
      }
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::only_prot_intensity() {

  double r, q, tmp, exponent, I0;
  double e_scattlen = nd.get_e_scatt_len();

  fill(intensity.begin(),intensity.end(),0);

  for( int i = 0; i < nq; i++ ) {
    for(int l = 0; l <= harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp = abs( beta.at( i, l, m ) );
        tmp *= tmp;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
      }
    }

    if( compute_scale ) {
      I0 = intensity[0] * e_scattlen * e_scattlen;

      double mean = 0;
      for( int i = 0; i < ni0; i++ ) {
        mean += rad[i][1]/ni0;
      }

      scale_factor = mean/I0;
    }
    compute_scale = false;
    intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor;
  }

}
//------------------------------------------------------------------------------

void BeadModeling::calc_intensity( vector<double> exp_q ) {

  double xr = fit.get_rough(); //4.609096 /* for only nanodisc benchmark */
  double background = fit.get_background(); //0.000346; /* for only nanodisc benchmark */
  double r, q, tmp, exponent, I0;
  double e_scattlen = nd.get_e_scatt_len();

  fill(intensity.begin(),intensity.end(),0);

  for( int i = 0; i < nq; i++ ) {
    q = exp_q[i];
    exponent = xr * q * xr * q;
    r = exp( - exponent / 2. );

    for(int l = 0; l <= harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp = abs( r * nd.get_alpha( i, l, m ) + beta.at( i, l, m ) );
        tmp *= tmp;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
      }
    }

    if( compute_scale ) {
      I0 = intensity[0] * e_scattlen * e_scattlen;

      double mean = 0;
      for( int i = 0; i < ni0; i++ ) {
        mean += rad[i][1]/ni0;
      }

      scale_factor = mean/I0;
    }

    compute_scale = false;
    intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor + background;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::calc_intensity_wcoil( vector<double> exp_q ) {

  double xr = fit.get_rough(); //4.609096 /* for only nanodisc benchmark */
  double background = fit.get_background(); //0.000346; /* for only nanodisc benchmark */
  double r, q, tmp, tmp2, exponent, I0, xf;
  double e_scattlen = nd.get_e_scatt_len();
  double rho = 0.14104;
  vector<double> coil( intensity.size() );

  fill(intensity.begin(),intensity.end(),0);
  fill(coil.begin(),coil.end(),0);

  for( int i = 0; i < nq; i++ ) {
    q = exp_q[i];
    exponent = xr * q * xr * q;
    r = exp( - exponent / 2. );
    xf = Rg * Rg * q * q;

    for(int l = 0; l <= harmonics_order; l++ ) {
      for(int m = 0; m <= l; m++ ) {
        tmp  = abs( r * nd.get_alpha( i, l, m ) + beta.at( i, l, m ) + nd.get_gamma(i, l, m) );
        tmp2 = abs( nd.get_gamma(i, l, m) );

        tmp  *= tmp;
        tmp2 *= tmp2;
        intensity[i] += ( (m > 0) + 1. ) * tmp;
        coil[i] += ( (m > 0) + 1. ) * tmp2;
      }
    }

    intensity[i] = intensity[i] - coil[i] + rho * rho * 2./(xf*xf) * (exp(-xf) - 1 + xf);

    if( compute_scale ) {
      I0 = intensity[0] * e_scattlen * e_scattlen;

      double mean = 0;
      for( int i = 0; i < ni0; i++ ) {
        mean += rad[i][1]/ni0;
      }

      scale_factor = mean/I0;
    }

    compute_scale = false;
    intensity[i] = intensity[i] * e_scattlen * e_scattlen * scale_factor + background;
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

    count1 = 0;
    count2 = 0;
    count3 = 0;

    for( unsigned int j = 0; j < nresidues; j++ ) {

      d = distances.at(i,j);

      //if( d < 10. ) ndist[ int(d) ] += 1. / nresidues;
      //if( d <= 10. ) ndist[ ceil(d) ] += 1. / nresidues;
      if( d < 12. ) ndist[ int(d) ] += 1. / nresidues;
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
  X = 0.;

  for( unsigned int i = 0; i < nq; i++ ) {
    tmp = intensity[i] - rad[i][1];
    err = rad[i][2];
    X += tmp * tmp / (err * err);
  }

  X /= (nq - 1);
}
//------------------------------------------------------------------------------

void BeadModeling::helix_ref_cmap() {

  double a, b, nturns, t, tmp;
  vector<vector<double> > helix( insertion, vector<double>(3) );

  nturns = insertion / 4.; /** number of helix turns */
  a      = 2.3; /** 6A is the typical width of a helix */
  b      = 0.9; /** 2b*pi is the helix pitch: a-helix pitch is 5.3A */ //1.5 * intra_seq.size() / ( 2. * nturns * M_PI);

  for( unsigned int i = 0; i < insertion; i++ ) {
    t = 2. * i * nturns * M_PI / insertion; /** arc-length of the helix */

    helix[i][0] = a * cos(t);
    helix[i][1] = a * sin(t);
    helix[i][2] = - b * t;
  }

  for( unsigned int i = 0; i < insertion; i++ ) {
    for( unsigned int j = 0; j < insertion; j++ ) {
      for( unsigned int k = 0; k < 3; k++ ) {
        tmp += (helix[i][k] - helix[j][k]) * (helix[i][k] - helix[j][k]);
      }
      ref_cmap[i][j] = sqrt( tmp );
      ref_mat_sum += tmp / (insertion * insertion);
    }
  }
}
//------------------------------------------------------------------------------

void BeadModeling::helix_cmap() {

  vector<vector<double> > helix( insertion, vector<double>(3) );
  vector<pair<double,int> > indices_z;
  double tmp;
  int n = 0;

  D = 0;

  for( unsigned int i = 0; i < nresidues; i++ ) {
    if( beads[i].type != 0 ) {
      indices_z.push_back( make_pair(beads[i].z,i) );
      n++;
    }
  }

  if( n > insertion ) n = insertion;

  if( n > 0 ) {

    sort( indices_z.begin(), indices_z.end() );

    for( unsigned int i = 0; i < n; i++ ) {
      helix[i][0] = beads[ indices_z[i].second ].x;
      helix[i][1] = beads[ indices_z[i].second ].y;
      helix[i][2] = beads[ indices_z[i].second ].z;
    }

    for( unsigned int i = 0; i < n; i++ ) {
      for( unsigned int j = 0; j < n; j++ ) {
        for( unsigned int k = 0; k < 3; k++ ) {
          tmp += (helix[i][k] - helix[j][k]) * (helix[i][k] - helix[j][k]);
        }
        D += ( ref_cmap[i][j] - sqrt( tmp ) ) * ( ref_cmap[i][j] - sqrt( tmp ) ) / (n*n);
      }
    }
  } else {
    D = 2. * ref_mat_sum;
  }

  D *= 1e-3;
}
//------------------------------------------------------------------------------

// void BeadModeling::type_penalty() {
//
//   double sum = 1. * ( nalkyl + nmethyl + nhead );
//   double tmp = sum - insertion;
//
//   T = T_strength * tmp * tmp;
//
// }
//------------------------------------------------------------------------------

// void BeadModeling::type_penalty() {
//
//   double sum = 1. * ( nalkyl + nmethyl + nhead );
//   double tmp = sum - insertion;
//
//   if( tmp < 0 ) {
//     T = T_strength * tmp * tmp;
//   } else {
//     T = 0.;
//   }
// }
//------------------------------------------------------------------------------

void BeadModeling::type_penalty() {

  double hbil    = nd.get_hlipid();
  double hmethyl = nd.get_hmethyl();
  double halkyl  = nd.get_hcore();
  double hhead   = hbil - hmethyl - halkyl;
  int tmp1, tmp2, tmp3, tmp4;

  tmp1 = nalkyl - (int)(insertion * halkyl/hbil);

  tmp2 = nmethyl - (int)(insertion * hmethyl/hbil);
  tmp3 = nhead - (int)(insertion * hhead/hbil);
  tmp4 = nalkyl + nmethyl + nhead - insertion;

  T = T_strength * ( tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3 + tmp4 * tmp4 );
}
//------------------------------------------------------------------------------

void BeadModeling::histogram_penalty() {

  double tmp1, tmp2, tmp3, tmp4;
  double nnum_len = nnum1_ref.size();
  double ndist_len = ndist_ref.size();
  H = 0;

  for( unsigned int i = 0; i < nnum_len; i++ ) {

    tmp1 = (nnum1_ref[i][1] - nnum1[i]) / nnum1_ref[i][2];
    tmp2 = (nnum2_ref[i][1] - nnum2[i]) / nnum2_ref[i][2];
    tmp3 = (nnum3_ref[i][1] - nnum3[i]) / nnum2_ref[i][2];

    if( i < ndist_len ) {
      tmp4 = (ndist_ref[i][1] - ndist[i]) / ndist_ref[i][2];
    } else {
      tmp4 = 0.;
    }

    H += ( 1. * ( tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3 ) + tmp4 * tmp4 );
    //H += tmp4 * tmp4;
  }

  H *= lambda;
}
//------------------------------------------------------------------------------

void BeadModeling::recursive_connect( int i, int s, int *pop ) {

  for( unsigned int j = 0; j < nresidues; j++ ) {
    if( distance(i,j) < conn_distance && i != j ) {
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

  int pop[nresidues];
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

  C = connect * fabs( log( (1. * nresidues) / max ) );
}
//------------------------------------------------------------------------------

void BeadModeling::compute_com() {

  com[0] = com[1] = com[2] = 0.;

  for( int i = 0; i < nresidues; i++ ) {
    com[0] += beads[i].x / nresidues;
    com[1] += beads[i].y / nresidues;
    com[2] += beads[i].z / nresidues;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::com_penalty() {

  double r0 = 0.1 * dmax;
  double tmp;

  compute_com();

  for( int k = 0; k < 3; k++ ) {
    tmp += com[k] * com[k];
  }

  tmp = sqrt(tmp) - r0;
  M = max( 0., tmp ) * max( 0., tmp );
}
//------------------------------------------------------------------------------

void BeadModeling::compact_penalty() {

  M = 0;
  double nc;

  for( int i = 0; i < nresidues; i++ ) {
    nc = 0.;
    for( int j = 0; j < nresidues; j++ ) {
      if( i != j && distances.at(i,j) <= 3.8 ) {
        nc += 1.;
      }
    }
    M += ( 1. - (exp(-0.5*nc) - exp(-6.)) ) / nresidues;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::penalty() {

  P = 0;
  T = 0;

  histogram_penalty();
  chi_squared();
  connect_penalty();
  //compact_penalty();
  if( with_nanodisc ) type_penalty();

  P += X + H + C + T; // + M;
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
//------------------------------------------------------------------------------

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

bool BeadModeling::msp_clash( unsigned const int i ) {

  double rmax = nd.get_radius_major();
  double rmin = nd.get_radius_minor();
  double hbelt = nd.get_hbelt();
  double wbelt = nd.get_wbelt();
  double fx    = fabs( beads[i].x );
  double fy    = fabs( beads[i].y );
  double fz    = fabs( beads[i].z );

  if( fz < hbelt/2. ) {
    //bead is clashing with MSP and the move cannot be accepted
    if( fx >= rmax && fx < rmax + wbelt && fy >= rmin && fy < rmin + wbelt ) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
//------------------------------------------------------------------------------

void BeadModeling::set_T0() {

  double P0 = P, P_av = 0.;
  int nincreases = 0;

  for( unsigned int l = 0; l < loops_per_pass; l++ ) {
      move_only_protein();

      if( P > P0 ) {
        P_av += -(P - P0);
        nincreases++;
      }

      reject_move();
  }

  P_av /= nincreases;
  T0 = P_av / std::log(0.7);
}
//------------------------------------------------------------------------------

// void BeadModeling::optimize_initial_position() {
//
//   //i and j are maJor and mInor, not the indices!
//   int rmj      = (int)( nd.get_radius_major() - dmax/4. );
//   int rmi      = (int)( nd.get_radius_minor() - dmax/4. );
//   double min_X = 1e20;
//
//   int k = 1;
//
//   for( int i = 0; i < rmj; i += 5 ) {
//     for( int j = 0; j < rmi; j += 5 ) {
//
//       beta.initialize(0);
//       initial_configuration( (double)(i), (double)(j) );
//
//       for( unsigned int i = 0; i < nresidues; i++ ) {
//         update_rho( i );
//       }
//
//       for( unsigned int j = 0; j < nq; j++ ) {
//         for( unsigned int i = 0; i < nresidues; i++ ) {
//           expand_sh( exp_q[j], j, i, 1, 0 );
//         }
//       }
//
//       calc_intensity( exp_q );
//       chi_squared();
//
//       if( X < min_X ) {
//         min_X = X;
//         opt_shift_x = (double)(i);
//         opt_shift_y = (double)(j);
//       }
//
//       cout << "# Initial optimization:    X^2 = " << min_X << "\r" << std::flush;
//     }
//   }
//   cout << endl;
// }
//------------------------------------------------------------------------------

void BeadModeling::optimize_initial_position() {

  //i and j are maJor and mInor, not the indices!
  int rmj      = (int)( nd.get_radius_major() - dmax/4. );
  int rmi      = (int)( nd.get_radius_minor() - dmax/4. );
  double min_X = 1e20;
  double iold = 0.;
  double jold = 0.;

  int k = 1;

  initial_configuration( 0., 0. );

  for( int i = 0; i < rmj; i += 5 ) {
    for( int j = 0; j < rmi; j += 5 ) {

      beta.initialize(0);

      for( unsigned int k = 0; k < nresidues; k++ ) {
        beads[k].x += i - iold;
        beads[k].y += j - jold;
      }

      for( unsigned int i = 0; i < nresidues; i++ ) {
        update_rho( i );
      }

      for( unsigned int j = 0; j < nq; j++ ) {
        for( unsigned int i = 0; i < nresidues; i++ ) {
          expand_sh( exp_q[j], j, i, 1, 0 );
        }
      }

      calc_intensity( exp_q );
      chi_squared();

      if( X < min_X ) {
        min_X = X;
        opt_shift_x = (double)(i);
        opt_shift_y = (double)(j);
      }

      cout << "# COM optimization:        X^2 = " << min_X << "\r" << std::flush;

      iold = i;
      jold = j;
    }
  }

  for( unsigned int i = 0; i < nresidues; i++ ) {
    beads[i].x += opt_shift_x - iold;
    beads[i].y += opt_shift_y - jold;
  }

  cout << endl;
}
//------------------------------------------------------------------------------

void BeadModeling::move_only_protein() {

  int s = 0, i, j;
  bool legal;
  save_old_config();

  double rmax, rmin, d2, z_ref;
  vector<double> vec(3);
  d2 = 3.8;//rng.in_range2( clash_distance, max_distance );

  do {

    s++;
    legal = true;
    if( s == 1 || s == 1001 ) { /*Try the same set of beads 1000 times*/
      do {
        s = 1;
        i = (int)( rng.in_range(0, nresidues) ); /*Pick a bead to be moved*/
        j = (int)( rng.in_range(0, nresidues) ); /*Pick another bead. n is to be placed in contact with m*/
      } while( i == j );
    }

    vec = rng.vector( d2 );
    beads[i].assign_position( beads[j].x + vec[0], beads[j].y + vec[1], beads[j].z + vec[2] );

    // if( distance(i,-1) > dmax/2. ) {
    //   legal = false;
    // }

    if( legal ) {
      legal = ! bead_clash( i );
    }

  } while( legal == false );

  for( unsigned int k = 0; k < nq; k++ ) {
    expand_sh( exp_q[k], k, i, -1, 0 ); //subtract the contribution of the previous position
    expand_sh( exp_q[k], k, i, 1, 0 );  //update beta with the new position
  }

  only_prot_intensity();
  distance_matrix();
  update_statistics();

  penalty();

}
//------------------------------------------------------------------------------

void BeadModeling::move( int l ) {

  int s = 0, i, j;
  bool legal;
  save_old_config();

  double rmax, rmin, d2, z_ref;
  vector<double> vec(3);
  //d2 = 3.8; //rng.in_range( clash_distance, max_distance );
  d2 = rng.in_range( clash_distance, max_distance );

  rmax = nd.get_radius_major() - 3.; //45.;//42.6;
  rmin = nd.get_radius_minor() - 3.; //32.;//29.0;

  do {

    s++;
    legal = true;
    if( s == 1 || s == 1001 ) { /*Try the same set of beads 1000 times*/
      do {
        s = 1;
        i = (int)( rng.in_range(0, nresidues) ); /*Pick a bead to be moved*/
        j = (int)( rng.in_range(0, nresidues) ); /*Pick another bead. n is to be placed in contact with m*/
      } while( i == j );
    }

    vec = rng.vector( d2 );
    beads[i].assign_position( beads[j].x + vec[0], beads[j].y + vec[1], beads[j].z + vec[2] );

    if( legal ) {
      legal = ! msp_clash( i );
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

  //if we have the disordered tail, account for it in the intensity
  if( n_dtail > 0 ) {
    calc_intensity_wcoil( exp_q );
  } else {
    calc_intensity( exp_q );
  }

  distance_matrix();
  update_statistics();

  penalty();
}
//------------------------------------------------------------------------------

void BeadModeling::SA_protein() {

  bool decreasing_p, metropolis, accept;
  int iterations = 1, rough_counter;
  double mean, sq_sum, stdev, scale_tmp, acc_ratio;
  bool relaxation = true;

  double rho_solvent = 1./3.;

  cout << endl;
  cout << "# PRELIMINARIES" << endl;
  cout << "# -------------" << endl;

  initial_configuration( 0, 0 );
  cout << "# Setup of initial configuration ..." << endl;

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

  only_prot_intensity();
  distance_matrix();
  update_statistics();
  penalty();

  T0 = X/t_ratio;

  cout << endl;
  cout << "# SIMULATED ANNEALING" << endl;
  cout << "# -------------------" << endl;
  cout << endl;

  B = T0; //effective temperature
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

        move_only_protein();

        decreasing_p = ( P < P_old );
        double tmp = rng.in_range(0.,1.);

        metropolis   = ( exp( - (P - P_old)/B ) > tmp );
        accept = ( decreasing_p || metropolis );

        if( !accept ) {
          reject_move();
        }

      } while( accept == false );

      penalty_file << iterations << "\t" << B << "\t" << X << "\t" << T << "\t" << H << "\t" << C << "\t" << P << endl;
      iterations++;
    }

    double mean = 0;
    for( int i = 0; i < ni0; i++ ) {
      mean += rad[i][1]/ni0;
    }
    scale_tmp = mean/intensity[0];
    scale_factor *= scale_tmp;

    acc_ratio = (1.*loops_per_pass)/attempts;

    cout << fixed << setprecision(2) << setfill('0');
    cout << setw(5) << "# Acceptance ratio:  " << acc_ratio << endl;
    cout << setw(5) << "# Temperature:       " << B << endl;
    cout << setw(5) << "# Chi squared:       " << X << endl;
    cout << setw(5) << "# Histogram penalty: " << H << endl;
    cout << setw(5) << "# Connect penalty:   " << C << endl;
    // cout << setw(5) << "# Compact penalty:   " << M << endl;
    cout << setw(5) << "# Total penalty:     " << P << endl;
    cout << setw(5) << "# I_exp[0]/I[0]:     " << setprecision(3) << scale_tmp << endl;
    cout << setw(5) << scientific << "# Scale factor:      " << scale_factor << endl;
    cout << endl;

    string pdb = outdir + "configurations/" + to_string(p) + ".pdb";
    string calc_intensity = outdir + "intensities/" + to_string(p) + ".dat";
    write_pdb( pdb );
    write_intensity( calc_intensity );
    // write_statistics( ndist, outdir + "ndist.dat" );
    // write_statistics( nnum1, outdir + "nnum1.dat" );
    // write_statistics( nnum2, outdir + "nnum2.dat" );
    // write_statistics( nnum3, outdir + "nnum3.dat" );

    B *= schedule;

    if( B < convergence_temp ) {
      cout << endl << "# CONVERGENCE!" << endl;
      exit(0);
    }
  }

  penalty_file.close();
}
//------------------------------------------------------------------------------

void BeadModeling::SA_nanodisc() {

  bool decreasing_p, metropolis, accept;
  int iterations = 1, rough_counter;
  double mean, sq_sum, stdev, scale_tmp;

  nd.nanodisc_form_factor( exp_q );
  if( n_dtail > 0 ) nd.gaussian_coil_form_factor( exp_q, Rg );

  optimize_initial_position();
  cout << "# Optimal sphere center:   [" << opt_shift_x << ", " << opt_shift_y << ", " << shift << "]" << endl;

  //initial_configuration( opt_shift_x, opt_shift_y );
  write_pdb( outdir + "configurations/initial.pdb" );

  nmethyl = 0;
  nalkyl = 0;
  nhead = 0;

  beta.initialize(0);

  //comment this for only nanodisc benchmark
  for( unsigned int i = 0; i < nresidues; i++ ) {
    update_rho( i );
  }

  for( unsigned int j = 0; j < nq; j++ ) {
    for( unsigned int i = 0; i < nresidues; i++ ) {
      expand_sh( exp_q[j], j, i, 1, 0 );
    }
  }

  if( n_dtail > 0 ) {
    calc_intensity_wcoil( exp_q );
  } else {
    calc_intensity( exp_q );
  }

  distance_matrix();
  update_statistics();

  cout << "# Background:              " << std::setprecision(2) << fit.get_background() << " (X^2 = " << fit.get_bck_chi2() << ")" << endl;
  cout << "# Scale factor (/1e15):    " << scale_factor/1.e15 << endl;

  cout << endl;
  cout << "# OPTIMIZATION" << endl;
  // cout << "# ------------" << endl;
  // cout << endl;

  cout << "# _______________________________________________________________________________________________" << endl;
  cout << "# Pass | Ac. R. |   Temp.   |     X     |     T     |     H     |     C     |     P     | Beads |" << endl;
  cout << "# -----------------------------------------------------------------------------------------------" << endl;

  penalty();

  B = X / t_ratio; //effective temperature

  cout << std::setw(5)  << -1 << "  |";
  cout << std::setw(6)  << std::fixed << std::setprecision(2) << "  ****  |";
  cout << std::setw(10)  << std::fixed << std::setprecision(3) << " ***** |";
  cout << std::setw(10) << std::fixed << std::setprecision(1) << X << " |";
  cout << std::setw(10) << std::fixed << std::setprecision(1) << T << " |";
  cout << std::setw(10) << std::fixed << std::setprecision(1) << H << " |";
  cout << std::setw(10) << std::fixed << std::setprecision(1) << C << " |";
  cout << std::setw(10) << std::fixed << std::setprecision(1) << P << " |";
  cout << std::setw(5)  << nhead + nalkyl + nmethyl << "  |" << endl;

  ofstream penalty_file;

  penalty_file.open( outdir + "penalty.dat" );
  penalty_file << "#Pass\tAccRatio\tTemperature\tChi2\tType\tHistogram\tConnect\tTotal\tBeads" << endl;

  for( unsigned int p = 0; p < npasses; p++ ) {

    int c = 0;
    int attempts = 0;

    //cout << "# PASS " << p << endl;
    for( unsigned int l = 0; l < loops_per_pass; l++ ) {

      do {

        attempts++;

        move( l );

        decreasing_p = ( P < P_old );
        double tmp = rng.in_range(0.,1.);
        metropolis   = ( exp( - (P - P_old)/B ) > tmp );
        accept = ( decreasing_p || metropolis );

        if( !accept ) {
          reject_move();
        }

      } while( accept == false );

      int diff = nalkyl + nmethyl + nhead - insertion;
      iterations++;
    }

    double mean = 0;
    for( int i = 0; i < ni0; i++ ) {
      mean += rad[i][1]/ni0;
    }
    scale_tmp = mean/intensity[0];
    scale_factor *= scale_tmp;

    cout << std::setw(5)  << p << "  |";
    cout << std::setw(6)  << std::fixed << std::setprecision(2) << (1.*loops_per_pass)/attempts << "  |";
    cout << std::setw(10)  << std::fixed << std::setprecision(3) << B << " |";
    cout << std::setw(10) << std::fixed << std::setprecision(1) << X << " |";
    cout << std::setw(10) << std::fixed << std::setprecision(1) << T << " |";
    cout << std::setw(10) << std::fixed << std::setprecision(1) << H << " |";
    cout << std::setw(10) << std::fixed << std::setprecision(1) << C << " |";
    cout << std::setw(10) << std::fixed << std::setprecision(1) << P << " |";
    cout << std::setw(5)  << nhead + nalkyl + nmethyl << "  |" << endl;

    penalty_file << p << "\t" << (1.*loops_per_pass)/attempts << "\t" << B << "\t" << X << "\t" << T << "\t" << H << "\t" << C << "\t" << P << "\t" <<  nhead + nalkyl + nmethyl << endl;

    if( p == (int)( 0.2 * npasses ) && C > 10 ) {
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
      cout << "# WARNING: troubles with keeping the protein connected." << endl;
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
    }
    if( p == (int)( 0.3 * npasses ) && C > 5 ) {
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
      cout << "# WARNING: the protein is likely disconnected from the nanodisc." << endl;
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
    }
    if( p > (int)( 0.5 * npasses ) && C > 5 ) {
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
      cout << "# WARNING: the protein is disconnected and at this stage it is likely irreversible. Consider restaring the simulation" << endl;
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
    }

    string pdb = outdir + "configurations/" + to_string(p) + ".pdb";
    string calc_intensity = outdir + "intensities/" + to_string(p) + ".dat";
    write_pdb( pdb );
    write_intensity( calc_intensity );

    if( B < convergence_temp && convergence_accr == 0 ) {
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
      cout << endl;
      cout << std::fixed << std::setprecision(3) << "# Convergence temperature " << convergence_temp << " has been reached." << endl;

      pdb = outdir + "/model.pdb";
      calc_intensity = outdir + "/best_fit.dat";
      write_pdb( pdb );
      write_intensity( calc_intensity );
      exit(0);

    } else if( convergence_accr > 0 && (1.*loops_per_pass)/attempts < convergence_accr ) {
      cout << "# -----------------------------------------------------------------------------------------------" << endl;
      cout << endl;
      cout << std::fixed << std::setprecision(3) << "# Convergence acceptance ratio " << convergence_accr << " has been reached." << endl;
      exit(0);
    }

    B *= schedule;
  }

  penalty_file.close();

  cout << "# -----------------------------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "# Maximum number of passes has been reached." << endl;

}
//------------------------------------------------------------------------------

BeadModeling::~BeadModeling() {
}
//------------------------------------------------------------------------------

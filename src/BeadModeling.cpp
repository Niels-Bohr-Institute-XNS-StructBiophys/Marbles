#include "BeadModeling.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ctime>
#define AV_RESIDUE_MASS 0.1 //average mass of a residue in kDa

using namespace std;

BeadModeling::BeadModeling( const string& filename ) {

  input_file = filename;
  sanity_check = false;
  sphere_generated = false;
  clash_distance = 1.8; //hardcoded because experimented
  sequence = "";
  shift = 50.;

  if( !sanity_check ) {
      load_input();
      nd.load_input( best_fit );
  }

  beads.resize( nresidues );

  sanity_check = true;
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

  ndist = load_matrix( ndist_file, 2 );
  nnum1 = load_matrix( nnum1_file, 2 );
  nnum2 = load_matrix( nnum2_file, 2 );
  nnum3 = load_matrix( nnum3_file, 2 );
}
//------------------------------------------------------------------------------

void BeadModeling::load_input() {

  ifstream file( input_file );
  string line;
  string d = " ";

  if( file.is_open() ) {

    rad_file       = parse_line( file, d );            //path to the experimental .rad file
    best_fit       = parse_line( file, d );            //path to the WillItFit results
    sequence_file  = parse_line( file, d );
    //nresidues      = stoi( parse_line( file, d ) );    //number of residues composing the protein

    // if( !nresidues ) {
    //   mass         = stoi( parse_line( file, d ) );    //molecular mass of the protein
    //   nresidues    = (int)( mass/AV_RESIDUE_MASS ); //number of residues computed from the molecular mass
    //
    //   cout << "# NOTE! You specified Residues: 0. The number of beads will be deduced from the molecular mass." << endl;
    //   cout << "# Using " << nresidues << " beads" << endl;
    // } else {
    //   cout << "# NOTE! You explicitely passed the number of residues: Mass parameter will be ignored." << endl;
    // }

    //nresidues = 495; //JUST FOR DEBUG!!!

    dmax           = stof( parse_line( file, d ) );   //dmax from fit
    npasses        = stoi( parse_line( file, d ) );   //total number of passes
    loops_per_pass = stoi( parse_line( file, d ) );   //number of performed loops per pass
    outdir         = parse_line( file, d );           //directory where results are saved
    lambda1        = stoi( parse_line( file, d ) );
    lambda2        = stoi( parse_line( file, d ) );
    connect        = stoi( parse_line( file, d ) );

    //cout << rad_file << "\t" << best_fit << "\t" << nresidues << "\t" << mass << endl;
    //cout << npasses << "\t" << loops_per_pass << "\t" << outdir << endl;
    //cout << lambda1 << "\t" << lambda2 << "\t" << connect << endl;

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

  cout << endl;
  cout << "# SUMMARY OF PARAMETERS" << endl;
  cout << "# ---------------------" << endl;
  cout << "# Number of beads: " << nresidues << endl;
  cout << "# Radius of initial sphere: " << 2. * dmax / 3. << endl;
  cout << "# Number of passes: " << npasses << endl;
  cout << "# Loops per pass: " << loops_per_pass << endl;
  cout << "# Storing results in: '" << outdir << "'" << endl;

  sanity_check = true;
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
  double x2 = x * x;

  double y = beads[i].y - beads[j].y;
  double y2 = y * y;

  double z = beads[i].z - beads[j].z;
  double z2 = z * z;

  return sqrt( x2 + y2 + z2 );
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

void BeadModeling::initial_configuration() {

  if( !sphere_generated ) {

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
            beads[i].assign_position( x, y, z + shift );

          } while( x*x +  y*y + z*z > r2 ); // condition that defines a sphere

          clash = bead_clash( i );

      } while( clash == true );
    }

    //for( unsigned int i = 0; i < nresidues; i++ ) {
    //  cout << i << " " << beads[i].x << " " << beads[i].y << " " << beads[i].z << " " << beads[i].v << " " << beads[i].rho << endl;
    //}


  } else {
    cout << "# NOTE! Skipping initial configuration because the the system is already set up.";
  }

}

void BeadModeling::test_flat() {

  int dim = rad.size();


  vector<double> exp_q( dim );
  for( int i = 0; i < dim; i++ ) {
    exp_q[i] = rad[i][0];
  }

  nd.nanodisc_form_factor( exp_q );
  update_rho();

}

// void BeadModeling::WritePDB()
// {
//     FILE * fil;
//     char* amin="rca";
//     fil=fopen("porcoddio.pdb","w");
//     for(int i=0;i<nresidues;i++){
//         fprintf(fil,"ATOM   %4d  CA  %s   %4d   %8.3lf%8.3lf%8.3lf  1.00 18.20           N\n",i+1,amin,i+1,beads[i].x, beads[i].y,beads[i].z);
//     }
//     fclose(fil);
// }

void BeadModeling::write_xyz() {

  ofstream pdb;

  pdb.open("test_file.xyz");
  if( pdb.is_open() ) {

    pdb << nresidues << endl << endl;

    for( unsigned int i = 0; i < nresidues; i++ ) {
      pdb << "CA " << beads[i].x << " " << beads[i].y << " " << beads[i].z << endl;
    }

  }
}

void BeadModeling::update_rho() {

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

  double a_endcaps               = radius_major * scale_endcaps;
  double a_endcaps_1             = 1. / a_endcaps;
  double b_endcaps               = radius_minor * scale_endcaps;
  double b_endcaps_1             = 1. / b_endcaps;
  double shift_endcaps           = - vertical_axis_ellipsoid / a_endcaps * sqrt( a_endcaps * a_endcaps - radius_major * radius_major );
  double c_endcaps_1             = 1. / vertical_axis_ellipsoid;
  double shift_z_core            = ( hcore / 2. + shift_endcaps ) * c_endcaps_1;
  double shift_z_lipid           = ( hlipid / 2. + shift_endcaps ) * c_endcaps_1;

  double x, y, z, fz, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp12;
  bool cnd1, cnd2, cnd3, cnd4;

  for( unsigned int i = 0; i < nresidues; i++ ) {

    x = beads[i].x;
    y = beads[i].y;
    z = beads[i].z;
    fz = fabs(z);

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

    if( fz < hmethyl * .5 ) {
      beads[i].type = 3;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_methyl;
    } else if( cnd1 || cnd2 || fz < hcore * .5 ) {
      beads[i].type = 2;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_alkyl;
    } else if( cnd3 || cnd4 || fz < hlipid * .5 ) {
      beads[i].type = 1;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_head;
    } else {
      beads[i].type = 0;
      beads[i].rho_modified = beads[i].rho - beads[i].v * cvprotein * rho_solvent;
    }

    cout << beads[i].rho_modified << endl;

  }
}

BeadModeling::~BeadModeling() {
}

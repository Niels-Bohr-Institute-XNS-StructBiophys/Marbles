#include "BeadModeling.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#define AV_RESIDUE_MASS 0.1 //average mass of a residue in kDa

using namespace std;

BeadModeling::BeadModeling( const string& filename ) {

  input_file = filename;
  sanity_check = false;
  sphere_generated = false;
  clash_distance = 1.8; //hardcoded because experimented
  sequence = "";
  //dmax = 60.;
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

float BeadModeling::distance( unsigned const int i, unsigned const int j ) {

  float x = beads[i].x - beads[j].x;
  float x2 = x * x;

  float y = beads[i].y - beads[j].y;
  float y2 = y * y;

  float z = beads[i].z - beads[j].z;
  float z2 = z * z;

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

    float x, y, z, r, r2;
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

BeadModeling::~BeadModeling() {
}

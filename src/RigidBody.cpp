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

#include "RigidBody.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <ctime>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <iomanip>
#include <numeric>
#include <algorithm>

using namespace std;

RigidBody::RigidBody( const string& data, const string& pdb, const string& out, int passes, double sched, double clash, double maxd, double tr ) {

  rad_file          = data;    //path to the experimental .rad file
  pdb_file          = pdb;     //path to the file with the protein sequence
  outdir            = out;     //directory where results are saved
  npasses           = passes;  //total number of passes

  clash_distance    = clash;   //clash distance
  max_distance      = maxd;    //maximum distance allowed between neighbours
  t_ratio           = tr;      //chi2 ratio to define initial temperature
  schedule          = sched;   //simulated annealing schedule

  convergence_temp  = 0.1;     //temperature at which simulation is considered convergent
  insertion         = 14;      //number of beads inserted into the nanodisc
  T_strength        = 5;       //strength of insetion penalty function
  compute_scale     = true;    //compute intensity scaling factor

  string mkdir = "mkdir " + outdir;
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "configurations/";
  system( mkdir.c_str() );

  mkdir = "mkdir " + outdir + "intensities/";
  system( mkdir.c_str() );

  load_rad(); //experimental file with SAXS or SANS data
  cout << "# File '" << rad_file << "' loaded " << endl;

  load_pdb(); //load protein configuration
  //cout << "# File '" << sequence_file << "' loaded" << endl;
  exit(-1);

  //nresidues = sequence.length();

  // load_rad(); //experimental file with SAXS or SANS data
  // cout << "# File '" << rad_file << "' loaded " << endl;
  //
  // cout << endl;
  // cout << "# SUMMARY OF PARAMETERS" << endl;
  // cout << "# ---------------------" << endl;
  // cout << "# Number of beads:          " << nresidues << endl;
  // cout << "# Radius of initial sphere: " << dmax / 2. << endl;
  // cout << "# Max number of passes:     " << npasses << endl;
  // cout << "# Loops per pass:           " << loops_per_pass << endl;
  // cout << "# Storing results in:      '" << outdir << "'" << endl;
  //
  // nq = rad.size();
  // nnnum = nnum1_ref.size();
  //
  // beads.resize( nresidues );
  // exp_q.resize( nq );
  // beta.resize_width( nq );
  // beta.initialize(0);
  // beta_old.resize_width( nq );
  // distances.resize( nresidues, nresidues );
  // distances_old.resize( nresidues, nresidues );
  // distances.initialize(0);
  // intensity.resize( nq );
  // intensity_old.resize( nq );
  //
  // for( unsigned int i = 0; i < nq; i++ ) {
  //   exp_q[i] = rad[i][0];
  // }
  //
  // logfile();
}
//------------------------------------------------------------------------------

void RigidBody::load_rad() {
  rad = load_matrix( rad_file, 3 );
}
//------------------------------------------------------------------------------

// void BeadModeling::logfile() {
//
//   ofstream log;
//   log.open(outdir + "parameters.log");
//
//   log << "# Number of beads:          " << nresidues << endl;
//   log << "# Radius of initial sphere: " << dmax / 2. << endl;
//   log << "# Number of passes:         " << npasses << endl;
//   log << "# Loops per pass:           " << loops_per_pass << endl;
//   log << "# Bead clash distance:      " << clash_distance << endl;
//   log << "# Maximum move distance:    " << max_distance << endl;
//   log << "# Connect distance:         " << conn_distance << endl;
//   log << "# Connectivity strength:    " << connect << endl;
//   log << "# Distributions strength:   " << lambda << endl;
//   log << "# Initial temperature:      " << "X2/10" << endl;
//   log << "# SA scheduling:            " << schedule << endl;
//   log << "# Storing results in:       " << outdir << endl;
//   log << "# Convergence at:           " << convergence_temp << endl;
//
//   log.close();
//
// }
//------------------------------------------------------------------------------

vector<double> RigidBody::atom_vol_scattlen( const string& atom_type ) {

  vector<double> vs(2);

  if( atom_type == "C" ) {
    vs[0] = 16.44;
    vs[1] =  6.;
  } else if( atom_type == "H" ) {
    vs[0] = 5.15;
    vs[1] = 1.;
  } else if( atom_type == "O" ) {
    vs[0] = 9.13;
    vs[1] = 8.;
  } else if( atom_type == "N" ) {
    vs[0] = 2.49;
    vs[1] = 7.;
  } else if( atom_type == "S" ) {
    vs[0] = 25.31;
    vs[1] = 16.;
  } else {
    vs[0] = 0.;
    vs[1] = 0.;
    cout << "Unrecognized atom type!" << endl;
  }

  return vs;
}
//------------------------------------------------------------------------------

void RigidBody::load_pdb() {

  ifstream file( pdb_file );

  if( file.is_open() ) { //open file
    string tmp; //will contain file line

    while( !file.eof() ) { //loop over the whole file
      getline( file, tmp ); //read file line

      if( tmp.rfind("ATOM", 0) == 0 ) { //if the word atom appears at the beginning of the string

        char cstr[tmp.size() + 1]; //allocate a char array
        vector<string> atom_tmp; //vector of string that will contain line tokens without spaces

        strcpy(cstr, tmp.c_str()); //convert line string in char
        char *token = strtok(cstr, " "); //use strtok to convert string into tokens

        while (token != NULL) { //loop until no more tokens
          string s(token); //convert token into string
          atom_tmp.push_back(s); //push string token in a vector
          token = strtok(NULL, " "); //next token
        }

        //now atom_tmp will look like:
        //ATOM 11639 HA PRO B 400 61.369 -12.358 100.981 1.00 37.60 H

        vector<string> atom;
        atom.push_back( atom_tmp[5] ); //residue number
        atom.push_back( atom_tmp[3] ); //residue type
        atom.push_back( atom_tmp[ atom_tmp.size() -1 ] ); //atom type
        atom.push_back( atom_tmp[6] ); //x
        atom.push_back( atom_tmp[7] ); //y
        atom.push_back( atom_tmp[8] ); //z

        pdb.push_back( atom );

      } else {
        continue;
      }
    }

    //compute centers of excess scattering length
    vector<double> vs;
    vector<double> cesl(3);
    string resn;

    for( int i = 0; i < pdb.size(); ) {

      double tot_esl = 0;
      resn = pdb[i][0];

       //center of excess scattering length
      fill(cesl.begin(), cesl.end(), 0);

      while( pdb[i][0] == resn ) {

        //cout << pdb[i][0] << " " << resn << " " << i << endl;

        vs = atom_vol_scattlen( pdb[i][2] );
        double esl = vs[1] - 1./3. * vs[0]; //excess scattering length

        for( int k = 0; k < 3; k++ ) {
          cesl[k] += stod(pdb[i][k+3]) * esl;
        }
        tot_esl += esl;

        i++;

        if( i >= pdb.size() ) {
          break;
        }
      }

      for( int k = 0; k < cesl.size(); k++ ) {
        cesl[k] /= tot_esl;
      }

      r.push_back( cesl );
    }

    // for( int i = 0; i < r.size(); i++ ) {
    //   cout << r[i][0] << " " << r[i][1] << " " << r[i][2] << endl;
    // }
    cout << "Number of residues: " << r.size() << endl;

  } else {
    cerr << "Cannot open " << pdb_file << endl;
  }

}
//------------------------------------------------------------------------------

RigidBody::~RigidBody() {
}
//------------------------------------------------------------------------------

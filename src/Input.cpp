#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Input.h"
#define AV_RESIDUE_MASS 0.11 //average mass of a residue in kDa

using namespace std;

string input::parse_line( ifstream& file ) {

  string str;

  getline( file, str ); //read line
  string delimiter = " "; //define delimiter
  string token = str.substr(0, str.find(delimiter)); //get the first substring before the delimeter
  str.erase(0, str.find(delimiter) + delimiter.length()); //delete it

  return str; //return substring, without delimeter and without parameter name
}
//------------------------------------------------------------------------------

vector<vector<float> > input::load_matrix( const string& filename, int ncols ) {

  ifstream file( filename );
  string line;
  int n = 0;
  vector<vector<float> > vec;

  if( file.is_open() ) {

    //get file length
    while( getline(file, line) )
      ++n;

    //go back to beginning of file
    file.clear();
    file.seekg(0, ios::beg);

    //load vector
    vec.resize( n, vector<float>(ncols) );
    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < ncols; j++ ) {
        file >> vec[i][j];
      }
    }

  } else {
    cerr << "Cannot open " << filename << endl;
    exit(-1);
  }

  file.close();

  return vec;

}
//------------------------------------------------------------------------------

void input::load_rad() {
  rad = load_matrix( rad_file, 3 );
}
//------------------------------------------------------------------------------

void input::load_statistics() {

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

input::input( const string& filename ) {

  ifstream file( filename );
  string line;

  if( file.is_open() ) {

    rad_file       = parse_line( file );            //path to the experimental .rad file
    best_fit       = parse_line( file );            //path to the WillItFit results
    nresidues      = stoi( parse_line( file ) );    //number of residues composing the protein

    if( !nresidues ) {
      mass         = stoi( parse_line( file ) );    //molecular mass of the protein
      nresidues    = (int)( mass/AV_RESIDUE_MASS ); //number of residues computed from the molecular mass

      cout << "# NOTE! You specified Residues: 0. The number of beads will be deduced from the molecular mass." << endl;
    } else {
      cout << "# NOTE! You explicitely passed the number of residues: Mass parameter will be ignored." << endl;
    }

    npasses        = stoi( parse_line( file ) );   //total number of passes
    loops_per_pass = stoi( parse_line( file ) );   //number of performed loops per pass
    outdir         = parse_line( file );           //directory where results are saved
    lambda1        = stoi( parse_line( file ) );
    lambda2        = stoi( parse_line( file ) );
    connect        = stoi( parse_line( file ) );

    //cout << rad_file << "\t" << best_fit << "\t" << nresidues << "\t" << mass << endl;
    //cout << npasses << "\t" << loops_per_pass << "\t" << outdir << endl;
    //cout << lambda1 << "\t" << lambda2 << "\t" << connect << endl;

  } else {
    cerr << "Cannot open " << filename << endl;
    exit(-1);
  }

  cout << "# '" << filename << "' loaded" << endl;

  load_rad(); //experimental file with SAXS or SANS data
  cout << "# '" << rad_file << "' loaded " << endl;

  load_statistics(); //statistics needed for penalty function
  cout << "# Statistics loaded" << endl;
}
//------------------------------------------------------------------------------

input::~input() {
}
//------------------------------------------------------------------------------

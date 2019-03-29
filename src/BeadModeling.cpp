#include "BeadModeling.h"
#define AV_RESIDUE_MASS 0.11 //average mass of a residue in kDa

using namespace std;

BeadModeling::BeadModeling( const std::string& filename ) {

  input_file = filename;
  sanity_check = false;

  if( !sanity_check ) {
      load_input();
      nd.load_input( best_fit );
  }

  sanity_check = true;

}

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
    nresidues      = stoi( parse_line( file, d ) );    //number of residues composing the protein

    if( !nresidues ) {
      mass         = stoi( parse_line( file, d ) );    //molecular mass of the protein
      nresidues    = (int)( mass/AV_RESIDUE_MASS ); //number of residues computed from the molecular mass

      cout << "# NOTE! You specified Residues: 0. The number of beads will be deduced from the molecular mass." << endl;
    } else {
      cout << "# NOTE! You explicitely passed the number of residues: Mass parameter will be ignored." << endl;
    }

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

  cout << "# '" << input_file << "' loaded" << endl;

  load_rad(); //experimental file with SAXS or SANS data
  cout << "# '" << rad_file << "' loaded " << endl;

  load_statistics(); //statistics needed for penalty function
  cout << "# Statistics loaded" << endl;

  sanity_check = true;

  //load_fit( &nd );
}
//------------------------------------------------------------------------------

BeadModeling::~BeadModeling() {
}

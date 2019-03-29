#include "Nanodisc.h"

using namespace std;

Nanodisc::Nanodisc() {
  //info about the file from which to load are not available at initialization time
  //therefore load_input() has to be called seperately
}

Nanodisc::~Nanodisc() {
}

void Nanodisc::load_input( const string& best_fit ) {

  ifstream file( best_fit );
  string line;
  string d = "=", d_ = " ";

  if( file.is_open() ) {

    skip_lines( file, 18 );
    hbelt   = stof( parse_line( file, d_ ) );
    nlipids = stof( parse_line( file, d_ ) );

    //cout << hbelt << " " << nlipids << endl;

    skip_lines( file, 1 );
    watheads = stof( parse_line( file, d_ ) );

    //cout << watheads << endl;

    skip_lines( file, 2 );
    xrough = stof( parse_line( file, d_ ) );
    cvbelt = stof( parse_line( file, d_ ) );
    cvlipids = stof( parse_line( file, d_ ) );
    cvmp = stof( parse_line( file, d_ ) );
    //cout << xrough << " " << cvbelt << " " << cvlipids << " " << cvmp << endl;

    skip_lines( file, 7 );
    cvwater = stof( parse_line( file, d_ ) );
    //cout << cvwater << endl;

    skip_lines( file, 1 );
    vertical_axis_endcaps = stof( parse_line( file, d_ ) );
    scale_endcaps = stof( parse_line( file, d_ ) );
    //cout << vertical_axis_endcaps << " " << scale_endcaps << endl;

    skip_lines( file, 31 );
    hlipid  = stof( parse_line( file, d ) );
    hcore   = stof( parse_line( file, d ) );
    hmethyl = stof( parse_line( file, d ) );

    skip_lines( file, 7 );
    radius_major = stof( parse_line( file, d ) );
    radius_minor = stof( parse_line( file, d ) );

    skip_lines( file, 2 );
    width_belt = stof( parse_line( file, d ) );

    cout << radius_major << " " << radius_minor << " " << width_belt << endl;

    skip_lines( file, 4 );
    rho_h2o = stof( parse_line( file, d_ ) );
    rho_d2o = stof( parse_line( file, d_ ) );
    rho_head = stof( parse_line( file, d_ ) );
    rho_alkyl = stof( parse_line( file, d_ ) );
    rho_methyl = stof( parse_line( file, d_ ) );
    rho_belt = stof( parse_line( file, d_ ) );

    skip_lines( file, 1 );
    rho_protein = stof( parse_line( file, d_ ) );

    cout << rho_h2o << " " << rho_d2o << " " << rho_head << " " << rho_alkyl << endl;
    cout << rho_methyl << " " << rho_belt << " " << rho_protein << endl;

  } else {
    cout << "problems" << endl;
  }

  file.close();
}

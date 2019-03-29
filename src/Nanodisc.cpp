#include "Nanodisc.h"
#include <math.h>

using namespace std;

Nanodisc::Nanodisc() {
  //info about the file from which to load are not available at initialization time
  //therefore load_input() has to be called seperately
}

Nanodisc::~Nanodisc() {
}

void Nanodisc::load_input( const string& best_fit ) {

  volume_tests = true;
  ifstream file( best_fit );
  string line;
  string d = "=", d_ = " ", d__ = ",";
  float tmp;

  if( file.is_open() ) {

    skip_lines( file, 18 );
    hbelt                 = stof( parse_line( file, d_ ) );
    nlipids               = stof( parse_line( file, d_ ) );

    skip_lines( file, 1 );
    wathead               = stof( parse_line( file, d_ ) );

    skip_lines( file, 2 );
    xrough                = stof( parse_line( file, d_ ) );
    cvbelt                = stof( parse_line( file, d_ ) );
    cvlipid               = stof( parse_line( file, d_ ) );
    cvprotein             = stof( parse_line( file, d_ ) );

    skip_lines( file, 7 );
    cvwater               = stof( parse_line( file, d_ ) );

    skip_lines( file, 1 );
    vertical_axis_endcaps = stof( parse_line( file, d_ ) );
    scale_endcaps         = stof( parse_line( file, d_ ) );

    skip_lines( file, 31 );
    hlipid                = stof( parse_line( file, d ) );
    hcore                 = stof( parse_line( file, d ) );
    hmethyl               = stof( parse_line( file, d ) );

    skip_lines( file, 7 );
    radius_major          = stof( parse_line( file, d ) );
    radius_minor          = stof( parse_line( file, d ) );

    skip_lines( file, 2 );
    width_belt            = stof( parse_line( file, d ) );

    skip_lines( file, 3 );
    rho_h2o               = stof( parse_double_delimiter( file, d_, d__ ) );
    rho_d2o               = stof( parse_double_delimiter( file, d_, d__ ) );
    rho_head              = stof( parse_double_delimiter( file, d_, d__ ) );
    rho_alkyl             = stof( parse_double_delimiter( file, d_, d__ ) );
    rho_methyl            = stof( parse_double_delimiter( file, d_, d__ ) );
    rho_belt              = stof( parse_double_delimiter( file, d_, d__ ) );

    skip_lines( file, 1 );
    rho_protein           = stof( parse_double_delimiter( file, d_, d__ ) );

    skip_lines( file, 2 );
    vh2o                  = stof( parse_double_delimiter( file, d_, d__ ) );
    vd2o                  = stof( parse_double_delimiter( file, d_, d__ ) );
    vhead                 = stof( parse_double_delimiter( file, d_, d__ ) );
    valkyl                = stof( parse_double_delimiter( file, d_, d__ ) );
    vmethyl               = stof( parse_double_delimiter( file, d_, d__ ) );
    vbelt                 = stof( parse_double_delimiter( file, d_, d__ ) );

    skip_lines( file, 1 );
    vprotein              = stof( parse_double_delimiter( file, d_, d__ ) );

    tmp = sqrt( scale_endcaps * scale_endcaps - 1. ) / scale_endcaps;
    vertical_axis_ellipsoid = vertical_axis_endcaps / ( 1. - tmp );

    //compute densities from the fit-corrected volumes
    rho_head    += wathead * rho_h2o;
    rho_head    /= ( vhead * cvlipid + wathead * vh2o );
    rho_alkyl   /= ( valkyl * cvlipid );
    rho_methyl  /= ( vmethyl * cvlipid );
    rho_belt    /= ( vbelt * cvbelt );
    rho_protein /= ( vprotein * cvprotein );
    rho_h2o     /= ( vh2o * cvwater );

    if( volume_tests ) {

      float vh1, vh2, va1, va2, vm1, vm2;

      vh1 = ( vhead * cvlipid + wathead * vh2o ) * nlipids;
      vh2 = ( hlipid - hcore ) * radius_minor * radius_major * M_PI;

      if( fabs( vh1 - vh2 ) > 0.5 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check failed: inconsistent volume of lipid heads." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }

      va1 = valkyl * cvlipid * nlipids;
      va2 = ( hcore - hmethyl ) * radius_minor * radius_major * M_PI;

      if( fabs( va1 - va2 ) > 60 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check failed: inconsistent volume of alkyl heads." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }

      vm1 = vmethyl * cvlipid * nlipids;
      vm2 = hmethyl * radius_minor * radius_major * M_PI;

      if( fabs( vm1 - vm2 ) > 0.1 ) {
        cout << "\n# ERROR IN PARSING" << endl;
        cout << "# Check failed: inconsistent volume of methyl groups." << endl;
        cout << "# This might indicate that some values have not been parsed correctly." << endl;
        cout << "# Set volume_tests to 'false' to ignore this and proceed anyway." << endl;
      }


    }

    //cout << hbelt << " " << nlipids << endl;
    //cout << watheads << endl;
    //cout << xrough << " " << cvbelt << " " << cvlipids << " " << cvmp << endl;
    //cout << cvwater << endl;
    //cout << vertical_axis_endcaps << " " << scale_endcaps << endl;
    //cout << radius_major << " " << radius_minor << " " << width_belt << endl;
    //cout << rho_h2o << " " << rho_d2o << " " << rho_head << " " << rho_alkyl << endl;
    //cout << rho_methyl << " " << rho_belt << " " << rho_protein << endl;
    //cout << v_h2o << " " << v_d2o << " " << v_head << " " << v_alkyl << endl;
    //cout << v_methyl << " " << v_belt << " " << v_protein << endl;

  } else {
    cout << "Cannot open '" << best_fit << "'" << endl;
  }

  file.close();
}

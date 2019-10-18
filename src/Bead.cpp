#include <iostream>
#include <vector>
#include <math.h>
#import "Bead.h"

using namespace std;

Bead::Bead() {
}
//------------------------------------------------------------------------------

Bead::~Bead() {
}
//------------------------------------------------------------------------------

void Bead::assign_position( double xx, double yy, double zz ) {
  x = xx;
  y = yy;
  z = zz;

  position_assigned = true;
}
//------------------------------------------------------------------------------

void Bead::assign_average( const string& sequence ) {

  vector<double> vs(2);
  double av_v = 0., av_s = 0.;

  for( unsigned int i = 0; i < sequence.size(); i++ ) {
    vs = volume_and_scattlen( string(1, sequence[i]) );
    av_v += vs[0];
    av_s += vs[1];
  }

  v   = av_v / sequence.size();
  rho = double( ceil( av_s / sequence.size() ) );
  res = "DUM";
}
//------------------------------------------------------------------------------

vector<double> Bead::volume_and_scattlen( const string& residue ) {

  vector<double> vs(2);

  if( residue == "M" || residue == "MET" ) {
      vs[0] = 163.296;
      vs[1] = 70.;
  } else if( residue == "A" || residue == "ALA" ) {
      vs[0] = 88.4593;
      vs[1] = 38.;
  } else if( residue == "L" || residue == "LEU" ) {
      vs[0] = 170.317;
      vs[1] = 62.;
  } else if( residue == "V" || residue == "VAL") {
      vs[0] = 143.031;
      vs[1] = 54.;
  } else if ( residue == "F" || residue == "PHE" ) {
      vs[0] = 210.133;
      vs[1] = 78.;
  } else if( residue == "Y" || residue == "TYR" ) {
      vs[0] = 219.449;
      vs[1] = 86.;
  } else if( residue == "G" || residue == "GLY" ) {
      vs[0] = 61.1736;
      vs[1] = 30.;
  } else if( residue == "T" || residue == "THR" ) {
      vs[0] = 125.061;
      vs[1] = 54.;
  } else if( residue == "H" || residue == "HIS" ) {
      vs[0] = 175.398;
      vs[1] = 76.;
  } else if( residue == "S" || residue == "SER" ) {
      vs[0] = 97.7757;
      vs[1] = 46.;
  } else if( residue == "K" || residue == "LYS" ) {
      vs[0] = 183.368;
      vs[1] = 71.;
  } else if( residue == "I" || residue == "ILE" ) {
      vs[0] = 170.317;
      vs[1] = 62.;
  } else if( residue == "P" || residue == "PRO" ) {
      vs[0] = 132.521;
      vs[1] = 52.;
  } else if( residue == "R" || residue == "ARG" ) {
      vs[0] = 188.449;
      vs[1] = 85.;
  } else if( residue == "D" || residue == "ASP" ) {
      vs[0] = 118.612;
      vs[1] = 59.;
  } else if( residue == "E" || residue == "GLU" ) {
      vs[0] = 145.898;
      vs[1] = 67.;
  } else if( residue == "N" || residue == "ASN" ) {
      vs[0] = 122.347;
      vs[1] = 60.;
  } else if( residue == "Q" || residue == "GLN" ) {
      vs[0] = 149.633;
      vs[1] = 68.;
  } else if( residue == "C" || residue == "CYS" ) {
      vs[0] = 108.725;
      vs[1] = 54.;
  } else if( residue == "W" || residue == "TRP" ) {
      vs[0] = 251.48;
      vs[1] = 98.;
  }

  return vs;
}
//------------------------------------------------------------------------------

void Bead::assign_volume_and_scattlen( const string& residue ) {

  if( residue == "M" || residue == "MET" ) {
      v   = 163.296; //160.03;
      rho = 70.; // 1.974e-11 / e_scatt_len;
      res = "MET";
  } else if( residue == "A" || residue == "ALA" ) {
      v   = 88.4593;//86.69;
      rho = 38.; //1.0716e-11 / e_scatt_len;
      res = "ALA";
  } else if( residue == "L" || residue == "LEU" ) {
      v   = 170.317;//166.91;
      rho = 62.; //1.7484e-11 / e_scatt_len;
      res = "LEU";
  } else if( residue == "V" || residue == "VAL") {
      v   = 143.031;//140.17;
      rho = 54.; //1.5228e-11 / e_scatt_len;
      res = "VAL";
  } else if ( residue == "F" || residue == "PHE" ) {
      v   = 210.133;//205.93;
      rho = 78.; //2.1996e-11 / e_scatt_len;
      res = "PHE";
  } else if( residue == "Y" || residue == "TYR" ) {
      v   = 219.449;//215.06;
      rho = 86.;//2.4252e-11 / e_scatt_len;
      res = "TYR";
  } else if( residue == "G" || residue == "GLY" ) {
      v   = 61.1736;//59.95;
      rho = 30.;//8.4600e-12 / e_scatt_len;
      res = "GLY";
  } else if( residue == "T" || residue == "THR" ) {
      v   = 125.061;//122.56;
      rho = 54.;//1.5228e-11 / e_scatt_len;
      res = "THR";
  } else if( residue == "H" || residue == "HIS" ) {
      v   = 175.398;//171.89;
      rho = 76.;//2.1432e-11 / e_scatt_len;
      res = "HIS";
  } else if( residue == "S" || residue == "SER" ) {
      v   = 97.7757;//95.82;
      rho = 46.;//1.2972e-11 / e_scatt_len;
      res = "SER";
  } else if( residue == "K" || residue == "LYS" ) {
      v   = 183.368;//179.7;
      rho = 71.;//2.0022e-11 / e_scatt_len;
      res = "LYS";
  } else if( residue == "I" || residue == "ILE" ) {
      v   = 170.317;//166.91;
      rho = 62.;//1.7484e-11 / e_scatt_len;
      res = "ILE";
  } else if( residue == "P" || residue == "PRO" ) {
      v   = 132.521;//129.87;
      rho = 52.;//1.4664e-11 / e_scatt_len;
      res = "PRO";
  } else if( residue == "R" || residue == "ARG" ) {
      v   = 188.449;//184.68;
      rho = 85.;//2.3970e-11 / e_scatt_len;
      res = "ARG";
  } else if( residue == "D" || residue == "ASP" ) {
      v   = 118.612;//116.24;
      rho = 59.;//1.6638e-11 / e_scatt_len;
      res = "ASP";
  } else if( residue == "E" || residue == "GLU" ) {
      v   = 145.898;//142.98;
      rho = 67.;//1.8894e-11 / e_scatt_len;
      res = "GLU";
  } else if( residue == "N" || residue == "ASN" ) {
      v   = 122.347;//119.9;
      rho = 60.;//1.6920e-11 / e_scatt_len;
      res = "ASN";
  } else if( residue == "Q" || residue == "GLN" ) {
      v   = 149.633;//146.64;
      rho = 68.;//1.9176e-11 / e_scatt_len;
      res = "GLN";
  } else if( residue == "C" || residue == "CYS" ) {
      v   = 108.725;//106.55;
      rho = 54.;//1.5228e-11 / e_scatt_len;
      res = "CYS";
  } else if( residue == "W" || residue == "TRP" ) {
      v   = 251.48;//246.45;
      rho = 98.;//2.7636e-11 / e_scatt_len;
      res = "TRP";
  } else { //non standard residue: using average values
      cout << "# NOTE! Unrecognized residue " << residue << ". Falling back to default bead." << endl;
      v   = 148.50;
      rho = 1.7921e-11 / e_scatt_len;
      res = "AVG";
  }

  selected = false;
}
//------------------------------------------------------------------------------

void Bead::save_old_config() {
  x_old    = x;
  y_old    = y;
  z_old    = z;
  rho_mold = rho_modified;
  type_old = type;
}
//------------------------------------------------------------------------------

void Bead::recover_old_config() {
  x = x_old;
  y = y_old;
  z = z_old;
  rho_modified = rho_mold;
  type = type_old;
}
//------------------------------------------------------------------------------

#include <iostream>
#import "Bead.h"

using namespace std;

Bead::Bead() {
}

Bead::~Bead() {
}

void Bead::assign_position( double xx, double yy, double zz ) {
  x = xx;
  y = yy;
  z = zz;

  position_assigned = true;
}

void Bead::assign_volume_and_scattlen( const string& residue ) {

  if( residue == "M" || residue == "MET" ) {
      v   = 160.03;
      rho = 1.974e-11 / e_scatt_len;
      res = "MET";
  } else if( residue == "A" || residue == "ALA" ) {
      v   = 86.69;
      rho = 1.0716e-11 / e_scatt_len;
      res = "ALA";
  } else if( residue == "L" || residue == "LEU" ) {
      v   = 166.91;
      rho = 1.7484e-11 / e_scatt_len;
      res = "LEU";
  } else if( residue == "V" || residue == "VAL") {
      v   = 140.17;
      rho = 1.5228e-11 / e_scatt_len;
      res = "VAL";
  } else if ( residue == "F" || residue == "PHE" ) {
      v   = 205.93;
      rho = 2.1996e-11 / e_scatt_len;
      res = "PHE";
  } else if( residue == "Y" || residue == "TYR" ) {
      v   = 215.06;
      rho = 2.4252e-11 / e_scatt_len;
      res = "TYR";
  } else if( residue == "G" || residue == "GLY" ) {
      v   = 59.95;
      rho = 8.4600e-12 / e_scatt_len;
      res = "GLY";
  } else if( residue == "T" || residue == "THR" ) {
      v   = 122.56;
      rho = 1.5228e-11 / e_scatt_len;
      res = "THR";
  } else if( residue == "H" || residue == "HIS" ) {
      v   = 171.89;
      rho = 2.1432e-11 / e_scatt_len;
      res = "HIS";
  } else if( residue == "S" || residue == "SER" ) {
      v   = 95.82;
      rho = 1.2972e-11 / e_scatt_len;
      res = "SER";
  } else if( residue == "K" || residue == "LYS" ) {
      v   = 179.7;
      rho = 2.0022e-11 / e_scatt_len;
      res = "LYS";
  } else if( residue == "I" || residue == "ILE" ) {
      v   = 166.91;
      rho = 1.7484e-11 / e_scatt_len;
      res = "ILE";
  } else if( residue == "P" || residue == "PRO" ) {
      v   = 129.87;
      rho = 1.4664e-11 / e_scatt_len;
      res = "PRO";
  } else if( residue == "R" || residue == "ARG" ) {
      v   = 184.68;
      rho = 2.3970e-11 / e_scatt_len;
      res = "ARG";
  } else if( residue == "D" || residue == "ASP" ) {
      v   = 116.24;
      rho = 1.6638e-11 / e_scatt_len;
      res = "ASP";
  } else if( residue == "E" || residue == "GLU" ) {
      v   = 142.98;
      rho = 1.8894e-11 / e_scatt_len;
      res = "GLU";
  } else if( residue == "N" || residue == "ASN" ) {
      v   = 119.9;
      rho = 1.6920e-11 / e_scatt_len;
      res = "ASN";
  } else if( residue == "Q" || residue == "GLN" ) {
      v   = 146.64;
      rho = 1.9176e-11 / e_scatt_len;
      res = "GLN";
  } else if( residue == "C" || residue == "CYS" ) {
      v   = 106.55;
      rho = 1.5228e-11 / e_scatt_len;
      res = "CYS";
  } else if( residue == "W" || residue == "TRP" ) {
      v   = 246.45;
      rho = 2.7636e-11 / e_scatt_len;
      res = "TRP";
  } else { //non standard residue: using average values
      cout << "# NOTE! Unrecognized residue " << residue << ". Falling back to default bead." << endl;
      v   = 148.50;
      rho = 1.7921e-11 / e_scatt_len;
      res = "AVG";
  }

  selected = false;
}

void Bead::save_old_config() {
  x_old    = x;
  y_old    = y;
  z_old    = z;
  rho_mold = rho_modified;
  type_old = type;
}

void Bead::recover_old_config() {
  x = x_old;
  y = y_old;
  z = z_old;
  rho_modified = rho_mold;
  type = type_old;
}

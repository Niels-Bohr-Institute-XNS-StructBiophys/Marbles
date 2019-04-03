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

  if( residue == "M" ) {
      v   = 160.03;
      rho = 1.974e-11;
  } else if( residue == "A" ) {
      v   = 86.69;
      rho = 1.0716e-11;
  } else if( residue == "L" ) {
      v   = 166.91;
      rho = 1.7484e-11;
  } else if( residue == "V" ) {
      v   = 140.17;
      rho = 1.5228e-11;
  } else if ( residue == "F" ) {
      v   = 205.93;
      rho = 2.1996e-11;
  } else if( residue == "Y" ) {
      v   = 215.06;
      rho = 2.4252e-11;
  } else if( residue == "G" ) {
      v   = 59.95;
      rho = 8.4600e-12;
  } else if( residue == "T" ) {
      v   = 122.56;
      rho = 1.5228e-11;
  } else if( residue == "H" ) {
      v   = 171.89;
      rho = 2.1432e-11;
  } else if( residue == "S" ) {
      v   = 95.82;
      rho = 1.2972e-11;
  } else if( residue == "K" ) {
      v   = 179.7;
      rho = 2.0022e-11;
  } else if( residue == "I" ) {
      v   = 166.91;
      rho = 1.7484e-11;
  } else if( residue == "P" ) {
      v   = 129.87;
      rho = 1.4664e-11;
  } else if( residue == "R" ) {
      v   = 184.68;
      rho = 2.3970e-11;
  } else if( residue == "D" ) {
      v   = 116.24;
      rho = 1.6638e-11;
  } else if( residue == "E" ) {
      v   = 142.98;
      rho = 1.8894e-11;
  } else if( residue == "N" ) {
      v   = 119.9;
      rho = 1.6920e-11;
  } else if( residue == "Q" ) {
      v   = 146.64;
      rho = 1.9176e-11;
  } else if( residue == "C" ) {
      v   = 106.55;
      rho = 1.5228e-11;
  } else if( residue == "W" ) {
      v   = 246.45;
      rho = 2.7636e-11;
  } else { //non standard residue: using average values
      v   = 148.50;
      rho = 1.7921e-11;
  }
}

#include <iostream>
#include <fstream>
#include <string>
#include "BeadModeling.h"
using namespace std;

/*! \mainpage Drawing Shapes
 *
 * This project helps user to draw shapes.
 * Currently two types of shapes can be drawn:
 * - \subpage drawingRectanglePage "How to draw rectangle?"
 *
 * - \subpage drawingCirclePage "How to draw circle?"
 *
 */

/*! \page drawingRectanglePage How to draw rectangle?
 *
 * Lorem ipsum dolor sit amet
 *
 */

/*! \page drawingCirclePage How to draw circle?
 *
 * This page is about how to draw a circle.
 * Following sections describe circle:
 * - \ref groupCircleDefinition "Definition of Circle"
 * - \ref groupCircleClass "Circle Class"
 */

int main() {

  cout << endl;
  cout << "#########################" << endl;
  cout << "# BEAD MODELING ROUTINE #" << endl;
  cout << "#########################" << endl;
  cout << endl;

  const string filename = "p450.input";
  BeadModeling BD = BeadModeling( filename );
  BD.initial_configuration();
  cout << "# Initial configuration set." << endl;

  BD.write_xyz( "test.xyz" );
  BD.test_flat();

  //BD.parse_FASTA();
  //BD.WritePDB();

  cout << endl;
  return 0;
}

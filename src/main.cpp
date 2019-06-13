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

  const string filename = "BSA.input";
  BeadModeling BD = BeadModeling( filename );

  //BD.generate_toy_model( "ub_CA_1.pdb", 0.1 );
  BD.SA_only_protein();

  //BD.test_flat();

  cout << endl;
  return 0;
}

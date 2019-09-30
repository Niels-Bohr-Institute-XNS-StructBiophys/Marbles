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

 int main( int argc, char *argv[] ) {

   cout << endl;
   cout << "#########################" << endl;
   cout << "# BEAD MODELING ROUTINE #" << endl;
   cout << "#########################" << endl;
   cout << endl;

   if( argc == 2 ) {
     string filename( argv[1] );

     BeadModeling BD = BeadModeling( filename );
     BD.SA_only_protein();

     cout << endl;

   } else if( argc > 2 ) {

     string sequence( argv[1] );
     string str_nano( argv[2] );

     bool w_nanodisc;

     if( str_nano == "True" ) {
       w_nanodisc = true;
     } else {
       w_nanodisc = false;
     }

     string data( argv[3] );
     string out( argv[4] );

     double dmax = stod( argv[5] );
     int npasses = stoi( argv[6] );
     int nloops  = stoi( argv[7] );
     double conn = stod( argv[8] );
     double neig = stod( argv[9] );
     double hstr = stod( argv[10] );
     double sche = stod( argv[11] );
     int tm_fact = stoi( argv[12] );
     double clas = stod( argv[13] );
     double maxd = stod( argv[14] );
     double cond = stod( argv[15] );

     // cout << "Sequence " << sequence << endl;
     // cout << "Data " << data << endl;
     // cout << "Out dir " << out << endl;
     // cout << "npasses " << npasses << endl;
     // cout << "nloops " << nloops << endl;
     // cout << "dmax " << dmax << endl;
     //
     // cout << "C " << conn << endl;
     // cout << "H " << neig << endl;
     // cout << "H0 " << hstr << endl;
     // cout << "schedule " << sche << endl;
     // cout << "clash " << clas << endl;
     // cout << "max dist " << maxd << endl;
     // cout << "conn dist " << cond << endl;
     // cout << "temp fact " << tm_fact << endl;

     BeadModeling BD = BeadModeling( sequence, data, out, npasses, nloops, dmax,
                                   conn, neig, hstr, sche, clas, maxd, cond, tm_fact );
     BD.SA_only_protein();

   }

  //const string filename = "/Users/simone/Desktop/NBI/BM/kinase/kinase.input";
  //BeadModeling BD = BeadModeling( filename );
  //BD.SA_only_protein();
  //cout << endl;

  return 0;
}

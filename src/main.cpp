/*******************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "RigidBody.h"
using namespace std;

 int main( int argc, char *argv[] ) {

   cout << endl;
   cout << "###########" << endl;
   cout << "# MARBLES #" << endl;
   cout << "###########" << endl;
   cout << endl;

   if( argc < 2 ) {

     cout << endl;
     cout << "Usage: python marbles.py --help" << endl;
     cout << endl;
     exit(-1);

   } else {

     string sequence( argv[1] );
     string str_nano( argv[2] );

     string data( argv[3] );
     string out( argv[4] );

     double dmax = stod( argv[5] );
     int npasses = stoi( argv[6] );
     int nloops  = stoi( argv[7] );
     double conn = stod( argv[8] );
     double neig = stod( argv[9] );
     double sche = stod( argv[10] );
     double tm_fact = stod( argv[11] );
     double clas = stod( argv[12] );
     double maxd = stod( argv[13] );
     double cond = stod( argv[14] );
     int qs_I0   = stoi( argv[18] );

     if( str_nano == "True" ) {

       string fit( argv[15] );
       int ins_res  = stoi( argv[16] );
       double t_str = stod( argv[17] );
       int n_dtail  = stoi( argv[19] );
       dhouble zs   = stoi( argv[20] );

       // cout << "Sequence  " << sequence << endl;
       // cout << "Data      " << data << endl;
       // cout << "Out dir   " << out << endl;
       // cout << "npasses   " << npasses << endl;
       // cout << "nloops    " << nloops << endl;
       // cout << "dmax      " << dmax << endl;
       // cout << "C         " << conn << endl;
       // cout << "H         " << neig << endl;
       // cout << "schedule  " << sche << endl;
       // cout << "clash     " << clas << endl;
       // cout << "max dist  " << maxd << endl;
       // cout << "conn dist " << cond << endl;
       // cout << "temp fact " << tm_fact << endl;
       // cout << "fit       " << fit << endl;
       // cout << "residues  " << ins_res << endl;
       // cout << "insertion " << t_str << endl;
       // exit(-1);

       BeadModeling BD = BeadModeling( sequence, data, fit, out, npasses, nloops, dmax,
                                     conn, neig, t_str, ins_res, sche, clas, maxd, cond, tm_fact, qs_I0, n_dtail, zs );
       BD.SA_nanodisc();

     } else {
       BeadModeling BD = BeadModeling( sequence, data, out, npasses, nloops, dmax,
                                     conn, neig, sche, clas, maxd, cond, tm_fact, qs_I0 );
       BD.SA_protein();
     }

   }

  return 0;
}

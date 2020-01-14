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
  cout << "##################################" << endl;
  cout << "#          MARBLES v.0.1         #" << endl;
  cout << "#                                #" << endl;
  cout << "# A software for the estimation  #" << endl;
  cout << "# of shapes of membrane proteins #" << endl;
  cout << "# inserted in membrane nanodiscs #" << endl;
  cout << "##################################" << endl;
  cout << endl;

  if( argc < 2 ) {

    cout << endl;
    cout << "Usage: python marbles.py --help" << endl;
    cout << endl;
    exit(-1);

  } else {

    string sequence( argv[1] );
    string data( argv[2] );
    string out( argv[3] );
    string fit( argv[14] );
    string sinfo( argv[25] );

    double dmax    = stod( argv[4]  );
    int npasses    = stoi( argv[5]  );
    int nloops     = stoi( argv[6]  );
    double conn    = stod( argv[7]  );
    double neig    = stod( argv[8]  );
    double sche    = stod( argv[9]  );
    double tm_fact = stod( argv[10] );
    double clas    = stod( argv[11] );
    double maxd    = stod( argv[12] );
    double cond    = stod( argv[13] );
    int ins_res    = stoi( argv[15] );
    double t_str   = stod( argv[16] );
    int qs_I0      = stoi( argv[17] );
    int n_dtail    = stoi( argv[18] );
    double zs      = stoi( argv[19] );
    int qs_b       = stoi( argv[20] );
    double convt   = stod( argv[21] );
    double convar  = stod( argv[22] );
    int istride    = stoi( argv[23] );
    int cstride    = stoi( argv[24] );

    // cout << "Sequence: " << sequence << endl;
    // cout << "Data:     " << data << endl;
    // cout << "Output    " << out << endl;
    // cout << "Dmax      " << dmax << endl;
    // cout << "Npasses   " << npasses << endl;
    // cout << "Nloops    " << nloops << endl;
    // cout << "Connect   " << conn << endl;
    // cout << "Neig      " << neig << endl;
    // cout << "Schedule  " << sche << endl;
    // cout << "Temp fact " << tm_fact << endl;
    // cout << "Clash     " << clas << endl;
    // cout << "Max dist  " << maxd << endl;
    // cout << "Conn dist " << cond << endl;
    // cout << "Insertion " << ins_res << endl;
    // cout << "Ins stren " << t_str << endl;
    // cout << "qs I0     " << qs_I0 << endl;
    // cout << "Tail      " << n_dtail << endl;
    // cout << "Shift     " << zs << endl;
    // cout << "qs B      " << qs_b << endl;
    // cout << "Conv temp " << convt << endl;
    // cout << "Conv acr  " << convar << endl;
    // cout << "Istride   " << istride << endl;
    // cout << "Cstride   " << cstride << endl;
    // exit(-1);

    BeadModeling BD = BeadModeling( sequence, data, fit, sinfo, out, npasses, nloops, dmax,
                                    conn, neig, t_str, ins_res, sche, clas, maxd, cond,
                                    tm_fact, qs_I0, n_dtail, zs, qs_b, convt, convar, istride, cstride );
    BD.SA_nanodisc();

  }

  return 0;
}

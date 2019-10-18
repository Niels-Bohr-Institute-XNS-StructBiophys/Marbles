#include <iostream>
#include <fstream>
#include <string>
#include "RigidBody.h"
using namespace std;

 int main( int argc, char *argv[] ) {

   cout << endl;
   cout << "#########################" << endl;
   cout << "# BEAD MODELING ROUTINE #" << endl;
   cout << "#########################" << endl;
   cout << endl;

   // RigidBody RB = RigidBody( "../rigid_body/kinase.dat", "../rigid_body/kin.pdb", "../rigid_body/tmp/", 100, 0.9, 1.8, 5.1, 10 );
   //
   // exit(-1);

   if( argc < 2 ) {

     cout << endl;
     cout << "Usage: python beads.py --help" << endl;
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
     int tm_fact = stoi( argv[11] );
     double clas = stod( argv[12] );
     double maxd = stod( argv[13] );
     double cond = stod( argv[14] );

     if( str_nano == "True" ) {

       string fit( argv[15] );
       int ins_res = stoi( argv[16] );
       double t_str = stod( argv[17] );

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
                                     conn, neig, t_str, ins_res, sche, clas, maxd, cond, tm_fact );
       BD.SA_nanodisc();

     } else {
       BeadModeling BD = BeadModeling( sequence, data, out, npasses, nloops, dmax,
                                     conn, neig, sche, clas, maxd, cond, tm_fact );
       BD.SA_protein();
     }

   }

  return 0;
}

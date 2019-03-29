#include <iostream>
#include <fstream>
#include <string>
#include "BeadModeling.h"
using namespace std;

int main() {

  cout << endl;
  cout << "#########################" << endl;
  cout << "# BEAD MODELING ROUTINE #" << endl;
  cout << "#########################" << endl;
  cout << endl;

  const string filename = "p450.input";
  BeadModeling BD = BeadModeling( filename );
  //input in = input( filename );

  cout << endl;
  return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include "Input.h"
using namespace std;

int main() {

  cout << endl;
  cout << "#########################" << endl;
  cout << "# BEAD MODELING ROUTINE #" << endl;
  cout << "#########################" << endl;
  cout << endl;

  const string filename = "p450.input";
  input in = input( filename );

  cout << endl;
  return 0;
}

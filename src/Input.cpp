#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Input.h"

using namespace std;

string Input::parse_line( ifstream& file, string delimiter ) {

  string str;

  getline( file, str ); //read line
  string token = str.substr(0, str.find(delimiter)); //get the first substring before the delimeter
  str.erase(0, str.find(delimiter) + delimiter.length()); //delete it

  return str; //return substring, without delimeter and without parameter name
}
//------------------------------------------------------------------------------

void Input::skip_lines( ifstream& file, int num ) {

  string str;
  for( int i = 0; i < num; i++ ) {
    getline( file, str );
  }

}
//------------------------------------------------------------------------------

vector<vector<float> > Input::load_matrix( const string& filename, int ncols ) {

  ifstream file( filename );
  string line;
  int n = 0;
  vector<vector<float> > vec;

  if( file.is_open() ) {

    //get file length
    while( getline(file, line) )
      ++n;

    //go back to beginning of file
    file.clear();
    file.seekg(0, ios::beg);

    //load vector
    vec.resize( n, vector<float>(ncols) );
    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < ncols; j++ ) {
        file >> vec[i][j];
      }
    }

  } else {
    cerr << "Cannot open " << filename << endl;
    exit(-1);
  }

  file.close();

  return vec;
}
//------------------------------------------------------------------------------

Input::Input() {
}
//------------------------------------------------------------------------------

Input::~Input() {
}
//------------------------------------------------------------------------------

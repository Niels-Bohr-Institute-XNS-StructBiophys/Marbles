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
#include <vector>
#include "Input.h"

using namespace std;

Input::Input() {
}
//------------------------------------------------------------------------------

string Input::parse_line( ifstream& file, string delimiter ) {

  /**
   * parses a numeric value from file line with single delimiter.
   *
   * Example
   * -------
   * from file_line = "CVLipid 50", the function parses 50 using " " as a delimiter
   *
   * Inputs
   * ------
   * ifstream& file:   address of a file stream
   * string delimiter: delimiter
   *
   * Returns
   * -------
   * string str: parsed output
   */

  string str;

  getline( file, str ); //read line
  string token = str.substr(0, str.find(delimiter)); //get the first substring before the delimeter
  str.erase(0, str.find(delimiter) + delimiter.length()); //delete it

  return str; //return substring, without delimeter and without parameter name
}
//------------------------------------------------------------------------------

string Input::parse_double_delimiter( ifstream& file, string delimiter1, string delimiter2 ) {

  /**
   * parses a numeric value from file line with double delimiter.
   *
   * Example
   * -------
   * from file_line = "!nu_X_h20 50 , comment", the function parses 50 using " " and "," as a delimiters
   *
   * Inputs
   * ------
   * ifstream& file:    address of a file stream
   * string delimiter1: first delimiter
   * string delimiter2: second delimiter
   *
   * Returns
   * -------
   * string token: parsed output
   */

  string str;

  getline( file, str ); //read line

  //first substring
  string token = str.substr(0, str.find(delimiter1)); //get the first substring before the delimeter
  str.erase(0, str.find(delimiter1) + delimiter1.length()); //delete it

  //second substring
  token = str.substr(0, str.find(delimiter2));

  return token; //return substring, without delimeter and without parameter name


}
//------------------------------------------------------------------------------

void Input::skip_lines( ifstream& file, int num ) {

  /**
   * skips a given number of lines while parsing a file
   *
   * Inputs
   * ------
   * ifstream& file: address of a file stream
   * int num:        number of lines to skip
   */

  string str;
  for( int i = 0; i < num; i++ ) {
    getline( file, str );
  }

}
//------------------------------------------------------------------------------

vector<vector<double> > Input::load_matrix( const string& filename, int ncols ) {

  /**
   * loads a matrix from file. Only the number of columns is needed, not the number of lines.
   *
   * Inputs
   * ------
   * const string& filename: path to a file
   * int ncols:              number of columns of the matrix
   *
   * Returns
   * -------
   * vector<vector<double> > vec: the loaded matrix
   */

  ifstream file( filename );
  string line;
  int n = 0;
  vector<vector<double> > vec;

  if( file.is_open() ) {

    //get file length
    while( getline(file, line) )
      ++n;

    //go back to beginning of file
    file.clear();
    file.seekg(0, ios::beg);

    //load vector
    vec.resize( n, vector<double>(ncols) );
    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < ncols; j++ ) {
        file >> vec[i][j];
      }
    }

  } else {
    cerr << "Cannot open '" << filename << "'" << endl;
    exit(-1);
  }

  file.close();

  return vec;
}
//------------------------------------------------------------------------------

vector<double> Input::load_vector_from_matrix( const string& filename, int col, int ncols ) {

  /**
   * loads a column from a matrix file. Only the number of the column is needed, not the number of lines.
   *
   * Inputs
   * ------
   * const string& filename: path to a file
   * int ncol:               the column that has to be loaded
   *
   * Returns
   * -------
   * vector<double> vec: the loaded vector
   */

   ifstream file( filename );
   string line;
   int n = 0;
   double tmp;
   //vector<vector<double> > tmp;
   vector<double> vec;

   if( file.is_open() ) {

     //get file length
     while( getline(file, line) )
       ++n;

     //go back to beginning of file
     file.clear();
     file.seekg(0, ios::beg);

     //load vector
     vec.resize( n );
     for( int i = 0; i < n; i++ ) {
       for( int j = 0; j < ncols; j++ ) {
         file >> tmp;
         if( j == col ) {
           vec[i] = tmp;
         }
       }
     }

   } else {
     cerr << "Cannot open '" << filename << "'" << endl;
     exit(-1);
   }

   file.close();

   return vec;
}

Input::~Input() {
}
//------------------------------------------------------------------------------

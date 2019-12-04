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
#include <vector>
#include <fstream>

/**
 * Class needed to interpret input files
 */

class Input {

  /* CLASSES */

  /* FLAGS */

  /* INPUT FILES */

  /* INFO VARIABLES */

  /* DETAILS OF THE CALCULATION */

  /* PRIVATE FUNCTIONS */

  public:
    Input();
    ~Input();

    /* PUBLIC UTILITIES */
    std::string parse_line( std::ifstream&, std::string );                          /** parses input line with single delimiter */
    std::string parse_double_delimiter( std::ifstream&, std::string, std::string ); /** parses input lines with double delimiters */
    std::vector<std::vector<double> > load_matrix( const std::string&, int );       /** loads a matrix from file */
    std::vector<double> load_vector_from_matrix( const std::string&, int, int );         /** loads a single column from a matrix file */
    void skip_lines( std::ifstream&, int );                                         /** skips a given number of lines when reading a file */

    /* GET FUNCTIONS */
};

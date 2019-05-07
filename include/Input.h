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

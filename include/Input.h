#include <iostream>
#include <vector>
#include <fstream>

class Input {

  /* Class needed to help interpreting input files */
  
  public:
    Input();
    ~Input();

    std::string parse_line( std::ifstream&, std::string );                   /* parses an input line */
    std::vector<std::vector<float> > load_matrix( const std::string&, int ); /* loads a matrix from file */
    void skip_lines( std::ifstream&, int );                                  /* skips a given number of lines when reading a file */
};

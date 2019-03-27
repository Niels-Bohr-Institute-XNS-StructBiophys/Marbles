#include <iostream>
#include <vector>
#include <fstream>
#include "Nanodisc.h"

class input {

  private:
    std::string rad_file; /* Experimental .rad file with SAXS data */
    std::string outdir;   /* directory where to store results */
    std::string best_fit; /* report file from WillItFit */

    unsigned int nresidues;      /* number of residues in the protein */
    unsigned int npasses;        /* number of passes to be executed */
    unsigned int loops_per_pass; /* number of loops per pass to be executed */

    float mass; /* molecular mass of the protein */
    float lambda1; /* TO BE CLEARED */
    float lambda2; /* TO BE CLEARED */
    float connect; /* TO BE CLEARED */

    std::vector<std::vector<float> > rad;   /* experimental SAXS value for different values of q */
    std::vector<std::vector<float> > ndist; /* histogram of number distances for selected ensemble */
    std::vector<std::vector<float> > nnum1; /* distribution of number neighbours for R = 5.3A */
    std::vector<std::vector<float> > nnum2; /* distribution of number neighbours for R = 6.8A */
    std::vector<std::vector<float> > nnum3; /* distribution of number neighbours for R = 8.3A */

    std::string parse_line( std::ifstream& );                                /* parses an input line */
    std::vector<std::vector<float> > load_matrix( const std::string&, int ); /* loads a matrix from file */
    void load_rad();                                                         /* loads the .rad experiment file */
    void load_statistics();                                                  /* loads the tabulated statistics files */
    void load_fit( nanodisc* );

  public:
    input( const std::string& );
    ~input();
};

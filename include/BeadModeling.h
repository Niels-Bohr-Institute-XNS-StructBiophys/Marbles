#include <iostream>
#include "Nanodisc.h"
#include "Bead.h"
#include "Random.h"

#define NH 17 //Order of harmonics
#define NTHETA ((NH+1)*2)
#define NPHI ((NH+1)*2)

/*
 * TODO:
 * find a more elegant way to allocate the 3D array
 * avoid having to specify the order of harmonics here and call it from Nanodisc
 * use Array2D for rad, ndist and nnum*
 */

class BeadModeling : public Input {

    private:

      /* CLASSES */
      Nanodisc nd;                 /** Nanodisc class for nanodisc handling */
      RandomNumbers rng;           /** RandomNumbers class fro random number generation */
      std::vector<Bead> beads;     /** Vector of Bead classes for protein beads handling */

      /* FLAGS */
      bool sanity_check;           /** Flag for sanity checks */
      bool sphere_generated;       /** Flag for avoiding regereating the initial sphere */

      /* INPUT FILES */
      std::string input_file;      /** Input file with run configurations */
      std::string rad_file;        /** Experimental .rad file with SAXS data */
      std::string outdir;          /** Directory where to store results */
      std::string best_fit;        /** Report file from WillItFit */
      std::string sequence_file;   /** FASTA sequence file of the protein */

      /* INFO VARIABLES */

      /* DETAILS OF THE CALCULATION */
      std::string sequence;        /** protein sequence */

      unsigned int nresidues;      /** Number of residues in the protein and, correspondigly, beads. */
      unsigned int npasses;        /** Number of Monte Carlo passes to be executed */
      unsigned int loops_per_pass; /** Number of loops per pass to be executed */

      double lambda1;              /** TO BE CLEARED */
      double lambda2;              /** TO BE CLEARED */
      double connect;              /** TO BE CLEARED */
      double dmax;                 /** maximum length detected from P(r) */
      double shift;                /** z shift of the initial sphere with respect to the nanodisc */
      double clash_distance;       /** distance below which a bead clash is called */

      const unsigned int harmonics_order = 17;
      const unsigned int ntheta = (harmonics_order + 1) * 2;
      const unsigned int nphi   = (harmonics_order + 1) * 2;

      std::vector<std::vector<double> > rad;   /* experimental SAXS value for different values of q */
      std::vector<std::vector<double> > ndist; /* histogram of number distances for selected ensemble */
      std::vector<std::vector<double> > nnum1; /* distribution of number neighbours for R = 5.3A */
      std::vector<std::vector<double> > nnum2; /* distribution of number neighbours for R = 6.8A */
      std::vector<std::vector<double> > nnum3; /* distribution of number neighbours for R = 8.3A */

      Array3D<std::complex<double>, 0, NH+1, NH+1> beta;

      std::vector<double> intensity;

      /* PRIVATE FUNCTIONS */
      void load_rad(); /* loads the .rad experiment file */
      void load_statistics(); /* loads the tabulated statistics files */
      void load_FASTA();
      void expand_sh( double, int, int );
      void calc_intensity( std::vector<double> );

      double distance( unsigned const int, unsigned const int ); /** measures the distance between beads **/

      bool bead_clash( unsigned const int ); /** checks wether the position of a bead clashes with another one **/

    public:
      BeadModeling( const std::string& );
      ~BeadModeling();

      /* PUBLIC UTILITIES */
      void load_input();
      void initial_configuration();
      void write_xyz();
      void test_flat();
      void update_rho();

      /* GET FUNCTIONS */

      //void WritePDB();
};

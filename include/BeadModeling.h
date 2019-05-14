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
      RandomNumbers rng;           /** RandomNumbers class for random number generation */
      std::vector<Bead> beads;     /** Vector of Bead classes for protein beads handling */

      /* FLAGS */
      bool sphere_generated;       /** Flag for avoiding regereating the initial sphere */
      bool init_type_penalty;      /** True if type_penalty is being called for the first time */
      bool init;                   /** True if penalty is called for the first time */

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
      unsigned int nalkyl;         /** Number of beads in the alkyl region of the nanodisc */
      unsigned int nmethyl;        /** Number of beads in the methyl region of the nanodisc */
      unsigned int nhead;          /** Number of beads in the head region of the nanodisc */
      unsigned int nalkyl_old;
      unsigned int nmethyl_old;
      unsigned int nhead_old;
      unsigned int insertion;      /** Number of residues that are required to be inserted in the nanodisc */
      unsigned int nq;             /** Length of the rad file, i.e. number of experimental q points */
      unsigned int nnnum;          /** Length of the nnum files */

      double lambda;              /** TO BE CLEARED */
      //double lambda2;              /** TO BE CLEARED */
      double connect;              /** TO BE CLEARED */
      double dmax;                 /** maximum length detected from P(r) */
      double shift;                /** z shift of the initial sphere with respect to the nanodisc */
      double clash_distance;       /** distance below which a bead clash is called */
      double X;                    /** chi squared of SAXS intensity */
      double T;                    /** value of the type penalty */
      double H;                    /** value of the histogram penalty */
      double C;                    /** value of the connect penalty */
      double P;                    /** value of the total penalty */
      double P_old;
      double B;
      double T_strength;           /** strength of the type penalty */
      double H_strength;           /** strength of the histogram penalty */

      const unsigned int harmonics_order = 17;
      const unsigned int ntheta = (harmonics_order + 1) * 2;
      const unsigned int nphi   = (harmonics_order + 1) * 2;

      std::vector<std::vector<double> > rad;   /* experimental SAXS value for different values of q */
      std::vector<double> exp_q;

      std::vector<double> ndist;
      std::vector<double> ndist_ref;
      std::vector<double> ndist_old;

      std::vector<double> nnum1;
      std::vector<double> nnum1_ref;
      std::vector<double> nnum1_old;

      std::vector<double> nnum2;
      std::vector<double> nnum2_ref;
      std::vector<double> nnum2_old;

      std::vector<double> nnum3;
      std::vector<double> nnum3_ref;
      std::vector<double> nnum3_old;

      Array3D<std::complex<double>, 0, NH+1, NH+1> beta;
      Array3D<std::complex<double>, 0, NH+1, NH+1> beta_old;

      Array2D<double, 1, 1> distances;
      Array2D<double, 1, 1> distances_old;

      std::vector<double> intensity;
      std::vector<double> intensity_old;

      /* PRIVATE FUNCTIONS */
      void load_rad(); /* loads the .rad experiment file */
      void load_statistics(); /* loads the tabulated statistics files */
      void load_FASTA();
      void expand_sh( double, int, int, int, int );
      void calc_intensity( std::vector<double> );
      void distance_matrix();
      void update_statistics();
      void recursive_connect( int, int, int* );
      void save_old_config();
      void move( int );
      void reject_move();
      void test_rho( int );

      //penalty functions
      void chi_squared();
      void type_penalty();
      void histogram_penalty();
      void connect_penalty();

      double distance( unsigned const int, unsigned const int ); /** measures the distance between beads **/
      double bead_distance( Bead, Bead );

      bool bead_clash( unsigned const int ); /** checks wether the position of a bead clashes with another one **/
      bool inside_ellipse( int, double, double );

    public:
      BeadModeling( const std::string& );
      ~BeadModeling();

      /* PUBLIC UTILITIES */
      void load_input();
      void initial_configuration();
      void write_xyz( const std::string& );
      void write_intensity( const std::string& );
      void test_flat();
      void update_rho( int );
      void penalty();

      /* GET FUNCTIONS */

      //void WritePDB();
};

#include <iostream>
#include "Bead.h"
#include "Random.h"
#include "Fit.h"

class BeadModeling : public Input {

    private:

      /* CLASSES */
      Nanodisc nd;                                                  /** class for nanodisc handling */
      RandomNumbers rng;                                            /** class for random number generation */
      Fit fit;                                                      /** class for fitting routines */
      std::vector<Bead> beads;                                      /** vector of Bead classes for protein beads handling */

      /* FLAGS */
      bool sphere_generated;                                        /** true if initial configuration has been generated */
      bool compute_scale;                                           /** true if intensity rescaling factor has to be computed */
      bool with_nanodisc;                                           /** true is simulated annealing is carried out in the presence of a nanodisc */

      /* INPUT FILES */
      std::string rad_file;                                         /** path to the file with experimenta SAXS data */
      std::string outdir;                                           /** path to directory where to store results */
      std::string best_fit;                                         /** path to result file from WillItFit */
      std::string sequence_file;                                    /** path to FASTA sequence file of the protein */
      std::string nano_model;                                       /** type of nanodisc model employed (with endcaps or flat) */

      /** DETAILS OF THE CALCULATION */
      std::string sequence;                                         /** protein sequence */

      unsigned int nresidues;                                       /** number of residues in the protein and, correspondigly, beads. */
      unsigned int npasses;                                         /** number of Monte Carlo passes to be executed */
      unsigned int loops_per_pass;                                  /** number of loops per pass to be executed */
      unsigned int nalkyl;                                          /** number of beads in the alkyl region of the nanodisc */
      unsigned int nmethyl;                                         /** number of beads in the methyl region of the nanodisc */
      unsigned int nhead;                                           /** number of beads in the head region of the nanodisc */
      unsigned int nalkyl_old;                                      /** number of beads in the alkyl region of the nanodisc in the previous loop */
      unsigned int nmethyl_old;                                     /** number of beads in the methyls region of the nanodisc in the previous loop */
      unsigned int nhead_old;                                       /** number of beads in the head region of the nanodisc in the previous loop */
      unsigned int insertion;                                       /** number of beads required to be inserted in the nanodisc */
      unsigned int nq;                                              /** length of the rad file, i.e. number of experimental q points */
      unsigned int nnnum;                                           /** length of the nnum files */
      static const unsigned int harmonics_order = 17;               /** number of spherical harmonics to be used in intensity estimation */
      static const unsigned int ntheta = (harmonics_order + 1) * 2; /** number of points to be used in theta integration */
      static const unsigned int nphi   = (harmonics_order + 1) * 2; /** number of points to be used in phi integration */

      double H0;                                                    /** strength of the histogram penalty function */
      double C0;                                                    /** strength of the connectivity penalty function */
      double T0;                                                    /** strength of the type penalty */
      double dmax;                                                  /** maximum length detected from P(r) */
      double shift;                                                 /** z shift of the initial sphere with respect to the nanodisc */
      double clash_distance;                                        /** distance below which a bead clash is called */
      double max_distance;                                          /** maximum distance allowed for bead MC move */
      double conn_distance;                                         /** maximum distance between beads to consider them connected */
      double t_ratio;                                               /** ratio of the chi squared to be used as initial temperature */
      double schedule;                                              /** simulated annealing schedule */
      double convergence_temp;                                      /** temperature at which convergence is achieved */
      double X;                                                     /** chi squared of SAXS intensity */
      double T;                                                     /** value of the type penalty */
      double H;                                                     /** value of the histogram penalty */
      double C;                                                     /** value of the connect penalty */
      double P;                                                     /** value of the total penalty */
      double P_old;                                                 /** value of the total penalty at previous loop */
      double B;                                                     /** effective temperature */
      double scale_factor;                                          /** intensity scale factor */

      std::vector<std::vector<double> > rad;                        /** experimental SAXS intensity for different values of q */
      std::vector<double> exp_q;                                    /** experimental q points */
      std::vector<double> ndist;                                    /** estimated histogram of neighbour distances between centers of scattering */
      std::vector<std::vector<double> > ndist_ref;                  /** referece histogram of neighbour distances between centers of scattering */
      std::vector<double> ndist_old;                                /** estimated histogram of neighbour distances between centers of scattering at previous loop */
      std::vector<double> nnum1;                                    /** estimated distribution of the number of neighbours within R = 5.3A */
      std::vector<std::vector<double> > nnum1_ref;                  /** reference distribution of the number of neighbours within R = 5.3A */
      std::vector<double> nnum1_old;                                /** estimated distribution of the number of neighbours within R = 5.3A at previous loop */
      std::vector<double> nnum2;                                    /** estimated distribution of the number of neighbours within R = 6.8A */
      std::vector<std::vector<double> > nnum2_ref;                  /** reference distribution of the number of neighbours within R = 6.8A */
      std::vector<double> nnum2_old;                                /** estimated distribution of the number of neighbours within R = 6.8A at previous loop */
      std::vector<double> nnum3;                                    /** estimated distribution of the number of neighbours within R = 8.3A */
      std::vector<std::vector<double> > nnum3_ref;                  /** reference distribution of the number of neighbours within R = 8.3A */
      std::vector<double> nnum3_old;                                /** estimated distribution of the number of neighbours within R = 8.3A at previous loop */

      Array3D<std::complex<double>, 0, harmonics_order + 1,
              harmonics_order + 1> beta;                            /** beads expansion coefficients on spherical harmonics basis */
      Array3D<std::complex<double>, 0, harmonics_order + 1,
              harmonics_order + 1> beta_old;                        /** beads expansion coefficients on spherical harmonics basis at previous loop */

      Array2D<double, 1, 1> distances;                              /** matrix of beads distances */
      Array2D<double, 1, 1> distances_old;                          /** matrix of beads distances at previous loop */

      std::vector<double> intensity;                                /** total estimated intensity */
      std::vector<double> intensity_old;                            /** total estimated intensity at previous loop */

      /* PRIVATE FUNCTIONS */
      void load_rad();                                              /** loads the experimental SAXS file */
      void load_statistics();                                       /** loads the tabulated statistics files */
      void load_FASTA();                                            /** load the protein sequence file */
      void logfile();                                               /** writes a summary log file with all the parameters used in the simulation */
      void expand_sh( double, int, int, int );                      /** computes the protein expansion coefficients */
      void calc_intensity( std::vector<double> );                   /** computes the overall SAXS intensity */
      void only_prot_intensity();                                   /** computes SAXS intensity in the absence of a nanodisc */
      void distance_matrix();                                       /** computes beads distance matrix */
      void update_statistics();                                     /** computes the estimated neighbour distributions */
      void recursive_connect( int, int, int* );                     /** computes the number of connected beads in the protein configuration */
      void save_old_config();                                       /** saves the current configuration before moving */
      void move( int );                                             /** Monte Carlo move in the prensence of nanodisc */
      void move_only_protein();                                     /** Monte Carlo move in the absence of nanodisc */
      void reject_move();                                           /** restores previous conformation when move is rejected */
      void chi_squared();                                           /** computes chi squared between estimated and experimental SAXS intensity */
      void type_penalty();                                          /** computes penalty enforcing the insertion of beads in the nanodisc */
      void histogram_penalty();                                     /** computes penalty enforcing the neighbouring distance distributions */
      void connect_penalty();                                       /** computes penalty enforcing protein connectedness */
      double distance( unsigned const int, unsigned const int );    /** measures the distance between beads */
      bool bead_clash( unsigned const int );                        /** returns true if the position of a bead clashes with another one **/
      bool inside_ellipse( int, double, double );                   /** returns true if a bead lies within the nanodisc ellipse */
      void initial_configuration();                                 /** generates initial spherical distribution of beads */
      void write_intensity( const std::string& );                   /** writes estimated intensity on file */
      void write_stat( std::vector<double>, const std::string& );   /** writes estimated neighbouring distance distributions */
      void write_pdb( const std::string& );                         /** writes pdb file with current protein shape */
      void update_rho( int );                                       /** updates excess scattering length of beads depending on their location */
      void penalty();                                               /** computes total penalty */

    public:
      BeadModeling( const std::string&, const std::string&,
                    const std::string&, int, int, double,
                    double, double, double, double, double,
                    double, double );                               /** class constructor in the absence of nanodisc */
      BeadModeling( const std::string&, const std::string&,
                    const std::string&, const std::string&,
                    int, int, double, double, double, double,
                    int, double, double, double, double, double );  /** class constructor in the presence of nanodisc */
      ~BeadModeling();                                              /** class destructor */

      /* PUBLIC UTILITIES */
      void SA_nanodisc();                                           /** runs simulated annealing in the presence of nanodisc */
      void SA_protein();                                            /** runs simulated annealing for solvated protein */

};

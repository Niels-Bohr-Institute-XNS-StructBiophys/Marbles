#include <iostream>
#include "Nanodisc.h"
#include "Bead.h"
#include "Random.h"

class BeadModeling : public Input {

    private:

      Nanodisc nd;
      RandomNumbers rng;
      std::vector<Bead> beads;

      bool sanity_check;
      bool sphere_generated;

      std::string input_file;
      std::string rad_file; /* Experimental .rad file with SAXS data */
      std::string outdir;   /* directory where to store results */
      std::string best_fit; /* report file from WillItFit */
      std::string sequence_file; /* FASTA sequence file of the protein */
      std::string sequence; /* protein sequence */

      unsigned int nresidues;      /* number of residues in the protein */
      unsigned int npasses;        /* number of passes to be executed */
      unsigned int loops_per_pass; /* number of loops per pass to be executed */

      //float mass; /* molecular mass of the protein */
      float lambda1; /* TO BE CLEARED */
      float lambda2; /* TO BE CLEARED */
      float connect; /* TO BE CLEARED */
      float dmax;
      float shift; /* z shift of the initial sphere with respect to the nanodisc */
      float clash_distance; /** distance below which a bead clash is called */

      std::vector<std::vector<float> > rad;   /* experimental SAXS value for different values of q */
      std::vector<std::vector<float> > ndist; /* histogram of number distances for selected ensemble */
      std::vector<std::vector<float> > nnum1; /* distribution of number neighbours for R = 5.3A */
      std::vector<std::vector<float> > nnum2; /* distribution of number neighbours for R = 6.8A */
      std::vector<std::vector<float> > nnum3; /* distribution of number neighbours for R = 8.3A */

      void load_rad(); /* loads the .rad experiment file */
      void load_statistics(); /* loads the tabulated statistics files */
      void load_FASTA();
      float distance( unsigned const int, unsigned const int ); /** measures the distance between beads **/
      float volume_and_scattering_dictionary( const std::string& );
      bool bead_clash( unsigned const int ); /** checks wether the position of a bead clashes with another one **/

    public:
      BeadModeling( const std::string& );
      ~BeadModeling();

      void load_input();
      void initial_configuration();
      void write_xyz();

      //void WritePDB();
};

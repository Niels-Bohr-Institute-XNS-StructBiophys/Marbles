#include <iostream>
#include "Input.h"

class Nanodisc : public Input {

  private:
    float hbelt;    /* height of the protein belt */
    float hlipid; /* height of the lipid bilayer */
    float hcore;    /* height of the hydrophobic bilayer */
    float hmethyl;  /* height of the methyl layer */

    float radius_major; /* major disc semiaxis */
    float radius_minor; /* minor disc semiaxis */
    float width_belt;   /* width of the protein belt */

    float nlipids;
    float watheads;
    float xrough ;
    float cvbelt;
    float cvlipids;
    float cvmp;
    float cvwater;

    float vertical_axis_endcaps;
    float scale_endcaps;

    float rho_h2o;
    float rho_d2o;
    float rho_head;
    float rho_alkyl;
    float rho_methyl;
    float rho_belt;
    float rho_protein;

  public:
    Nanodisc();
    ~Nanodisc();

    void load_input( const std::string& );
};

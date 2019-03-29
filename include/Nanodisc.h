#include <iostream>
#include "Input.h"

class Nanodisc : public Input {

  private:
    bool volume_tests;
    
    float hbelt;    /* height of the protein belt */
    float hlipid; /* height of the lipid bilayer */
    float hcore;    /* height of the hydrophobic bilayer */
    float hmethyl;  /* height of the methyl layer */

    float radius_major; /* major disc semiaxis */
    float radius_minor; /* minor disc semiaxis */
    float width_belt;   /* width of the protein belt */

    float nlipids;
    float wathead;
    float xrough ;
    float cvbelt;
    float cvlipid;
    float cvprotein;
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

    float vh2o;
    float vd2o;
    float vhead;
    float valkyl;
    float vmethyl;
    float vbelt;
    float vprotein;

    float vertical_axis_ellipsoid;

  public:
    Nanodisc();
    ~Nanodisc();

    void load_input( const std::string& );

};

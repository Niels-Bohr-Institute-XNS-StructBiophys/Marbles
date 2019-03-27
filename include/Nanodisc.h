#include <iostream>

class nanodisc {

  private:
    float hbelt;  /* height of the protein belt */
    float hlipid; /* height of the lipid bilayer */
    float hcore;  /* height of the hydrophobic bilayer */
    float hmethyl; /* height of the methyl layer */

    float radius_major; /* major disc semiaxis */
    float radius_minor; /* minor disc semiaxis */
    float width_belt;   /* width of the protein belt */


  public:
    nanodisc();
    ~nanodisc();
};

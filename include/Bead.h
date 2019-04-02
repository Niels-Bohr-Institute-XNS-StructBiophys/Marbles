#include <iostream>

class Bead {

  private:

  public:
    Bead();
    ~Bead();

    float x; /* position along x */
    float y; /* position along y */
    float z; /* position along z */
    float rho; /* excess scattering length */
    float v; /* volume */
    float nn; /* number of neighbors */

    bool position_assigned;

    void assign_position( float, float, float );
    void assign_volume_and_scattlen( const std::string& );

};

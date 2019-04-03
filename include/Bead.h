#include <iostream>

class Bead {

  private:

  public:
    Bead();
    ~Bead();

    double x; /* position along x */
    double y; /* position along y */
    double z; /* position along z */
    double rho; /* excess scattering length */
    double v; /* volume */
    double nn; /* number of neighbors */

    bool position_assigned;

    void assign_position( double, double, double );
    void assign_volume_and_scattlen( const std::string& );

};

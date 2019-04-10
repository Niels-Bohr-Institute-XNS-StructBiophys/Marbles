#include <iostream>

class Bead {

  private:
    double e_scatt_len = 2.82e-13;

  public:
    Bead();
    ~Bead();

    double x; /* position along x */
    double y; /* position along y */
    double z; /* position along z */
    double rho; /* excess scattering length, associated via the corresponding residue */
    double rho_modified; /* excess scattering length that changes depending on the position */
    double v; /* volume */
    double nn; /* number of neighbors */

    int type;

    bool position_assigned;

    void assign_position( double, double, double );
    void assign_volume_and_scattlen( const std::string& );

};

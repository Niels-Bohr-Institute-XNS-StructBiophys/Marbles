#include <iostream>

class Bead {

  private:
    double e_scatt_len = 2.82e-13;

  public:
    Bead();
    ~Bead();

    double x;            /* position along x */
    double y;            /* position along y */
    double z;            /* position along z */
    double rho;          /* excess scattering length, associated via the corresponding residue */
    double rho_modified; /* excess scattering length that changes depending on the position */
    double v;            /* volume */
    double nn;           /* number of neighbors */
    std::string res;     /* residue type (for pruning) */

    //info on the configuration after the attempted MC move
    double x_old;
    double y_old;
    double z_old;
    double rho_mold;
    unsigned int type_old;

    unsigned int type;   /* location of the bead inside the nanodisc */
    unsigned int burn;

    bool position_assigned;
    bool selected;

    void assign_position( double, double, double );
    void update_position( double, double, double );
    void assign_volume_and_scattlen( const std::string& );
    std::vector<double> volume_and_scattlen( const std::string& );
    void assign_average( const std::string& );
    void save_old_config();
    void recover_old_config();
    void accept(); /** assign new configuration */

};

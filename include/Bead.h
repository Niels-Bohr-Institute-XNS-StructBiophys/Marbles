/*******************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

#include <iostream>
#include <vector>

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

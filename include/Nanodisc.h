#include <iostream>
#include <complex>
#include "Input.h"
#include "Tools.h"

#define NH 17 //Order of harmonics
#define NTHETA ((NH+1)*2)
#define NPHI ((NH+1)*2)

class Nanodisc : public Input {

  private:
    bool volume_tests;

    double e_scatt_len = 2.82e-13; /* scattering length of an electron in cm */

    double hbelt;    /* height of the protein belt */
    double hlipid; /* height of the lipid bilayer */
    double hcore;    /* height of the hydrophobic bilayer */
    double hmethyl;  /* height of the methyl layer */

    double radius_major; /* major disc semiaxis */
    double radius_minor; /* minor disc semiaxis */
    double width_belt;   /* width of the protein belt */

    double nlipids;
    double wathead;
    double xrough ;
    double cvbelt;
    double cvlipid;
    double cvprotein;
    double cvwater;

    double vertical_axis_endcaps;
    double scale_endcaps;

    double rho_h2o;
    double rho_d2o;
    double rho_head;
    double rho_alkyl;
    double rho_methyl;
    double rho_belt;
    double rho_protein;

    double vh2o;
    double vd2o;
    double vhead;
    double valkyl;
    double vmethyl;
    double vbelt;
    double vprotein;

    double vertical_axis_ellipsoid;

    const unsigned int harmonics_order = 17;
    const unsigned int ntheta = (harmonics_order + 1) * 2;
    const unsigned int nphi   = (harmonics_order + 1) * 2;

    Array3D<double, 0, NTHETA, NPHI> F;
    Array3D<std::complex<double>, 0, NH+1, NH+1> alpha;

    void flat_disc_form_factor( double, double, double, double, double, int );
    double PsiEllipticCylinderWithEndcaps(double, double, double, double, double, double, double, double);
    void disc_w_endcaps_form_factor( double, double, double, double, double, int );
    //double expand_sh( int );
    double expand_sh2( int );

  public:
    Nanodisc();
    ~Nanodisc();

    void load_input( const std::string& );
    void nanodisc_form_factor( std::vector<double> );

    //get functions
    double get_radius_major();
    double get_radius_minor();
    double get_scale_endcaps();
    double get_vertical_axis_ellipsoid();
    double get_rho_solvent();
    double get_hcore();
    double get_hlipid();
    double get_hmethyl();
    double get_rho_methyl();
    double get_rho_alkyl();
    double get_rho_head();
    double get_cvprotein();
    double get_xrough();
    std::complex<double> get_alpha( int, int, int );
    int get_harmonics_order();

};

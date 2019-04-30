#include <iostream>
#include <complex>
#include "Input.h"
#include "Tools.h"

#define NH 17 //Order of harmonics
#define NTHETA ((NH+1)*2)
#define NPHI ((NH+1)*2)

/*
 * TODO
 * Find a more elegant solution for the allocation of 3D vectors
 */

class Nanodisc : public Input {

  private:

    /* CLASSES */

    /* FLAGS */
    bool volume_tests;                                     /** Consistency check for the quantities loaded from WillItFit */

    /*INPUT FILES */

    /* INFO VARIABLES */
    //To be obtained from WillItFit output */
    double e_scatt_len = 2.82e-13;  /** scattering length of an electron in cm */
    double hbelt;                   /** height of the protein belt */
    double hlipid;                  /** height of the lipid bilayer */
    double hcore;                   /** height of the hydrophobic bilayer */
    double hmethyl;                 /** height of the methyl layer */
    double radius_major;            /** major semiaxis of the disc */
    double radius_minor;            /** minor semiaxis of the disc */
    double width_belt;              /** width of the protein belt */
    double nlipids;                 /** number of lipids in the nanodisc */
    double wathead;                 /** number of water molecules at the lipid polar heads */
    double xrough ;                 /** roughness coefficient */
    double cvbelt;                  /** correction factor for the protein belt volume */
    double cvlipid;                 /** correction factor for the lipids volume */
    double cvprotein;               /** correction factor for the membrane protein volume */
    double cvwater;                 /** correction factor for the water volume */
    double vertical_axis_endcaps;   /** vertical axis of the nanodisc endcaps */
    double scale_endcaps;           /** scale factor of the nanodisc endcaps */
    double vertical_axis_ellipsoid; /** vertical axis of the ellipsoid */
    double rho_h2o;                 /** scatterling length of water */
    double rho_d2o;                 /** scattering length of D2O */
    double rho_head;                /** scattering length of lipid heads */
    double rho_alkyl;               /** scattering length of the alkyl chains */
    double rho_methyl;              /** scattering length of the methyl groups */
    double rho_belt;                /** scattering length of the belt protein */
    double rho_protein;             /** scattering length of the membrane protein */
    double vh2o;                    /** volume of water */
    double vd2o;                    /** volume of D2O */
    double vhead;                   /** volume of lipid heads */
    double valkyl;                  /** volume of the alkyl chains */
    double vmethyl;                 /** volume of the methyl groups */
    double vbelt;                   /** volume of the protein belt */
    double vprotein;                /** volume of the membrane protein */

    /* DETAILS OF THE CALCULATION */
    const unsigned int harmonics_order = 17;               /** Number of spherical harmonics used in the calculations */
    const unsigned int ntheta = (harmonics_order + 1) * 2; /** Number of points in theta for the angular integration */
    const unsigned int nphi   = (harmonics_order + 1) * 2; /** Number of points in phi for the angular integration */

    Array3D<double, 0, NTHETA, NPHI> F;                    /** temporarily stores the values of the form factor during calculation */
    Array3D<std::complex<double>, 0, NH+1, NH+1> alpha;    /** form factor of the nanodisc expanded on the basis of spherical harmonics */

    /* PRIVATE FUNCTIONS */
    void flat_disc_form_factor( double, double, double, double, double, int );      /** computes the form factor of a flat disc */
    void disc_w_endcaps_form_factor( double, double, double, double, double, int ); /** computes the form factor of a disc with elliptic endcaps */

    double PsiEllipticCylinderWithEndcaps( double, double, double, double, double, double, double, double ); /** will get rid of this in flat version! */
    double expand_sh( int );                                                        /** expands the form factor in the basis of spherical harmonics */

  public:
    Nanodisc();  /** class constructor */
    ~Nanodisc(); /** class destructor */

    /* PUBLIC UTILITIES */
    void load_input( const std::string& );            /** reads the output of WillItFit to obtain info on the nanodisc */
    void nanodisc_form_factor( std::vector<double> ); /** computes the form factor of the nanodisc */

    /* GET FUNCTIONS */
    // Return the current value of private members
    double get_radius_major();                        /** returns the major semiaxis of the nanodisc */
    double get_radius_minor();                        /** returns the minor semiaxis of the nanodisc */
    double get_scale_endcaps();                       /** returns the scale factor of endcaps */
    double get_vertical_axis_ellipsoid();             /** returns the vertical axis of the ellipsoid */
    double get_rho_solvent();                         /** returns the scattering length of the solvent */
    double get_hcore();                               /** returns the height of the hydrophobic bilayer */
    double get_hlipid();                              /** returns the height of the lipid bilayer */
    double get_hmethyl();                             /** returns the height of the methyl groups */
    double get_rho_methyl();                          /** returns the scattering length of methyl groups */
    double get_rho_alkyl();                           /** returns the scattering length of the alkyl chains */
    double get_rho_head();                            /** returns the scattering length of the lipid heads */
    double get_cvprotein();                           /** returns the correction factor for the membrane protein volume */
    double get_xrough();                              /** returns the roughness coefficient */
    double get_e_scatt_len();                         /** returns electron scattering length in cm */

    std::complex<double> get_alpha( int, int, int );  /** returns the value of the expanded form factor at position i,l,m */

    int get_harmonics_order();                        /** returns the number of spherical harmonics used in the calculation */

};

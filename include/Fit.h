#include <vector>
#include <nlopt.hpp>
#include <complex>
#include "Nanodisc.h"

#define NH 17 //Order of harmonics

/*
 * Class for the fit of quantities
 */

typedef struct {
  unsigned int len;
  std::vector<double> exper;
  std::vector<double> err;
} bck_data;

typedef struct {
  Array3D<std::complex<double>, 0, NH+1, NH+1> alpha;
  Array3D<std::complex<double>, 0, NH+1, NH+1> beta;
  double e_scattlen;
  double scale_factor;
  double background;
  unsigned int nq;
  unsigned int harmonics_order;
  std::vector<double> exp_q;
  std::vector<double> ref_intensity;
  std::vector<double> err;
} int_data;

class Fit {
    public:

        Fit();
        ~Fit();

        void fit_background( std::vector<std::vector<double> >, unsigned int );
        void fit_intensity( std::vector<std::complex<double> >, std::vector<std::complex<double> >, std::vector<std::vector<double> >, double, unsigned int );
        double get_background();
        double get_rough();
        double get_bck_chi2();
        double get_rough_chi2();
        void set_default_roughness( double );

    private:
        bool setup;
        bck_data datab;
        int_data datai;
        std::vector<double> y_bck; //result of background optimization
        std::vector<double> y_int; //result of background optimization
        std::vector<double> y_int_old; //result of background optimization
        double chi2_bck; //chi^2 obtained from the background optimization
        double chi2_int; //chi^2 obtained from the background optimization
        double chi2_int_old;


        static double chi2_background( const std::vector<double>&, std::vector<double>&, void* ); //without static doesn't work!
        static double chi2_intensity( const std::vector<double>&, std::vector<double>&, void* );
};

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
} int_data;

class Fit {
    public:

        Fit();
        ~Fit();

        void fit_background( std::vector<std::vector<double> >, unsigned int );
        void fit_intensity( std::vector<std::complex<double> >, std::vector<std::complex<double> >, std::vector<double>, std::vector<double>, double, double, double, unsigned int, unsigned int );

    private:
        bool setup;
        bck_data datab;
        int_data datai;
        std::vector<double> y_bck; //result of background optimization
        std::vector<double> y_int; //result of background optimization
        double chi2_bck; //chi^2 obtained from the background optimization
        double chi2_int; //chi^2 obtained from the background optimization

        static double chi2_background( const std::vector<double>&, std::vector<double>&, void* ); //without static doesn't work!
        static double chi2_intensity( const std::vector<double>&, std::vector<double>&, void* );
};

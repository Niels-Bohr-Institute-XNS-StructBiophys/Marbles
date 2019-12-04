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
        //double chi2_int_old;


        static double chi2_background( const std::vector<double>&, std::vector<double>&, void* ); //without static doesn't work!
        static double chi2_intensity( const std::vector<double>&, std::vector<double>&, void* );
};

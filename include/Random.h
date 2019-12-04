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

#include <gsl/gsl_rng.h>
#include <vector>

/*
 * Class for the computation of random numbers
 * uses the functions provided by gsl_rng.h
 */
class RandomNumbers {
    public:
        /*
         * Public functions provided:
         * constructor
         * destructor
         * uniform random number
         * gaussian random number
         */
        RandomNumbers(); /** class constructor: allocates specific gsl arrays for the calculation */
        ~RandomNumbers(); /** class destructor: deallocates gsl arrays */
        double uniform(); /** returns a random number uniformly distributed in [0,1] */
        double in_range( double, double ); /** return a random number uniformly distributed in [a,b] **/
        double in_range2( double, double );
        double gaussian( double ); /** returns a random number gaussian distributed with variance tau */
        std::vector<double> vector( double );
        std::vector<double> vector3( double );
        std::vector<double> vector2( double, double );


    private:
        /*
         * Private informations provided:
         * two gsl_arrays allocated in the constructor
         */
        const gsl_rng_type *T; /** holds static information about each type of generator */
        gsl_rng *r; /** describes an instance of a generator created from a given gsl_rng_type */
};

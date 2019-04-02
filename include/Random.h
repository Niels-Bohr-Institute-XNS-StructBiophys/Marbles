#include <gsl/gsl_rng.h>

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
        float uniform(); /** returns a random number uniformly distributed in [0,1] */
        float in_range( float, float ); /** return a random number uniformly distributed in [a,b] **/
        double in_range2( float, float );
        float gaussian( float ); /** returns a random number gaussian distributed with variance tau */

    private:
        /*
         * Private informations provided:
         * two gsl_arrays allocated in the constructor
         */
        const gsl_rng_type *T; /** holds static information about each type of generator */
        gsl_rng *r; /** describes an instance of a generator created from a given gsl_rng_type */

        bool ready; /** true if the gaussian random generator has already a number ready **/
};

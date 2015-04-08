/*
 * This program calculates \log( \gamma ( 1-s ) )
 * where s = sigma + i * t.
 */

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
using namespace std;

int main() {
  // determines number of bits the answer will have
  int mpfr_bits = 300;
  // initialize mpfr variables
  mpfr_t sigma, t, arg;
  mpfr_init(sigma,mpfr_bits);
  mpfr_init(t,mpfr_bits);
  mpfr_init(arg,mpfr_bits);
  mpfr_set_d(sigma, 0.5, GMP_RNDN); // sigma
  mpfr_set_ui(t, 1000000, GMP_RNDN); // t
  mpfr_set_ui(arg,0,GMP_RNDN);

  argument(t,arg);

  // clear mprf variables
  mpfr_clear(sigma);
  mpfr_clear(t);
  mpfr_clear(arg);
  return 0;
}

/*
 * Calculates Arg( \log( 1/4 + i * t/2 ) )
 */
void argument(mpfr_t t, mpfr_t arg) {
  mpfr_div_ui(arg,t,8,GMP_RNDN);
  mpfr_atan(arg,arg,GMP_RNDN);
}


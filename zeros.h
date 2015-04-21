#include <iostream>
#include <gmp.h>
#include <mpfr.h>
using namespace std;

void logGamma(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits);
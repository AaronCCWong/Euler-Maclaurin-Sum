/*
 * Calculates \theta (t) = Argument of \pi^{-it/2} * \Gamma( \frac{1}{4} + i*\frac{t}{2} )
 */
#include "zeros.h"

void argumentPi(mpfr_t t, mpfr_t argPi, int mpfr_bits);

int main()
{
	// determines number of bits the answer will have
	int mpfr_bits = 300;
	// initialize mpfr variables
	mpfr_t sigma, t, epsilon, argPi, logGammaReal, logGammaImag, theta;
	mpfr_init2(sigma, mpfr_bits);
	mpfr_init2(t, mpfr_bits);
	mpfr_init2(epsilon, mpfr_bits);
	mpfr_init2(argPi, mpfr_bits);
	mpfr_init2(logGammaReal, mpfr_bits);
	mpfr_init2(logGammaImag, mpfr_bits);
	mpfr_init2(theta, mpfr_bits);
	mpfr_set_d(sigma, 0.5, GMP_RNDN); // sigma
	mpfr_set_ui(t, 1000000, GMP_RNDN); // t
	mpfr_set_d(epsilon, 1E-80, GMP_RNDN); // epsilon
	mpfr_set_d(argPi, 0, GMP_RNDN);
	mpfr_set_d(logGammaReal, 0, GMP_RNDN);
	mpfr_set_d(logGammaImag, 0, GMP_RNDN);
	mpfr_set_d(theta, 0, GMP_RNDN);

	argumentPi(t, argPi, mpfr_bits);
	logGamma(sigma, t, epsilon, logGammaReal, logGammaImag, mpfr_bits);
	mpfr_add(theta, logGammaImag, argPi, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(sigma);
	mpfr_clear(t);
	mpfr_clear(epsilon);
	mpfr_clear(argPi);
	mpfr_clear(logGammaReal);
	mpfr_clear(logGammaImag);
	mpfr_clear(theta);

	return 0;
}

/*
 * Calculates the argument of \pi^{-it/2}:
 *            \pi^{-it/2} = e ^ { -i * \frac{t * log (\pi)}{2} }
 * where the argument is given by
 *              \frac{t * log (\pi)}{2}.
 */
void argumentPi(mpfr_t t, mpfr_t argPi, int mpfr_bits) {
	mpfr_t pi_holder;
	mpfr_init2(pi_holder, mpfr_bits);
	mpfr_const_pi(pi_holder, GMP_RNDN);
    
    mpfr_log(argPi, pi_holder, GMP_RNDN);
    mpfr_mul(argPi, argPi, t, GMP_RNDN);
    mpfr_div_ui(argPi, argPi, 2, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(pi_holder);
}
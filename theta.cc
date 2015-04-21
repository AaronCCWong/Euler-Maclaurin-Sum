/*
 * Calculates \theta (t) = Argument of \pi^{-it/2} * \Gamma( \frac{1}{4} + i*\frac{t}{2} )
 */
#include <iostream>
#include <gmp.h>
#include <mpfr.h>

using namespace std;

int main()
{
	// determines number of bits the answer will have
	int mpfr_bits = 300;
	// initialize mpfr variables
	mpfr_t sigma, t, epsilon, argPi;
	mpfr_init(sigma, mpfr_bits);
	mpfr_init(t, mpfr_bits);
	mpfr_init(epsilon, mpfr_bits);
	mpfr_init(argPi, mpfr_bits);
	mpfr_set_d(sigma, 0.5, GMP_RNDN); // sigma
	mpfr_set_ui(t, 1000000, GMP_RNDN); // t
	mpfr_set_d(epsilon, 1E-80, GMP_RNDN); // epsilon
	mpfr_set_d(argPi, 0, GMP_RNDN); // epsilon

	// clear mpfr variables
	mpfr_clear(sigma);
	mpfr_clear(t);
	mpfr_clear(epsilon);
	mpfr_clear(argPi)

	return 0;
}

/*
 * Calculates the argument of \pi^{-it/2} using the formula:
 *				ArcTan[ -Sin( t/2 * log\pi ) / Cos( t/2 * log\pi ) ]
 */
void argumentPi(mpfr_t t, mpfr_t argPi, int mpfr_bits) {
	mpfr_t pi, counter1;
	mpfr_init(counter1, mpfr_bits);
	mpfr_init(pi, GMP_RNDN);
	mpfr_const_pi(pi, GMP_RNDN);
	mpfr_set_ui(counter1, 0, GMP_RNDN);
	
	// calculates cos( t/2 * log\pi )
	mpfr_div_ui(counter1, t, 2, GMP_RNDN);
	mpfr_log(argPi, pi, GMP_RNDN);
	mpfr_mul(counter1, counter1, argPi, GMP_RNDN);
	mpfr_cos(counter1, counter1, GMP_RNDN);
	// calculates -sin( t/2 * log\pi )
	mpfr_mul(argPi, argPi, t, GMP_RNDN);
	mpfr_div_ui(argPi, argPi, 2, GMP_RNDN);
	mpfr_sin(argPi, argPi, GMP_RNDN);
	mpfr_mul_si(argPi, argPi, -1, GMP_RNDN);

	mpfr_div(counter1, argPi, counter1, GMP_RNDN);
	mpfr_atan(argPi, counter1, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(pi);
	mpfr_clear(counter1)
}
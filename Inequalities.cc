#include <iostream>
#include <gmp.h>
#include <mpfr.h>
using namespace std;

int inequalityOfL1(mpfr_t sigma, mpfr_t t, double epsilon, int mpfr_bits);
int inequalityOfN(mpfr_t sigma, mpfr_t t, int L1, int mpfr_bits);

int inequality() {
	int mpfr_bits = 100;
	double epsilon = .00001;
	mpfr_t sigma, t;
	mpfr_init2(sigma, mpfr_bits);
	mpfr_init2(t, mpfr_bits);
	mpfr_set_d(sigma, 0.5, GMP_RNDN);
	mpfr_set_d(t, 10000000, GMP_RNDN);

	int L1 = inequalityOfL1(sigma, t, epsilon, mpfr_bits);
	int N = inequalityOfN(sigma, t, L1, mpfr_bits);
	cout << L1 << "--------" << N;

	mpfr_clear(sigma);
	mpfr_clear(t);
	return 0;
}

int inequalityOfL1(mpfr_t sigma, mpfr_t t, double epsilon, int mpfr_bits) {
	int L1(1); // L1 that will be returned to be used as endSum of Euler-Maclaurin remainder
	// initialize mpfr variables
	mpfr_t Lmpfr;
	mpfr_t counter, counter2, counter3, epsmpfr;
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_init2(counter3, mpfr_bits);
	mpfr_init2(Lmpfr, mpfr_bits);
	mpfr_init2(epsmpfr, mpfr_bits);
	mpfr_set_si(Lmpfr, 1, GMP_RNDN);
	mpfr_set_si(counter, 0, GMP_RNDN);
	mpfr_set_d(epsmpfr, epsilon, GMP_RNDN);
	/*
		When s = (0,0) and L1 = 1, the term on the right hand side of the inequality is
		equal to 11.51 so we assume this is the smallest value as to where the while loop
		starts. We set counter = 12 and L1 = 1. We then increment L1 until the left hand
		side of the inequality is greater than the right hand side.
	*/
	mpfr_set_si(counter2, 12, GMP_RNDN);
	mpfr_set_si(counter3, 0, GMP_RNDN);

	while (mpfr_cmp(counter, counter2) < 0) {
		mpfr_add_si(Lmpfr, Lmpfr, 1, GMP_RNDN);
		/*
		2 * L1 - 1
		*/
		mpfr_mul_si(counter, Lmpfr, 2, GMP_RNDN);
		mpfr_sub_ui(counter, counter, 1, GMP_RNDN);
		/*
		0.5 * log|s + 2*L1 - 1| - log( epsilon )

		*/
		mpfr_mul_si(counter2, Lmpfr, 2, GMP_RNDN);
		mpfr_sub_ui(counter2, counter2, 1, GMP_RNDN);
		mpfr_add(counter2, counter2, sigma, GMP_RNDN);
		mpfr_pow_ui(counter2, counter2, 2, GMP_RNDN); // (sigma + 2L1 -1)^2
		mpfr_mul(counter3, t, t, GMP_RNDN); // t^2
		mpfr_add(counter2, counter2, counter3, GMP_RNDN);
		mpfr_sqrt(counter2, counter2, GMP_RNDN); // sqrt((sigma + 2L1 -1)^2 + t^2)
		mpfr_log(counter2, counter2, GMP_RNDN); // log(sqrt((sigma + 2L1 -1)^2 + t^2))
		mpfr_mul_d(counter2, counter2, 0.5, GMP_RNDN);
		mpfr_log(counter3, epsmpfr, GMP_RNDN);
		mpfr_sub(counter2, counter2, counter3, GMP_RNDN);
	}
	L1 = mpfr_get_si(Lmpfr, GMP_RNDN); // convert Lmpfr to an integer

	// clear mpfr variables
	mpfr_clear(Lmpfr);
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(counter3);
	return L1;
}

int inequalityOfN(mpfr_t sigma, mpfr_t t, int L1, int mpfr_bits) {
	int N; // to be returned to be used as the endSum of the total sum
	// initialize mpfr variables
	mpfr_t counter, counter2, one;
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_init2(one, mpfr_bits);
	mpfr_set_si(counter, L1, GMP_RNDN);
	mpfr_set_si(counter2, 0, GMP_RNDN);
	mpfr_set_si(one, 1, GMP_RNDN);
	
	/*
		e|s + 2*L1 - 1| / 2*pi
	*/
	mpfr_mul_si(counter, counter, 2, GMP_RNDN);
	mpfr_add(counter, counter, sigma, GMP_RNDN);
	mpfr_sub_ui(counter, counter, 1, GMP_RNDN);
	mpfr_pow_ui(counter, counter, 2, GMP_RNDN);
	mpfr_pow_ui(counter2, t, 2, GMP_RNDN);
	mpfr_add(counter, counter, counter2, GMP_RNDN);
	mpfr_sqrt(counter, counter, GMP_RNDN); // sqrt((sigma + 2L1 -1)^2 + t^2)
	mpfr_exp(counter2, one, GMP_RNDN); // e
	mpfr_mul(counter, counter, counter2, GMP_RNDN);
	mpfr_const_pi(counter2, GMP_RNDN); // pi
	mpfr_mul_si(counter2, counter2, 2, GMP_RNDN);
	mpfr_div(counter, counter, counter2, GMP_RNDN);

	// gets the ceiling of counter and sets it equal to N
	mpfr_ceil(counter, counter);
	N = mpfr_get_si(counter, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(one);
	return N;
}
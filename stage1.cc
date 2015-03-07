#include "stage1.h"

int stage1() {
    return 0;
}

// Calculates the terms in the partial sum
void summand (mpfr_t sigma, mpfr_t t, mpz_t n, mpfr_t rsum, mpfr_t isum, int mpfr_bits) {
    mpfr_t var1, var2, var3, nmpfr;
    mpfr_init2(var1,mpfr_bits);
    mpfr_init2(var2,mpfr_bits);
    mpfr_init2(var3,mpfr_bits);
    mpfr_init2(nmpfr,mpfr_bits);
	mpfr_set_d(var1, 0, GMP_RNDN);
	mpfr_set_d(var2, 0, GMP_RNDN);
	mpfr_set_d(var3, 0, GMP_RNDN);
    
    // Calculates the real part of the summand
    mpfr_mul_si (var1, sigma, -1.0, GMP_RNDN); // -sigma
    mpfr_set_z(nmpfr,n,GMP_RNDN);
    mpfr_log (var2, nmpfr, GMP_RNDN); // log(n)
    mpfr_mul (var1, var1, var2,GMP_RNDN); // -sigma*log(n)
    mpfr_exp (var1, var1, GMP_RNDN); // exp(-sigma*log(n))
    mpfr_mul (var3, t, var2, GMP_RNDN); // t*log(n)
    mpfr_cos (var2, var3, GMP_RNDN); // cos(t*log(n))
    mpfr_mul (rsum, var1, var2, GMP_RNDN); // exp(-sigma*log(n))*cos(t*log(n))

    // Calculates the imaginary part of the summand
    mpfr_mul_si (var1, var1, -1.0, GMP_RNDN); // -exp(-sigma*log(n))
    mpfr_sin (var2, var3, GMP_RNDN); // sin(t*log(n))
    mpfr_mul (isum, var1, var2, GMP_RNDN); // -exp(-sigma*log(n))*sin(t*log(n))
    
    mpfr_clear(var1);
    mpfr_clear(var2);
    mpfr_clear(var3);
}

// Calculates the partial sum using MPFR
void partial_sum_mpfr(mpfr_t sigma, mpfr_t t, int M1, mpfr_t rresult, mpfr_t iresult, int mpfr_bits) {
    mpz_t n, end;
    mpfr_t rsum, isum;
    
    mpz_init (n);
    mpz_init (end);
    mpfr_init2 (rsum,mpfr_bits);
    mpfr_init2 (isum,mpfr_bits);

    mpz_set_ui(n,1);
    mpz_set_ui(end,M1-1);
	mpfr_set_d(rsum, 0, GMP_RNDN);
	mpfr_set_d(isum, 0, GMP_RNDN);
	mpfr_set_d(rresult, 0, GMP_RNDN);
	mpfr_set_d(iresult, 0, GMP_RNDN);

    while (mpz_cmp(end, n) >= 0) { // while M1-1>n for n starting at 1
        summand (sigma, t, n, rsum, isum,mpfr_bits); // n^(-s) where s=sigma+i*t
		mpfr_add(rresult, rresult, rsum, GMP_RNDN);
		mpfr_add(iresult, iresult, isum, GMP_RNDN);
        mpz_add_ui (n, n, 1);
    }

    mpz_clear(n);
    mpz_clear(end);
    mpfr_clear(rsum);
    mpfr_clear(isum);
}
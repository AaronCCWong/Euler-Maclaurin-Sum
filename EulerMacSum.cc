//
// Calculates the Euler-Maclaurin sum using stage 1.cc
//

#include "stage1.h"
#include "Misc.h"

using namespace std;

int EMSum()
{
	/*
	int precision(100);
	double epsilon(.00001);
	// initialize and set mpfr variables
	mpfr_t sigma, t, rresult, iresult;
	mpfr_init2(sigma,precision);
	mpfr_init2(t,precision);
	mpfr_init2(rresult,precision);
	mpfr_init2(iresult,precision);
	mpfr_set_d(sigma, 0.5, GMP_RNDN);
	mpfr_set_d(t, 1000, GMP_RNDN);
	mpfr_set_d(rresult, 0, GMP_RNDN);
	mpfr_set_d(iresult , 0, GMP_RNDN);

	// initialize and set mpz variables
	mpz_t N, M, L;
	mpz_init(N);
	mpz_init(M);
	mpz_init(L);
	mpz_set_ui(N, 1);
	Inequalities(sigma, t, M, L, epsilon, precision);

	// Calculates \sum_{n=N}^{N+M-1} n^(-sigma+it)
	partial_sum_mpfr(sigma, t, N, M, precision, rresult, iresult);
    
	// calculates the remainder terms
	EMsum(sigma, t, M, L, rresult, iresult, precision);

	cout << endl;
	cout << "The real parts are:";
	cout << endl;
	mpfr_out_str(stdout, 10, precision, rresult, GMP_RNDN);

	cout << endl;
	cout << "The imaginary parts are:";
	cout << endl;
	mpfr_out_str(stdout, 10, precision, iresult, GMP_RNDN);
    cout << endl;

	// clear mpfr variables
	mpfr_clear(sigma);
	mpfr_clear(t);
	mpfr_clear(rresult);
	mpfr_clear(iresult);

	// clear mpz variables
	mpz_clear(N);
	mpz_clear(M);
	mpz_clear(L);
	*/

	return 0;
}

// Calculates the Euler-Maclaurin sum remainder from n=1 to N-1
void EMsum(mpfr_t sigma, mpfr_t t, int M1, int L1, mpfr_t rresult, mpfr_t iresult, int precision) {
	mpfr_t Mmpfr, Lmpfr, l, l1, counter, token, token2, token3, token4, token5, token6, token7, token8, token9, endprod;
	mpfr_t bnum, bden, bernoulli, token10, token11;
	mpfr_init2(endprod, precision);
	mpfr_init2(Mmpfr, precision);
	mpfr_init2(Lmpfr, precision);
	mpfr_init2(l, precision);
	mpfr_init2(l1, precision);
	mpfr_init2(bnum, precision);
	mpfr_init2(bden, precision);
	mpfr_init2(counter, precision);
	mpfr_init2(token,precision);
	mpfr_init2(token2,precision);
	mpfr_init2(token3,precision);
	mpfr_init2(token4,precision);
	mpfr_init2(token5,precision);
	mpfr_init2(token6,precision);
	mpfr_init2(token7,precision);
	mpfr_init2(token8,precision);
	mpfr_init2(token9,precision);
	mpfr_init2(token10, precision);
	mpfr_init2(token11, precision);
	mpfr_init2(bernoulli, precision);
	mpfr_set_d(counter, 0, GMP_RNDN);
	mpfr_set_d(token, 0, GMP_RNDN);
	mpfr_set_d(token2, 0, GMP_RNDN);
	mpfr_set_d(token3, 0, GMP_RNDN);
	mpfr_set_d(token4, 0, GMP_RNDN);
	mpfr_set_d(token5, 0, GMP_RNDN);
	mpfr_set_d(token6, 0, GMP_RNDN);
	mpfr_set_d(token7, 0, GMP_RNDN);
	mpfr_set_d(token8, 0, GMP_RNDN);
	mpfr_set_d(token9, 0, GMP_RNDN);
	mpfr_set_d(token10, 0, GMP_RNDN);
	mpfr_set_d(token11, 0, GMP_RNDN);
	mpfr_set_d(bernoulli, 0, GMP_RNDN);
	mpfr_set_si(l, 1, GMP_RNDN);
	mpfr_set_ui(Mmpfr, M1, GMP_RNDN);
    mpfr_add_ui(Mmpfr, Mmpfr, 1, GMP_RNDN);
	mpfr_set_ui(Lmpfr, L1, GMP_RNDN);
	mpfr_set_si(l1, 1, GMP_RNDN);


	// (M^(-s))/2
	mpfr_log(token, Mmpfr, GMP_RNDN); // log(M)
	mpfr_mul(token2, token, sigma, GMP_RNDN); // sigma*log(M)
	mpfr_mul_d(token2, token2, -1, GMP_RNDN); // -sigma*log(M)
	mpfr_exp(token2, token2, GMP_RNDN); // exp(-sigma*log(M))
	mpfr_mul_d(token2, token2, 0.5, GMP_RNDN); // (1/2)*exp(-sigma*log(M))
	mpfr_mul(token, token, t, GMP_RNDN); // t*log(M)
	mpfr_cos(token4, token, GMP_RNDN); // cos(t*log(M))
	mpfr_mul(token3, token2, token4, GMP_RNDN); // (1/2)*exp(-sigma*log(M))*cos(t*log(M))

	mpfr_add(rresult, rresult, token3, GMP_RNDN);

	mpfr_mul_d(token2, token2, -1, GMP_RNDN); // -(1/2)*exp(-sigma*log(M))
	mpfr_sin(token4, token, GMP_RNDN); // sin(t*log(M))
	mpfr_mul(token4, token2, token4, GMP_RNDN); // -(1/2)*exp(-sigma*log(M))*sin(t*log(M))

	mpfr_add(iresult, iresult, token4, GMP_RNDN);
    

	// (M^(1-s))/(s-1)
    mpfr_mul_ui(token, Mmpfr, 2, GMP_RNDN);
    mpfr_sub_ui(token2, sigma, 1,GMP_RNDN);
    mpfr_mul(token, token, token2, GMP_RNDN);
    
    mpfr_mul(token2, token2, token2, GMP_RNDN);
    mpfr_mul(token6, t, t, GMP_RNDN);
    mpfr_add(token5, token2, token6, GMP_RNDN);
    mpfr_div(token, token, token5, GMP_RNDN);
    
    mpfr_mul(token8, Mmpfr, t, GMP_RNDN);
    mpfr_mul_ui(token8, token8, 2, GMP_RNDN);
    mpfr_mul_si(token8, token8, -1, GMP_RNDN);
    mpfr_div(token2, token8, token5, GMP_RNDN);
    
    mpfr_mul(token6, token, token3, GMP_RNDN);
    mpfr_mul(token7, token2, token4, GMP_RNDN);
    mpfr_sub(token6, token6, token7, GMP_RNDN);
    
    mpfr_mul(token5, token, token4, GMP_RNDN);
    mpfr_mul(token7, token2, token3, GMP_RNDN);
    mpfr_add(token5, token5, token7, GMP_RNDN);
    
    mpfr_add(rresult, rresult, token6, GMP_RNDN);
    mpfr_add(iresult, iresult, token5, GMP_RNDN);
    
	// \sum_{l1=1}^{L} T_{l1, M}(s)
	mpfr_set_si(token7, 0, GMP_RNDN); // real part of our sum
	mpfr_set_si(token8, 0, GMP_RNDN); // imaginary part of our sum

	// Constant part of sum
	int arr_index = mpfr_get_ui(l1, GMP_RNDN);
	int arr_index1 = arr_index - 1;

	mpfr_set_z(bnum, numerator_arr[arr_index1].get_mpz_t(), GMP_RNDN);
	mpfr_set_z(bden, denominator_arr[arr_index1].get_mpz_t(), GMP_RNDN);

	mpfr_div(counter, bnum, bden, GMP_RNDN); // Bernoulli number B_{2*l1}
	int darr_index = 2 * arr_index;
	mpfr_fac_ui(token9, darr_index, GMP_RNDN);
	mpfr_div(bernoulli, counter, token9, GMP_RNDN); // B_{2*l1}/(2*l1!)
	while (mpfr_cmp(Lmpfr, l1) >= 0) {
		// M^{-s}
		mpfr_mul_si(token9, sigma, -1, GMP_RNDN);
		mpfr_log(token5, Mmpfr, GMP_RNDN);
		mpfr_mul(token9, token9, token5, GMP_RNDN);
		mpfr_exp(token9, token9, GMP_RNDN); // e^{-sigma*logM}
		mpfr_mul(token6, t, token5, GMP_RNDN);
		mpfr_cos(token6, token6, GMP_RNDN);
		mpfr_mul(token6, token9, token6, GMP_RNDN); // e^{-sigma*logM}cos(tlogM)
		mpfr_mul(token5, t, token5, GMP_RNDN);
		mpfr_sin(token5, token5, GMP_RNDN);
		mpfr_mul_si(token5, token5, -1, GMP_RNDN);
		mpfr_mul(token9, token9, token5, GMP_RNDN); // -e^{-sigma*logM}sin(tlogM)

		mpfr_mul(counter, bernoulli, token6, GMP_RNDN); // real part
		mpfr_mul(token9, bernoulli, token9, GMP_RNDN); // imaginary part

		// endprod = 2 * l1 - 2
		mpfr_mul_ui(endprod, l1, 2, GMP_RNDN);
		mpfr_sub_ui(endprod, endprod, 2, GMP_RNDN);

		mpfr_div(token, sigma, Mmpfr, GMP_RNDN); // real part of our product
		mpfr_div(token2, t, Mmpfr, GMP_RNDN); // imaginary part of our product
        
        
		// \prod_{l=0}^{2*l1-2} (s+l)/M
		while (mpfr_cmp(endprod, l) >= 0) {
			mpfr_add(token3, sigma, l, GMP_RNDN);
			mpfr_div(token3, token3, Mmpfr, GMP_RNDN); // (sigma+l)/M
			mpfr_div(token4, t, Mmpfr, GMP_RNDN); // t/M
			// (a+bi)*(c+di)=(ac-bd)+(ad+bc)i
			mpfr_mul(token5, token, token3, GMP_RNDN);
			mpfr_mul(token6, token2, token4, GMP_RNDN);
			mpfr_mul(token10, token, token4, GMP_RNDN);
			mpfr_mul(token11, token2, token3, GMP_RNDN);
			mpfr_sub(token, token5, token6, GMP_RNDN);
			mpfr_add(token2, token10, token11, GMP_RNDN);

			mpfr_add_ui(l, l, 1, GMP_RNDN); // increment l by 1
		}
        mpfr_set_si(l, 1, GMP_RNDN);
        
		mpfr_mul(token5, counter, token, GMP_RNDN);
		mpfr_mul(token6, token9, token2, GMP_RNDN);
		mpfr_sub(token5, token5, token6, GMP_RNDN);
		mpfr_add(token7, token7, token5, GMP_RNDN);

		mpfr_mul(token5, counter, token2, GMP_RNDN);
		mpfr_mul(token6, token9, token, GMP_RNDN);
		mpfr_add(token5, token5, token6, GMP_RNDN);
		mpfr_add(token8, token8, token5, GMP_RNDN);

		mpfr_add_ui(l1, l1, 1, GMP_RNDN); // increment l1 by 1
	}

	mpfr_add(rresult, rresult, token7, GMP_RNDN);
	mpfr_add(iresult, iresult, token8, GMP_RNDN);
    
	mpfr_clear(token);
	mpfr_clear(token2);
	mpfr_clear(token3);
	mpfr_clear(token4);
	mpfr_clear(token5);
	mpfr_clear(token6);
	mpfr_clear(token7);
	mpfr_clear(token8);
	mpfr_clear(token9);
	mpfr_clear(token10);
	mpfr_clear(token11);
	mpfr_clear(Mmpfr);
	mpfr_clear(l);
	mpfr_clear(Lmpfr);
	mpfr_clear(l1);
	mpfr_clear(bernoulli);
}
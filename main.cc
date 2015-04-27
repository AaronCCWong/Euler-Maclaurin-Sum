/* 
    program to compute values of partial dirichlet sum

              \sum_{n=N}^{N+M-1} n^(-sigma+it)

    to within eps_rnd in roundoff error.
*/

#include "main.h"
#include <time.h>

// chooses the number of bits to use in stage1 and stage 2
int numOfBits(mpfr_t sigma, mpfr_t epsilon, int N, int M1, int mpfr_bits);

int main(int argc, char * argv[]) {
    clock_t time;
    
    string sigma_string = "";
    string t_string = "";
    string epsilon_string = "";

    mpfr_t t, sigma, epsilon; // s = sigma + it
    double td,sigmad; // double versions of sigma and t
    double Nd,Md; // double versions of N and M
    bool stage1 = 0;
	bool stage2 = 0;

    while(1) {
        enum {SIGMA, T, EPSILON, STAGE1, STAGE2};

        static struct option long_options[] =
            {
                {"sigma", required_argument, 0, SIGMA},
                {"t", required_argument, 0, T},
                {"epsilon", required_argument, 0, EPSILON},
                {"stage1",no_argument, 0, STAGE1},
				{"stage2", no_argument, 0, STAGE2},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        int cc = getopt_long(argc, argv, "", long_options, &option_index);
        if (cc == -1)  break;

        switch(long_options[option_index].val) {
            case SIGMA:
                sigma_string = optarg;
                sigmad = (double)atof(optarg);
                break;
            case T:
                t_string = optarg;
                td = (double)atof(optarg);
                break;
            case EPSILON:
                epsilon_string = optarg;
                break;
            case STAGE1:
                stage1 = 1;
                break;
			case STAGE2:
				stage2 = 1;
				break;
        }
    }
    time = clock();
    
	// precision = ceiling(log(eps_rnd/(2*M*max(N^(-sigma),(N+M-1)^(-sigma)))))
	// double mm1 = max(pow(Nd, -sigmad), pow(Nd + Md - 1, -sigmad));
	// double mm2 = (abs(td) + abs(sigmad))*log(Nd + Md);
    
	// precision = ceil(-log(eps_rnd / (2 * Md*mm1*mm2)) / log(2.));
    int mpfr_bits = 1000;

    mpfr_init2(t, mpfr_bits);
	mpfr_init2(sigma, mpfr_bits);
    mpfr_init2(epsilon, mpfr_bits);

	// set mpz and mpfr variable to supplied values
	mpfr_set_str(sigma, sigma_string.c_str(), 10, GMP_RNDN);
	mpfr_set_str(t, t_string.c_str(), 10, GMP_RNDN);
    mpfr_set_str(epsilon, epsilon_string.c_str(), 10, GMP_RNDN);

	// sum_{n=N}^{N+M-1} n^(-s) = r_result + i*i_result
	mpfr_t rresult, iresult;
	mpfr_init2(rresult, mpfr_bits);
	mpfr_init2(iresult, mpfr_bits);
    
    mpfr_set_d(rresult, 0, GMP_RNDN);
    mpfr_set_d(iresult, 0, GMP_RNDN);

    int L1 = inequalityOfL1(sigma, t, epsilon, mpfr_bits);
    int N = inequalityOfN(sigma, t, L1, mpfr_bits);
    int M1 = endStage1(sigma, t, N, mpfr_bits);
	// mpfr_bits = numOfBits(sigma, epsilon, N, M1, mpfr_bits);

    if (stage1) {
        partial_sum_mpfr(sigma, t, M1, rresult, iresult, mpfr_bits);
		EMsum(sigma, t, N, L1, rresult, iresult, mpfr_bits);
    } 
	if (stage2) {
		partialSum2MPFR(sigma, t, M1, N, rresult, iresult, epsilon, mpfr_bits);
	}
	
    cout << endl;
    mpfr_out_str(stdout, 10, mpfr_bits, rresult, GMP_RNDN);
    cout << "            ";
    mpfr_out_str(stdout, 10, mpfr_bits, iresult, GMP_RNDN);
    cout << endl;

    //clear mpfr variables
    mpfr_clear(t);
    mpfr_clear(sigma);
    mpfr_clear(epsilon);
    
    time = clock() - time;
    cout << "Program took " << ((float)time)/CLOCKS_PER_SEC << " seconds to run." << endl;

    return 0;
}

// chooses when to stop stage 1 and begin stage 2
int endStage1(mpfr_t sigma, mpfr_t t, int N, int mpfr_bits) {
    int M1 = 0;
    mpfr_t counter, counter2;
    mpfr_init2(counter, mpfr_bits);
    mpfr_init2(counter2, mpfr_bits);
    mpfr_set(counter, sigma, GMP_RNDN);
    mpfr_set(counter2, t, GMP_RNDN);
    mpfr_mul(counter, counter, counter, GMP_RNDN);

    mpfr_mul(counter2, counter2, counter2, GMP_RNDN);
    mpfr_add(counter, counter, counter2, GMP_RNDN);
    mpfr_sqrt(counter, counter, GMP_RNDN);
    mpfr_root(counter, counter, 3, GMP_RNDN);
    mpfr_add_ui(counter, counter, 100, GMP_RNDN);
    mpfr_ceil(counter, counter);
    
    int M = mpfr_get_ui(counter, GMP_RNDN);
    M1 = min(N, M);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);

    return M1;
}

int numOfBits(mpfr_t sigma, mpfr_t epsilon, int N, int M1, int mpfr_bits) {
	/*
		We use the following formula for the number of bits to use:
		             _                                                                _ 
		mpfr_bits = | log( (epsilon/4) / ( 2 * M1 * max( N^{-sigma}, M1^{-sigma} ) ) ) |
	*/
	mpfr_t bits, negSigma, counter, counter2;
	mpfr_init2(bits, mpfr_bits);
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_init2(negSigma, mpfr_bits);
	mpfr_set(bits, epsilon, GMP_RNDN);
	mpfr_set_ui(counter, 0, GMP_RNDN);
	mpfr_set_ui(counter2, 0, GMP_RNDN);
	mpfr_set(negSigma, sigma, GMP_RNDN);

	mpfr_mul_si(negSigma, sigma, -1, GMP_RNDN); // -sigma
	mpfr_div_ui(bits, bits, 4, GMP_RNDN);
	mpfr_ui_pow(counter, N, negSigma, GMP_RNDN);
	mpfr_ui_pow(counter2, M1, negSigma, GMP_RNDN);
	mpfr_max(counter, counter, counter2, GMP_RNDN); // max( N^{-sigma}, M1^{-sigma} )
	mpfr_mul_ui(counter, counter, 2, GMP_RNDN);
	mpfr_mul_ui(counter, counter, M1, GMP_RNDN);
	mpfr_div(bits, bits, counter, GMP_RNDN);
	mpfr_ceil(bits, bits); // mpfr_bits

	mpfr_bits = mpfr_get_si(bits, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(negSigma);
	mpfr_clear(bits);

	return mpfr_bits;
}
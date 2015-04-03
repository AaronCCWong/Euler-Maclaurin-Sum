#include "main.h"

int numOfPartitions(int M1, int N);
int endPartition(int N, int v, int M1);
int numOfTerms(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, int v, int M1, int N, int mpfr_bits);
// Calculates the coefficients of the real part of our taylor expansion
void coefficientRealCalculator(mpfr_t sigma, int v, mpfr_t *coefficientImagValues,
														int numOfTerms, int mpfr_bits);
// Calculates the coefficients of the imaginary part of our taylor expansion
void coefficientImagCalculator(mpfr_t t, int v, mpfr_t *coefficientImagValues, 
														int numOfTerms, int mpfr_bits);
// Calculates the Taylor series
void taylorSeries(mpfr_t sigma, mpfr_t t, int k, int numOfTerms,
					mpfr_t realTaylorSeries, mpfr_t imagTaylorSeries, 
					mpfr_t *coefficientRealValues, mpfr_t *coefficientImagValues, int mpfr_bits);
// Calculates the terms of the partial sum
void sumTerms(mpfr_t sigma, mpfr_t t, mpfr_t realTaylorSeries, mpfr_t imagTaylorSeries,
				mpfr_t realTerm, mpfr_t imagTerm, int mpfr_bits);

int stage2() {
	int mpfr_bits, M1, N;
	mpfr_bits = 300;
	M1 = 201;
	N = 432628;

	// initialize mpfr_t variables
	mpfr_t sigma, t, epsilon, rresult, iresult, realTaylorSeries, imagTaylorSeries;
	mpfr_init2(sigma, mpfr_bits);
	mpfr_init2(t, mpfr_bits);
	mpfr_init2(epsilon, mpfr_bits);
	mpfr_init2(rresult, mpfr_bits);
	mpfr_init2(iresult, mpfr_bits);
	mpfr_init2(realTaylorSeries, mpfr_bits);
	mpfr_init2(imagTaylorSeries, mpfr_bits);
	mpfr_set_d(sigma, 0.5, GMP_RNDN);
	mpfr_set_d(t, 1000000, GMP_RNDN);
	mpfr_set_d(epsilon, 1E-80, GMP_RNDN);
	mpfr_set_ui(rresult, 0, GMP_RNDN);
	mpfr_set_ui(iresult, 0, GMP_RNDN);
	mpfr_set_ui(realTaylorSeries, 0, GMP_RNDN);
	mpfr_set_ui(imagTaylorSeries, 0, GMP_RNDN);

	partialSum2MPFR(sigma, t, M1, N, rresult, iresult, epsilon, mpfr_bits);

	// clear mpfr_t variables
	mpfr_clear(sigma);
	mpfr_clear(t);
	mpfr_clear(epsilon);
	mpfr_clear(rresult);
	mpfr_clear(iresult);
	return 0;
}

int numOfPartitions(int M1, int N) {
	int R = ceil(2 * M1*log(N / M1) + 1);
	return R;
}

int numOfTerms(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, int v, int M1, int N, int mpfr_bits) {
	int R = numOfPartitions(M1, N);
    //cout << "R=" << R << endl;
	int K = endPartition(N, v, M1) - v + 1;
	// initialize mpfr variables
	mpfr_t counter, counter2, counter3;
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_init2(counter3, mpfr_bits);
	mpfr_set(counter, epsilon, GMP_RNDN);
	mpfr_set(counter2, sigma, GMP_RNDN);
	mpfr_set(counter3, t, GMP_RNDN);
	/*
		numOfTerms = ceil( log( epsilon/(2*R)/(2*(|s|+1)) ) / log( K/v ) )
	*/
	mpfr_div_ui(counter, counter, 2, GMP_RNDN);
	mpfr_div_ui(counter, counter, R, GMP_RNDN);
	mpfr_div_ui(counter, counter, 2, GMP_RNDN); // epsilon/(2*R)/2
	mpfr_pow_ui(counter2, counter2, 2, GMP_RNDN);
	mpfr_pow_ui(counter3, counter3, 2, GMP_RNDN);
	mpfr_add(counter2, counter2, counter3, GMP_RNDN); // sigma^2 + t^2
	mpfr_sqrt(counter2, counter2, GMP_RNDN); // sqrt( sigma^2 + t^2 )
	mpfr_add_ui(counter2, counter2, 1, GMP_RNDN);
	mpfr_div(counter, counter, counter2, GMP_RNDN);
	mpfr_log(counter, counter, GMP_RNDN); // log(epsilon / (2 * R) / (2 * (| s | +1)))
  
	mpfr_set_si(counter2, K, GMP_RNDN);
	mpfr_div_ui(counter2, counter2, v, GMP_RNDN);
	mpfr_log(counter2, counter2, GMP_RNDN); // log( K/v )
	mpfr_div(counter, counter, counter2, GMP_RNDN);
	mpfr_ceil(counter, counter);
	int numberOfTerms = mpfr_get_si(counter, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(counter3);
	return numberOfTerms;
}

int endPartition(int N, int v, int M1) {
	/*
		k = min( floor( v/M1 + 1 ), N - v )
	*/
	int K = fmin(floor(v / M1 + 1), N - v);
	return v + K - 1;
}

// Calculates the coefficients of the real part of our taylor expansion
void coefficientRealCalculator(mpfr_t sigma, int v, mpfr_t *coefficientRealValues,
													int numOfTerms, int mpfr_bits) {
	mpfr_t Vmpfr;
	mpfr_init2(Vmpfr, mpfr_bits);
	mpfr_set_ui(Vmpfr, v, GMP_RNDN);
	for (int i = 0; i < numOfTerms; i++) {
		mpfr_set(coefficientRealValues[i], sigma, GMP_RNDN);
		mpfr_div_si(coefficientRealValues[i], coefficientRealValues[i], 
											pow(-1,i+1)*(i+1), GMP_RNDN);
		mpfr_div(coefficientRealValues[i], coefficientRealValues[i], Vmpfr, GMP_RNDN);
		mpfr_mul_ui(Vmpfr, Vmpfr, v, GMP_RNDN);
	}

	mpfr_clear(Vmpfr);
}

// Calculates the coefficients of the imaginary part of our taylor expansion
void coefficientImagCalculator(mpfr_t t, int v, mpfr_t *coefficientImagValues, 
													int numOfTerms, int mpfr_bits) {
	mpfr_t pi, one, Vmpfr;
	/*
		alpha, beta, gamma will be used to help us calculate the imaginary part
		of our taylor expansion
	*/
	mpfr_t alpha, beta, gamma;
	mpfr_init2(alpha, mpfr_bits);
	mpfr_init2(beta, mpfr_bits);
	mpfr_init2(gamma, mpfr_bits);
	mpfr_init2(pi, mpfr_bits);
	mpfr_init2(one, mpfr_bits);
	mpfr_set(alpha, t, GMP_RNDN);
	mpfr_set_d(beta, 0, GMP_RNDN);
	mpfr_set_d(gamma, 0, GMP_RNDN);
	mpfr_set_ui(one, 1, GMP_RNDN);
	mpfr_const_pi(pi, GMP_RNDN); // pi
	/*
		alpha = t/(2*pi*v) mod 1
		beta = -alpha/(2*v) mod 1
		gamma = -2*beta/(3*v) mod 1
	*/
	mpfr_div(alpha, alpha, pi, GMP_RNDN);
	mpfr_div_ui(alpha, alpha, 2, GMP_RNDN);
	mpfr_div_ui(alpha, alpha, v, GMP_RNDN);
	mpfr_div_ui(beta, alpha, v, GMP_RNDN);
	mpfr_div_si(beta, beta, -2, GMP_RNDN);
	mpfr_fmod(alpha, alpha, one, GMP_RNDN); // alpha mod 1
	mpfr_mul_si(gamma, beta, -2, GMP_RNDN);
	mpfr_div_ui(gamma, gamma, 3, GMP_RNDN);
	mpfr_div_ui(gamma, gamma, v, GMP_RNDN);
	mpfr_fmod(beta, beta, one, GMP_RNDN); // beta mod 1
	if (mpfr_cmp_ui(beta, 0) < 0) {
		mpfr_add_ui(beta, beta, 1, GMP_RNDN);
	}
	mpfr_t temp; // temporary variable to help us with our calculations
	mpfr_init2(temp, mpfr_bits);
	mpfr_set(temp, gamma, GMP_RNDN); // we want to use the unmodded gamma value
	mpfr_mul_ui(temp, temp, 3, GMP_RNDN);
	mpfr_fmod(gamma, gamma, one, GMP_RNDN); // gamma mod 1
	// Store the values of the coefficients into the array
	mpfr_set(coefficientImagValues[0], alpha, GMP_RNDN);
	mpfr_mul_si(coefficientImagValues[0], coefficientImagValues[0], -1, GMP_RNDN);
	mpfr_set(coefficientImagValues[1], beta, GMP_RNDN);
	mpfr_mul_si(coefficientImagValues[1], coefficientImagValues[1], -1, GMP_RNDN);
	mpfr_set(coefficientImagValues[2], gamma, GMP_RNDN);
	mpfr_mul_si(coefficientImagValues[2], coefficientImagValues[2], -1, GMP_RNDN);
	/*
		We stored the previous values manually. The rest of the coefficients will
		use gamma for calculations.
	*/
	mpfr_init2(Vmpfr, mpfr_bits);
	mpfr_set_ui(Vmpfr, v, GMP_RNDN);
	for (int i = 3; i < numOfTerms; i++) {
		mpfr_div_si(coefficientImagValues[i], temp, 
											pow(-1,i+1) * (i+1), GMP_RNDN);
		mpfr_div(coefficientImagValues[i], coefficientImagValues[i], Vmpfr, GMP_RNDN);
		mpfr_mul_ui(Vmpfr, Vmpfr, v, GMP_RNDN);
	}

	// clear mpfr variables
	mpfr_clear(alpha);
	mpfr_clear(beta);
	mpfr_clear(gamma);
	mpfr_clear(one);
	mpfr_clear(pi);
	mpfr_clear(temp);
	mpfr_clear(Vmpfr);
}

void taylorSeries(mpfr_t sigma, mpfr_t t, int k, int numOfTerms,
					mpfr_t realTaylorSeries, mpfr_t imagTaylorSeries, 
					mpfr_t *coefficientRealValues, mpfr_t *coefficientImagValues, int mpfr_bits) {
	// initialize mpfr variables
	mpfr_t counter, counter2;
	mpfr_t pi;
	mpfr_init2(pi, mpfr_bits);
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_set_ui(counter, 0, GMP_RNDN);
	mpfr_set_ui(counter2, 0, GMP_RNDN);
	mpfr_const_pi(pi, GMP_RNDN); // pi
	
	mpfr_t Kmpfr;
	mpfr_init2(Kmpfr, mpfr_bits);
	mpfr_set_ui(Kmpfr, k, GMP_RNDN);
	//for (int i = 0; i < 15; i++) {
	for (int i = 0; i < numOfTerms; i++) {
		mpfr_mul(counter, coefficientRealValues[i], Kmpfr, GMP_RNDN);
		mpfr_mul(counter2, coefficientImagValues[i], Kmpfr, GMP_RNDN);
		mpfr_mul_ui(Kmpfr, Kmpfr, k, GMP_RNDN);
		/*
			Real part of the taylor series
			-\sigma * \sum_{i=1}^{numOfTerms} (-1)^{i-1} * k^{i} / (i * v^{i})
		*/
		mpfr_add(realTaylorSeries, realTaylorSeries, counter, GMP_RNDN);
		/*
			Imaginary part of the taylor series
			-t * \sum_{i=1}^{numOfTerms} (-1)^{i-1} * k^{i} / (i * v^{i})
		*/
		mpfr_add(imagTaylorSeries, imagTaylorSeries, counter2, GMP_RNDN);
	}
	/*
		We multiply the imaginary part of the taylor series expansion
		by 2*pi because we divided its coefficients by 2*pi in the calculations
	*/
	mpfr_mul(imagTaylorSeries, imagTaylorSeries, pi, GMP_RNDN);
	mpfr_mul_ui(imagTaylorSeries, imagTaylorSeries, 2, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(pi);
	mpfr_clear(Kmpfr);
}

/*
	Note that all of the places where there are comments for log( 1+ K/v ) we are actually using 
	the taylor approximation of log( 1 + K/v ) calculated to a predetermined number of terms
*/
void sumTerms(mpfr_t sigma, mpfr_t t, mpfr_t realTaylorSeries, 
					mpfr_t imagTaylorSeries, mpfr_t realTerm, mpfr_t imagTerm, int mpfr_bits) {
	// initialize mpfr variables
	mpfr_t counter, counter2;
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_set_d(counter, 0, GMP_RNDN);
	mpfr_set_d(counter2, 0, GMP_RNDN);
	/*
		Real part of the term to be summed by partialSum2MPFR:
		exp( -sigma * log( 1 + K/v)  ) * cos( -t*log( 1 + K/v ) ) 
	*/
	mpfr_cos(counter, imagTaylorSeries, GMP_RNDN);
	mpfr_exp(counter2, realTaylorSeries, GMP_RNDN);
	mpfr_mul(counter, counter, counter2, GMP_RNDN);
	mpfr_add(realTerm, realTerm, counter, GMP_RNDN);
	/*
		Imaginary part of the term to be summed by partialSum2MPFR:
		exp( -sigma * log( 1 + K/v)  ) * sin( -t*log( 1 + K/v ) )
	*/
	mpfr_sin(counter, imagTaylorSeries, GMP_RNDN);
	mpfr_mul(counter, counter, counter2, GMP_RNDN);
	mpfr_add(imagTerm, imagTerm, counter, GMP_RNDN);

	// clear mpfr variables
	mpfr_clear(counter);
	mpfr_clear(counter2);
}

void partialSum2MPFR(mpfr_t sigma, mpfr_t t, int M1, int N, mpfr_t rresult, 
									mpfr_t iresult, mpfr_t epsilon, int mpfr_bits) {
	int v = M1; // initial starting point of stage2
	// initialize counters for calculations
	mpfr_t counter, counter2, counter3, counter4, counter5, counter6;
	mpfr_init2(counter, mpfr_bits);
	mpfr_init2(counter2, mpfr_bits);
	mpfr_init2(counter3, mpfr_bits);
	mpfr_init2(counter4, mpfr_bits);
	mpfr_init2(counter5, mpfr_bits);
	mpfr_init2(counter6, mpfr_bits);
	mpfr_set_d(counter, 0, GMP_RNDN);
	mpfr_set_d(counter2, 0, GMP_RNDN);
	mpfr_set_d(counter5, 0, GMP_RNDN);
	mpfr_set_d(counter6, 0, GMP_RNDN);
	// initialize real part of taylor series
	mpfr_t realTaylorSeries;
	mpfr_init2(realTaylorSeries, mpfr_bits);
	mpfr_set_d(realTaylorSeries, 0, GMP_RNDN);
	// initialize imaginary part of taylor series
	mpfr_t imagTaylorSeries;
	mpfr_init2(imagTaylorSeries, mpfr_bits);
	mpfr_set_d(imagTaylorSeries, 0, GMP_RNDN);
	// initialize real part of term to sum
	mpfr_t realTerm;
	mpfr_init2(realTerm, mpfr_bits);
	// initialize imaginary part of taylor series
	mpfr_t imagTerm;
	mpfr_init2(imagTerm, mpfr_bits);

    // initialize array for mpfr variables
	mpfr_t *coefficientRealValues;
	mpfr_t *coefficientImagValues;
	coefficientRealValues = new mpfr_t[100];
	coefficientImagValues = new mpfr_t[100];
	for (int i = 0; i < 100; i++) {
		mpfr_init2(coefficientRealValues[i], mpfr_bits);
		mpfr_init2(coefficientImagValues[i], mpfr_bits);
	}
	while (v < N) {
		int numberOfTerms = numOfTerms(sigma, t, epsilon, v, M1, N, mpfr_bits);
		mpfr_set_ui(counter, 0, GMP_RNDN);
		mpfr_set_ui(counter2, 0, GMP_RNDN);
		mpfr_set_ui(counter3, v, GMP_RNDN);
		mpfr_set_ui(counter4, v, GMP_RNDN);
		/*
			The factored out real part
			exp( -sigma * log( v ) ) * cos( -t*log( v ) )
		*/
		mpfr_log(counter3, counter3, GMP_RNDN);
		mpfr_mul_si(counter3, counter3, -1, GMP_RNDN);
		mpfr_mul(counter3, counter3, t, GMP_RNDN);
		mpfr_cos(counter5, counter3, GMP_RNDN); // cos( -t*log( v ) )

		mpfr_log(counter4, counter4, GMP_RNDN);
		mpfr_mul(counter4, counter4, sigma, GMP_RNDN);
		mpfr_mul_si(counter4, counter4, -1, GMP_RNDN);
		mpfr_exp(counter6, counter4, GMP_RNDN); // exp( -sigma * log( v ) )
		mpfr_mul(counter4, counter6, counter5, GMP_RNDN);
		/*
			The factored out imaginary part
			exp( -sigma * log( v ) ) * sin( -t*log*( v ) )
		*/
		mpfr_sin(counter5, counter3, GMP_RNDN); // sin( -t*log*( v ) )
		mpfr_mul(counter5, counter5, counter6, GMP_RNDN);
		
		coefficientRealCalculator(sigma, v, coefficientRealValues, numberOfTerms, mpfr_bits);
		coefficientImagCalculator(t, v, coefficientImagValues, numberOfTerms, mpfr_bits);
		/*
		coefficientRealCalculator(sigma, v, coefficientRealValues,
			15, mpfr_bits);
		coefficientImagCalculator(t, v, coefficientImagValues,
			15, mpfr_bits);*/
		/*
			The sum itself
			\sum_{k=0}^{K-1} [ exp( -sigma * log( 1 + k/v)  ) * cos( -t*log( 1 + k/v ) ) +
			i * exp( -sigma * log( 1 + k/v)  ) * sin( -t*log( 1 + k/v ) ) ]
		*/
		int K = endPartition(N, v, M1) - v + 1;
		for (int k = 0; k <= K - 1; k++) {
			mpfr_set_ui(realTaylorSeries, 0, GMP_RNDN);
			mpfr_set_ui(imagTaylorSeries, 0, GMP_RNDN);
			mpfr_set_ui(realTerm, 0, GMP_RNDN);
			mpfr_set_ui(imagTerm, 0, GMP_RNDN);
			taylorSeries(sigma, t, k, numberOfTerms, realTaylorSeries, imagTaylorSeries,
											coefficientRealValues, coefficientImagValues, mpfr_bits);
            
            sumTerms(sigma, t, realTaylorSeries, imagTaylorSeries, realTerm, imagTerm, mpfr_bits);
			mpfr_add(counter, counter, realTerm, GMP_RNDN);
			mpfr_add(counter2, counter2, imagTerm, GMP_RNDN);
		}
		/*
			We need to multiply the factored out part by the sum calculated above. Since
			both numbers are complex numbers we multiply as follows:
			(counter4 + i*counter5)(counter + i*counter2) = 
						(counter*counter4 - counter2*counter5) + i*(counter2*counter4 + counter*counter5)
			The real part here will be added to our variable rresult and the imaginary part will be added
			to iresult.
		*/
		// We start with the real part of the above mentioned mulplication
		mpfr_mul(counter3, counter, counter4, GMP_RNDN);
		mpfr_mul(counter6, counter2, counter5, GMP_RNDN);
		mpfr_sub(counter3, counter3, counter6, GMP_RNDN);
		mpfr_add(rresult, rresult, counter3, GMP_RNDN);


		// Imaginary part of the above mentioned multiplcation
		mpfr_mul(counter3, counter2, counter4, GMP_RNDN);
		mpfr_mul(counter6, counter, counter5, GMP_RNDN);
		mpfr_add(counter3, counter3, counter6, GMP_RNDN);
		mpfr_add(iresult, iresult, counter3, GMP_RNDN);

		// starting value of next partition
		v = endPartition(N, v, M1) + 1;
	}
    // clear array of mpfr variables
	for (int i = 0; i < 25; i++) {
		mpfr_clear(coefficientRealValues[i]);
		mpfr_clear(coefficientImagValues[i]);
	}
	delete[] coefficientRealValues;
	delete[] coefficientImagValues;
	
	// clear mpfr variables
	mpfr_clear(realTaylorSeries);
	mpfr_clear(imagTaylorSeries);
	mpfr_clear(realTerm);
	mpfr_clear(imagTerm);
	mpfr_clear(counter);
	mpfr_clear(counter2);
	mpfr_clear(counter3);
	mpfr_clear(counter4);
	mpfr_clear(counter5);
	mpfr_clear(counter6);
}
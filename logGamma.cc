/*
 * This program calculates \log( \gamma ( 1-s ) )
 * where s = sigma + i * t. We use the asymptotic 
 * version of Stirling's formula to estimate.
 * We assume sigma > 0.
 */
#include "misc.h"
#include "zeros.h"

int endSum(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, int mpfr_bits);
void argument(mpfr_t sigma, mpfr_t t, mpfr_t arg);
void powerOfS(mpfr_t sigma, mpfr_t t, int power, mpfr_t realPower, mpfr_t imagPower, int mpfr_bits);
void firstTerm(mpfr_t sigma, mpfr_t t, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits);
void secondTerm(mpfr_t sigma, mpfr_t t, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits);
void thirdTerm(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits);

int main() {
    // determines number of bits the answer will have
    int mpfr_bits = 300;
    // initialize mpfr variables
    mpfr_t sigma, t, epsilon, logGammaReal, logGammaImag;
    mpfr_init2(sigma,mpfr_bits);
    mpfr_init2(t,mpfr_bits);
    mpfr_init2(epsilon, mpfr_bits);
    mpfr_init2(logGammaReal, mpfr_bits);
    mpfr_init2(logGammaImag, mpfr_bits);
    mpfr_set_d(sigma, 1/2, GMP_RNDN); // sigma
    mpfr_set_ui(t, 100, GMP_RNDN); // t
    mpfr_set_d(epsilon, 1E-80, GMP_RNDN); // epsilon
    mpfr_set_ui(logGammaReal, 0, GMP_RNDN);
    mpfr_set_ui(logGammaImag, 0, GMP_RNDN);
    
    //cout << endSum(sigma, t, epsilon, mpfr_bits) << endl;
    //logGamma(sigma, t, epsilon, logGammaReal, logGammaImag, mpfr_bits);
    
    mpfr_t realPower, imagPower;
    mpfr_init2(realPower, mpfr_bits);
    mpfr_init2(imagPower, mpfr_bits);
    mpfr_set_ui(realPower, 0, GMP_RNDN); // t
    mpfr_set_ui(imagPower, 0, GMP_RNDN); // t
    //powerOfS(sigma,t,10,realPower, imagPower, mpfr_bits);
     thirdTerm(sigma, t, epsilon, logGammaReal, logGammaImag, mpfr_bits);
    /*cout << endl;
    mpfr_out_str(stdout, 10, mpfr_bits, logGammaReal, GMP_RNDN);
    cout << endl << "******************" << endl;
    mpfr_out_str(stdout, 10, mpfr_bits, logGammaImag, GMP_RNDN);
    cout << endl;*/
    thirdTerm(sigma, t, epsilon, logGammaReal, logGammaImag, mpfr_bits);
    // clear mprf variables
    mpfr_clear(sigma);
    mpfr_clear(t);
	mpfr_clear(epsilon);
    mpfr_clear(logGammaReal);
    mpfr_clear(logGammaImag);
    
    return 0;
}

/*
 * Calculates the number of terms needed in the sum to bound the error term by epsilon
 */
int endSum(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, int mpfr_bits) {
    int m = 0;
    int arrIndex;
    mpfr_t remainder, bnum, bden, bernoulli;
    mpfr_t counter1,counter2;
    mpfr_init2(remainder, mpfr_bits);
    mpfr_init2(bnum, mpfr_bits);
    mpfr_init2(bden, mpfr_bits);
    mpfr_init2(bernoulli, mpfr_bits);
    mpfr_init2(counter1, mpfr_bits);
    mpfr_init2(counter2, mpfr_bits);
    mpfr_set_ui(remainder, 1000, GMP_RNDN); // set to a number clearly larger than epsilon
    mpfr_set_ui(counter1, 0, GMP_RNDN);
    mpfr_set_ui(counter2, 0, GMP_RNDN);
    mpfr_set_ui(bnum, 0, GMP_RNDN);
    mpfr_set_ui(bden, 0, GMP_RNDN);
    mpfr_set_ui(bernoulli, 0, GMP_RNDN);
    while (mpfr_cmp(remainder, epsilon) > 0) {
        m++;
        arrIndex = m;
        mpfr_set_z(bnum, numerator_arr[arrIndex].get_mpz_t(), GMP_RNDN);
        mpfr_set_z(bden, denominator_arr[arrIndex].get_mpz_t(), GMP_RNDN);
        mpfr_div(bernoulli, bnum, bden, GMP_RNDN); // B_{2k}
        mpfr_abs(bernoulli, bernoulli, GMP_RNDN);
        mpfr_div_ui(bernoulli, bernoulli, (2*m)+1, GMP_RNDN);
        mpfr_div_ui(bernoulli, bernoulli, (2*m)+2, GMP_RNDN);
        mpfr_mul(counter1, sigma, sigma, GMP_RNDN);
        mpfr_mul(counter2, t, t, GMP_RNDN);
        mpfr_add(counter1, counter1, counter2, GMP_RNDN);
        mpfr_sqrt(counter1, counter1, GMP_RNDN);
        mpfr_pow_ui(counter1, counter1, (2*m)+1, GMP_RNDN);
        mpfr_div(bernoulli, bernoulli, counter1, GMP_RNDN);
        argument(sigma, t, counter1);
        mpfr_div_ui(counter1, counter1, 2, GMP_RNDN);
        mpfr_cos(counter1, counter1, GMP_RNDN);
        mpfr_pow_ui(counter1, counter1, (2*m)+2, GMP_RNDN);
        mpfr_div(remainder, bernoulli, counter1, GMP_RNDN);
    }
    
    // clear mpfr variables
    mpfr_clear(remainder);
    
    return m;
}

/*
 * Calculates the argument of s
 */
void argument(mpfr_t sigma, mpfr_t t, mpfr_t arg) {
    mpfr_div(arg, t, sigma, GMP_RNDN);
    mpfr_atan(arg, arg, GMP_RNDN);
}

/*
 * Calculates the first term of the expansion:
 * (s - \frac{1}{2}) * \log( s )
 */
void firstTerm(mpfr_t sigma, mpfr_t t, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits) {
    mpfr_t realPart, imagPart;
    mpfr_t counter1,counter2, counter3;
    mpfr_init2(realPart, mpfr_bits);
    mpfr_init2(imagPart, mpfr_bits);
    mpfr_init2(counter1, mpfr_bits);
    mpfr_init2(counter2, mpfr_bits);
    mpfr_init2(counter3, mpfr_bits);
    mpfr_set(realPart, sigma, GMP_RNDN);
    mpfr_set(imagPart, t, GMP_RNDN);
    mpfr_set_ui(counter1, 0, GMP_RNDN);
    mpfr_set_ui(counter2, 0, GMP_RNDN);
    mpfr_set_ui(counter3, 0, GMP_RNDN);
    /*
     * First we calculate log(s)
     */
    mpfr_mul(counter1, sigma, sigma, GMP_RNDN);
    mpfr_mul(imagPart, imagPart, t, GMP_RNDN);
    mpfr_add(counter1, counter1, imagPart, GMP_RNDN);
    mpfr_sqrt(counter1, counter1, GMP_RNDN);
    mpfr_log(counter1, counter1, GMP_RNDN); // real part of log(s)
    argument(sigma, t, imagPart); // imaginary part of log(s)
    /*
     * Calculates: s - \frac{1}{2}
     */
    mpfr_sub_d(realPart, realPart, 0.5, GMP_RNDN);
    /*
     * Complex multiplication of the above parts:
     * (s - \frac{1}{2}) * log(s)
     */
    // real part
    mpfr_mul(counter2, realPart, counter1, GMP_RNDN);
    mpfr_mul(counter3, t, imagPart, GMP_RNDN);
    mpfr_sub(counter2, counter2, counter3, GMP_RNDN);
    // imag part
    mpfr_mul(counter1, counter1, t, GMP_RNDN);
    mpfr_mul(counter3, realPart, imagPart, GMP_RNDN);
    mpfr_add(counter3, counter3, counter1, GMP_RNDN);
    // set result to real and imaginary parts
    mpfr_add(logGammaReal, logGammaReal, counter2, GMP_RNDN);
    mpfr_add(logGammaImag, logGammaImag, counter3, GMP_RNDN);
    
    // clear mpfr variables
    mpfr_clear(realPart);
    mpfr_clear(imagPart);
    mpfr_clear(counter1);
    mpfr_clear(counter2);
    mpfr_clear(counter3);
}

/*
 * Calculates: -s + \frac{1}{2} * \log( 2 * pi )
 */
void secondTerm(mpfr_t sigma, mpfr_t t, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits) {
    mpfr_t pi;
    mpfr_init2(pi, mpfr_bits);
    mpfr_const_pi(pi, GMP_RNDN); // pi
    // updates logGamma
    mpfr_sub(logGammaReal, logGammaReal, sigma, GMP_RNDN);
    mpfr_sub(logGammaImag, logGammaImag, t, GMP_RNDN);
    mpfr_mul_ui(pi,pi,2,GMP_RNDN);
    mpfr_log(pi, pi, GMP_RNDN);
    mpfr_div_ui(pi, pi, 2, GMP_RNDN);
    // updates logGamma
    mpfr_add(logGammaReal, logGammaReal, pi, GMP_RNDN);
    
    // clear mpfr variables
    mpfr_clear(pi);
}

/*
 * Calcualtes s ^ power
 */
void powerOfS(mpfr_t sigma, mpfr_t t, int power, mpfr_t realPower, mpfr_t imagPower, int mpfr_bits) {
    mpfr_t counter1, counter2, counter3;
    mpfr_init2(counter1, mpfr_bits);
    mpfr_init2(counter2, mpfr_bits);
    mpfr_init2(counter3, mpfr_bits);
    mpfr_set_ui(counter1, 0, GMP_RNDN);
    mpfr_set_ui(counter2, 0, GMP_RNDN);
    mpfr_set_ui(counter3, 0, GMP_RNDN);
    mpfr_set(realPower, sigma, GMP_RNDN);
    mpfr_set(imagPower, t, GMP_RNDN);
    for (int i = 1; i < power; i++) {
        // real part
        mpfr_mul(counter1, realPower, sigma, GMP_RNDN);
        mpfr_mul(counter2, imagPower, t, GMP_RNDN);
        mpfr_sub(counter1, counter1, counter2, GMP_RNDN);
        // imag part
        mpfr_mul(counter2, imagPower, sigma, GMP_RNDN);
        mpfr_mul(counter3, realPower, t, GMP_RNDN);
        mpfr_add(counter2, counter2, counter3, GMP_RNDN);
        // updates realPower and imagPower
        mpfr_set(realPower, counter1, GMP_RNDN);
        mpfr_set(imagPower, counter2, GMP_RNDN);
    }
    
    // clear mpfr variables
    mpfr_clear(counter1);
    mpfr_clear(counter2);
    mpfr_clear(counter3);
}

/*
 * Calculates:
 * \sum_{k = 1}^{m} s^{1-2k} * (2k)^{-1} * (2k-1)^{-1} * B_{2k}
 */
void thirdTerm(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits) {
    int m = endSum(sigma, t, epsilon, mpfr_bits);
    cout << m << endl;
    int arrIndex;
    mpfr_t realPowerOfS, imagPowerOfS, bnum, bden, bernoulli;
    mpfr_t counter1, counter2;
    mpfr_init2(realPowerOfS, mpfr_bits);
    mpfr_init2(imagPowerOfS, mpfr_bits);
    mpfr_init2(bnum, mpfr_bits);
    mpfr_init2(bden, mpfr_bits);
    mpfr_init2(bernoulli, mpfr_bits);
    mpfr_init2(counter1, mpfr_bits);
    mpfr_init2(counter2, mpfr_bits);
    mpfr_set_ui(realPowerOfS, 0, GMP_RNDN);
    mpfr_set_ui(imagPowerOfS, 0, GMP_RNDN);
    mpfr_set_ui(counter1, 0, GMP_RNDN);
    mpfr_set_ui(counter2, 0, GMP_RNDN);
    /*
     * After calling powerOfS, we multiply by the conjugate to break up the real and imaginary
     * parts of the sum.
     */
    for (int k = 1; k <= m; k++) {
        arrIndex = k - 1;
        mpfr_set_z(bnum, numerator_arr[arrIndex].get_mpz_t(), GMP_RNDN);
        mpfr_set_z(bden, denominator_arr[arrIndex].get_mpz_t(), GMP_RNDN);
        mpfr_div(bernoulli, bnum, bden, GMP_RNDN); // B_{2k}
        mpfr_out_str(stdout, 10, mpfr_bits, bernoulli, GMP_RNDN);
        cout << endl << "******************" << endl;

        mpfr_div_ui(bernoulli, bernoulli, 2*k, GMP_RNDN);
        mpfr_div_ui(bernoulli, bernoulli, (2*k)-1, GMP_RNDN);
        powerOfS(sigma, t, (2*k)-1, realPowerOfS, imagPowerOfS, mpfr_bits);
        // real part of sum
        mpfr_mul(counter1, realPowerOfS, realPowerOfS, GMP_RNDN);
        mpfr_mul(counter2, imagPowerOfS, imagPowerOfS, GMP_RNDN);
        mpfr_add(counter2, counter1, counter2, GMP_RNDN);
        mpfr_mul(counter1, bernoulli, realPowerOfS, GMP_RNDN);
        mpfr_div(counter1, counter1, counter2, GMP_RNDN);
        
        // updates real part of logGamma
        mpfr_add(logGammaReal, logGammaReal, counter1, GMP_RNDN);
        
        
        
        // imag part of sum
        mpfr_mul(counter1, bernoulli, imagPowerOfS, GMP_RNDN);
        mpfr_div(counter1, counter1, counter2, GMP_RNDN);
        // updates imaginary part of logGamma
        mpfr_sub(logGammaImag, logGammaImag, counter1, GMP_RNDN);
    }
    
    
    mpfr_out_str(stdout, 10, mpfr_bits, logGammaImag, GMP_RNDN);
    cout << endl;
    
    // clear mpfr variables
    mpfr_clear(realPowerOfS);
    mpfr_clear(imagPowerOfS);
    mpfr_clear(counter1);
    mpfr_clear(counter2);
}

void logGamma(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t logGammaReal, mpfr_t logGammaImag, int mpfr_bits) {
	firstTerm(sigma, t, logGammaReal, logGammaImag, mpfr_bits);
	secondTerm(sigma, t, logGammaReal, logGammaImag, mpfr_bits);
	thirdTerm(sigma, t, epsilon, logGammaReal, logGammaImag, mpfr_bits);
}
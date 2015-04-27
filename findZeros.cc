#include "main.h"
#include "zeros.h"

void valueCalc(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t theta,
          mpfr_t zetaReal, mpfr_t zetaImag, mpfr_t value, int mpfr_bits);
void zeroSearch(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t theta,
           mpfr_t zetaReal, mpfr_t zetaImag, mpfr_t beginInterval,
           mpfr_t endInterval, mpfr_t zero, int mpfr_bits);

int main(int argc, char * argv[])
{
  int mpfr_bits = 300;
  mpfr_t sigma, t, value, zero, epsilon, zetaReal, zetaImag;
  mpfr_t theta, beginInterval, endInterval;
  mpfr_init2(sigma, mpfr_bits);
  mpfr_init2(t, mpfr_bits);
  mpfr_init2(theta, mpfr_bits);
  mpfr_init2(beginInterval, mpfr_bits);
  mpfr_init2(endInterval, mpfr_bits);
  mpfr_init2(value, mpfr_bits);
  mpfr_init2(zero, mpfr_bits);
    mpfr_init2(zetaReal, mpfr_bits);
    mpfr_init2(zetaImag, mpfr_bits);
     mpfr_init2(epsilon, mpfr_bits);
  mpfr_set_d(theta, 0, GMP_RNDN);
  mpfr_set_d(value, 0, GMP_RNDN);
  mpfr_set_d(zero, 0, GMP_RNDN);
    mpfr_set_d(zetaReal, 0, GMP_RNDN);
    mpfr_set_d(zetaImag, 0, GMP_RNDN);
       mpfr_set_d(epsilon, 0, GMP_RNDN);

  string begin_string = "";
  string end_string = "";
  string sigma_string = "";
  string t_string = "";
  while(1) {
       enum {sigma, t, beginInterval, endInterval};

       static struct option long_options[] =
           {
               {"sigma", required_argument, 0, beginInterval},
               {"t", required_argument, 0, beginInterval},
               {"begin", required_argument, 0, beginInterval},
               {"end", required_argument, 0, endInterval},
               {0, 0, 0, 0}
           };

       int option_index = 0;
       int cc = getopt_long(argc, argv, "", long_options, &option_index);
       if (cc == -1)  break;

       switch(long_options[option_index].val) {
           case beginInterval:
               begin_string = optarg;
               break;
           case endInterval:
               end_string = optarg;
               break;
       }
   }
   mpfr_set_str(sigma, sigma_string.c_str(), 10, GMP_RNDN);
   mpfr_set_str(t, t_string.c_str(), 10, GMP_RNDN);
   mpfr_set_str(beginInterval, begin_string.c_str(), 10, GMP_RNDN);
   mpfr_set_str(endInterval, end_string.c_str(), 10, GMP_RNDN);

   valueCalc(sigma, beginInterval, epsilon, theta, zetaReal, zetaImag, value, mpfr_bits);
   if (mpfr_cmp_ui(value, 0) > 0) {
     valueCalc(sigma, endInterval, epsilon, theta, zetaReal, zetaImag, value, mpfr_bits);
     if (mpfr_cmp_ui(value, 0) < 0) {

     } else {
       cout << "No sign change in interval." << endl;
     }
   } else if (mpfr_cmp_ui(value, 0) < 0) {
     valueCalc(sigma, endInterval, epsilon, theta, zetaReal, zetaImag, value, mpfr_bits);
     if (mpfr_cmp_ui(value, 0) > 0) {

     } else {
       cout << "No sign change in interval." << endl;
     }
   }

   // clear mpfr variables
   mpfr_clear(sigma);
   mpfr_clear(t);
   mpfr_clear(theta);
   mpfr_clear(beginInterval);
   mpfr_clear(endInterval);
   mpfr_clear(value);
   return 0;
}

void valueCalc(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t theta,
                mpfr_t zetaReal, mpfr_t zetaImag, mpfr_t value, int mpfr_bits) {
  mpfr_t valueImag;
  mpfr_t counter1, counter2, counter3;
  mpfr_init2(valueImag, mpfr_bits);
  mpfr_init2(counter1, mpfr_bits);
  mpfr_init2(counter2, mpfr_bits);
  mpfr_init2(counter3, mpfr_bits);
  mpfr_set_d(counter1, 0, GMP_RNDN);
  mpfr_set_d(counter2, 0, GMP_RNDN);
  mpfr_set_d(counter3, 0, GMP_RNDN);
  mpfr_set_d(valueImag, 0, GMP_RNDN);

  calcTheta(sigma, t, epsilon, theta, mpfr_bits);
  mpfr_cos(value, theta, GMP_RNDN);
  mpfr_sin(valueImag, theta, GMP_RNDN);

  int L1 = inequalityOfL1(sigma, t, epsilon, mpfr_bits);
  int N = inequalityOfN(sigma, t, L1, mpfr_bits);
  int M1 = endStage1(sigma, t, N, mpfr_bits);
  partial_sum_mpfr(sigma, t, M1, zetaReal, zetaImag, mpfr_bits);
  partialSum2MPFR(sigma, t, M1, N, zetaReal, zetaImag, epsilon, mpfr_bits);
    EMsum(sigma, t, N, L1, zetaReal, zetaImag,mpfr_bits);

  mpfr_mul(counter1, value, zetaReal, GMP_RNDN);
  mpfr_mul(counter2, valueImag, zetaImag, GMP_RNDN);
  mpfr_sub(counter3, counter1, counter2, GMP_RNDN);

  mpfr_mul(counter1, value, zetaImag, GMP_RNDN);
  mpfr_mul(counter2, valueImag, zetaReal, GMP_RNDN);
  mpfr_add(valueImag, counter1, counter2, GMP_RNDN);
  mpfr_set(value, counter3, GMP_RNDN);

  // clear mpfr variables
  mpfr_clear(valueImag);
  mpfr_clear(counter1);
  mpfr_clear(counter2);
  mpfr_clear(counter3);
}

void zeroSearch(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, mpfr_t theta,
                mpfr_t zetaReal, mpfr_t zetaImag, mpfr_t beginInterval,
                          mpfr_t endInterval, mpfr_t zero, int mpfr_bits) {
  mpfr_t intervalHalf, counter1, counter2;
  mpfr_init2(intervalHalf, mpfr_bits);
  mpfr_init2(counter1, mpfr_bits);
  mpfr_init2(counter2, mpfr_bits);
  mpfr_set_d(intervalHalf, 0, GMP_RNDN);
  mpfr_set_d(counter1, 0, GMP_RNDN);
  mpfr_set_d(counter2, 0, GMP_RNDN);
  mpfr_add(intervalHalf, beginInterval, endInterval, GMP_RNDN);
  mpfr_div_ui(intervalHalf, intervalHalf, 2, GMP_RNDN);

  int L1 = inequalityOfL1(sigma, t, epsilon, mpfr_bits);
  int N = inequalityOfN(sigma, t, L1, mpfr_bits);
  int M1 = endStage1(sigma, t, N, mpfr_bits);
  partial_sum_mpfr(sigma, t, M1, zetaReal, zetaImag, mpfr_bits);
  partialSum2MPFR(sigma, t, M1, N, zetaReal, zetaImag, epsilon, mpfr_bits);
    EMsum(sigma, t, N, L1, zetaReal, zetaImag,mpfr_bits);
    
    bool flag = true;
  while(flag) {
    valueCalc(sigma, intervalHalf, epsilon, theta, zetaReal, zetaImag, counter1, mpfr_bits);
    valueCalc(sigma, beginInterval, epsilon, theta, zetaReal, zetaImag, counter2, mpfr_bits);
    if (mpfr_cmp_ui(counter2, 0) > 0 && mpfr_cmp_ui(counter1, 0) > 0) {
      mpfr_set(beginInterval, intervalHalf, GMP_RNDN);
    } else if (mpfr_cmp_ui(counter2, 0) > 0 && mpfr_cmp_ui(counter1, 0) < 0) {
      mpfr_set(endInterval, intervalHalf, GMP_RNDN);
    } else if (mpfr_cmp_ui(counter2, 0) < 0 && mpfr_cmp_ui(counter1, 0) > 0) {
      mpfr_set(endInterval, intervalHalf, GMP_RNDN);
    } else if (mpfr_cmp_ui(counter2, 0) < 0 && mpfr_cmp_ui(counter1, 0) < 0) {
      mpfr_set(beginInterval, intervalHalf, GMP_RNDN);
    } else if (mpfr_cmp_ui(counter2, 0) == 0) {
      flag = false;
    }
  }
}

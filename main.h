#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <cmath>
#include <gmp.h>
#include <mpfr.h>

using namespace std;

// Call to stage 1
int endStage1(mpfr_t sigma, mpfr_t t, int N, int mpfr_bits);
void partial_sum_mpfr(mpfr_t sigma, mpfr_t t, int M1, mpfr_t rresult, mpfr_t iresult, int mpfr_bits);

// Call to stage 2
void partialSum2MPFR(mpfr_t sigma, mpfr_t t, int M1, int N, mpfr_t rresult, mpfr_t iresult,
						mpfr_t epsilon, int mpfr_bits);

// Inequalities
int inequalityOfN(mpfr_t sigma, mpfr_t t, int L1, int mpfr_bits);
int inequalityOfL1(mpfr_t sigma, mpfr_t t, mpfr_t epsilon, int mpfr_bits);

// Calculates the remainder terms of the Euler-Maclaurin sum
void EMsum(mpfr_t sigma, mpfr_t t, int N, int L1, mpfr_t rresult, mpfr_t iresult, int mpfr_bits);
# Euler Maclaurin Sum
Calculates the Riemann Zeta Function with arbitrary precison where s = sigma + i * t.

Use 

g++ main.cc EulerMacSum.cc Inequalities.cc stage1.cc stage2.cc -o main -lgmp -lmpfr

to compile in terminal. Use

./main --sigma 0.5 --t 1E6 --epsilon 1E-20 --stage1 --stage2

to run where values for sigma, t, and epsilon can be any value. The smaller the value of epsilon, the more accurate the result. These commands will change when the program is made more efficient and userfriendly. 

The next update will be to implement a calculation of mpfr_bits that depends on user supplied values of s and epsilon. This will also include a separate pre-defined function to calculate the number of bits used for modded variables in stage2.

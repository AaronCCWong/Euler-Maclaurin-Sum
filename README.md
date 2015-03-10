# Euler Maclaurin Sum
Calculates the Riemann Zeta Function with arbitrary precison where s = sigma + i * t.

Use 

g++ main.cc EulerMacSum.cc Inequalities.cc stage1.cc stage2.cc -o main -lgmp -lmpfr

to compile in terminal. Use

./main --sigma 0.5 --t 1E6 --epsilon 1E-20 --stage1 --stage2

to run where values for sigma, t, and epsilon can be any value. 
These commands will change when the program is made more efficient and userfriendly. 

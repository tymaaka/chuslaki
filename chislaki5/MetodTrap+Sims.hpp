//
//  MetodTrap+Sims.hpp
//  chislaki5
//
//  Created by Liza on 13.12.2023.
//

#ifndef MetodTrap_Sims_hpp
#define MetodTrap_Sims_hpp

#include <stdio.h>

#endif /* MetodTrap_Sims_hpp */

double F(double x);
double FindMetodTrap(double a, double b, double h, int n, double h2, int n2);
double FindMetodSimpsona(double a, double b, double h, int n);
double FindSumSimps(double a, double b, double h, int n);
double FindSumTrap(double a, double b, double h, int n);
double FindMetodSimpsona2( double(&f31)(double, double), double a1, double b1, double c1, double d1, int n1, int m1);
double FindSumSimps2( double a1, double b1, double c1, double d1,  int n1,  int m1);
double F31(double x, double y);

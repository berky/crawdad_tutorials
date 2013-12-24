#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>

double masses[] = {
  0.000000, // Ghost
  1.007825, // H
  4.002603, // He
  6.015123, // Li
  9.012182, // Be
  10.012937, // B
  12.000000, // C
  14.003074, // N
  15.994915, // O
  18.998403, // F
  19.992440 // Ne
};

double amu2kg = 1.6605402e-27;
double bohr2m = 0.529177249e-10;
double planck = 6.6260755e-34;
// double planck = 6.62606957e-34;
double pi = acos(-1.0);
double planckbar = planck/(2*pi);
double speed_of_light = 299792458;

double rot_constant = planck/(8*pi*pi*speed_of_light);

#endif /* CONSTANTS_HPP */

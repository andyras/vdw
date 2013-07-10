#ifndef CONSTANTS
#define CONSTANTS

const double h = 6.6260755e-34; // Planck's constant --- Joule.second
const double me = 9.1093897e-31; // electron rest mass --- kilograms
const double e0 = 8.854187816e-12; // permittivity of free vacuum --- coulomb^2.joule^-1.meter^-1
const double e = 1.60217733e-19; // proton charge --- coulomb
const double c = 2.99792458e8; // speed of light --- meter.second^-1
const double pi = 3.141592653589793; // unitless
const double ang_m = 1e10; // angstroms per meter

const double a = e * e / ( 2.0 * h * e0 * c ); // fine structure constant --- unitless
const double a0 = ( h / ( 2.0 * pi )) / ( me * c * a  ) * ang_m ; // bohr radius --- angstrom

#endif

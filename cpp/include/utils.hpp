// utils.hpp

/* A very basis header file that includes the standard C++ headers so that they
do not need to be individually added to all the other user-defined header files
*/

#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cfloat>
#include <ctime>
#include <complex>
#include <vector>

////////////////////////////////////////////////////////
// Simple, but non-standard functions and definitions //
////////////////////////////////////////////////////////

// Shorten complex type declaration
typedef std::complex<double> Complex;

// DEFINE imaginary number i as ii
const std::complex<double> I(0., 1.);

// Shorten vector declarations
typedef std::vector<double> Vector;
typedef std::vector<Complex> ComplexVector;

typedef std::vector<int> List;
typedef std::vector<ComplexVector> ComplexMatrix;
typedef std::vector<Vector> RealMatrix;
typedef std::vector<ComplexMatrix> ComplexTensor;
typedef std::vector<RealMatrix> RealTensor;

Vector downsample_base_2(Vector v, int sampleRate);
Vector unwrap_phase(const Vector& phase);

#endif

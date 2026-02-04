/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef Hmod_UtilVectors
#define Hmod_UtilVectors

#include "Global.h"
#include "util/errors.h"
#include "util/checkpoint_containers.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace OM {
namespace util {

/** Various utilities acting on vectors. */
namespace vectors {

  ///@brief Basic operations on std::vector
  //@{
  /// Scale all elements of a vector by a in-situ
  inline void scale (vector<double>& vec, double a) {
    for(size_t i = 0; i < vec.size(); ++i)
      vec[i] *= a;
  }
  
  /// Return sum of all elements
  inline double sum (const vector<double>& vec) {
    double r = 0.0;
    for(size_t i = 0; i < vec.size(); ++i)
      r += vec[i];
    return r;
  }
  
  /// Return sum of all elements
 //double sum (const gsl_vector *vec);
  
  /// Return arithmetric mean
  inline double mean (const vector<double>& vec){
    return sum(vec) / vec.size();
  }
  
  /// Add one vector into another (x += y)
  inline void addTo (vector<double>& x, vector<double>& y){
    assert( x.size() == y.size() );
    for( size_t i=0; i<x.size(); ++i ){
      x[i] += y[i];
    }
  }
  //@}
  
  
  ///@brief Comparissons on std::vector
  //@{
  /** Return true if, approximately, a == b
   *
   * In detail: true when (fabs(a-b) <= max(fabs(a),fabs(b)) * lim_fact). */
  inline bool approxEqual (const double a, const double b, const double lim_fact = 1e-6) {
    bool aE = (fabs(a-b) <= max(fabs(a),fabs(b)) * lim_fact);
#ifndef NDEBUG
    if (!aE){
        cerr<<"not approx equal: "<<a<<", "<<b<<endl;
    }
#endif
    return aE;
}
  
  /// Returns true if vec1 and vec2 have equal length and all elements are
  /// approximately equal (see approxEqual for double parameters).
  inline bool approxEqual (const vector<double>& vec1, const vector<double>& vec2, const double lim_fact) {
    if (vec1.size() != vec2.size())
      return false;
    for(size_t i = 0; i < vec1.size(); ++i) {
      if (!approxEqual (vec1[i], vec2[i], lim_fact))
        return false;
    }
    return true;
  }
  //@}
  
  /** Calculate what is essentially a discrete Fourier transform of log values.
   * 
   * Encoding is slightly different:
   *  FC[0] is the real component of the first component of the frequency
   *  domain over T, while FC[2*n-1], FC[2*n] are the real and complex
   * components of component n over T. (I.e. the first imaginary component,
   * which is always zero, is not stored, and all components are divided by T).
   * 
   * @param iArray Input values; a vector of size T.
   * @param FC Output values. Should already be allocated, with length
   *    2*N-1 for 1 ≤ N ≤ T. */
  inline void logDFT(const vector<double>& iArray, vector<double>& FC) {
    if (FC.size() > 2*iArray.size() || mod_nn(FC.size(), 2) == 0)
        throw TRACED_EXCEPTION_DEFAULT("Require DFT FC.size() to be 2*n-1 for n<=iArray.size()");
    
    size_t T = iArray.size();
    size_t N = (FC.size() + 1) / 2;
    
    double w = 2.0 * M_PI / T;
    FC.assign(FC.size(), 0.0);
    
    for(size_t t = 0; t < T; ++t) {
        double val = log( iArray[t] );
        FC[0] += val;
        for(size_t n = 1; n < N; ++n ){
            FC[2*n-1] += val * cos(w*n*t);
            FC[2*n  ] += val * sin(w*n*t);
        }
    }
    scale(FC, 1.0 / T);
  }

  /** Calculate Fourier Series coefficients for log values; slightly different
   * method to (log) DFT.
   * 
   * @param iArray Integrals of inputs on curve.
   * @param FC Fourier Series coefficients.
   */
  inline void logFourierCoefficients(const vector<double>& iArray, vector<double>& FC) {
    logDFT(iArray, FC);
    for(size_t i=1; i<FC.size(); ++i)
        FC[i] *= 2;
  }

  /** The inverse of logDFT (or an approximation, when N&lt;T or
   * tArray.size() ≠ T). Result may also be rotated.
   * 
   * (This was called calcExpFourierSeries, and does essentially the same
   * thing.)
   *
   * @param tArray Array to fill with exponated values from Fourier series.
   * Length should already be set. Need not have the same length as the
   * array used to calculate FC.
   * @param FC Fourier coefficients (a0, a1,b1, a2,b2, ...); can be any length
   * so long as it is odd.
   * @param rAngle Angle to rotate generated series by in radians: [0,2π] */
  inline void expIDFT(std::vector< double >& tArray, const std::vector< double >& FC, double rAngle ){
    if (mod_nn(FC.size(), 2) == 0)
        throw TRACED_EXCEPTION_DEFAULT("The number of Fourier coefficents should be odd.");
    
    size_t T2 = tArray.size();
    size_t N = (FC.size() + 1) / 2;
    
    double w = 2.0 * M_PI / T2;
    
    // Calculate inverse discrete Fourier transform
    // TODO: This may not interpolate sensibly. See for example
    // https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
    for( size_t t = 0; t < T2; t++ ){
        double temp = FC[0];
        double wt = w*t - rAngle;
        for(size_t n = 1; n < N; ++n) {
            temp += FC[2*n-1]*cos(n*wt) + FC[2*n]*sin(n*wt);
        }
        tArray[t] = exp(temp);
    }
  }
}

/// Utility to print a vector (operator must be in namespace)
template<class T>
ostream& operator<< (ostream& out, vector<T> vec) {
  out << '[';
  if (vec.size())
    out << vec[0];
  for(size_t i = 1; i < vec.size(); ++i)
    out << ", " << vec[i];
  out << ']';
  return out;
}

} }
#endif

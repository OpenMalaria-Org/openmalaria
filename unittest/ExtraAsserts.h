/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
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

/** This module adds some extra *TS_ASSERT_* macros for CxxTest, and a test
 * suite to check these work as intended.
 * 
 * General methodology: use TS_* for most asserts; ETS_* when a failed
 * assertion should stop the test (prevent invalid defreferences, etc.).
 * #define a new test rather than use [E]TSM_* asserts.
 * 
 * TS_ASSERT_IS_NAN (x):
 * Simply asserts that x != x.
 * 
 * TS_ASSERT_APPROX (x,y), TS_ASSERT_APPROX_TOL (x,y, rel, abs):
 * Check that doubles x and y are approximately equal, optionally specifying
 * relative and absolute tolerances (both default to 1e-7).
 * See ExtraAsserts:approx for details of how this works.
 * 
 * TS_ASSERT_VECTOR_APPROX (x,y):, TS_ASSERT_VECTOR_APPROX_TO (x,y, rel, abs):
 * Approximate equality test for std::vector<T>, gsl_vector or gsl_matrix.
 * (The gsl_* types are only supported if the appropriate gsl headers are
 * imported first. This means tests including them should appear before this
 * file on the cxxtestget.p? command line if this test suite is generated).
 * First asserts the dimensions are the same, then recursively checks values
 * with TS_ASSERT_APPROX.
 * Caveat: error messages aren't reported in the nicest way (they could be with
 * a little extra code).
 * 
 * ETS_* versions are included, but not [E]TSM_* versions since extending the
 * asserts like this gives nicer results than using functions to run tests,
 * passing messages. */
#ifndef Hmod_ExtraAsserts
#define Hmod_ExtraAsserts

#include "Global.h"
#include <cxxtest/TestSuite.h>
#include <limits>
#include <cmath>

// Note: cxxtest already defines those macros, use a trailing underscore to avoid redefinition
// Check for NaN:
#   define TS_ASSERT_IS_NAN_(x) _TS_ASSERT_DIFFERS(__FILE__,__LINE__,x,x)
#   define TSM_ASSERT_IS_NAN_(m,x) _TSM_ASSERT_DIFFERS(__FILE__,__LINE__,m,x,x)


// TS_ASSERT_APPROX
#define DEF_REL_PRECISION 1e-7
#define DEF_ABS_PRECISION 1e-7

#   define ___ETS_ASSERT_APPROX(f,l,x,y,r,a) ExtraAsserts::assert_approx( (f), (l), (x), (y), (r), (a) )
#   define ___TS_ASSERT_APPROX(f,l,x,y,r,a) { _TS_TRY { ___ETS_ASSERT_APPROX(f,l,x,y,r,a); } __TS_CATCH(f,l) }

// Explicit tolerance:
#   define ETS_ASSERT_APPROX_TOL(x,y,r,a) ___ETS_ASSERT_APPROX(__FILE__,__LINE__,x,y,r,a)
#   define TS_ASSERT_APPROX_TOL(x,y,r,a) ___TS_ASSERT_APPROX(__FILE__,__LINE__,x,y,r,a)

// No explicit accuracy:
#   define ETS_ASSERT_APPROX(x,y) ___ETS_ASSERT_APPROX(__FILE__,__LINE__,x,y, DEF_REL_PRECISION, DEF_ABS_PRECISION)
#   define TS_ASSERT_APPROX(x,y) ___TS_ASSERT_APPROX(__FILE__,__LINE__,x,y, DEF_REL_PRECISION, DEF_ABS_PRECISION)


// TS_ASSERT_VECTOR_APPROX
// Actually supports std::vector<T>, gsl_vector and gsl_matrix
#   define ___ETS_ASSERT_VECTOR_APPROX(f,l,x,y,r,a) ExtraAsserts::doAssertVector( (f), (l), #x, (x), #y, (y), (r), (a) )
#   define ___TS_ASSERT_VECTOR_APPROX(f,l,x,y,r,a) { _TS_TRY { ___ETS_ASSERT_VECTOR_APPROX(f,l,x,y,r,a); } __TS_CATCH(f,l) }

// Explicit tolerance:
#   define ETS_ASSERT_VECTOR_APPROX_TOL(x,y,r,a) ___ETS_ASSERT_VECTOR_APPROX(__FILE__,__LINE__,x,y,r,a)
#   define TS_ASSERT_VECTOR_APPROX_TOL(x,y,r,a) ___TS_ASSERT_VECTOR_APPROX(__FILE__,__LINE__,x,y,r,a)

// No explicit accuracy:
#   define _ETS_ASSERT_VECTOR_APPROX(f,l,x,y) ___ETS_ASSERT_VECTOR_APPROX(f,l,x,y, DEF_REL_PRECISION, DEF_ABS_PRECISION)
#   define _TS_ASSERT_VECTOR_APPROX(f,l,x,y) ___TS_ASSERT_VECTOR_APPROX(f,l,x,y, DEF_REL_PRECISION, DEF_ABS_PRECISION)

#   define ETS_ASSERT_VECTOR_APPROX(x,y) _ETS_ASSERT_VECTOR_APPROX(__FILE__,__LINE__,x,y)
#   define TS_ASSERT_VECTOR_APPROX(x,y) _TS_ASSERT_VECTOR_APPROX(__FILE__,__LINE__,x,y)


/** Functions to check approximate equality and some other IEEE 754 stuff. */
namespace ExtraAsserts {
  /** Calculate the delta to which x and y should be equal in an approximate
   * equality test.
   * 
   * May return NaN, since in any case (x <= NaN) will evaluate false.
   * May not return inf, since (x <= inf) should not pass. */
  double tolerance (double x, double y, double relPrecision = 1e-7, double absPrecision = 1e-7) {
    double tol = relPrecision * max(fabs(x), fabs(y));
    if (tol > (numeric_limits<double>::max)())
      return numeric_limits<double>::quiet_NaN();
    if (tol < absPrecision)
      return absPrecision;
    return tol;
  }
  
  /** Wrapper function. Used because implementation as a macro involves re-evaluation of x and y. */
  void assert_approx (const char *file, unsigned line, double x, double y, double r, double a) {
    const double rel = fabs(x / y) - 1.0;
    if (rel < r && rel > -r) {
      // relative difference small enough; skip delta check
      return;
    }
    // check absolute difference via the provided routine:
    ___ETS_ASSERT_DELTA(file,line,x,y,a,"not approx equal");
  }
  
  /** Basic approximate equality test for doubles, using relative precision.
   *
   * Should work the same as TS_ASSERT_APPROX (x,y) when precision is not
   * specified, hence writing the test like this and not: (fabs(x,y) < d).
   * 
   * Check x and y are approximately equal. Return true if:
   *   x equals y to at least log10(relPrecision) significant figures
   *   or at least log10(absPrecision) decimal places.
   * This should work for small and large values, when one is zero, and when
   * either is infinite or an NaN. */
  bool approx (double x, double y, double relPrecision = 1e-7, double absPrecision = 1e-7) {
    double d = tolerance (x,y, relPrecision, absPrecision);
    return ((y >= x - d) && (y <= x + d));
  }
  
  template<class T>
    void doAssertVector (const char *file, unsigned line,
                         const char *xExpr, vector<T> x,
                         const char *yExpr, vector<T> y,
			 T relP, T absP)
  {
    // This first test must throw on failure, so invalid array indices are not accessed:
    _ETS_ASSERT_EQUALS (file, line, x.size(), y.size());
    for (size_t i = 0; i < x.size(); ++i)
      ___TS_ASSERT_APPROX (file, line, x[i], y[i], relP, absP);
  }
  
# ifdef __GSL_VECTOR_H__
  template<class T>
    void doAssertVector (const char *file, unsigned line,
                         const char *xExpr, gsl_vector *x,
                         const char *yExpr, gsl_vector *y,
			 T relP, T absP)
  {
    _ETS_ASSERT_EQUALS (file, line, x->size, y->size);
    for (size_t i = 0; i < x->size; ++i)
      ___TS_ASSERT_APPROX (file, line, gsl_vector_get(x,i), gsl_vector_get(y,i), relP, absP);
  }
# endif
  
# ifdef __GSL_MATRIX_H__
  template<class T>
    void doAssertVector (const char *file, unsigned line,
                         const char *xExpr, gsl_matrix *x,
                         const char *yExpr, gsl_matrix *y,
			 T relP, T absP)
  {
    _ETS_ASSERT_EQUALS (file, line, x->size1, y->size1);
    _ETS_ASSERT_EQUALS (file, line, x->size2, y->size2);
    for (size_t i = 0; i < x->size1; ++i)
      for (size_t j = 0; j < x->size2; ++j)
	___TS_ASSERT_APPROX (file, line, gsl_matrix_get(x,i,j), gsl_matrix_get(y,i,j), relP, absP);
  }
# endif
  
  class ExtraAssertsSuite : public CxxTest::TestSuite {
  public:
    ExtraAssertsSuite () :
      NaN (numeric_limits<double>::quiet_NaN()),
      Inf (numeric_limits<double>::infinity())
    {}
    
    /// Check some IEEE 754 compliance
    void testIEEE754 () {
      TS_ASSERT_DIFFERS (NaN, NaN);
      TS_ASSERT_DIFFERS (NaN, Inf);
      TS_ASSERT_DIFFERS (NaN, -Inf);
      TS_ASSERT_DIFFERS (Inf, -Inf);
      
      TS_ASSERT_EQUALS (Inf, Inf);
      TS_ASSERT_EQUALS (-Inf, -Inf);
      
      TS_ASSERT_IS_NAN_ (Inf - Inf);
    }
    
    /// Check ExtraAsserts::approx works as expected.
    void testApproxEq () {
      // It might be preferable to use TS_ASSERT_APPROX() and
      // TS_ASSERT_THROWS (ETS_ASSERT_APPROX(), except) to be closer to what's
      // usually used, and get rid of the approx() function.
      
      TS_ASSERT (!approx (NaN, 1))
      TS_ASSERT (!approx (NaN, 0))
      TS_ASSERT (!approx (NaN, NaN))
      TS_ASSERT (!approx (NaN, Inf))
      
      // These 2 tests pass without explicitly checking for inf in tolerance
      // because an NaN is produced in the test:
      TS_ASSERT (!approx (Inf, 1))
      TS_ASSERT (!approx (Inf, 0))
      // However, these 2 fail in this case, because (inf <= inf):
      TS_ASSERT (!approx (1, Inf))
      TS_ASSERT (!approx (0, Inf))
      // This is contrary to IEEE 754 rules, however preferred for unit testing:
      TS_ASSERT (!approx (Inf, Inf))
      TS_ASSERT (!approx (Inf, -Inf))
      
      TS_ASSERT (!approx (0, 1e-6))
      TS_ASSERT ( approx (0, 1e-7))
      
      TS_ASSERT (!approx (1, 0))
      TS_ASSERT (!approx (1, 0.9999998))
      TS_ASSERT ( approx (1, 0.99999995))
      TS_ASSERT (!approx (10000000, 9999998))
      TS_ASSERT ( approx (10000000, 9999999.5))
      // these are considered equal because of absolute precision limitation rather than relative:
      TS_ASSERT ( approx (0.0000001, 0.00000005))
      // this is roughly on the verge of what isn't considered equal:
      TS_ASSERT (!approx (0.0000001, 0.0000003))
      // if we only want to test relative precision:
      TS_ASSERT (!approx (0.0000001, 0.00000009999998,  1e-7, 0))
      TS_ASSERT ( approx (0.0000001, 0.000000099999995, 1e-7, 0))
    }
    
  protected:
    double NaN, Inf;
  };
}

// The convertor script is too simplistic to do this:
#define ExtraAssertsSuite ExtraAsserts::ExtraAssertsSuite

#endif

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

#ifndef Hmod_CheckpointSuite
#define Hmod_CheckpointSuite

#include <cxxtest/TestSuite.h>
#include "util/checkpoint.h"
#include <sstream>
#include <limits>
#include <climits>
#include <iomanip>

using namespace OM::util::checkpoint;

class CheckpointSuite : public CxxTest::TestSuite
{
public:
    CheckpointSuite () {
	test = &testObj;
    }
    
    void setUp () {
	testObj = orig;
    }
    void tearDown () {
    }
    
    void testCheckpointing () {
	stream << setprecision(20);
	ostream& os (stream);
	istream& is (stream);
	
	(*test) & os;
	test->clear ();
	(*test) & is;
	orig.assert_equals (*test);
    }
    
    struct TestObject {
	TestObject () : x(-23263) {}
	virtual ~TestObject () {}
	
	virtual void clear () {
	    x = 0;
	}
	
	virtual void checkpoint (istream& stream) {
	    x & stream;
	}
	virtual void checkpoint (ostream& stream) {
	    x & stream;
	}
	
	template<class S>
	void operator& (S& stream) {
	    checkpoint (stream);
	}
	
	virtual void assert_equals (TestObject& that) {
	    TS_ASSERT_EQUALS (x, that.x);
	}
	
	int x;
    };
    struct DerivedObject : TestObject {
	DerivedObject () :
	  y(3.4422e7),
	  b(true), c(-57),
	  s(-2843), l(LONG_MIN), ll(LONG_MIN),
	  uc(250), us(USHRT_MAX), ui(UINT_MAX), ul(ULONG_MAX), ull(0x1000000000ull),
	  f((numeric_limits<float>::min)()), ld((numeric_limits<long double>::max)()), n(numeric_limits<double>::quiet_NaN())
	{}
	
	virtual void clear () {
	    x = 0;
	    y = 0;
	    c = 0;
	    s = 0;
	    l = 0;
	    ll = 0;
	    uc = 0;
	    us = 0;
	    ui = 0;
	    ul = 0;
	    ull = 0;
	    f = 0;
	    ld = 0;
	    n = 0;
	}
	
	virtual void checkpoint (istream& stream) {
	    x & stream;
	    y & stream;
	    b & stream;
	    c & stream;
	    s & stream;
	    l & stream;
	    ll & stream;
	    uc & stream;
	    us & stream;
	    ui & stream;
	    ul & stream;
	    ull & stream;
	    f & stream;
	    ld & stream;
	    n & stream;
	}
	virtual void checkpoint (ostream& stream) {
	    x & stream;
	    y & stream;
	    b & stream;
	    c & stream;
	    s & stream;
	    l & stream;
	    ll & stream;
	    uc & stream;
	    us & stream;
	    ui & stream;
	    ul & stream;
	    ull & stream;
	    f & stream;
	    ld & stream;
	    n & stream;
	}
	
	virtual void assert_equals (TestObject& p) {
	    DerivedObject& that = dynamic_cast<DerivedObject&> (p);
	    TS_ASSERT_EQUALS (x, that.x);
	    TS_ASSERT_EQUALS (y, that.y);
	    TS_ASSERT_EQUALS (b, that.b);
	    TS_ASSERT_EQUALS (c, that.c);
	    TS_ASSERT_EQUALS (s, that.s);
	    TS_ASSERT_EQUALS (l, that.l);
	    TS_ASSERT_EQUALS (ll, that.ll);
	    TS_ASSERT_EQUALS (uc, that.uc);
	    TS_ASSERT_EQUALS (us, that.us);
	    TS_ASSERT_EQUALS (ui, that.ui);
	    TS_ASSERT_EQUALS (ul, that.ul);
	    TS_ASSERT_EQUALS (ull, that.ull);
	    TS_ASSERT (f != 0.0f);
	    TS_ASSERT_DELTA (f, that.f, (numeric_limits<float>::min)());
	    TS_ASSERT_EQUALS (f, that.f);
	    TS_ASSERT (ld == ld);	// not NaN
	    TS_ASSERT (ld != numeric_limits<double>::infinity());
	    TS_ASSERT_DELTA (ld, that.ld, numeric_limits<long double>::epsilon());
	    TS_ASSERT_EQUALS (ld, that.ld);
	    TS_ASSERT (n != that.n);	// is an NaN
	}
	
	double y;
	bool b;
	signed char c;
	short s;
	long l;
	long long ll;
	unsigned char uc;
	unsigned short us;
	unsigned int ui;
	unsigned long ul;
	unsigned long long ull;
	float f;
	long double ld;
	double n;
    };
    DerivedObject orig, testObj;
    TestObject* test;
    std::stringstream stream;
};

#endif

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

#ifndef OM_util_checkpoint_containers
#define OM_util_checkpoint_containers

#ifndef Hmod_Global
#error "Please include Global.h first."
// otherwise "using ..." declaration in Global.h won't work
#endif

/** Provides some extra functions. See checkpoint.h. */
namespace OM {
namespace util {
namespace checkpoint {
    ///@brief Operator& for smart pointers
    //@{
    template<class T>
    void operator& (unique_ptr<T> & x, ostream& stream) {
        (*x) & stream;
    }

    template<class T>
    void operator& (unique_ptr<T>& x, istream& stream) {
        (*x) & stream;
    }
    //@}
    
    ///@brief Operator& for stl containers
    //@{
    template<class U, class V, class S>
    void operator& (pair<U,V> x, S& stream) {
        x.first & stream;
        x.second & stream;
    }
    
    template<class T>
    void operator& (vector<T>& x, ostream& stream) {
        x.size() & stream;
        for (T& y : x) {
            y & stream;
        }
    }

    template<class T>
    void operator& (vector<T>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        for (T& y : x) {
            y & stream;
        }
    }

    /// Version of above taking an element to initialize each element from.
    template<class T>
    void checkpoint (vector<T>& x, istream& stream, T templateInstance) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l, templateInstance);
        for (T& y : x) {
            y & stream;
        }
    }
    
    template<class T>
    void operator& (list<T> x, ostream& stream) {
        x.size() & stream;
        for (T& y : x) {
            y & stream;
        }
    }
    template<class T>
    void operator& (list<T>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        for (T& y : x) {
            y & stream;
        }
    }

    template<class K, class V>
    void operator& (std::map<K,V> const& m, std::ostream& stream) {
        size_t n = m.size();
        n & stream; // write count
        for (auto const& kv : m) {
            kv.first  & stream;
            kv.second & stream;
        }
    }

    template<class K, class V>
    void operator& (std::map<K,V>& m, std::istream& stream) {
        size_t n;
        n & stream; // read count
        validateListSize(n);
        m.clear();
        for (size_t i = 0; i < n; ++i) {
            K key{stream};
            V val;
            val & stream;
            m.emplace(std::move(key), std::move(val));
        }
    }
    
    void operator& (const set<interventions::ComponentId>&x, ostream& stream);
    void operator& (set<interventions::ComponentId>& x, istream& stream);
    
    void operator& (const map<string,double>& x, ostream& stream);
    void operator& (map<string, double >& x, istream& stream);
    
    void operator& (const map<double,double>& x, ostream& stream);
    void operator& (map<double, double>& x, istream& stream);
    
    void operator& (const map<interventions::ComponentId,SimTime>& x, ostream& stream);
    void operator& (map<interventions::ComponentId,SimTime>& x, istream& stream);
    
    void operator& (const multimap<double,double>& x, ostream& stream);
    void operator& (multimap<double, double>& x, istream& stream);
    //@}
    
} } }   // end of namespaces
#endif

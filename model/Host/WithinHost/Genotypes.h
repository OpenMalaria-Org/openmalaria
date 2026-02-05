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

#ifndef Hmod_WH_Genotypes
#define Hmod_WH_Genotypes

#include "Global.h"
#include "util/random.h"
#include <iostream>
#include <set>

namespace OM { namespace WithinHost {

using util::LocalRng;

/** Represents infection genotypes. */
class Genotypes {
public:
    /// Represent a combination of alleles, each from a different locus
    struct Genotype{
        Genotype( uint32_t allele, double iFreq, double fit, bool hrp2_del ):
            init_freq(iFreq), fitness(fit), hrp2_deficient(hrp2_del)
        {
            alleles.insert(allele);
        }
        Genotype cross( const Genotype& that )const;
        std::set<uint32_t> alleles;      // set of codes of all alleles
        double init_freq;
        double fitness;
        bool hrp2_deficient;
    };
    
    /** Initialise with a single genotype. */
    static void initSingle();
    
    /** Initialise from XML data. Call this before other static methods are
     * used (from PK/PD code). */
    static void init( const scnXml::Scenario& scenario );
    
    /** Switch to whichever mode has been enabled for the main simulation. */
    static void preMainSimInit();
    
    /** Map a locus name and allele name to an allele code.
     * 
     * Note that two alleles from different loci will always have different
     * codes.
     * 
     * Returns max value when no match is found. */
    static uint32_t findAlleleCode( const string& locus, const string& allele );
    
    /** Get a reference to the list of all genotypes. */
    static const vector<Genotype>& getGenotypes();
    
    /** Sample the genotype using the configured approach.
     * 
     * @param genotype_weights When in tracking mode, this vector gives the
     *  weights of each genotype for use in sampling. Total need not be one.
     *  Also, passing a zero-length vector is a signal to use initial
     *  frequencies in sampling. */
    static uint32_t sampleGenotype( LocalRng& rng, std::vector<double>& genotype_weights );
    
    /** Get the number of genotypes. Functions like sampleGenotype use values
     * from 0 to one less than this. */
    inline static size_t N(){ return N_genotypes; }
    
    /** Get the initial frequency of some genotype. */
    static double initialFreq( size_t genotype );
    
    static void staticCheckpoint( std::ostream& stream );
    static void staticCheckpoint( std::istream& stream );
    
private:
    static size_t N_genotypes;
};

}
}
#endif

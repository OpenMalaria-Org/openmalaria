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

#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/vectors.h"
#include "util/CommandLine.h"
#include "schema/scenario.h"

#include <iomanip>

namespace OM {
namespace WithinHost {

namespace GT /*for genotype impl details*/{
// ———  Model constants (after init)  ———
// keys are cumulative probabilities; last entry should equal 1; values are genotype codes
map<double,uint32_t> cum_initial_freqs;

// we give each allele of each loci a unique code
map<string, map<string, uint32_t> > alleleCodes;
uint32_t nextAlleleCode = 0;

// Represent a set of loci: all possible combinations of alleles
// This is actually just machinery to calculate the list of all genotypes
struct LocusSet{
    LocusSet( const scnXml::ParasiteLocus& elt_l ){
        alleles.reserve( elt_l.getAllele().size() );
        double cum_p = 0.0;
        for( const scnXml::ParasiteAllele& elt_a : elt_l.getAllele() ){
            uint32_t alleleCode = (nextAlleleCode++);
            alleleCodes[elt_l.getName()][elt_a.getName()] = alleleCode;
            cum_p += elt_a.getInitialFrequency();
            alleles.push_back( Genotypes::Genotype( alleleCode, elt_a.getInitialFrequency(), elt_a.getFitness(), elt_a.getHrp2_deletion() ) );
        }
        if( cum_p < 0.999 || cum_p > 1.001 )
            throw util::xml_scenario_error("expected sum of initial probabilities of alleles to be 1, but for the " 
                + to_string(alleles.size()) + " alleles under locus " + string(elt_l.getName()) + " this is " + to_string(cum_p));

        alleles[0].init_freq += 1.0 - cum_p;     // account for any small errors by adjusting the first frequency
    }
    void include( LocusSet& that ){
        vector<Genotypes::Genotype> newAlleles;
        newAlleles.reserve( alleles.size() * that.alleles.size() );
        for( size_t i = 0; i < alleles.size(); i += 1 ){
            for( size_t j = 0; j < that.alleles.size(); j += 1 ){
                // note that we generate new codes each time; we just waste the old ones
                newAlleles.push_back( alleles[i].cross(that.alleles[j]) );
            }
        }
        alleles.swap(newAlleles);
    }
    
    vector<Genotypes::Genotype> alleles;
};

vector<Genotypes::Genotype> genotypes;
// ———  Model variables  ———
enum SampleMode{
    SAMPLE_FIRST,      // always choose first genotype (essentially the off switch)
    SAMPLE_INITIAL,    // sample from initial probabilities
    SAMPLE_TRACKING    // sample from tracked success at genotype level (no recombination)
};
// Mode to use now (until switched) and from the start of the intervention period.
SampleMode current_mode = SAMPLE_FIRST, interv_mode = SAMPLE_FIRST;
}
size_t Genotypes::N_genotypes = 1;

Genotypes::Genotype Genotypes::Genotype::cross( const Genotype& that )const{
    Genotype result(*that.alleles.begin(),
                    init_freq*that.init_freq,
                    fitness*that.fitness,
                    hrp2_deficient || that.hrp2_deficient);
    result.alleles.insert( alleles.begin(), alleles.end() );
    result.alleles.insert( that.alleles.begin(), that.alleles.end() );
    return result;
}

void Genotypes::initSingle()
{
    // no specification implies there is a single genotype
    GT::genotypes.assign( 1, Genotypes::Genotype(
        0 /*allele code*/, 1.0/*frequency*/, 1.0/*fitness*/, false /*hrp2 deficiency*/) );
    N_genotypes = 1;
}

// utility function: does a vector contain an element?
template<class T>
bool contains(const vector<T>& vec, const T& x){
    for( auto it = vec.begin(), end = vec.end(); it != end; ++it ){
        if( *it == x ) return true;
    }
    return false;
}

void Genotypes::init( const scnXml::Scenario& scenario ){
    if( scenario.getParasiteGenetics().present() ){
        const scnXml::ParasiteGenetics& genetics =
            scenario.getParasiteGenetics().get();
        
        GT::current_mode = GT::SAMPLE_INITIAL;      // turn on sampling
        if( genetics.getSamplingMode() == "initial" ){
            GT::interv_mode = GT::SAMPLE_INITIAL;
        }else if( genetics.getSamplingMode() == "tracking" ){
            GT::interv_mode = GT::SAMPLE_TRACKING;
            if( !scenario.getEntomology().getVector().present() ){
                throw util::xml_scenario_error( "incompatibility; either use "
                    "entomology/vector (not nonVector) or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
            if( scenario.getEntomology().getMode() != "dynamic" ){
                throw util::xml_scenario_error( "incompatibility; either set "
                    "entomology/mode to \"dynamic\" (not \"forced\") or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
        }else{
            throw util::xml_scenario_error( "parasiteGenetics/samplingMode: expected \"initial\" or \"tracking\"" );
        }
        
        // Build the list of all allele combinations by iterating over loci (plural of locus):
        GT::LocusSet loci( genetics.getLocus()[0] );
        for( size_t i = 1; i < genetics.getLocus().size(); i += 1 ){
            GT::LocusSet newLocus(genetics.getLocus()[i]);
            loci.include( newLocus );
        }
        GT::genotypes.swap( loci.alleles );
        N_genotypes = GT::genotypes.size();
        
        double cum_p = 0.0;
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            cum_p += GT::genotypes[i].init_freq;
            uint32_t genotype_id = static_cast<uint32_t>(i);
            assert( genotype_id == i );
            GT::cum_initial_freqs.insert( make_pair(cum_p, genotype_id) );
        }
        
        // Test cum_p is approx. 1.0 in case the input tree is wrong.
        if (cum_p < 0.999 || cum_p > 1.001)
            throw util::xml_scenario_error ("decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cum_p));

        // last cum_p might be slightless less than 1 due to arithmetic errors; add a failsafe:
        GT::cum_initial_freqs[1.0] = GT::genotypes.size() - 1;
    }else{
        initSingle();
        GT::cum_initial_freqs[1.0] = GT::genotypes.size() - 1;
    }
    
    if( util::CommandLine::option( util::CommandLine::PRINT_GENOTYPES ) ){
        // reorganise GT::alleleCodes so that we can look up codes, not names
        vector<pair<string,string> > allele_codes( GT::cum_initial_freqs.size() );
        for( auto i = GT::alleleCodes.begin(), iend = GT::alleleCodes.end(); i != iend; ++i ) {
            const string locus = i->first;
            for( auto j = i->second.begin(),
                jend = i->second.end(); j != jend; ++j )
            {
                uint32_t code = j->second;
                const string allele = j->first;
                assert( code < allele_codes.size() );
                allele_codes[code] = make_pair( locus, allele );
            }
        }
        
        // determine our columns
        map<string,uint32_t> longest;      // longest name in column; key is locus
        for( auto i = GT::genotypes.begin(),
            iend = GT::genotypes.end(); i != iend; ++i )
        {
            if( longest.size() == 0 ){
                for( auto j = i->alleles.begin();
                    j != i->alleles.end(); ++j ){
                    const string& locus = allele_codes[*j].first;

                    longest[locus] = locus.length();  // locus name is included in column
                }
            }else assert( longest.size() == i->alleles.size() );
            
            for( auto j = i->alleles.begin();
                j != i->alleles.end(); ++j ){
                auto it = longest.find(allele_codes[*j].first);
                assert( it != longest.end() );
                uint32_t len_allele = allele_codes[*j].second.length();
                if( len_allele > it->second ) it->second = len_allele;
            }
        }
        
        // find original loci order
        vector<string> loci;
        loci.reserve( longest.size() );
        for( size_t i = 0; loci.size() < longest.size(); ++i ){
            assert( i < allele_codes.size() );
            const string& locus = allele_codes[i].first;
            if( !contains(loci, locus) ){
                loci.push_back( locus );
            }
        }
        
//         string fm = "|%1$-12d|%|14t|%2$-12d|";
        cout << endl;
        stringstream fmt;
        fmt << "|%8d|";
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            fmt << "%" << longest[*it] << "s|";
        fmt << "%9.3f|%7.3f|";
        
        // Table header:
        cout << "|" << std::setw(8) << "Genotype"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            cout << "|" << std::setw(longest[*it]) << *it;
        cout << "|" << std::setw(9) << "init freq" << "|" << std::setw(7) << "fitness" << "|" << endl;
        
        cout << "|" << std::setw(8) << "--------"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
        {
            if(longest[*it] > 0)
                cout << "|" << "---------";
            else
                cout << "|";
        }

        cout << "|" << std::setw(9) << "---------" << "|" << std::setw(7) << "-------" << "|" << endl;

        cout << std::fixed << setprecision(3);
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            const Genotype& genotype = GT::genotypes[i];
            map<string,string> locus_allele;    // to find allele for each locus
            for( auto a = genotype.alleles.begin(); a != genotype.alleles.end(); ++a ){
                assert( *a < allele_codes.size() );
                locus_allele[allele_codes[*a].first] = allele_codes[*a].second;
            }

            cout << "|" << std::setw(8) << (i*100000);
            for( auto it = loci.begin(); it !=loci.end(); ++it ){
                auto la = locus_allele.find( *it );
                assert( la != locus_allele.end() );
                cout << "|" << std::setw(longest[*it]) << la->second << "|";

            }
            cout << std::setw(9) << genotype.init_freq << "|" << std::setw(7) << genotype.fitness << "|" << endl;
        }
    }
}

void Genotypes::preMainSimInit(){
    GT::current_mode = GT::interv_mode;
}

uint32_t Genotypes::findAlleleCode(const string& locus, const string& allele){
    auto it = GT::alleleCodes.find(locus);
    if( it == GT::alleleCodes.end() ) return numeric_limits<uint32_t>::max();
    auto it2 = it->second.find( allele );
    if( it2 == it->second.end() ) return numeric_limits<uint32_t>::max();
    return it2->second;
}

const vector< Genotypes::Genotype >& Genotypes::getGenotypes(){
    return GT::genotypes;
}

uint32_t Genotypes::sampleGenotype( LocalRng& rng, vector<double>& genotype_weights ){
    if( GT::current_mode == GT::SAMPLE_FIRST ){
        return 0;       // always the first genotype code
    }else if( GT::current_mode == GT::SAMPLE_INITIAL
            || genotype_weights.size() == 0 )
    {
        double sample = rng.uniform_01();
        auto it = GT::cum_initial_freqs.upper_bound( sample );
        assert( it != GT::cum_initial_freqs.end() );
        return it->second;
    }else{
        assert( GT::current_mode == GT::SAMPLE_TRACKING );
        assert( genotype_weights.size() == N_genotypes );
        double weight_sum = util::vectors::sum( genotype_weights );
        assert( weight_sum >= 0.0 && weight_sum < 1e5 );        // possible loss of precision or other error
        double sample = rng.uniform_01() * weight_sum;
        double cum = 0.0;
        for( size_t g = 0; g < N_genotypes; ++g ){
            cum += genotype_weights[g];
            if( sample < cum ) return g;
        }
        return 0;       // just to be safe (could happen if weight_sum == 0.0)
    }
}

double Genotypes::initialFreq( size_t genotype ){
    if( GT::genotypes.size() == 0 ){
        assert( genotype == 0 );
        assert( N_genotypes == 1 );
        return 1.0;     // only 1 genotype, thus initial frequency is 100%
    }else{
        return GT::genotypes[genotype].init_freq;
    }
}


// ———  checkpointing  ———

void Genotypes::staticCheckpoint( ostream& stream ){
    int t = GT::current_mode;
    t & stream;
}
void Genotypes::staticCheckpoint( istream& stream ){
    int t;
    t & stream;
    GT::current_mode = static_cast<GT::SampleMode>(t);
    assert( GT::current_mode == GT::SAMPLE_FIRST ||
        GT::current_mode == GT::SAMPLE_INITIAL ||
        GT::current_mode == GT::SAMPLE_TRACKING );
}

}
}

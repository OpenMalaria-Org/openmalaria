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

#include "Host/WithinHost/WHVivax.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/Infection.h"
#include "Host/Human.h"
#include "mon/reporting.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/checkpoint_containers.h"
#include "util/UnitParse.h"
#include "util/CommandLine.h"
#include "util/sampler.h"
#include <schema/scenario.h>
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

struct HypnozoiteReleaseDistribution {
    HypnozoiteReleaseDistribution(std::unique_ptr<util::LognormalSampler> sampler, double latentRelapse) :
        sampler(std::move(sampler)),
        latentRelapse(latentRelapse)
        {}
    
    /// Sample the time until next release
    SimTime sampleReleaseDelay(LocalRng& rng) const {
        double liverStageMaximumDays = 16.0*30.0; // maximum of about 16 months in liver stage 
        double delay = numeric_limits<double>::quiet_NaN();       // in days
        int count = 0;
        int maxcount = 1e3;
        
        do{
            delay = sampler->sample(rng);
            count += 1;
            
            if( count >= maxcount ){
                throw util::xml_scenario_error( "<vivax><hypnozoiteRelease>  [random delay calculation causes probably an indefinite loop]:\n The hypnozoite release distribution seems off, sigma of secondRelease could be too high. We except the hypnozoite to reside a maximum of 16 months in the liver stage. Sigma choose well, dear padawan." );
            }
        } while( delay > liverStageMaximumDays || delay < 0.0  );
        
        assert( delay >= 0 && delay < liverStageMaximumDays );
        return sim::roundToTSFromDays( delay + latentRelapse );
    }
    
private:
    std::unique_ptr<util::LognormalSampler> sampler;
    double latentRelapse;   // days
};

std::unique_ptr<HypnozoiteReleaseDistribution> createHypnozoiteReleaseDistribution(const scnXml::HypnozoiteReleaseDistribution& elt)
{
    return std::make_unique<HypnozoiteReleaseDistribution>(createSampler(elt),  elt.getLatentRelapse());
}

// ———  parameters  ———

// Set from the parameters block, or a default value is assumed.
SimTime latentP = sim::never();       // attribute on <parameters> element.

// Set from <vivax .../> element:
double probBloodStageInfectiousToMosq = numeric_limits<double>::signaling_NaN();
int maxNumberHypnozoites = -1;
double baseNumberHypnozoites = numeric_limits<double>::signaling_NaN();
std::unique_ptr<HypnozoiteReleaseDistribution> latentRelapse1st, latentRelapse2nd;
double pSecondRelease = numeric_limits<double>::signaling_NaN();
SimTime bloodStageProtectionLatency = sim::never();
WeibullSampler bloodStageLength;    // units: days
double pPrimaryA = numeric_limits<double>::signaling_NaN(),
    pPrimaryB = numeric_limits<double>::signaling_NaN();
double pRelapseOneA = numeric_limits<double>::signaling_NaN(),
    pRelapseOneB = numeric_limits<double>::signaling_NaN();
double pRelapseTwoA = numeric_limits<double>::signaling_NaN(),
    pRelapseTwoB = numeric_limits<double>::signaling_NaN();
double pEventIsSevere = numeric_limits<double>::signaling_NaN();
std::string vivaxClinOption = " ";

// Set from healthSystem element:
bool ignoreNoPQ = false;
double pHetNoPQ = numeric_limits<double>::signaling_NaN();
double pReceivePQ = numeric_limits<double>::signaling_NaN();
double effectivenessPQ = numeric_limits<double>::signaling_NaN();

#ifdef WHVivaxSamples
WHVivax *sampleHost = 0;
VivaxBrood *sampleBrood = 0;
#endif


// ———  individual models  ———


// number of hypnozoites per brood:
map<double,int> nHypnozoitesProbMap;
void initNHypnozoites(){
    assert(baseNumberHypnozoites <= 1 && baseNumberHypnozoites >= 0);
    double total = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n )
        total += pow( baseNumberHypnozoites, n );
    
    double cumP = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n ){
        cumP += pow( baseNumberHypnozoites, n ) / total;
        // pair n with the cumulative probability of sampling n:
        nHypnozoitesProbMap[cumP] = n;
    }
}


int sampleNHypnozoites(LocalRng& rng){
    double x = rng.uniform_01();
    // upper_bound finds the first key (cumulative probability) greater than x:
    auto it = nHypnozoitesProbMap.upper_bound( x );
    assert( it != nHypnozoitesProbMap.end() );  // i.e. that we did find a real key
    return it->second;  // corresponding n
}


// time to hypnozoite release after initial release:
SimTime sampleReleaseDelay(LocalRng& rng){
    bool isFirstRelease = true;
    if( pSecondRelease == 1.0 ||
        (pSecondRelease > 0.0 && rng.bernoulli(pSecondRelease)) ){
        isFirstRelease = false;
    }
    
    if (isFirstRelease){
        return latentRelapse1st->sampleReleaseDelay(rng);
    } else {
        return latentRelapse2nd->sampleReleaseDelay(rng);
    }
}


// ———  per-brood code  ———

VivaxBrood::VivaxBrood( LocalRng& rng, int origin, WHVivax *host ) :
        primaryHasStarted( false ),
        relapseHasStarted( false ),
		relapsebHasStarted( false),
        hadEvent( false ),
        hadRelapse( false ),
		origin(origin)
 {
    set<SimTime> releases;     // used to initialise releaseDates; a set is better to use now but a vector later
    
    // primary blood stage plus hypnozoites (relapses)
    releases.insert( sim::nowOrTs0() + latentP );
    int numberHypnozoites = sampleNHypnozoites(rng);
    
    for( int i = 0; i < numberHypnozoites; ){
        SimTime randomReleaseDelay = sampleReleaseDelay(rng);
        SimTime timeToRelease = sim::nowOrTs0() + latentP + randomReleaseDelay;
        bool inserted = releases.insert( timeToRelease ).second;
        if( inserted ) ++i;     // successful
        // else: sample clash with an existing release date, so resample
    }

    // Copy times to the vector, backwards (smallest last):
    releaseDates.insert( releaseDates.end(), releases.rbegin(), releases.rend() );
    
#ifdef WHVivaxSamples
    if( sampleHost == host && sampleBrood == 0 ){
        sampleBrood = this;
        cout << "New sample brood";
        for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it )
            cout << '\t' << *it;
        cout << endl;
    }
#endif
}
VivaxBrood::~VivaxBrood(){
#ifdef WHVivaxSamples
    if( sampleBrood == this ){
        sampleBrood = 0;
        cout << "Brood terminated" << endl;
    }
#endif
}

void VivaxBrood::checkpoint( ostream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    relapseHasStarted & stream;
	relapsebHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
	origin & stream;
}
VivaxBrood::VivaxBrood( istream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    relapseHasStarted & stream;
	relapsebHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
	origin & stream;
}


VivaxBrood::UpdResult VivaxBrood::update(LocalRng& rng){
    if( bloodStageClearDate == sim::ts0() ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage infections for a short time.
    }
    
    UpdResult result;
    while( releaseDates.size() > 0 && releaseDates.back() == sim::ts0() ){
        releaseDates.pop_back();
        
#ifdef WHVivaxSamples
        if( sampleBrood == this ){
            cout << "Time\t" << sim::ts0();
            for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it )
                cout << '\t' << *it;
            cout << endl;
        }
#endif
        
        // an existing or recently terminated blood stage from the same brood
        // protects against a newly released hypnozoite for relapse classification B 
        if( (vivaxClinOption=="B1j" || vivaxClinOption=="B2j") && (bloodStageClearDate + bloodStageProtectionLatency >= sim::ts0()) ) continue;
        		
        if( !relapsebHasStarted && relapseHasStarted ){
            relapsebHasStarted = true;
            result.newRelapsebBS = true;
        }		
        if( !relapseHasStarted && primaryHasStarted ){
            relapseHasStarted = true;
            result.newRelapseBS = true;
        }
        if( !primaryHasStarted ){
            primaryHasStarted = true;
            result.newPrimaryBS = true;
        }
        result.newBS = true;
        
        double lengthDays = bloodStageLength.sample(rng);
        bloodStageClearDate = sim::ts0() + sim::roundToTSFromDays( lengthDays );
        // Assume gametocytes emerge at the same time (they mature quickly and
        // we have little data, thus assume coincedence of start)
    }
    
    result.isFinished = (!isPatent()) && (releaseDates.size() == 0);
    return result;
}

void VivaxBrood::treatmentBS(){
    // Blood stage treatment: clear both asexual and sexual parasites from the
    // blood. NOTE: we assume infections removed via treatment do not leave
    // protective immunity since the patient was unable to self-clear.
    bloodStageClearDate = sim::never();
}

void VivaxBrood::treatmentLS(){
    releaseDates.clear();       // 100% clearance
    
    /* partial clearance code, in case of need:
    vector<SimTime> survivingZoites;
    survivingZoites.reserve( releaseDates.size() );   // maximum size we need
    for( auto it = releaseDates.begin(); it != releaseDates.end(); ++it ){
        if( !rng.bernoulli( pClearEachHypnozoite ) ){
            survivingZoites.push_back( *it );    // copy            
        }
    }
    swap( releaseDates, survivingZoites );
    */
}


// ———  per-host code  ———

WHVivax::WHVivax( LocalRng& rng, double comorbidityFactor ) :
    cumPrimInf(0),
	treatExpiryLiver(0), treatExpiryBlood(0),
    pEvent( numeric_limits<double>::quiet_NaN() ),
    pFirstRelapseEvent( numeric_limits<double>::quiet_NaN() ),
    pSevere( 0.0 )
{
    if( comorbidityFactor != 1.0 ){
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
    }
#ifdef WHVivaxSamples
    if( sampleHost == 0 ){
        sampleHost = this;
        cout << "New host" << endl;
    }
#endif
    noPQ = ( pHetNoPQ > 0.0 && rng.bernoulli(pHetNoPQ) );
}

WHVivax::~WHVivax(){
#ifdef WHVivaxSamples
    if( this == sampleHost ){
        sampleHost = 0;
        cout << "Host terminates" << endl;
    }
#endif
}



double WHVivax::probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l)const{
    assert( WithinHost::Genotypes::N() == 1 );
    for(auto inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probBloodStageInfectiousToMosq;
        }
    }
    return 0;   // no gametocytes
}


bool WHVivax::summarize(Host::Human& human) const{
    if( infections.size() == 0 ) return false;  // no infections: not patent, nothing to report
    mon::reportStatMHI( mon::MHR_INFECTED_HOSTS, human, 1 );
    bool patentHost = false;
    // (patent) infections are reported by genotype, even though we don't have
    // genotype in this model
    mon::reportStatMHGI( mon::MHR_INFECTIONS, human, 0, infections.size() );
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        if (inf->isPatent()){
            mon::reportStatMHGI( mon::MHR_PATENT_INFECTIONS, human, 0, 1 );
            patentHost = true;
        }
    }
    if( patentHost ) mon::reportStatMHI( mon::MHR_PATENT_HOSTS, human, 1 );
    return patentHost;
}

void WHVivax::importInfection(LocalRng& rng){
    // this means one new liver stage infection, which can result in multiple blood stages
    infections.push_back( VivaxBrood( rng, WithinHost::InfectionOrigin::Imported, this ) );
}

void WHVivax::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    pSevere = 0.0;
    
    // create new infections, letting the constructor do the initialisation work:
    for( int i = 0; i < nNewInfs_i; ++i )
        infections.push_back( VivaxBrood( rng, InfectionOrigin::Introduced, this ) );

    for( int i = 0; i < nNewInfs_l; ++i )
        infections.push_back( VivaxBrood( rng, InfectionOrigin::Indigenous, this ) );

	
    
    // update infections
    // NOTE: currently no BSV model
    morbidity = Pathogenesis::NONE;
    uint32_t oldCumInf = cumPrimInf;
    bool treatmentLiver = treatExpiryLiver > sim::ts0();
    bool treatmentBlood = treatExpiryBlood > sim::ts0();
    double matImmClin = 1 - (0.90 * exp(-2.53*ageInYears));
    //double matImmClin = 1 - (1 * exp(-(0.639*8*ageInYears)/2.53));
    auto inf = infections.begin();
    while( inf != infections.end() ){
        if( treatmentLiver ) inf->treatmentLS();
        if( treatmentBlood ) inf->treatmentBS();        // clearance due to treatment; no protection against reemergence
        VivaxBrood::UpdResult result = inf->update(rng);
        if( result.newPrimaryBS ) cumPrimInf += 1;
        
        if( result.newBS ){
            // Sample for each new blood stage infection: the chance of some clinical event.
            // model variant: no illness from relapses possible unless there was illness from the primary infection
            
            bool clinicalEvent = false;
            if( result.newPrimaryBS ){
                // Blood stage is primary. oldCumInf wasn't updated yet.
                //double pPrimaryInfEvent = matImmClin * pPrimaryA * pPrimaryB / (pPrimaryB+oldCumInf);
				double pPrimaryInfEvent = matImmClin * pPrimaryA * exp(-pPrimaryB * oldCumInf);
                clinicalEvent = rng.bernoulli( pPrimaryInfEvent );
                inf->setHadEvent( clinicalEvent );
					
			} else if ( result.newRelapseBS ){
                    // Blood stage is a relapse. oldCumInf wasn't updated yet. Subtract 1 from oldCumInf not
                    // to count the current brood in the number of cumulative primary infections 
					
                    double pFirstRelapseEvent = matImmClin * pRelapseOneA * exp(-pRelapseOneB * (oldCumInf-1));
                    clinicalEvent = rng.bernoulli( pFirstRelapseEvent );     
                    inf->setHadRelapse( clinicalEvent );
      
            } else if ( result.newRelapsebBS ){
                             
				    if (vivaxClinOption=="A1j" || vivaxClinOption=="A2j"){			 
                       double pFirstRelapseEvent = matImmClin * pRelapseOneA * exp(-pRelapseOneB * (oldCumInf-1));
                       clinicalEvent = rng.bernoulli( pFirstRelapseEvent );     
                       inf->setHadRelapse( clinicalEvent );
				    }
				    if (vivaxClinOption=="B1j" || vivaxClinOption=="B2j"){			 
                       double pSecondRelapseEvent = matImmClin * pRelapseTwoA * exp(-pRelapseTwoB * (oldCumInf-1));
                       clinicalEvent = rng.bernoulli( pSecondRelapseEvent );
					   inf->setHadRelapse( clinicalEvent );
			        }
				                    
            } else {
                     double pSecondRelapseEvent = matImmClin * pRelapseTwoA * exp(-pRelapseTwoB * (oldCumInf-1));
                     clinicalEvent = rng.bernoulli( pSecondRelapseEvent );
                    }
             

            if( clinicalEvent ){
                pSevere = pSevere + (1.0 - pSevere) * pEventIsSevere;
                if( rng.bernoulli( pEventIsSevere ) )
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_SEVERE );
                else
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_MALARIA );
            }
        }
        
        if( result.isFinished ) inf = infections.erase( inf );
        else ++inf;
    }
    
    
    
    //NOTE: currently we don't model co-infection or indirect deaths
    if( morbidity == Pathogenesis::NONE ){
        morbidity = Pathogenesis::PathogenesisModel::sampleNMF( rng, ageInYears );
    }
}

bool WHVivax::diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const{
    //TODO(monitoring): this shouldn't ignore the diagnostic (especially since
    // it should always return true if diagnostic.density=0)
    for(auto inf = infections.begin(); inf != infections.end(); ++inf) {
        if (inf->isPatent())
            return true;        // at least one patent infection
    }
    return false;
}

Pathogenesis::StatePair WHVivax::determineMorbidity( Host::Human& human, double ageYears, bool ){
    mon::reportStatMHF( mon::MHF_EXPECTED_SEVERE, human, pSevere );
    Pathogenesis::StatePair result;     // no indirect mortality in the vivax model
    result.state = morbidity;
    return result;
}

void WHVivax::clearImmunity(){
    throw TRACED_EXCEPTION_DEFAULT( "vivax model does not include immune suppression" );
}

void WHVivax::treatment( Host::Human& human, TreatmentId treatId ){
    const Treatments& treat = Treatments::select( treatId );
    treatSimple( human, treat.liverEffect(), treat.bloodEffect() );
    
    // triggered intervention deployments:
    treat.deploy( human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
}
void WHVivax::optionalPqTreatment( Host::Human& human ){
    // PQ clears liver stages. We don't worry about the effect of PQ on
    // gametocytes, because these are always cleared by blood-stage drugs with
    // Vivax, and PQ is not given without BS drugs. NOTE: this ignores drug failure.
    if (pReceivePQ > 0.0 && (ignoreNoPQ || !noPQ) && human.rng.bernoulli(pReceivePQ)){
        if( human.rng.bernoulli(effectivenessPQ) ){
            for( auto it = infections.begin(); it != infections.end(); ++it ){
                it->treatmentLS();
            }
        }
        mon::reportEventMHI( mon::MHT_LS_TREATMENTS, human, 1 );
    }
}
bool WHVivax::treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    //TODO: this should be implemented properly (allowing effects on next
    // update instead of now)
    
    // liver-stage treatment is only via "LiverStageDrug" option, if at all
    if( timeLiver != sim::zero() ){
        if( pReceivePQ > 0.0 ){
            // This is only really disallowed to prevent simultaneous treatment through both methods
            throw util::xml_scenario_error("simple treatment for vivax liver "
            "stages is incompatible with case-management pUseUncomplicated "
            "(liverStageDrug) option; it is suggested to use the former over the latter");
        }
        if( (ignoreNoPQ || !noPQ) && (effectivenessPQ == 1.0 || human.rng.bernoulli(effectivenessPQ)) ){
            if( timeLiver >= sim::zero() ){
                treatExpiryLiver = max( int(treatExpiryLiver), sim::nowOrTs1() + timeLiver );
            }else{
                for( auto it = infections.begin(); it != infections.end(); ++it ){
                    it->treatmentLS();
                }
            }
        }
        mon::reportEventMHI( mon::MHT_LS_TREATMENTS, human, 1 );
    }
    
    // there probably will be blood-stage treatment
    if( timeBlood != sim::zero() ){
        if( timeBlood < sim::zero() ){
            // legacy mode: retroactive clearance
            for( auto it = infections.begin(); it != infections.end(); ++it ){
                it->treatmentBS();
            }
        }else{
            treatExpiryBlood = max( int(treatExpiryBlood), sim::nowOrTs1() + timeBlood );
        }
        return true;    // blood stage treatment
    }
    return false;    // no blood stage treatment
}


// ———  boring stuff: checkpointing and set-up  ———

void WHVivax::checkpoint(istream& stream){
    WHInterface::checkpoint(stream);
    size_t len;
    len & stream;
    for( size_t i = 0; i < len; ++i ){
        infections.push_back( VivaxBrood( stream ) );
    }
    noPQ & stream;
    int morbidity_i;
    morbidity_i & stream;
    morbidity = static_cast<Pathogenesis::State>( morbidity_i );
    cumPrimInf & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}
void WHVivax::checkpoint(ostream& stream){
    WHInterface::checkpoint(stream);
    infections.size() & stream;
    for( auto it = infections.begin(); it != infections.end(); ++it ){
        it->checkpoint( stream );
    }
    noPQ & stream;
    static_cast<int>( morbidity ) & stream;
    cumPrimInf & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}

char const*const not_impl = "feature not available in Vivax model";
void WHVivax::treatPkPd(size_t schedule, size_t dosages, double age, double delay_d){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getTotalDensity() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_h() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_Y() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

void WHVivax::init( const OM::Parameters& parameters, const scnXml::Model& model ){

    // Default value to use if no <parameters> element is present in XML.
    latentP = UnitParse::readShortDuration("15d", UnitParse::DAYS);
    try
    {
        if (model.getParameters().present())
        {
            latentP = UnitParse::readShortDuration(
                model.getParameters().get().getLatentp(), UnitParse::NONE);
        }
    }
    catch( const util::format_error& e )
    {
            throw util::xml_scenario_error( string("model/parameters/latentP: ").append(e.message()) );
    }

    if( !model.getVivax().present() )
        throw util::xml_scenario_error( "no vivax model description in scenario XML" );
    const scnXml::Vivax& elt = model.getVivax().get();
    probBloodStageInfectiousToMosq = elt.getProbBloodStageInfectiousToMosq().getValue();
    maxNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getMax();
    baseNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getBase();
    const scnXml::HypnozoiteRelease& hr = elt.getHypnozoiteRelease();
    latentRelapse1st = createHypnozoiteReleaseDistribution(hr.getFirstReleaseDays());
    if(hr.getSecondReleaseDays().present()){
        latentRelapse2nd = createHypnozoiteReleaseDistribution(hr.getSecondReleaseDays().get());
        pSecondRelease = hr.getPSecondRelease();
        assert( pSecondRelease >= 0 && pSecondRelease <= 1 );
    } // else pSecondRelease is NaN and other values don't get used
    bloodStageProtectionLatency = sim::roundToTSFromDays( elt.getBloodStageProtectionLatency().getValue() );
    bloodStageLength.setParams(elt.getBloodStageLengthDays());
    
    const scnXml::ClinicalEvents& ce = elt.getClinicalEvents();
    pPrimaryA = ce.getPPrimaryInfection().getA();
    pPrimaryB = ce.getPPrimaryInfection().getB();
    pRelapseOneA = ce.getPRelapseOne().getA();
    pRelapseOneB = ce.getPRelapseOne().getB();
    pRelapseTwoA = ce.getPRelapseTwoPlus().getA();
    pRelapseTwoB = ce.getPRelapseTwoPlus().getB();
    pEventIsSevere = ce.getPEventIsSevere().getValue();
    vivaxClinOption = ce.getVivaxClinOption();
    
	// Define accepted values for vivaxClinOption
    const std::set<std::string> acceptedOptions = {"A1j", "A2j", "B1j", "B2j"};

   // Check if vivaxClinOption is valid
    if (acceptedOptions.find(vivaxClinOption) == acceptedOptions.end()) {
        throw util::xml_scenario_error("Invalid vivaxClinOption: " + vivaxClinOption +
                                     ". Accepted values are: A1j, A2j, B1j, B2j.");
    }
		
    initNHypnozoites();
    Pathogenesis::PathogenesisModel::init( parameters, model.getClinical(), true );
}
void WHVivax::setHSParameters(const scnXml::LiverStageDrug* elt){
    double oldPHetNoPQ = pHetNoPQ;
    if( elt == 0 ){
        ignoreNoPQ = false;
        pHetNoPQ = 0.0;
        pReceivePQ = 0.0;
        effectivenessPQ = 1.0;  // sensible default: does not affect simple liver stage treatment option
    }else{
        ignoreNoPQ = elt->getIgnoreCannotReceive().present() ?
            elt->getIgnoreCannotReceive().get().getValue() : false;
        pHetNoPQ = elt->getPHumanCannotReceive().getValue();
        pReceivePQ = elt->getPUseUncomplicated().present() ?
            elt->getPUseUncomplicated().get().getValue() :  0.0;
        if( pReceivePQ > 0.0 ){
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: pUseUncomplicated is deprecated; it is "
                    "suggested to use the liver stage simple treatment option instead." << endl;
            }
        }
        effectivenessPQ = elt->getEffectivenessOnUse().getValue();
    }
    
    if( !std::isnan(oldPHetNoPQ) && oldPHetNoPQ != pHetNoPQ ){
        throw util::xml_scenario_error( "changeHS cannot change pHumanCannotReceive value" );
    }
}

}
}

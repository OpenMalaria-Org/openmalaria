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

#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/Drug/LSTMDrugOneComp.h"
#include "PkPd/LSTMModel.h"
#include "PkPd/LSTMMedicate.h"
#include "PkPd/LSTMTreatments.h"
#include "mon/reporting.h"
#include "util/checkpoint_containers.h"
#include "util/errors.h"

#include "schema/scenario.h"

#include <cassert>

namespace OM { namespace PkPd {

// ———  static  ———

void LSTMModel::init( const scnXml::Scenario& scenario ){
    if (scenario.getPharmacology().present()) {
        LSTMDrugType::init(scenario.getPharmacology().get().getDrugs());
        LSTMTreatments::init(scenario.getPharmacology().get().getTreatments());
    }
}

// ———  non-static set up / tear down functions  ———

void LSTMModel::checkpoint (istream& stream) {
    // We need an RNG to call createInstance, but override the values later
    // anyway. Hacky solution: create a local RNG.
    LocalRng rng(0, 0);
    
    size_t numDrugs;	// type must be same as m_drugs.size()
    numDrugs & stream;
    util::checkpoint::validateListSize (numDrugs);
    for(size_t i=0; i<numDrugs; ++i) {
        size_t index;
        index & stream;
        m_drugs.push_back( LSTMDrugType::createInstance(rng, index) );
        m_drugs.back() & stream;
    }
    medicateQueue & stream;
}

void LSTMModel::checkpoint (ostream& stream) {
    m_drugs.size() & stream;
    for( auto& drug : m_drugs ){
        drug->getIndex() & stream;
        drug & stream;
    }
    medicateQueue & stream;
}


// ———  non-static simulation time functions  ———

void LSTMModel::prescribe(size_t schedule, size_t dosage, double age, double body_mass, double delay_d){
    DosageTable& table = dosages[dosage];
    double key = table.useMass ? body_mass : age;
    double doseMult = table.getMultiplier( key );
    for( MedicateData& medicateData : schedules[schedule].medications ){
        MedicateData data = medicateData.multiplied(doseMult);
        data.time += delay_d;
        medicateQueue.push_back( std::move(data) );
    }
}

void LSTMModel::medicate(LocalRng& rng){
    if( medicateQueue.empty() ) return;
    
    // Process pending medications (in interal queue) and apply/update:
    auto iter = medicateQueue.begin();
    while( iter != medicateQueue.end() ){
        if( iter->time < 1.0 ){   // Medicate medications to be prescribed starting at the next time-step
            // This function could be inlined, except for uses in testing:
            medicateDrug (rng, iter->drug, iter->qty, iter->time);
            iter = medicateQueue.erase (iter);
        }else{  // and decrement treatment seeking delay for the rest
            iter->time -= 1.0;
            ++iter;
        }
    }
}

void LSTMModel::medicateDrug(LocalRng& rng, size_t typeIndex, double qty, double time) {
    //TODO: might be a little faster if m_drugs was pre-allocated with a slot for each drug type, using a null pointer
    for( auto& drug : m_drugs ){
        if (drug->getIndex() == typeIndex){
            drug->medicate (time, qty);
            return;
        }
    }
    // No match, so insert one:
    m_drugs.push_back( LSTMDrugType::createInstance(rng, typeIndex) );
    (*m_drugs.back()).medicate (time, qty);
}

double LSTMModel::getDrugConc (size_t drug_index) const{
    double c = 0.0;
    double d = 0.0;
    for( auto& drug : m_drugs ){
        d = drug->getConcentration(drug_index);
        if(d < 0.0){
            // FIXME: @tph-thuering, 2015-08-11:
            // This is trying to circumvent an issue we have with negative metabolite concentration values in the implementation of the conversion model. see #53
            cerr << "ERROR! concentration is lower than zero, setting it to zero for now: " << d << " drug_index: " << drug_index << endl;
            d = 0.0;
        }
        c += d;
    }
    assert(c >= 0.0);
    return c;
}

double LSTMModel::getDrugFactor (LocalRng& rng, WithinHost::CommonInfection *inf, double body_mass) const{
    double factor = 1.0; //no effect
    
    for( auto drug = m_drugs.begin(), end = m_drugs.end(); drug != end; ++drug ){
        double drugFactor = (*drug)->calculateDrugFactor(rng, inf, body_mass);
        factor *= drugFactor;
    }
    return factor;
}

void LSTMModel::decayDrugs (double body_mass) {
    // Update concentrations for each drug.
    // TODO: previously we removed drugs with negligible concentration here. What now, just set concentration to 0?
    for( auto& drug : m_drugs ){
        drug->updateConcentration(body_mass);
    }
}

void LSTMModel::summarize(const Host::Human& human) const{
    const vector<size_t> &drugsInUse( LSTMDrugType::getDrugsInUse() );
    for( size_t index : drugsInUse ){
        for( auto& drug : m_drugs ){
            double conc = drug->getConcentration(index);
            if( conc > 0.0 ){
                mon::reportStatMHPI( mon::MHR_HOSTS_POS_DRUG_CONC, human, index, 1 );
                mon::reportStatMHPF( mon::MHF_LOG_DRUG_CONC, human, index, log(conc) );
            }
        }
    }
}

} }

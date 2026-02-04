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

#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/util.h"

#include <sstream>
#include <iostream>

namespace OM { namespace util {
    std::bitset<NUM_OPTIONS> ModelOptions::options;
    
    // Utility: converts option strings to codes and back
    class OptionCodeMap {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> codeMap;
        set<string> ignoreOptions;
	
    public:
	OptionCodeMap () {
// 	    codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	    codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	    // codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	    codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	    codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	    codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	    codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	    // 	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
	    codeMap["DUMMY_WITHIN_HOST_MODEL"] = DUMMY_WITHIN_HOST_MODEL;
	    codeMap["PREDETERMINED_EPISODES"] = PREDETERMINED_EPISODES;
	    codeMap["NON_MALARIA_FEVERS"] = NON_MALARIA_FEVERS;
        ignoreOptions.insert("INCLUDES_PK_PD");
	    codeMap["CLINICAL_EVENT_SCHEDULER"] = CLINICAL_EVENT_SCHEDULER;
	    codeMap["MUELLER_PRESENTATION_MODEL"] = MUELLER_PRESENTATION_MODEL;
	    codeMap["TRANS_HET"] = TRANS_HET;
	    codeMap["COMORB_HET"] = COMORB_HET;
	    codeMap["TREAT_HET"] = TREAT_HET;
	    codeMap["COMORB_TRANS_HET"] = COMORB_TRANS_HET;
	    codeMap["TRANS_TREAT_HET"] = TRANS_TREAT_HET;
	    codeMap["COMORB_TREAT_HET"] = COMORB_TREAT_HET;
	    codeMap["TRIPLE_HET"] = TRIPLE_HET;
	    codeMap["EMPIRICAL_WITHIN_HOST_MODEL"] = EMPIRICAL_WITHIN_HOST_MODEL;
        codeMap["MOLINEAUX_WITHIN_HOST_MODEL"] = MOLINEAUX_WITHIN_HOST_MODEL;
        codeMap["PENNY_WITHIN_HOST_MODEL"] = PENNY_WITHIN_HOST_MODEL;
        codeMap["GARKI_DENSITY_BIAS"] = GARKI_DENSITY_BIAS;
//      codeMap["IPTI_SP_MODEL"] = IPTI_SP_MODEL;
//      codeMap["REPORT_ONLY_AT_RISK"] = REPORT_ONLY_AT_RISK; 
	    codeMap["MEAN_DURATION_GAMMA"] = MEAN_DURATION_GAMMA;
	    codeMap["FIRST_LOCAL_MAXIMUM_GAMMA"] = FIRST_LOCAL_MAXIMUM_GAMMA;
	    codeMap["IMMUNE_THRESHOLD_GAMMA"] = IMMUNE_THRESHOLD_GAMMA;
	    codeMap["UPDATE_DENSITY_GAMMA"] = UPDATE_DENSITY_GAMMA;
	    codeMap["PARASITE_REPLICATION_GAMMA"] = PARASITE_REPLICATION_GAMMA;
        codeMap["VECTOR_LIFE_CYCLE_MODEL"] = VECTOR_LIFE_CYCLE_MODEL;
        codeMap["VECTOR_SIMPLE_MPD_MODEL"] = VECTOR_SIMPLE_MPD_MODEL;
        codeMap["MOLINEAUX_PAIRWISE_SAMPLE"] = MOLINEAUX_PAIRWISE_SAMPLE;
        ignoreOptions.insert("PROPHYLACTIC_DRUG_ACTION_MODEL");
        codeMap["VIVAX_SIMPLE_MODEL"] = VIVAX_SIMPLE_MODEL;
        codeMap["INDIRECT_MORTALITY_FIX"] = INDIRECT_MORTALITY_FIX;
        codeMap["VACCINE_GENOTYPE"] = VACCINE_GENOTYPE;
        codeMap["CFR_PF_USE_HOSPITAL"] = CFR_PF_USE_HOSPITAL;
        codeMap["HEALTH_SYSTEM_MEMORY_FIX"] = HEALTH_SYSTEM_MEMORY_FIX;
	}
	
	OptionCodes operator[] (const string s) {
	    auto codeIt = codeMap.find (s);
	    if( codeIt != codeMap.end() ) return codeIt->second;
            if( ignoreOptions.count(s) ){
                if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
                    cerr << "Deprecation warning: model option " << s << " is no longer used" << endl;
                }
                return IGNORE;
            }
            
            ostringstream msg;
            if( s == "PENALISATION_EPISODES" || s == "ATTENUATION_ASEXUAL_DENSITY" ){
                msg << "Please use schema 31 or earlier to use option "
                    << s << "; it is not available in later versions.";
            }else if( s == "IPTI_SP_MODEL" ){
                msg << "The IPT model is no longer available. Use MDA instead.";
            }else if( s == "REPORT_ONLY_AT_RISK" ){
                msg << "Option " << s << " has been replaced by <SurveyOptions onlyNewEpisode=\"true\">.";
            }else{
                msg << "Unrecognised model option: " << s;
            }
            throw xml_scenario_error(msg.str());
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const OptionCodes code) {
	    for(auto codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw TRACED_EXCEPTION_DEFAULT ("toString called with unknown code");	// this is a code error
	}
    };

    std::bitset<NUM_OPTIONS> ModelOptions::getBaseModelOptions()
    {
        // The "base" model is very simple: *no* model options are turned on.
        bitset<NUM_OPTIONS> allOptionsOff;

        return allOptionsOff;
    }

    std::bitset<NUM_OPTIONS> ModelOptions::getLegacyDefaultModelOptions()
    {
        bitset<NUM_OPTIONS> defaultOptSet;

        // Kept for legacy reasons.  That is, when an input XML does not specify a model name
        // and does not set these model options to false, they are turned on *for backwards
        // compatibility*.
        defaultOptSet.set (MAX_DENS_CORRECTION);
        defaultOptSet.set (INNATE_MAX_DENS);
        defaultOptSet.set (INDIRECT_MORTALITY_FIX);
        defaultOptSet.set (HEALTH_SYSTEM_MEMORY_FIX);

        return defaultOptSet;
    }

    void ModelOptions::init(scnXml::Model::ModelOptionsOptional& optionsElt,
                            util::ModelNameProvider mnp)
    {
    OptionCodeMap codeMap;

	// Set options to legacy defaults, then override any given in the XML file:
	options = getLegacyDefaultModelOptions();

    const util::ModelNames namedModelToUse = mnp.GetModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            // This completely discards any legacy model options set above.
            options = getBaseModelOptions();
        }
        else
        {
            throw util::xml_scenario_error(
                "No pre-set model options are available for the specified model name.");
        }
    }

    // Apply any user-specified overrides.
    if (optionsElt.present())
    {
	    const scnXml::OptionSet::OptionSequence& optSeq = optionsElt.get().getOption();
	    for(auto it = optSeq.begin(); it != optSeq.end(); ++it) {
                OptionCodes opt = codeMap[it->getName()];
                if( opt != IGNORE ) options[opt] = it->getValue();
	    }
    }
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_OPTIONS)) {
        const std::bitset<NUM_OPTIONS> defaultOptSet = getLegacyDefaultModelOptions();
	    cout << "Non-default model options:";
	    for(int i = 0; i < NUM_OPTIONS; ++i) {
		if (options[i] != defaultOptSet[i])
		    cout << "\t" << codeMap.toString(OptionCodes(i)) << "=" << options[i];
	    }
	    cout << endl;
	}
	
	// Test for incompatible options
	
	// An incompatibility triangle, listing options incompatible with each option
	// Doesn't check required versions
	// Originally from "description of variables for interface" excel sheet
	bitset<NUM_OPTIONS> incompatibilities[NUM_OPTIONS];	// all default to 0
	
	incompatibilities[NEGATIVE_BINOMIAL_MASS_ACTION]
	    .set(LOGNORMAL_MASS_ACTION)
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[LOGNORMAL_MASS_ACTION]
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	
	// Note: MAX_DENS_CORRECTION is irrelevant when using new
	// within-host models, but we don't mark it incompatible so that we can
	// leave MAX_DENS_CORRECTION on by default.
	incompatibilities[DUMMY_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(EMPIRICAL_WITHIN_HOST_MODEL)
	    .set(MOLINEAUX_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);

	incompatibilities[EMPIRICAL_WITHIN_HOST_MODEL]
            .set(MOLINEAUX_WITHIN_HOST_MODEL)
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	    
	incompatibilities[MOLINEAUX_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	incompatibilities[PENNY_WITHIN_HOST_MODEL]
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA);
	    
	incompatibilities[NON_MALARIA_FEVERS]
	    .set(MUELLER_PRESENTATION_MODEL);
	
	incompatibilities[TRANS_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TREAT_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TRANS_HET]
	    .set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TRANS_TREAT_HET]
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TREAT_HET]
	    .set(TRIPLE_HET);
        
        incompatibilities[MOLINEAUX_PAIRWISE_SAMPLE]
            .set(FIRST_LOCAL_MAXIMUM_GAMMA)
            .set(MEAN_DURATION_GAMMA);
        
        
        incompatibilities[VIVAX_SIMPLE_MODEL]
            .set( DUMMY_WITHIN_HOST_MODEL )
            .set( EMPIRICAL_WITHIN_HOST_MODEL )
            .set( MOLINEAUX_WITHIN_HOST_MODEL )
            .set( PENNY_WITHIN_HOST_MODEL );
        
	for(size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if (options [i] && (options & incompatibilities[i]).any()) {
		ostringstream msg;
		msg << "Incompatible model options: " << codeMap.toString(OptionCodes(i)) << "=" << options[i]
			<< " is incompatible with flags:";
		bitset<NUM_OPTIONS> incompat = (options & incompatibilities[i]);
		for(int j = 0; j < NUM_OPTIONS; ++j) {
		    if (incompat[j])
			msg << "\t" << codeMap.toString(OptionCodes(j)) << "=" << options[j];
		}
		throw xml_scenario_error (msg.str());
	    }
	}
	
        // Required options (above table can't check these):
        if (options[INNATE_MAX_DENS] && !options[MAX_DENS_CORRECTION])
            throw xml_scenario_error ("INNATE_MAX_DENS requires MAX_DENS_CORRECTION");
        if( !options[MOLINEAUX_WITHIN_HOST_MODEL] && (
            options[FIRST_LOCAL_MAXIMUM_GAMMA] ||
            options[MEAN_DURATION_GAMMA] ||
            options[PARASITE_REPLICATION_GAMMA] ) )
            throw xml_scenario_error( "Molineaux model option used without MOLINEAUX_WITHIN_HOST_MODEL option" );
        if( !options[PENNY_WITHIN_HOST_MODEL] && (
            options[IMMUNE_THRESHOLD_GAMMA] ||
            options[UPDATE_DENSITY_GAMMA] ) )
            throw xml_scenario_error( "Penny model option used without PENNY_WITHIN_HOST_MODEL option" );
        
        if( sim::oneTS() == sim::fromDays(5) ){
            // 5 day TS is okay; some tests specific to this TS:
            bitset<NUM_OPTIONS> require1DayTS;
            require1DayTS
                .set( CLINICAL_EVENT_SCHEDULER );
            
            for(size_t i = 0; i < NUM_OPTIONS; ++i) {
                if (options [i] && require1DayTS[i]) {
                    ostringstream msg;
                    msg << "Model option " << codeMap.toString(OptionCodes(i)) << " is only compatible with a 1-day time step.";
                    throw xml_scenario_error (msg.str());
                }
            }
        }else if( sim::oneTS() == sim::fromDays(1) ){
            // 1 day TS is also okay
        }else{
            ostringstream msg;
            msg << "Time step set to " << sim::oneTS() << " days but only 1 and 5 days are supported.";
            throw xml_scenario_error (msg.str());
        }
    }
    
} }

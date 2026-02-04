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

#include "mon/Continuous.h"
#include "mon/info.h"   // lineEnd
#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "schema/monitoring.h"

#include <vector>
#include <map>
#include <fstream>
#include <gzstream/gzstream.h>

namespace OM {
    namespace mon {
        using util::xml_scenario_error;

        /// File we output to
        string cts_filename;

        /// This is used to output some statistics in a tab-deliminated-value file.
        /// (It used to be csv, but German Excel can't open csv directly.)
        fstream ctsOStream;

        /* Record last position in file (as position minus start), for checkpointing.
         * Don't use a streampos directly, because I'm not convinced we can save and
         * reload a streampos and use on a new file. */
        streamoff streamOff;
        streampos streamStart;

        map<string, tuple<string, std::function<void(Population&, ostream&)>> > registered;

        // List that we report.
        vector< tuple<string, std::function<void(Population&, ostream&)>> > toReport;

        SimTime ctsPeriod = sim::zero();

        bool duringInit = false;

        ContinuousType Continuous;

        ContinuousType::~ContinuousType() {}

        /* Initialise: enable outputs registered and requested in XML.
         * Search for Continuous::registerCallback to see outputs available. */
        void ContinuousType::init (const scnXml::Monitoring& monitoring, bool isCheckpoint) {
            const scnXml::Monitoring::ContinuousOptional& ctsOpt = monitoring.getContinuous();
            if( ctsOpt.present() == false ) {
                ctsPeriod = sim::zero();
                return;
            }
            try{
                //NOTE: if changing XSD, this should not have a default unit:
                ctsPeriod = UnitParse::readShortDuration( ctsOpt.get().getPeriod(), UnitParse::STEPS );
                if( ctsPeriod < sim::oneTS() )
                    throw util::format_error("must be >= 1 time step");
            }catch( const util::format_error& e ){
                throw xml_scenario_error( string("monitoring/continuous/period: ").append(e.message()) );
            }

            if( ctsOpt.get().getDuringInit().present() )
                duringInit = ctsOpt.get().getDuringInit().get();

            cts_filename = util::CommandLine::getCtsoutName();

            ctsOStream.width (0);

            if( isCheckpoint ){
                scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
                for(scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
                    auto reg_it = registered.find( it->getName() );
                    if( reg_it == registered.end() )
                        throw xml_scenario_error("monitoring.continuous: no output " + string(it->getName()));
                    if( it->getValue() ){
                        toReport.push_back( reg_it->second );
                    }
                }

                // When loading a check-point, we resume reporting to this file.
                // Use "ate" mode and seek to desired pos.
                ctsOStream.open (cts_filename.c_str(), ios::binary|ios::ate|ios::in|ios::out );
                if( ctsOStream.fail() )
                    throw util::checkpoint_error ("Continuous: resume error (no file)");
                ctsOStream.seekp( 0, ios_base::beg );
                streamStart = ctsOStream.tellp();
                // we set position later, in staticCheckpoint
            }else{
                ctsOStream.open( cts_filename.c_str(), ios::binary|ios::out );
                streamStart = ctsOStream.tellp();
                ctsOStream << "##\t##" << endl;	// live-graph needs a deliminator specifier when it's not a comma

                if( duringInit )
                    ctsOStream << "simulation time\t";
                ctsOStream << "timestep";   //TODO: change to days or remove or leave?
                scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
                for(scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
                    auto reg_it = registered.find( it->getName() );
                    if( reg_it == registered.end() )
                        throw xml_scenario_error("monitoring.continuous: no output " + string(it->getName()));
                    if( it->getValue() ){
                        ctsOStream << std::get<0>(reg_it->second);
                        toReport.push_back( reg_it->second );
                    }
                }
                ctsOStream << mon::lineEnd << flush;
                streamOff = ctsOStream.tellp() - streamStart;
            }
        }

        void ContinuousType::checkpoint (ostream& stream){
            if( ctsPeriod == sim::zero() )
                return;	// output disabled

            streamOff & stream;
        }
        void ContinuousType::checkpoint (istream& stream){
            if( ctsPeriod == sim::zero() )
                return;	// output disabled

            /* We attempt to resume output correctly after a reload by recording
            * the last position, and relocating there.
            * 
            * (Keeping results in memory until end of sim would be another,
            * slightly safer, option, but loses real-time output.) */
            streamOff & stream;
            // We skip back to the last write-point, so anything written after the
            // last checkpoint will be repeated:
            ctsOStream.seekp( streamOff, ios_base::beg );

            if( ctsOStream.fail() )
                throw util::checkpoint_error ("Continuous: resume error (bad pos/file)");
        }

        void ContinuousType::registerCallback (string optName, string titles, function<void(ostream&)> f){
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population&, ostream& ostream){ f(ostream); };
            registered[optName] = {titles, _f};
        }

        void ContinuousType::registerCallback (string optName, string titles, function<void(const vector<Host::Human> &, ostream&)> f){
            assert(registered.count(optName) == 0); // name clash/registered twice?
            function<void(const Population&, ostream&)> _f = [f](const Population &p, ostream& ostream){ f(p.humans, ostream); };
            registered[optName] = {titles, _f};
        }

        void ContinuousType::registerCallback (string optName, string titles, function<void(Population &, ostream&)> f){
            assert(registered.count(optName) == 0); // name clash/registered twice?
            registered[optName] = {titles, f};
        }

        void ContinuousType::update (Population &population){
            if( ctsPeriod == sim::zero() )
                return;	// output disabled
            if( !duringInit ){
                if( sim::intervTime() < sim::zero()
                    || mod_nn(sim::intervTime(), ctsPeriod) != sim::zero() )
                    return;
            } else {
                if( mod_nn(sim::now(), ctsPeriod) != sim::zero() )
                    return;
                ctsOStream << sim::inSteps(sim::now()) << '\t';
            }

            if( duringInit && sim::intervTime() < sim::zero() ){
                ctsOStream << "nan";
            }else{
                // NOTE: we could switch this to output dates, but (1) it would be
                // breaking change and (2) it may be harder to use.
                ctsOStream << sim::inSteps(sim::intervTime());
            }
            for( size_t i = 0; i < toReport.size(); ++i )
                std::get<1>(toReport[i])( population, ctsOStream );

            // We must flush often to avoid temporarily outputting partial lines
            // (resulting in incorrect real-time graphs).
            ctsOStream << mon::lineEnd << flush;

            streamOff = ctsOStream.tellp() - streamStart;
        }
    }
}

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

#ifndef Hmod_util_CommandLine
#define Hmod_util_CommandLine

#include "Global.h"
#include <string>
#include <set>
#include <bitset>
#include <limits>
using namespace std;

namespace OM { namespace util {
    /// Command line options and processing
	class CommandLine {
	public:
	/// Boolean command-line options
		enum Options {
	    /** Outputs non-default "ModelOptions" values in a human-readable form. */
			PRINT_MODEL_OPTIONS = 0,
	    /// Verbose output
			VERBOSE,
		/// Show progress
			PROGRESS,
	    /// Forces checkpointing just before starting the main phase.
			CHECKPOINT,
	    /// Exit after writing checkpoint
			CHECKPOINT_STOP,
	    /** Do initialisation and error checks, but don't run simulation. */
			SKIP_SIMULATION,
            /** Compress output.txt file. */
			COMPRESS_OUTPUT,
	    /** Print the annual EIR. */
			PRINT_ANNUAL_EIR,
            /** Outputs samples from the active interpolation methods of all
             * age-group data suitible for graphing. */
			SAMPLE_INTERPOLATIONS,
            /** Show details of vector-parameter fitting.
             * 
             * The fitting methods used aren't guaranteed to work. If they don't, this output should help work out why. */
			DEBUG_VECTOR_FITTING,
            /** Print out details about interventions. */
			PRINT_INTERVENTIONS,
            /** Warn on use of deprecated features; that is recommend the use
             * of replacement features. */
			DEPRECATION_WARNINGS,
            /** Print times of all surveys. */
			PRINT_SURVEY_TIMES,
			PRINT_GENOTYPES,
			NUM_OPTIONS
		};

	/** Return true if given option (from CommandLine::Options) is active. */
		static inline bool option(size_t code) {
			return options.test(code);
		}

	/** If path is relative, prepend it with clResourcePath. */
		static string lookupResource (const string& path);

        /** Get the name of the output file. */
		static inline string getOutputName (){
			return outputName;
		}

    /** Get the name of the ctsout file. */
		static inline string getCtsoutName (){
			return ctsoutName;
		}

     /** Get the name of the checkpoint file. */
		static inline string getCheckpointName (){
			return checkpointFileName;
		}

	/** Looks through all command line options.
	*
	* @returns The name of the scenario XML file to use.
	*
	* Throws cmd_exit in the case a help message is printed. Help
	* is printed to cout.
	* 
	* In other cases command-line parameters cause variables to be set in Global
	* to achieve the desired result. */
		static string parse (int argc, char* argv[]);

	/** @brief Checkpointing.
	*
	* Not really required; more to confirm things are expected. */
		static void staticCheckpoint (istream& stream);
	static void staticCheckpoint (ostream& stream);	///< ditto
	
private:
	// Static parameters: don't require checkpointing (set by parse())
	static bitset<NUM_OPTIONS> options;	// boolean options
	static string resourcePath;
	
	//Output filename (for main output file "output.txt")
	static string outputName;
	static string ctsoutName;
	static string checkpointFileName;
};
} }
#endif

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

#include "Global.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/DocumentLoader.h"
/* if you get compile errors like "version.h not found", run CMake first */
#include "util/version.h"

#include <sstream>
#include <iostream>
#include <cassert>

namespace OM { namespace util {
	bitset<CommandLine::NUM_OPTIONS> CommandLine::options;
	string CommandLine::resourcePath;
	string CommandLine::outputName;
	string CommandLine::ctsoutName;
	string CommandLine::checkpointFileName;

	string parseNextArg (int argc, char* argv[], int& i) {
		++i;
		if (i >= argc)
			throw cmd_exception ("Expected an argument following the last option");
		return string(argv[i]);
	}

	string CommandLine::lookupResource (const string& path) {
		string ret;
		if (path.size() >= 1 && path[0] == '/') {
		    // UNIX absolute path
		} else if (path.size() >= 3 && path[1] == ':'
			&& (path[2] == '\\' || path[2] == '/')) {
		    // Windows absolute path.. at least probably
		} else {	// relative
			ret = resourcePath;
		}
		ret.append (path);
		return ret;
	}

	string CommandLine::parse (int argc, char* argv[]) {
		bool cloHelp = false, cloVersion = false, cloError = false;
		string scenarioFile = "";
		outputName = "";
		ctsoutName = "";
#	ifdef OM_STREAM_VALIDATOR
		string sVFile;
#	endif

		// Set DEPRECATION_WARNINGS to ON by default
		options.set (DEPRECATION_WARNINGS);

	/* Simple command line parser. Seems to work fine.
	* If an extension is wanted, http://tclap.sourceforge.net/ looks good. */
		for(int i = 1; i < argc; ++i) {
			string clo = argv[i];

	    // starts "--"
			if (clo.size() >= 2 && *clo.data() == '-' && clo.data()[1] == '-') {
				clo.assign (clo, 2, clo.size()-2);
				if (clo == "verbose") {
					options.set (VERBOSE);
				} else if (clo == "progress") {
					options.set (PROGRESS);
				} else if (clo == "resource-path") {
					if (resourcePath.size())
						throw cmd_exception ("--resource-path (or -p) may only be given once");
					resourcePath = parseNextArg (argc, argv, i).append ("/");
				} else if (clo == "scenario") {
					if (scenarioFile != ""){
						throw cmd_exception ("--scenario argument may only be given once");
					}
					scenarioFile = parseNextArg (argc, argv, i);
				} else if (clo == "output") {
					if (outputName != ""){
						throw cmd_exception ("--output argument may only be given once");
					}
					outputName = parseNextArg (argc, argv, i);
				} else if (clo == "compress-output") {
					options.set (COMPRESS_OUTPUT);
				} else if (clo == "ctsout") {
					if (ctsoutName != ""){
						throw cmd_exception ("--ctsout argument may only be given once");
					}
					ctsoutName = parseNextArg (argc, argv, i);
				} else if (clo == "name") {
					if (ctsoutName != "" || outputName != "" || scenarioFile != ""){
						throw cmd_exception ("--name may not be used along with --scenario, --output or --ctsout");
					}
					string name = parseNextArg (argc, argv, i);
					(scenarioFile = "scenario").append(name).append(".xml");
					(outputName = "output").append(name).append(".txt");
					(ctsoutName = "ctsout").append(name).append(".txt");
				} else if (clo == "validate-only") {
					options.set (SKIP_SIMULATION);
				} else if (clo == "no-deprecation-warnings") {
					options.reset (DEPRECATION_WARNINGS);
				} else if (clo == "print-model") {
					options.set (PRINT_MODEL_OPTIONS);
					options.set (SKIP_SIMULATION);
				} else if (clo == "print-EIR") {
					options.set (PRINT_ANNUAL_EIR);
					options.set (SKIP_SIMULATION);
				} else if (clo == "print-interventions") {
					options.set (PRINT_INTERVENTIONS);
					options.set (SKIP_SIMULATION);
				} else if (clo == "print-survey-times") {
					options.set (PRINT_SURVEY_TIMES);
					options.set (SKIP_SIMULATION);
				} else if (clo == "print-genotypes") {
					options.set (PRINT_GENOTYPES);
					options.set (SKIP_SIMULATION);
				} else if (clo == "sample-interpolations") {
					options.set (SAMPLE_INTERPOLATIONS);
					options.set (SKIP_SIMULATION);
				} else if (clo == "checkpoint") {
					options.set (CHECKPOINT);
				} else if (clo == "checkpoint-file") {
					if (checkpointFileName != ""){
						throw cmd_exception ("--checkpoint-file argument may only be given once");
					}
					options.set (CHECKPOINT);
					checkpointFileName = parseNextArg (argc, argv, i);
				} else if (clo == "checkpoint-stop") {
					options.set (CHECKPOINT);
					options.set (CHECKPOINT_STOP);
				} else if (clo == "debug-vector-fitting") {
					options.set (DEBUG_VECTOR_FITTING);
#	ifdef OM_STREAM_VALIDATOR
				} else if (clo == "stream-validator") {
					if (sVFile.size())
						throw cmd_exception ("--stream-validator may only be given once");
					sVFile = parseNextArg (argc, argv, i);
#	endif
				} else if (clo == "version") {
					cloVersion = true;
				} else if (clo == "help") {
					cloHelp = true;
				} else {
					cerr << "Unrecognised command-line option: --" << clo << endl;
					cloError = true;
				}
	    } else if (clo.size() >= 1 && *clo.data() == '-') {	// single - (not --)
	    	for(size_t j = 1; j < clo.size(); ++j) {
	    		if (clo[j] == 'V') {
	    			options.set (VERBOSE);
	    		} else if (clo[j] == 'p') {
	    			if (j + 1 != clo.size())
	    				throw cmd_exception ("a path must be given as next argument after -p");
	    			if (resourcePath.size())
	    				throw cmd_exception ("--resource-path (or -p) may only be given once");
	    			resourcePath = parseNextArg (argc, argv, i).append ("/");
	    		} else if (clo[j] == 'm') {
	    			options.set (PRINT_MODEL_OPTIONS);
	    			options.set (SKIP_SIMULATION);
	    		} else if (clo[j] == 's') {
	    			if (scenarioFile != ""){
	    				throw cmd_exception ("-s argument may only be given once");
	    			}
	    			scenarioFile = parseNextArg (argc, argv, i);
	    		} else if (clo[j] == 'o') {
	    			if (outputName != ""){
	    				throw cmd_exception ("-o argument may only be given once");
	    			}
	    			outputName = parseNextArg (argc, argv, i);
	    		} else if (clo[j] == 'n') {
	    			if (ctsoutName != "" || outputName != "" || scenarioFile != ""){
	    				throw cmd_exception ("--name may not be used along with --scenario, --output or --ctsout");
	    			}
	    			string name = parseNextArg (argc, argv, i);
	    			(scenarioFile = "scenario").append(name).append(".xml");
	    			(outputName = "output").append(name).append(".txt");
	    			(ctsoutName = "ctsout").append(name).append(".txt");
	    		} else if (clo[j] == 'c') {
	    			options.set (CHECKPOINT);
	    		} else if (clo[j] == 'v') {
	    			cloVersion = true;
	    		} else if (clo[j] == 'z') {
	    			options.set (COMPRESS_OUTPUT);
	    		} else if (clo[j] == 'h') {
	    			cloHelp = true;
	    		} else {
	    			cerr << "Unrecognised command-line option: -" << clo[j] << endl;
	    			cloError = true;
	    		}
	    	}
	    } else {
	    	cerr << "Unexpected parameter: " << clo << endl << endl;
	    	cloError = true;
	    }
	}
	
	if (cloVersion || cloHelp){
		cerr<<"OpenMalaria simulator of malaria epidemiology and control."<<endl<< endl
		<<"For more information, see https://github.com/SwissTPH/openmalaria/wiki"<<endl<<endl
		<<"\tschema version: \t"   <<SCHEMA_VERSION<<endl
		<<"\tprogram version:\t" << util::semantic_version <<endl<<endl
		<<"\tCopyright © 2005-2026 Swiss Tropical and Public Health Institute"<<endl
		<<"\tCopyright © 2005-2015 Liverpool School Of Tropical Medicine"<<endl
		<<"\tCopyright © 2020-2026 University of Basel"<<endl
		<<"\tCopyright © 2025-2026 The Kids Research Institute Australia"<<endl<<endl
		<<"OpenMalaria comes with ABSOLUTELY NO WARRANTY. This is free software, and you"<<endl<<endl
		<<"are welcome to redistribute it under certain conditions. See the file COPYING"<<endl
		<<"or http://www.gnu.org/licenses/gpl-2.0.html for details of warranty or terms of"<<endl
		<<"redistribution."<<endl<<endl;
	}

	if (cloHelp || cloError) {
		cerr << "Usage: " << argv[0] << " [options]" << endl << endl
		<< "Options:"<<endl
		<< " -V --verbose		Verbose output" <<endl
		<< " -p --resource-path	Path to look up input resources with relative URLs (defaults to"<<endl
		<< "			working directory). Not used for output files."<<endl
		<< " -s --scenario file.xml	Uses file.xml as the scenario. If not given, scenario.xml is used." << endl
		<< "			If path is relative (doesn't start '/'), --resource-path is used."<<endl
		<< " -o --output file.txt	Uses file.txt as output file name. If not given, output.txt is used." << endl
		<< "    --ctsout file.txt	Uses file.txt as ctsout file name. If not given, ctsout.txt is used." << endl
		<< " -n --name NAME		Equivalent to --scenario scenarioNAME.xml --output outputNAME.txt \\"<<endl
		<< "			--ctsout ctsoutNAME.txt" <<endl
		<< " -z --compress-output	Compress output with gzip (writes output.txt.gz)." << endl
		<< "    --validate-only	Initialise and validate scenario, but don't run simulation." << endl
		<< "    --no-deprecation-warnings" << endl
		<< "			OpenMalaria warn about the use of features deemed error-prone and where" << endl
		<< "			more flexible alternatives are available. Use this option to silence it." << endl
		<< endl
		<< "Debugging options:"<<endl
		<< " -m --print-model	Print all model options with a non-default value and exit." << endl
		<< "    --print-EIR		Print the annual EIR (of each species in vector mode) and exit." << endl
		<< "    --print-interventions" << endl
		<< "			Print intervention deployment details and exit." << endl
		<< "    --print-survey-times" << endl
		<< "			Print out the times of all surveys and exit." << endl
		<< "    --print-genotypes" << endl
		<< "                        Print out genotype ids and exit." << endl
		<< "    --sample-interpolations" <<endl
		<< "			Output samples of all used age-group data according to active"<<endl
		<< "			interpolation method and exit."<<endl
		<< " -c --checkpoint	Write a checkpoint just before starting the main phase."<<endl
		<< "			This may be used to skip redundant computation when multiple"<<endl
		<< "			simulations differ only during the intervention phase."<<endl
		<< "    --checkpoint-file file	Checkpoint as above. Uses file as checkpoint file name. If not given, checkpoint is used." << endl
		<< "    --checkpoint-stop	Checkpoint as above, then stop immediately afterwards. Can be used with --checkpoint-file."<<endl
		<< "    --debug-vector-fitting"<<endl
		<< "			Show details of vector-parameter fitting. The fitting methods used" <<endl
		<< "			aren't guaranteed to work. If they don't, this output should help"<<endl
		<< "			work out why."<<endl
#	ifdef OM_STREAM_VALIDATOR
		<< "    --stream-validator PATH" <<endl
		<< "			Use StreamValidator to validate against reference file PATH." <<endl
		<< "			(note: PATH must be absolute or relative to resource path)." <<endl
#	endif
		<< " -v --version           Display the current schema version of OpenMalaria." << endl
		<< " -h --help              Print this message." << endl<<endl
		;
		if( cloError )
			throw cmd_exception( "bad argument" );
	}
	if( cloVersion || cloHelp ){
		throw cmd_exception("Printed help",Error::None);
	}
	
#	ifdef OM_STREAM_VALIDATOR
	if( sVFile.size() )
		StreamValidator.loadStream( sVFile );
#	endif
	
	if (scenarioFile == ""){
		scenarioFile = "scenario.xml";
	}
	if (outputName == ""){
		outputName = "output.txt";
	}
	if (ctsoutName == ""){
		ctsoutName = "ctsout.txt";
	}

	return lookupResource(scenarioFile);
}

    /* These check parameters are as expected. They only really serve to make
     * sure important command-line parameters didn't change (and only in DEBUG mode)! */
void CommandLine::staticCheckpoint (istream& stream) {
	string tOpt;
	string tResPath;
	tOpt & stream;
	tResPath & stream;
	assert (tOpt == options.to_string());
	assert (tResPath == resourcePath);
}
void CommandLine::staticCheckpoint (ostream& stream) {
	options.to_string() & stream;
	resourcePath & stream;
}

} }

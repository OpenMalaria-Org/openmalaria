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

#include <iostream>
#include <fstream>
#include <gzstream/gzstream.h>

#include "Global.h"
#include "mon/Continuous.h"
#include "mon/management.h"
#include "interventions/InterventionManager.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"

#include "checkpoint.h"

namespace OM
{

/** @brief checkpointing functions
*
* readCheckpoint/writeCheckpoint prepare to read/write the file,
* and read/write read and write the actual data. */
//@{
int readCheckpointNum (const string &checkpointFileName)
{
    ifstream checkpointFile;
    checkpointFile.open(checkpointFileName, fstream::in);
    int checkpointNum=0;
    checkpointFile >> checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
        throw util::checkpoint_error ("error reading from file \"checkpoint\"");
    return checkpointNum;
}

void checkpoint (istream& stream, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission)
{
    try
    {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        Host::NeonatalMortality::staticCheckpoint (stream);
        Host::InfectionIncidenceModel::staticCheckpoint (stream);
        Clinical::InfantMortality::staticCheckpoint(stream);
        WithinHost::Genotypes::staticCheckpoint(stream);
        mon::Continuous & stream;
        mon::checkpoint( stream );
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        sim::s_interv & stream;
        endTime & stream;
        estEndTime & stream;
        transmission & stream;
        population.checkpoint(stream);
        interventions::InterventionManager::checkpoint(stream);
        interventions::InterventionManager::loadFromCheckpoint(population.humans, transmission);
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::s_t0 & stream;
        sim::s_t1 & stream;
        util::master_RNG.checkpoint(stream);
    } catch (const util::checkpoint_error& e) { // append " (pos X of Y bytes)"
        ostringstream pos;
        pos<<" (pos "<<stream.tellg()<<" of ";
        stream.ignore (numeric_limits<streamsize>::max()-1);    // skip to end of file
        pos<<stream.tellg()<<" bytes)";
        throw util::checkpoint_error( e.what() + pos.str() );
    }
    
    
    stream.ignore (numeric_limits<streamsize>::max()-1);        // skip to end of file
    if (stream.gcount () != 0) {
        ostringstream msg;
        msg << "Checkpointing file has " << stream.gcount() << " bytes remaining." << endl;
        throw util::checkpoint_error (msg.str());
    } else if (stream.fail())
        throw util::checkpoint_error ("stream read error");
}

void checkpoint (ostream& stream, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission) {
    util::checkpoint::header (stream);
    if (!stream.good())
        throw util::checkpoint_error ("Unable to write to file");

    util::CommandLine::staticCheckpoint (stream);
    Host::NeonatalMortality::staticCheckpoint (stream);
    Host::InfectionIncidenceModel::staticCheckpoint (stream);
    Clinical::InfantMortality::staticCheckpoint(stream);
    WithinHost::Genotypes::staticCheckpoint(stream);
    mon::Continuous & stream;
    mon::checkpoint( stream );
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator & stream;
# endif
    
    sim::s_interv & stream;
    endTime & stream;
    estEndTime & stream;
    transmission & stream;
    population.checkpoint(stream);
    interventions::InterventionManager::checkpoint( stream );
    
    sim::s_t0 & stream;
    sim::s_t1 & stream;
    util::master_RNG.checkpoint(stream);
    
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

void writeCheckpoint(const bool startedFromCheckpoint, const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission)
{
    // We alternate between two checkpoints, in case program is closed while writing.
    const int NUM_CHECKPOINTS = 2;
    
    int oldCheckpointNum = 0, checkpointNum = 0;
    if (startedFromCheckpoint)
    {
        oldCheckpointNum = readCheckpointNum(checkpointFileName);
        checkpointNum = mod_nn(oldCheckpointNum + 1, NUM_CHECKPOINTS); // Get next checkpoint number:
    }
    
    {   // Open the next checkpoint file for writing:
        ostringstream name;
        name << checkpointFileName << checkpointNum << ".gz";
        ogzstream out(name.str().c_str(), ios::out | ios::binary);
        checkpoint (out, endTime, estEndTime, population, transmission);
        out.close();
    }
    
    {   // Indicate which is the latest checkpoint file.
        ofstream checkpointFile;
        checkpointFile.open(checkpointFileName,ios::out);
        checkpointFile << checkpointNum;
        checkpointFile.close();
        if (!checkpointFile)
            throw util::checkpoint_error ("error writing to file \"checkpoint\"");
    }

    // Truncate the old checkpoint to save disk space, when it existed
    if( oldCheckpointNum != checkpointNum ){
        ostringstream name;
        name << checkpointFileName << oldCheckpointNum << ".gz";
        ofstream out(name.str().c_str(), ios::out | ios::binary);
        out.close();
    }
}

void readCheckpoint(const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, Transmission::TransmissionModel &transmission)
{
    int checkpointNum = readCheckpointNum(checkpointFileName);
    
    // Open the latest file
    ostringstream name;
    name << checkpointFileName << checkpointNum << ".gz";
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    //Note: gzstreams are considered "good" when file not open!
    if ( !( in.good() && in.rdbuf()->is_open() ) )
        throw util::checkpoint_error ("Unable to read file");
    checkpoint (in, endTime, estEndTime, population, transmission);
    in.close();
  
    cerr << sim::inSteps(sim::now()) << "t loaded checkpoint" << endl;
}

}
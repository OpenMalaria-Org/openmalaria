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

#ifndef Hmod_WH_Path_Submodels
#define Hmod_WH_Path_Submodels

#include "Global.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"

using namespace std;

namespace OM { namespace WithinHost { namespace Pathogenesis {

/// Müller presentation model.
class MuellerPathogenesis : public PathogenesisModel {
public:
    MuellerPathogenesis(double cF) :
        PathogenesisModel(cF) {}
    ~MuellerPathogenesis() {}
    
    virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);

    /// Read parameters
    static void init( const Parameters& parameters );
};

/// Pyrogenic threshold presentation model.
class PyrogenPathogenesis : public PathogenesisModel {
protected:
    /// Critical density for fever (clinical episodes)
    double _pyrogenThres;
    /// Determine the current pyrogenic threshold.
    virtual void updatePyrogenThres(double totalDensity);

public:
    PyrogenPathogenesis(double cF);
    virtual ~PyrogenPathogenesis() {}
    virtual void summarize (const Host::Human& human);
    virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);
    
    /// Read parameters from XML
    static void init( const OM::Parameters& parameters );
    
protected:
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
};

/// Predetermined episodes presentation model.
class PredetPathogenesis : public PyrogenPathogenesis {
public:
    PredetPathogenesis (double cF) :
        PyrogenPathogenesis(cF) {}
    ~PredetPathogenesis() {}
    
    virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);
};

} } }
#endif

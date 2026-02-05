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

#ifndef OM_INTERVENTIONS_INTERFACES
#define OM_INTERVENTIONS_INTERFACES

#include "Global.h"
#include "mon/reporting.h"
#include <schema/interventions.h>
#include <limits>

namespace scnXml{ class DeploymentBase; }

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {

/** Enumeration of all available components, in the order that these should be
 * deployed in within a single intervention. */
namespace Component { enum Type {
    RECRUIT_ONLY,     // selection for a sub-population (without other effects)
    SCREEN,     // screening, e.g. as part of MSAT
    TREAT_SIMPLE,       // treatment using the simple model (e.g. for MDA)
    TREAT_PKPD, // treatment using the PK/PD model (e.g. MDA)
    CM_DT,      // execute a "case management" decision tree (e.g. MSAT or more complicated things)
    PEV,        // pre-erythrocytic vaccine
    BSV,        // blood-stage vaccine
    TBV,        // transmission-blocking vaccine
    IPT,        // intermittent preventative treatment
    ITN,        // insecticide treated net
    IRS,        // indoor residual spraying
    GVI,        // generic vector intervention
    CLEAR_IMMUNITY,     // reset accumulated immunity to zero
}; }

/** Specifies limits on the number of existing doses when deciding whether to
 * vaccinate a human. */
struct VaccineLimits{
    VaccineLimits() : minPrevDoses( 0 ), maxCumDoses( numeric_limits<uint32_t>::max() ) {}
    void set( const scnXml::DeploymentBase& );
    uint32_t minPrevDoses, maxCumDoses;
};

/** Essentially just an integer, used as a vector index.
 * 
 * This class wraps the integer to guard against unintended conversions to or
 * from an integer (which could be mis-use). */
struct ComponentId{
    explicit inline ComponentId( size_t id ) : id(id) {}
    explicit inline ComponentId( istream& stream ){ id & stream; }
    inline void operator& (istream& stream) { id & stream; }
    inline void operator& (ostream& stream) const{ id & stream; }
    inline bool operator== (const ComponentId that) const{ return id == that.id; }
    inline bool operator< (const ComponentId that) const{ return id < that.id; }
    
    /// Special value indicating the whole population
    static ComponentId wholePop() {
        return ComponentId{ numeric_limits<uint32_t>::max() };
    }
    
    size_t id;
};

/** A description of one component of a human intervention.
 * 
 * Note that one component can have several "actions", but that deployment and
 * decay of these "actions" is usually related.
 * 
 * This is a base class. */
class HumanInterventionComponent {
public:
    virtual ~HumanInterventionComponent() {}
    
    /** Deploy the component to a pre-selected human.
     * 
     * This may be called between updates (usual intervention deployment time)
     * or from case management. It should use time sim::nowOrTs1().
     * 
     * @param human Individual receiving the intervention
     * @param method Channel of deployment (mass, continuous)
     */
    virtual void deploy( Host::Human& human, mon::Deploy::Method method,
        VaccineLimits vaccLimits ) const =0;
    
    /** Get the component identifier. */
    inline ComponentId id()const{ return m_id; }
    
    /** Get the duration. */
    inline SimTime duration()const{ return m_duration; }
    
    /** Returns the appropriate descriptor from the Component::Type enum.
     * 
     * This is only used a small number of times during setup, so doesn't need
     * to be fast. */
    virtual Component::Type componentType() const=0;
    
    virtual void print_details( std::ostream& out )const =0;
    
    // Only for use by InterventionManager:
    inline void setExpireAfter( SimTime duration ){
        m_duration = duration;
    }
    
protected:
    /** Construct (from a derived class).
     * 
     * @param id Component identifier; used as an identifier for cumulative
     *  deployment as well as to match human-specific components to general
     *  parameters (i.e. objects of the class extending this one).
     */
    explicit HumanInterventionComponent(ComponentId id) : m_id(id) {}
    
private:
    /** Don't copy (this may be possible but shouldn't be needed). */
    HumanInterventionComponent( const HumanInterventionComponent& ) = delete;
    HumanInterventionComponent& operator= ( const HumanInterventionComponent& ) = delete;
    
    ComponentId m_id;
    SimTime m_duration = sim::never();
};

/** A description of a human intervention (as a list of components). */
class HumanIntervention {
public:
    //NOTE: it would be preferable to use named constructors rather than rely
    // on selection of the correct overloaded constructor here, but doing so
    // would result in an extra copy without C++11's move semantics.
    
    /** Create from a list of XML elements: list of intervention components
     * with unconditional deployment (for triggered deployments). **/
    explicit HumanIntervention(
        const xsd::cxx::tree::sequence<scnXml::Component>& componentList );
    /** Create from a list of XML elements: list of intervention components
     * and list of conditions of deployment. **/
    explicit HumanIntervention(
        const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
        const xsd::cxx::tree::sequence<scnXml::Condition>& conditionList );
    /** Create from a list of XML elements: list of intervention components
     * with unconditional deployment (for triggered deployments). **/
    explicit HumanIntervention(
        const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList );
    
    /** Deploy all components to a pre-selected human. */
    void deploy( Host::Human& human, mon::Deploy::Method method,
        VaccineLimits vaccLimits ) const;
    
    void print_details( std::ostream& out )const;
    
protected:
    // List of pointers to components. Does not manage memory (InterventionManager::humanComponents does that).
    vector<const HumanInterventionComponent*> components;
    // List of conditions of deployment. All must be satisfied to deploy.
    vector<size_t> conditions;
};

//TODO: this class is likely more complicated than it needs to be (e.g. SubList
// could just store ComponentId instead of inheriting HumanIntervention)
class TriggeredDeployments {
public:
    TriggeredDeployments( const scnXml::TriggeredDeployments& elt );
    
    void deploy( Host::Human& human, mon::Deploy::Method method,
            VaccineLimits vaccLimits )const;
    
private:
    struct SubList : protected HumanIntervention{
        SubList( const scnXml::TriggeredDeployments::DeployType& elt );
        
        void deploy( Host::Human& human, mon::Deploy::Method method,
            VaccineLimits vaccLimits )const;
        
        // Deployment restrictions:
        SimTime minAge, maxAge;
        double coverage;
    };
    
    vector<SubList> lists;
};

} }
#endif

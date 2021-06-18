/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "KLKOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

KLKOdeSystem::KLKOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(5)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<KLKOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Corneodesmosome concentration for this cell
     * 1 - KLK enzyme concentration for this cell
     * 2 - LEKTI concentration for this cell
     * 3 - Substrate Complex concentration for this cell
     * 4 - Inhibitor complex concentration for this cell
     */

    // These will be overwritten by the initial conditions defined in the Srn Model
    SetDefaultInitialCondition(0, 1.0); 
    SetDefaultInitialCondition(1, 0.0); 
    SetDefaultInitialCondition(2, 0.0); 
    SetDefaultInitialCondition(3, 0.0); 
    SetDefaultInitialCondition(4, 1.0); 

    // Add an element for each required parameter
    this->mParameters.push_back(0.0); // Add element for z location
    this->mParameters.push_back(0.0); // Add element for substrate concentration (s0)
    this->mParameters.push_back(0.0); // Add element for enzyme concentration (eT)
    this->mParameters.push_back(0.0); // Add element for inhibitor concentration (iT)

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

KLKOdeSystem::~KLKOdeSystem()
{
}

double KLKOdeSystem::pH(double z)
{
    // Not valid for z < 0
    assert(z>=0.0);
    // When z gets above 1, we keep it constant
    z = std::min(z,1.0);
    // pH fit from the Ohman and Vahlquist (1994) paper (doi: 10. 2340/0001555574375379)
    double pH_val = 6.8482 - 0.3765*z - 5.1663*z*z + 3.1792*z*z*z;
    return(pH_val);
}

double KLKOdeSystem::k_p1(double pH_val)
{
    // k_{+1} is calculated from the Caubet et al. (2004) paper (doi: 10.1074/jbc.M607567200)
    // This is currently a constant, but we leave the option open for this to change in the future
    double kp1 = 1.49e8;
    return(kp1);
}

double KLKOdeSystem::k_m1(double pH_val)
{
    // k_{-1} is calculated from the Caubet et al. (2004) paper (doi: 10.1074/jbc.M607567200)
    // This is currently a constant, but we leave the option open for this to change in the future
    double km1 = 0.0;
    return(km1);
}

double KLKOdeSystem::k_2(double pH_val)
{
    // k_2 is calculated from the Caubet et al. (2004) paper (doi: 10.1074/jbc.M607567200)
    // This is currently a constant, but we leave the option open for this to change in the future
    double k2 = 6.87e3;
    return(k2);
}

double KLKOdeSystem::k_p3(double pH_val)
{
    // Formula fit from data in Deraison (2007) et al. paper (doi: 10.1091/mbc.e07-02-0124)
    double kp3 = (15.6*pH_val-58.5)*1.0e7;
    return(kp3); 
}

double KLKOdeSystem::k_m3(double pH_val)
{
    // Formula fit from data in Deraison (2007) et al. paper (doi: 10.1091/mbc.e07-02-0124)
    double km3 = 6.90e6 * std::exp(-3.0*pH_val);
    return(km3); 
}


void KLKOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // 0: s/desmosome, 1: e/protease, 2: i/inhibitor, 3: cs/substrate complex, 4: ci/inhibitor complex
    double s = rY[0]; // s = substrate = desmosome
    double e = rY[1]; // e = enzyme = protease
    double i = rY[2]; // i = inhibitor = LEKTI
    double cs = rY[3]; // cs 
    double ci = rY[4]; // ci
    double z_loc = this->GetParameter("z_location");

    // Check all the values
    assert(s > -1.0e-3 && s < (1.0+1.0e-3));
    assert(e > -1.0e-3);
    assert(i > -1.0e-3);
    assert(cs > -1.0e-3);
    assert(ci > -1.0e-3);
    assert((e+cs+ci) < (1.0+1.0e-3));

    // If the z location is negative, we are not in the corneum, so nothing happens
    if ( z_loc < 0.0)
    {
	    rDY = std::vector<double>(rDY.size(),0.0);
        return;
    }

    // Calculate the pH and rate parameters
    double pH_val = pH(z_loc);
    double kp1 = k_p1(pH_val);
    double km1 = k_m1(pH_val);
    double k2 = k_2(pH_val);
    double kp3 = k_p3(pH_val);
    double km3 = k_m3(pH_val);

    // Get the reactant concentrations
    double eT = this->GetParameter("eT");    
    double s0 = this->GetParameter("s0");
    double iT = this->GetParameter("iT");
    assert(s0 > 0);
    assert(i < (iT/s0+1.0e-3));

    // The next lines define the ODE system 
    // 0: s/desmosome, 1: e/protease, 2: i/inhibitor, 3: cs/substrate complex, 4: ci/inhibitor complex
    rDY[0] = -kp1*e*eT*s + km1*eT*cs/s0; // dsdt
    rDY[1] = -kp1*e*s*s0+(k2+km1)*cs-kp3*s0*e*i+km3*ci; //dedt
    rDY[2] = -kp3*e*i*eT + km3*ci*eT/s0; //didt
    rDY[3] = kp1*e*s*s0-(k2+km1)*cs; // dcsdt
    rDY[4] = kp3*e*i*s0-km3*ci; //dcidt
}

template<>
void CellwiseOdeSystemInformation<KLKOdeSystem>::Initialise()
{
    // 0: s/desmosome, 1: e/protease, 2: i/inhibitor, 3: cs/substrate complex, 4: ci/inhibitor complex
    // Note: initial conditions will be overwritten in initialisation of ODE system
    this->mVariableNames.push_back("CND_substrate");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); 

    this->mVariableNames.push_back("klk_enzyme");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("LEKTI_inhibitor");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("substrate_complex");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);   
    
    this->mVariableNames.push_back("inhibitor_complex");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);   

    // Parameters
    this->mParameterNames.push_back("z_location");
    this->mParameterUnits.push_back("CD");
    this->mParameterNames.push_back("s0");
    this->mParameterUnits.push_back("M");
    this->mParameterNames.push_back("eT");
    this->mParameterUnits.push_back("M");
    this->mParameterNames.push_back("iT");
    this->mParameterUnits.push_back("M");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(KLKOdeSystem)

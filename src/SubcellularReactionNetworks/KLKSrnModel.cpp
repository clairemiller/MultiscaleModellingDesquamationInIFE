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

#include "KLKSrnModel.hpp"

KLKSrnModel::KLKSrnModel(double s0, double eT, double iT, boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(5, pOdeSolver), ms0(s0), meT(eT), miT(iT)
{
    assert(s0 > 0.0);
    assert(eT > 0.0); 
    assert(iT >= 0.0);

    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<KLKSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<KLKSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());

}

KLKSrnModel::KLKSrnModel(const KLKSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    ms0 = rModel.ms0;
    meT = rModel.meT;
    miT = rModel.miT;
    
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new KLKOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));

    // Set the reactant concentrations in the ODE system
    assert(mpOdeSystem);
    mpOdeSystem->SetParameter("s0", ms0);
    mpOdeSystem->SetParameter("eT", meT);
    mpOdeSystem->SetParameter("iT", miT);
}

AbstractSrnModel* KLKSrnModel::CreateSrnModel()
{
    return new KLKSrnModel(*this);
}

double KLKSrnModel::GetInitialSubstrateConcentration() const
{
    return(ms0);
}

double KLKSrnModel::GetTotalEnzymeConcentration() const
{
    return(meT);
}

double KLKSrnModel::GetTotalInhibitorConcentration() const
{
    return(miT);
}


void KLKSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateZLocation();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();

    // Update the cell data
    UpdateAdhesiveProteinLevel();
    UpdateKLKLevel();
    UpdateLEKTILevel();
}

void KLKSrnModel::Initialise()
{
    // Calculate the initial conditions before initialising the ODE system
    double e_fs = std::max( (meT-miT)/meT, 0.0 ); // The enzyme not in complex with inhibitor
    double i0 = std::max( (miT-meT)/ms0, 0.0 );
    // The modified initial condition for enzyme for numerical reasons
    double k2 = 6.87e3;
    double kp1 = 1.49e8;
    double e0 = k2*e_fs/(kp1*ms0+k2);

    assert(e0 >= 0.0 && e0 <= e_fs);
    assert(e_fs <= (1.0+1.0e-8));

    // 0: s/desmosome, 1: e/protease, 2: i/inhibitor, 3: cs/substrate complex, 4: ci/inhibitor complex
    std::vector<double> initial_conditions;
    initial_conditions.push_back(1.0); // s
    initial_conditions.push_back(e0); // e
    initial_conditions.push_back(i0); // i
    initial_conditions.push_back(e_fs-e0); // cs
    initial_conditions.push_back(1.0-e_fs); // ci
    SetInitialConditions(initial_conditions);

    // Initialise the ODE system, which will also set the initial conditions
    AbstractOdeSrnModel::Initialise(new KLKOdeSystem);

    // Set the reactant concentrations
    assert(mpOdeSystem);
    mpOdeSystem->SetParameter("s0", ms0);
    mpOdeSystem->SetParameter("eT", meT);
    mpOdeSystem->SetParameter("iT", miT);
    
    // Initialise the cell with cell data
    UpdateAdhesiveProteinLevel();
    UpdateKLKLevel();
    UpdateLEKTILevel();
}

void KLKSrnModel::UpdateZLocation()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double z_loc = mpCell->GetCellData()->GetItem("ZLocation");
    mpOdeSystem->SetParameter("z_location", z_loc);
}

void KLKSrnModel::UpdateAdhesiveProteinLevel()
{
    // Update the cell data with the new desmosome value
    double desmosome = GetCNDSubstrate();
    mpCell->GetCellData()->SetItem("AdhesiveProteinLevel",desmosome);
}

void KLKSrnModel::UpdateKLKLevel()
{
    // Update the cell data with the amount of klk (for output)
    double klk = GetKLKEnzyme();
    mpCell->GetCellData()->SetItem("KLKLevel", klk);
}

void KLKSrnModel::UpdateLEKTILevel()
{
    // Update the cell data with the amount of LEKTI (for output)
    double lekti = GetLEKTIInhibitor();
    mpCell->GetCellData()->SetItem("LEKTILevel",lekti);
}

double KLKSrnModel::GetCNDSubstrate()
{
    assert(mpOdeSystem != nullptr);
    double desmosome = mpOdeSystem->rGetStateVariables()[0];
    return desmosome;
}

double KLKSrnModel::GetKLKEnzyme()
{
    assert(mpOdeSystem != nullptr);
    double conc_klk = mpOdeSystem->rGetStateVariables()[1];
    return(conc_klk);
}

double KLKSrnModel::GetLEKTIInhibitor()
{
    assert(mpOdeSystem != nullptr);
    double con_lekti = mpOdeSystem->rGetStateVariables()[2];
    return(con_lekti);
}

void KLKSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(KLKSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(KLKSrnModel)

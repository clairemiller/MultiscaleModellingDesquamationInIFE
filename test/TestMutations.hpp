#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>

// To run in parallel
// #include "PetscSetupAndFinalize.hpp"
// When run in serial
#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

// Core chaste classes
#include "CellBasedEventHandler.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"

// My classes
#include "TopOfTissueTrackingModifier.hpp"
#include "TopOfTissueForce.hpp"
#include "CellAgeAtDeathWriter.hpp"
#include "DetachedCellKillerWithWriter.hpp"
#include "TissueHeightWriter.hpp"

// Fill file
#include "FillTissueFunctionsProject3.hpp"



class TestMutations : public AbstractCellBasedTestSuite
{
public:
    void TestDifferentMutationCounts() throw(Exception)
    {
        for (unsigned i = 0; i < 5; i++)
        {
            unsigned mut_count = 25*i;
            RunSimulationWithMutations(mut_count);
        }        
    }

    void RunSimulationWithMutations(unsigned mut_count) throw(Exception)
    {
        assert(mut_count >= 0);
        assert(mut_count <= 100);

        double fill_length = 15.0*24.0;
        std::stringstream fill_dir_stream;
        fill_dir_stream << "Project3WithAdaptiveHeight/TestWithInitialRestriction/Cured/NDiseased";

        fill_dir_stream << ZeroFill("",mut_count,3);
        fill_dir_stream << "/Seed";

        RunFillTissue(fill_dir_stream.str(),fill_length,mut_count,0.0);

        // Load the simulation results
        unsigned seed = GetSeed();
        std::string input_dir = ZeroFill(fill_dir_stream.str(),seed);
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_dir, fill_length);

        // Set up output folder
        std::stringstream output_folder_stream;
        output_folder_stream << "Project3WithAdaptiveHeight/Cured/NDiseased";
	    output_folder_stream << ZeroFill("",mut_count,3);
	    output_folder_stream << "/Seed";
        std::string output_folder = ZeroFill(output_folder_stream.str(),seed);
        p_simulator->SetOutputDirectory(output_folder);

        // Set up the reporting
        CellBasedEventHandler::Reset();

        // The timing for the sim   
        double sim_length = 60.0*24.0;
        unsigned output_frequency = 120*24;
        p_simulator->SetSamplingTimestepMultiple(output_frequency);
        p_simulator->SetEndTime(fill_length + sim_length);

        // Add the top of tissue tracking
        MAKE_PTR(TopOfTissueTrackingModifier<3>, p_top_modifier);
        p_simulator->AddSimulationModifier(p_top_modifier);

        // Add the force at the top
        c_vector<double,3> force = zero_vector<double>(3);
        force[2] = 0.5;
        MAKE_PTR_ARGS(TopOfTissueForce<3>, p_top_force, (force,60.0,1.5));
        p_simulator->AddForce(p_top_force);

        // Get a reference to the cell population 
        NodeBasedCellPopulation<3>& r_population = dynamic_cast<NodeBasedCellPopulation<3>& >(p_simulator->rGetCellPopulation());

        // Add the cell height and death point writers
        r_population.AddCellPopulationCountWriter<TissueHeightWriter>();
        boost::shared_ptr<CellAgeAtDeathWriter<3,3> > p_cell_writer(new CellAgeAtDeathWriter<3,3>());
        r_population.AddPopulationWriter(p_cell_writer);

        // Add the cell killer
        MAKE_PTR_ARGS(DetachedCellKillerWithWriter<3>, p_detached_killer, (&r_population,0.7,p_cell_writer));
        p_simulator->AddCellKiller(p_detached_killer);

        // Run solver
        p_simulator->Solve();

        // Reporting
        PRINT_VARIABLE(output_folder_stream.str());
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        // Save the simulation at this time point
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);

        // Delete pointer and reset the singletons
        delete p_simulator;
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Destroy();
    }
};

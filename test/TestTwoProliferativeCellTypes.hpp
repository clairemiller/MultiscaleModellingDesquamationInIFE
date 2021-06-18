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



class TestTwoProliferativeCellTypes : public AbstractCellBasedTestSuite
{
public:
    void TestArithmeticMean() throw(Exception)
    {
        unsigned i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-index");
        assert(i<5);

        std::stringstream output_folder_stream;
        output_folder_stream << "ArithmeticMean/CycleA";
        
        double base_cycle_length = 15.0;
        double cycle_diff = ((double) i) + 1.0;
        double cycleA = base_cycle_length-cycle_diff;
        double cycleB = base_cycle_length + cycle_diff;
        output_folder_stream << (unsigned) cycleA << "hr_CycleB";
        output_folder_stream << (unsigned) cycleB << "hr";

        RunSimulationWithTwoProliferativePopulations(cycleA,cycleB, output_folder_stream.str());
    }

    void TestHarmonicMean() throw(Exception)
    {
        unsigned i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-index");
        assert(i<5);

        double base_cycle_length = 15.0;
        double deltaT = ((double) i) + 1.0;
        double cycleA = 0.5*(base_cycle_length-deltaT) + 0.5*std::sqrt(std::pow(deltaT-base_cycle_length,2.0) + 2*base_cycle_length*deltaT);
        double cycleB = cycleA + deltaT;

        std::stringstream output_folder_stream;
        output_folder_stream << "HarmonicMean/";
        output_folder_stream << "DeltaT" << (unsigned) deltaT << "hr_";
        output_folder_stream << "CycleA" << (unsigned) (cycleA*100) << "e-2hr_";
        output_folder_stream << "CycleB" << (unsigned) (cycleB*100) << "e-2hr/";

        RunSimulationWithTwoProliferativePopulations(cycleA,cycleB, output_folder_stream.str());
    }


    void RunSimulationWithTwoProliferativePopulations(double cycleA,double cycleB, std::string folder_ext)
    {
        // Run the fill simulation first
        //-----------------------------------------------------------------------------------------

        // Generate the output folder for the fill
        std::stringstream fill_dir_stream;
        fill_dir_stream << "Project3WithAdaptiveHeight/TestWithInitialRestriction/" << folder_ext;
        fill_dir_stream << "/Seed";

        double fill_length = 15.0*24.0;

        // Get the simulation info from the function and reseed
        unsigned seed = GetSeed();
        std::string fill_dir = ZeroFill(fill_dir_stream.str(),seed);
        RandomNumberGenerator::Instance()->Reseed(seed);

        // Set up the reporting
        CellBasedEventHandler::Reset();

        // Set up the nodes in standard setup
        unsigned width = 10;
        const unsigned n_cells = width*width;
        std::vector<Node<3>*> nodes(n_cells);
        for ( unsigned i=0; i < width; i++ )
        {
            for (unsigned j = 0; j < width; j++)
            {
                unsigned id = j*width + i;
                nodes[id] = new Node<3>( id,false, (double)(i),(double)(j),0.0);
            }
        }

        // Construct the mesh
        std::vector<double> periodic_widths(2, width);
        PeriodicNdNodesOnlyMesh<3> mesh(periodic_widths,true,true,false);
        mesh.ConstructNodesWithoutMesh(nodes,2.0);

    // Create a vector of cell cycles
    RandomNumberGenerator* rand_gen = RandomNumberGenerator::Instance();
    std::vector<double> cycle_vec(n_cells,cycleA);
    unsigned cycleB_count = 0;
    for (unsigned i = 0; i < (n_cells/2); i++)
    {
        unsigned i_cycleB = (unsigned) (rand_gen->ranf()*n_cells);
        if ( cycle_vec[i_cycleB] == cycleA ) {
            cycle_vec[i_cycleB] = cycleB;
            cycleB_count++;
        }
        else {
            // This cell has already been allocated as cycleB, need to choose a different cell
            i--;
        }
    }
    assert(cycleB_count == 50);

    // Create the cells
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    for (unsigned i=0; i<(n_cells); i++)
    {
        // Create a wild type state
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        // Set up the cell cycle model
        // Note: preset durations are S=5, G2=4, and M=1. 
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel;
        p_cell_cycle_model->SetDimension(3);
        p_cell_cycle_model->SetMaxTransitGenerations(0);
        // We need to do (cycle_length-11) as the uniform cell cycle model add [0,4] to each cell cycle
        // which means the cell cycle is really G1+S+G2+M+2 (S+G2+M+2=5+3+1+2=11) in length
        assert(cycle_vec[i]>11.0);
        p_cell_cycle_model->SetSDuration(5.0);
        p_cell_cycle_model->SetG2Duration(3.0);
        p_cell_cycle_model->SetStemCellG1Duration(cycle_vec[i]-11.0);
        // Set up the srn model
        // Input: s0, eT, iT
        KLKSrnModel* p_srn_model = new KLKSrnModel(10.0e-6,0.1e-9,0.1e-9);

        // Create the cell and assign stem cell type and birth time
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
        p_cell->SetBirthTime(birth_time);

        // Add to the cell vector
        cells.push_back(p_cell);
    }

    // Create cell population
    NodeBasedCellPopulation<3> cell_population(mesh, cells);
    cell_population.SetAbsoluteMovementThreshold(1.5);

    // Set up simulator
    OffLatticeSimulation<3> simulator(cell_population);
    simulator.SetOutputDirectory(fill_dir);
    simulator.SetSamplingTimestepMultiple(24.0*120.0);
    simulator.SetEndTime(fill_length);

    // Add the required modifier for the SRN model
    // (input: start height and expected height)
    MAKE_PTR_ARGS(KLKOdeModelCellModifierWithAdaptiveHeight<3>, p_klk_modifier, (4.0));
    simulator.AddSimulationModifier(p_klk_modifier);

    // Add the output of the SRN model details
    boost::shared_ptr<CellDataItemWriter<3,3> > p_loc_writer(new CellDataItemWriter<3,3>("ZLocation"));
    cell_population.AddCellWriter(p_loc_writer);
    boost::shared_ptr<CellDataItemWriter<3,3> > p_adh_writer(new CellDataItemWriter<3,3>("AdhesiveProteinLevel"));
    cell_population.AddCellWriter(p_adh_writer);
    boost::shared_ptr<CellDataItemWriter<3,3> > p_klk_writer(new CellDataItemWriter<3,3>("KLKLevel"));
    cell_population.AddCellWriter(p_klk_writer);
    boost::shared_ptr<CellDataItemWriter<3,3> > p_lekti_writer(new CellDataItemWriter<3,3>("LEKTILevel"));
    cell_population.AddCellWriter(p_lekti_writer);

    // Add adhesive and repulsive forces
    MAKE_PTR_ARGS(KLKDrivenPalssonAdhesionForce<3>, p_adh_force, (1.0));
    p_adh_force->SetMeinekeDivisionRestingSpringLength(0.001);
    simulator.AddForce(p_adh_force);
    MAKE_PTR( RepulsionForce<3>, p_rep_force);
    p_rep_force->SetMeinekeDivisionRestingSpringLength(0.001);
    simulator.AddForce(p_rep_force);

    // Add the sloughing at the top
    c_vector<double,3> pt = zero_vector<double>(3);
    c_vector<double,3> nml = zero_vector<double>(3);
    pt[2] = 10.0;
    nml[2] = 1.0;
    MAKE_PTR_ARGS(PlaneBasedCellKiller<3>,p_killer,(&cell_population,pt,nml));
    simulator.AddCellKiller(p_killer);

    // Add the bottom boundary
    MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>, p_base_bc, (&cell_population));
    p_base_bc->SetUseJiggledBottomCells(true);
    simulator.AddCellPopulationBoundaryCondition(p_base_bc);
    MAKE_PTR(UndulatingBaseMembraneAdhesionForce<3>, p_base_force);
    simulator.AddForce(p_base_force);

    // Add the division direction and rotational force
    c_vector<double,3> div_vec = zero_vector<double>(3);
    div_vec[2] = 0.001;
    MAKE_PTR_ARGS(RotationalDivisionForce<3>,p_rot_force, (10.0));
    simulator.AddForce(p_rot_force);
    MAKE_PTR_ARGS(FixedDirectionCentreBasedDivisionRule<3>,p_div_rule,(div_vec));
    cell_population.SetCentreBasedDivisionRule(p_div_rule);

    // Add the modifier to change stem cells to differentiated at a certain height
    MAKE_PTR_ARGS(HeightDependentDifferentiationCellModifier<3>, p_diff_modifier, (2.0));
    simulator.AddSimulationModifier(p_diff_modifier);

    // Add the output of the mutation states, velocity, and divisions
    cell_population.AddCellWriter<CellMutationStatesWriter>();
    cell_population.AddPopulationWriter<NodeVelocityWriter>();
    cell_population.AddCellWriter<CellAgesWriter>();
    simulator.SetOutputDivisionLocations(true);

    // Add the pinning
    MAKE_PTR_ARGS(VerticallyFixedStemCellBoundaryCondition<3>, p_sc_bc, (&cell_population));
    simulator.AddCellPopulationBoundaryCondition(p_sc_bc);

    // Run solver
    simulator.Solve();

    // Remove the sloughing 
    simulator.RemoveAllCellKillers();

    // Output run time data
    PRINT_VARIABLE(seed);
    PRINT_VARIABLE(cycleA);
    PRINT_VARIABLE(cycleB);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);



    // Then run the actual simulation
    //-----------------------------------------------------------------------------------------

    // Set up output folder
    std::stringstream output_folder_stream;
    output_folder_stream << "Project3WithAdaptiveHeight/TwoProliferativePopulations/" << folder_ext;
	output_folder_stream << "/Seed";
    std::string output_folder = ZeroFill(output_folder_stream.str(),seed);
    simulator.SetOutputDirectory(output_folder);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The timing for the simulation  
    double sim_length = 60.0*24.0;
    unsigned output_frequency = 120*24;
    simulator.SetSamplingTimestepMultiple(output_frequency);
    simulator.SetEndTime(fill_length + sim_length);

    // Add the top of tissue tracking
    MAKE_PTR(TopOfTissueTrackingModifier<3>, p_top_modifier);
    simulator.AddSimulationModifier(p_top_modifier);

    // Add the force at the top
    c_vector<double,3> force = zero_vector<double>(3);
    force[2] = 0.5;
    MAKE_PTR_ARGS(TopOfTissueForce<3>, p_top_force, (force,60.0,1.5));
    simulator.AddForce(p_top_force);

    // Add the cell height and death point writers
    cell_population.AddCellPopulationCountWriter<TissueHeightWriter>();
    boost::shared_ptr<CellAgeAtDeathWriter<3,3> > p_cell_writer(new CellAgeAtDeathWriter<3,3>());
    cell_population.AddPopulationWriter(p_cell_writer);

    // Add the cell killer
    MAKE_PTR_ARGS(DetachedCellKillerWithWriter<3>, p_detached_killer, (&cell_population,0.7,p_cell_writer));
    simulator.AddCellKiller(p_detached_killer);

    // Run solver
    simulator.Solve();

    // Reporting
    PRINT_VARIABLE(output_folder_stream.str());
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the simulation at this time point
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);

    // Delete pointer and reset the singletons
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);
    RandomNumberGenerator::Destroy();
}
};

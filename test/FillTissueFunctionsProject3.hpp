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
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "RepulsionForce.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CellDataItemWriter.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellAgesWriter.hpp"

// My classes
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "KLKSrnModel.hpp"
#include "KLKDrivenPalssonAdhesionForce.hpp"
#include "KLKOdeModelCellModifierWithAdaptiveHeight.hpp"
#include "FlatBaseMembraneBoundaryCondition.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"
#include "RotationalDivisionForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "HeightDependentDifferentiationCellModifier.hpp"
#include "VerticallyFixedStemCellBoundaryCondition.hpp"

/* Functions used for directory naming/access in the fill and subsequent simulations
*-------------------------------------------------------------------------------------
*/

// Zero fill
std::string ZeroFill(std::string dir, unsigned seed, unsigned n_digits = 2)
{
    std::stringstream dir_stream;
    dir_stream << dir;
    for ( unsigned d = (n_digits-1); d > 0; d-- )
    {
        if ( seed/std::pow(10,d) < 1 )
        {
            dir_stream << "0";
        }
    }
    dir_stream << seed;
    return (dir_stream.str());
}

// Get the seed
unsigned GetSeed()
{
    unsigned seed = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed");
    return seed;
}

/* The fill simulation function
*-------------------------------------------------------------------------------------
*/
void RunFillTissue(std::string output_dir, double fill_length, unsigned n_mutated=0, double iT_mut=0.0) throw(Exception)
{
    // Get the simulation info from the function and reseed
    unsigned seed = GetSeed();
    output_dir = ZeroFill(output_dir,seed);
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

    // Create a vector of iT
    double iT_healthy = 0.1e-9;
    RandomNumberGenerator* rand_gen = RandomNumberGenerator::Instance();
    std::vector<double> iT_vec(n_cells,iT_healthy);
    unsigned mut_count = 0;
    if (n_mutated == 100) {
        iT_vec = std::vector<double>(n_cells,iT_mut);
    }
    else {
        for (unsigned i = 0; i < n_mutated; i++)
        {
            unsigned i_mut = (unsigned) (rand_gen->ranf()*100.0);
            if ( iT_vec[i_mut] == iT_healthy ) {
                iT_vec[i_mut] = iT_mut;
                mut_count++;
            }
            else {
                // This cell has already been allocated as mutated, need to choose a different cell
                i--;
            }
        }
        assert(mut_count == n_mutated);
    }

    // Create the cells
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    for (unsigned i=0; i<(n_cells); i++)
    {
        // Create a wild type state
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        // Set up the cell cycle model
        // Note: durations are S=5, G2=4, and M=1. So G1=5 gives a total of 15 hours
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel;
        p_cell_cycle_model->SetDimension(3);
        p_cell_cycle_model->SetMaxTransitGenerations(0);
        p_cell_cycle_model->SetStemCellG1Duration(3); 
        
        // get the amount of enzyme
        double iT = iT_vec[i];
        // Set up the srn model
        // Input: s0, eT, iT
        KLKSrnModel* p_srn_model = new KLKSrnModel(10.0e-6,0.1e-9,iT);

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
    simulator.SetOutputDirectory(output_dir);
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
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);
}

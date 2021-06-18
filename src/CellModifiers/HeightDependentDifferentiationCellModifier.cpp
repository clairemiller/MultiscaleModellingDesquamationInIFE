
#include "HeightDependentDifferentiationCellModifier.hpp"

template<unsigned DIM>
// Default constructor
HeightDependentDifferentiationCellModifier<DIM>::HeightDependentDifferentiationCellModifier(double maxHeight)
: mMaxHeight(maxHeight)
{}

template<unsigned DIM>
void HeightDependentDifferentiationCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void HeightDependentDifferentiationCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}
	
// Override UpdatAtEndOfTimeStep method
template<unsigned DIM>
void HeightDependentDifferentiationCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();

	// Get the simulation time
	//double simTime = SimulationTime :: Instance()->GetTime();

	// If the cell is below a certain height, it becomes a differentiated
	MAKE_PTR(DifferentiatedCellProliferativeType,p_diff_type);
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{
		bool diff_cell = (*cell_iter)->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>();
		if (!diff_cell)
		{
			c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
			double cell_height = cell_loc[DIM-1];
			if (cell_height > mMaxHeight)
			{
				(*cell_iter)->SetCellProliferativeType(p_diff_type);
			}
		}
	}
}

// Get member variables functions
template<unsigned DIM>
double HeightDependentDifferentiationCellModifier<DIM>::GetMaxHeight() const
{
	return mMaxHeight;
}

// Set member functions
template<unsigned DIM>
void HeightDependentDifferentiationCellModifier<DIM>::SetMaxHeight(double height)
{
	mMaxHeight = height;
}

// Override output method
template<unsigned DIM>
void HeightDependentDifferentiationCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<MaxHeight>" << mMaxHeight << "</MaxHeight>\n";
	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template class HeightDependentDifferentiationCellModifier<1>;
template class HeightDependentDifferentiationCellModifier<2>;
template class HeightDependentDifferentiationCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HeightDependentDifferentiationCellModifier)


#include "TopOfTissueForce.hpp"

// Default constructor
template<unsigned DIM>
TopOfTissueForce<DIM>::TopOfTissueForce(c_vector<double,DIM> appliedForce, double minAge, double minHeight, unsigned vertical)
 		: AbstractForce<DIM,DIM>(), mAppliedForce(appliedForce), mMinHeight(minHeight), mVertical(vertical), mMinAge(minAge)
{
	assert(mVertical < DIM);
} 

// Overridden AddForceContribution method
template<unsigned DIM>
void TopOfTissueForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	// Check that it is a node based cell population
	NodeBasedCellPopulation<DIM>* pNodeBasedPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
	if ( pNodeBasedPopulation == NULL )
	{
		EXCEPTION("Currently only valid for node based population.\n");
	}

    // Loop over cell population and determine the force
	for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter!= rCellPopulation.End(); ++cell_iter )
	{
		// Check if this cell qualifies for the force
		if ( (*cell_iter)->GetMutationState()->template IsType<TopOfTissueCellMutationState>() )
		{
			if (rCellPopulation.GetLocationOfCellCentre(*cell_iter)[mVertical] > mMinHeight && (*cell_iter)->GetAge() > mMinAge)
			{
				// Find the index of the cell
				unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
				// Add force
				rCellPopulation.GetNode(index)->AddAppliedForceContribution(mAppliedForce);
			}
		}
	}
}

template<unsigned DIM>
void TopOfTissueForce<DIM>::SetAppliedForce(c_vector<double, DIM> appliedForce)
{
	mAppliedForce = appliedForce;
}

template<unsigned DIM>
c_vector<double,DIM> TopOfTissueForce<DIM>::GetAppliedForce() const
{
	return mAppliedForce;
}

template<unsigned DIM>
double TopOfTissueForce<DIM>::GetMinimumHeight() const
{
	return mMinHeight;
}

template<unsigned DIM>
double TopOfTissueForce<DIM>::GetMinimumAge() const
{
	return mMinAge;
}


// Overridden OutputForceParameters method
template<unsigned DIM>
void TopOfTissueForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	for ( unsigned i = 0; i < DIM; i++ )
	{
		*rParamsFile << "\t\t\t<AppliedForceDim" << i << ">" << mAppliedForce[i] << "</AppliedForceDim" << i << ">\n";
	}
	*rParamsFile << "\t\t\t<MinimumHeight>" << mMinHeight << "</MinimumHeight>\n";
	*rParamsFile << "\t\t\t<MinimumAge>" << mMinAge << "</MinimumAge>\n";
	AbstractForce<DIM,DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class TopOfTissueForce<1>;
template class TopOfTissueForce<2>;
template class TopOfTissueForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TopOfTissueForce)
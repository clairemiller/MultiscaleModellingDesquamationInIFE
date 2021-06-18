
#include "UndulatingBaseMembraneAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
UndulatingBaseMembraneAdhesionForce<DIM>::UndulatingBaseMembraneAdhesionForce(unsigned verticalDirection)
		 : AbstractForce<DIM,DIM>(), mVert(verticalDirection), mInteractionDistance(1.5), mAlphaProlifCell(50.0), mAlphaDiffCell(0.0)
{
	// Alpha value of 50 (500 micro-m) comes from; Li et al. 2013, Skin Stem Cell Hypotheses and Lon Term Clone Survival
	if ( DIM == 1 )
	{
		EXCEPTION("Undulating membrane force non-sensical for 1 dimension");
	}
}

// Function to extract the membrane point from CellData
template<unsigned DIM>
c_vector<double,DIM> UndulatingBaseMembraneAdhesionForce<DIM>::GetVectorFromCellData(CellPtr p_cell)
{
	c_vector<double,DIM> vec;
	switch(DIM) {
		case 1: NEVER_REACHED;
		case 3:
			vec[2] = p_cell->GetCellData()->GetItem("attachment_point_z");
		case 2:
			vec[1] = p_cell->GetCellData()->GetItem("attachment_point_y");
			vec[0] = p_cell->GetCellData()->GetItem("attachment_point_x");
	}
	return vec;
}


// Overridden AddForceContribution method
template<unsigned DIM>
void UndulatingBaseMembraneAdhesionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
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
        
		// Find the index of the cell
		unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		
		// If it is on the membrane or below, it will be fixed in the boundary condition
		// If it is above the membrane, it experiences adhesion back downwards
		
		if ( (*cell_iter)->GetMutationState()->template IsType<AttachedCellMutationState>() )
		{
			// Get the cell location
        	c_vector<double,DIM> c_loc = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

			// Get the appropriate membrane location
			c_vector<double,DIM> attachment_point = GetVectorFromCellData(*cell_iter);

			// Check it is above the membrane (if not, dealt with in boundary condition)
			TS_ASSERT_DELTA(attachment_point[mVert],0.0,1e-6); // Remove once sure this is working
			if ( c_loc[mVert] > attachment_point[mVert] )
			{
				// Calculate the distance
				c_vector<double,DIM> force_dir_vec = rCellPopulation.rGetMesh().GetVectorFromAtoB(c_loc,attachment_point);
				double scalar_dist = norm_2(force_dir_vec);
				// If it is less than the interaction distance, it is attracted to its attachment point
				if ( scalar_dist < mInteractionDistance )
				{
					force_dir_vec = force_dir_vec / scalar_dist;
					// Equation taken from: Li. et al. 2013. DOI: 10.1038/srep01904
					double lambda = 7.0;
					double c1 = 1.0 / sqrt(2.0*lambda);
					double bij = scalar_dist*2.0;
					double scalar_force = -1.0 * ( (bij+c1)*exp(-lambda*pow(bij+c1,2.0)) - 
						c1*exp(-lambda*(bij*bij + c1*c1)) );
					// Check for factor if 
					boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = (*cell_iter)->GetCellProliferativeType();
					if ( p_cell_type->IsType<StemCellProliferativeType>() || p_cell_type->IsType<TransitCellProliferativeType>() )
					{
						scalar_force = mAlphaProlifCell*scalar_force;
						// Assert with a small tolerance
						assert(scalar_force > (-1.0e-9));
					}
					else
					{
						scalar_force = mAlphaDiffCell*scalar_force;
						assert( p_cell_type->IsType<DifferentiatedCellProliferativeType>() );
					}
					// Convert back to vector force
					c_vector<double,DIM> force = force_dir_vec * scalar_force;
					// Assert with a small tolerance
					assert(force[mVert] < 1.0e-9);
					// Add force
					rCellPopulation.GetNode(index)->AddAppliedForceContribution(force); 

					// Add the force magnitude to cell data
					(*cell_iter)->GetCellData()->SetItem("membrane_force_magnitude",scalar_force);
				}
				// Otherwise remove the attachment point
				else
				{
					(*cell_iter)->SetMutationState(rCellPopulation.GetCellPropertyRegistry()->template Get<WildTypeCellMutationState>());
					(*cell_iter)->GetCellData()->SetItem("membrane_force_magnitude",-1.0);
				}
			}
		}
	}
}

template<unsigned DIM>
double UndulatingBaseMembraneAdhesionForce<DIM>::GetAlphaProlifCell() const
{
	return mAlphaProlifCell;

}

template<unsigned DIM>
double UndulatingBaseMembraneAdhesionForce<DIM>::GetAlphaDiffCell() const
{
	return mAlphaDiffCell;
}

template<unsigned DIM>
unsigned UndulatingBaseMembraneAdhesionForce<DIM>::GetVerticalDirection() const
{
	return mVert;
}

template<unsigned DIM>
void UndulatingBaseMembraneAdhesionForce<DIM>::SetAdhesionParameters(double alpha_prolif, double alpha_diff)
{
	mAlphaProlifCell = alpha_prolif;
	mAlphaDiffCell = alpha_diff;
	assert(alpha_prolif >= 0.0 && alpha_diff >= 0.0);
}


// Overridden OutputForceParameters method
template<unsigned DIM>
void UndulatingBaseMembraneAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<VerticalDirection>" << mVert << "</VerticalDirection>\n";
	*rParamsFile << "\t\t\t<InteractionDistance>" << mInteractionDistance << "</InteractionDistance>\n";
	*rParamsFile << "\t\t\t<AlphaProlifCell>" << mAlphaProlifCell << "</AlphaProlifCell>\n";
	*rParamsFile << "\t\t\t<AlphaDiffCell>" << mAlphaDiffCell << "</AlphaDiffCell>\n";

	AbstractForce<DIM,DIM>::OutputForceParameters(rParamsFile);
}

// Overridden WriteDataToVisualisedSetupFile method
template<unsigned DIM>
void UndulatingBaseMembraneAdhesionForce<DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
	AbstractForce<DIM,DIM>::WriteDataToVisualizerSetupFile(pVizSetupFile);
}

// Explicit instantiation
template class UndulatingBaseMembraneAdhesionForce<1>;
template class UndulatingBaseMembraneAdhesionForce<2>;
template class UndulatingBaseMembraneAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(UndulatingBaseMembraneAdhesionForce)

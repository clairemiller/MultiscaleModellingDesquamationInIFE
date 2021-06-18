
#include "AbstractUndulatingBaseMembraneBoundaryCondition.hpp"

// Default constructor
template<unsigned DIM>
AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::AbstractUndulatingBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection)
		 : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation), mVert(verticalDirection), mUseJiggledBottomCells(false), mIsInitialised(false)
{
	assert(verticalDirection < DIM);

	if ( DIM == 1 )
	{
		EXCEPTION("Undulating membrane boundary condition non-sensical for 1 dimension");
	}
}

template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::InitialisePopulation()
{
	// Check the cell data is not stored
	typename AbstractCellPopulation<DIM>::Iterator iter_cell_0 = this->mpCellPopulation->Begin();
	std::vector<std::string> v_keys = iter_cell_0->GetCellData()->GetKeys();
	if (std::find(v_keys.begin(), v_keys.end(), "attachment_point_x") != v_keys.end())
	{
		EXCEPTION("Base membrane boundary condition initialisation boolean is false, but attachment point found in cell data.");
	}

	// Initialise cell vectors to max double
	c_vector<double,DIM> vec_maxdbl = scalar_vector<double>(DIM, DBL_MAX);
	for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin(); cell_iter!= this->mpCellPopulation->End(); ++cell_iter )
	{
		AddVectorToCellData(*cell_iter,vec_maxdbl);
	}

	mIsInitialised = true;
}

// Overridden ImposeBoundaryCondition() method
template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
	// Check it is initialised
	if (!mIsInitialised)
	{
		InitialisePopulation();
	}

	// Check that it is a node based cell population
	NodeBasedCellPopulation<DIM>* pNodeBasedPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);
	if ( pNodeBasedPopulation == NULL )
	{
		EXCEPTION("Currently only valid for node based population.\n");
	}

    // Loop over cell population and determine the force
	for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = pNodeBasedPopulation->Begin(); cell_iter!= pNodeBasedPopulation->End(); ++cell_iter )
	{
        // Get the cell location
        c_vector<double,DIM> c_loc = pNodeBasedPopulation->GetLocationOfCellCentre(*cell_iter);

        // Calculate the membrane location vertically above/below
		double f_c_loc = BaseShapeFunction(c_loc);
		
		// If the point is on the membrane, assign a new attachment point
		if ( c_loc[mVert] == f_c_loc )
		{
			// Assign the attachment point and mutation
			(*cell_iter)->SetMutationState(pNodeBasedPopulation->GetCellPropertyRegistry()->template Get<AttachedCellMutationState>());
			AddVectorToCellData(*cell_iter,c_loc);
		}
		// If the point is below, move to the membrane and assign new attachment point
		else if ( c_loc[mVert] < f_c_loc )
        {
			// Find the index of the cell
			unsigned index = pNodeBasedPopulation->GetLocationIndexUsingCell(*cell_iter);
			// Calculate the appropriate membrane location
			c_vector<double,DIM> attachment_point = DetermineClosestMembraneLocation(c_loc);
			// Move the force back to the attachment point
			pNodeBasedPopulation->GetNode(index)->rGetModifiableLocation() = attachment_point;
			// Add a 'jiggle' if desired (in the direction the cell is pushed back)
			if (mUseJiggledBottomCells)
			{
				double jiggle_magnitude = 0.05*RandomNumberGenerator::Instance()->ranf();
				(pNodeBasedPopulation->GetNode(index)->rGetModifiableLocation()[DIM-1]) += jiggle_magnitude;
			}
			// Assign the attachment point and mutation
			(*cell_iter)->SetMutationState(pNodeBasedPopulation->GetCellPropertyRegistry()->template Get<AttachedCellMutationState>());
			AddVectorToCellData(*cell_iter,attachment_point);
		}
	}
}
   

// Function for determing the force direction
template<unsigned DIM>
c_vector<double,DIM> AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::DetermineClosestMembraneLocation(c_vector<double,DIM> p)
{
	/* 
	 * Let p be a point with coordinates (x0,y0,z0), and z1 be value of z if p were raised vertically to the membrane
	 * Let pm be the closest point to p on the membrane, with coordinates xm, ym, and zm
	 * We firstly assume that we can approximate the surface between (x0,y0,z1) and pm with the tangent surface at (x0,y0,z1)
	 * This means we can approximate the normal at pm with the normal, n, at (x0,y0,z1). Therefore n = (-dz/dx,-dz/dy,1)
	 * Then we know that the closest point on the membrane is given by: pm = p + k*n
	 * Therefore: zm = z0 + k, ym = y0 - k*dz/dy, and xm = x0 - k*dz/dy
	 * The tangent surface at (x0,y0,z1), and therefore pm, is given by: z-z1 = (dz/dx)*(x-x0) + (dz/dy)*(y-y0)
	 * Substituting we find that: k = (z1-z0) / ( (dz/dx)^2 + (dz/dy)^2 + 1 ), giving us xm and ym, then we calculate zm
	 * In 2D this simplifies to k = (y1-y0)/( (dy/dx)^2 + 1 ) and therefore xm = x0 + (dy/dx)*(y0-y1) / ( (dy/dx)^2 + 1 )
	*/

	// Get the derivatives and points
	c_vector<double,DIM> dz = CalculateDerivativesAtPoint(p);
	double z1 = BaseShapeFunction(p);

	c_vector<double,DIM> pm = zero_vector<double>(DIM);
	
	// Calculate the parameter k (noting that the denominator is the square of the norm of the derivatives)
	double norm_dz = norm_2(dz);
	double k =  ( z1 - p[mVert] ) / (norm_dz*norm_dz);
	// Now calculate the pm coordinates for the horizontal directions
	for ( unsigned i = 0; i < DIM; i++ )
	{
		// xm = x0 - k*(dz/dx)
		if (i != mVert)
		{
			pm[i] = p[i] - k*dz[i]; 
		}
	}

	// Then add in the z value on the membrane
	pm[mVert] = BaseShapeFunction(pm);
	
	return(pm);
}

// Function to add the membrane point to CellData
template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::AddVectorToCellData(CellPtr p_cell, c_vector<double,DIM> vec)
{
	switch(DIM) {
		case 1:
			NEVER_REACHED;
		case 3:
			p_cell->GetCellData()->SetItem("attachment_point_z",vec[2]);
		case 2:
			p_cell->GetCellData()->SetItem("attachment_point_y",vec[1]);
			p_cell->GetCellData()->SetItem("attachment_point_x",vec[0]);
	}
}

template<unsigned DIM>
unsigned AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::GetVerticalDirection() const
{
	return mVert;
}

template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::SetUseJiggledBottomCells(bool useJiggledBottomCells)
{
	mUseJiggledBottomCells = useJiggledBottomCells;
}

template<unsigned DIM>
bool AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::GetUseJiggledBottomCells() const
{
	return(mUseJiggledBottomCells);
}

template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::SetIsInitialised(bool isInitialised)
{
	mIsInitialised = isInitialised;
}


// Overridden OutputCellPopulationBoundaryConditionParameters method
template<unsigned DIM>
void AbstractUndulatingBaseMembraneBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<VerticalDirection>" << mVert << "</VerticalDirection>\n";
	*rParamsFile << "\t\t\t<UseJiggledBottomCells>" << mUseJiggledBottomCells << "</UseJiggledBottomCells>\n";

	AbstractCellPopulationBoundaryCondition<DIM,DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractUndulatingBaseMembraneBoundaryCondition<1>;
template class AbstractUndulatingBaseMembraneBoundaryCondition<2>;
template class AbstractUndulatingBaseMembraneBoundaryCondition<3>;

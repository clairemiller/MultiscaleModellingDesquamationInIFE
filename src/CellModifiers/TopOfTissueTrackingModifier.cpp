/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "TopOfTissueTrackingModifier.hpp"

template<unsigned DIM>
TopOfTissueTrackingModifier<DIM>::TopOfTissueTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(), mGridSize(1.0)
{
}

template<unsigned DIM>
TopOfTissueTrackingModifier<DIM>::~TopOfTissueTrackingModifier()
{
}

template<unsigned DIM>
void TopOfTissueTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TopOfTissueTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TopOfTissueTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Get the domain
    c_vector<unsigned, DIM-1> domainSize = scalar_vector<double>(DIM-1,0);
    for ( unsigned i = 0; i < DIM-1; i++ )
    {
        domainSize[i] = (unsigned) std::max( std::ceil(rCellPopulation.GetWidth(i)/mGridSize), 1.0 );
    }

    // Check we are a node based cell populatoin
    NodeBasedCellPopulation<DIM>* pPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if ( ! (bool) pPopulation )
    {
        EXCEPTION("Neighbour tracking modifier only implemented for node based cell populations.\n");
    }

    // Iterate over cell population and store the highest locations
    // Top cells contains the cell index as the key, then the value is the pair of cell id and height
    typename std::map<std::vector<unsigned>, std::pair<unsigned, double> > top_cells;
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
        double cell_height = cell_loc[DIM-1];
        cell_loc/=mGridSize;
        std::vector<unsigned> cell_index(cell_loc.begin(), cell_loc.end()-1);

        assert(cell_index.size()==DIM-1);
        typename std::map<std::vector<unsigned>, std::pair<unsigned, double> >::iterator it_topcell;
        it_topcell = top_cells.find(cell_index);
        if ( it_topcell != top_cells.end() )
        {
            double current_height = (it_topcell->second).second;
            if (cell_height > current_height)
            {
                unsigned cell_id = rCellPopulation.GetLocationIndexUsingCell(*cell_it);
                (it_topcell->second) = std::make_pair(cell_id,cell_height);
            }
        }
        else
        {
            unsigned cell_id = rCellPopulation.GetLocationIndexUsingCell(*cell_it);
            top_cells[cell_index] = std::make_pair(cell_id,cell_height);
        }

        // Remove any previous assignment of the force
        if ( (*cell_it)->GetMutationState()->template IsType<TopOfTissueCellMutationState>() )
        {
            (*cell_it)->SetMutationState(rCellPopulation.GetCellPropertyRegistry()->template Get<WildTypeCellMutationState>());
        }
    }

    // Now iterate over the highest cells and assign mutation state
    typename std::map<std::vector<unsigned>, std::pair<unsigned, double> >::iterator it_topcell = top_cells.begin();
    for (; it_topcell != top_cells.end(); it_topcell++)
    {
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex( (it_topcell->second).first );
        // Check it is not an attached cell
        if ( !(p_cell->GetMutationState()->IsType<AttachedCellMutationState>()) )
        {
            p_cell->SetMutationState(rCellPopulation.GetCellPropertyRegistry()->template Get<TopOfTissueCellMutationState>());
        }
    }
}

template<unsigned DIM>
double TopOfTissueTrackingModifier<DIM>::GetGridSize() const
{
    return(mGridSize);
}

template<unsigned DIM>
void TopOfTissueTrackingModifier<DIM>::SetGridSize(double gridSize)
{
    assert(gridSize > 0.0);
    mGridSize = gridSize;
}

template<unsigned DIM>
void TopOfTissueTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GridSize>" << mGridSize << "</GridSize>\n";
    // Now call parent method
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class TopOfTissueTrackingModifier<1>;
template class TopOfTissueTrackingModifier<2>;
template class TopOfTissueTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TopOfTissueTrackingModifier)


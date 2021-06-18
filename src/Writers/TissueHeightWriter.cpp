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

#include "TissueHeightWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::TissueHeightWriter()
    : AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>("tissueheight.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Order: vector[0] is wild type, vector[1] is the top of tissue type
    c_vector<unsigned,2> n_cells = scalar_vector<unsigned>(2,0);
    c_vector<double,2> max_height = scalar_vector<double>(2,0.0);
    c_vector<double,2> mean_height = scalar_vector<double>(2,0.0);

    // Iterate over cell population
    for (typename NodeBasedCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
            cell_iter != pCellPopulation->End();
            ++cell_iter)
    {
        double height = (pCellPopulation->GetLocationOfCellCentre(*cell_iter))[SPACE_DIM-1];
        boost::shared_ptr<AbstractCellMutationState> cell_mut_state = (*cell_iter)->GetMutationState();
        unsigned cell_type;
        if ( cell_mut_state->IsType<WildTypeCellMutationState>() )
        {
            cell_type = 0;
        }
        else if ( cell_mut_state->IsType<TopOfTissueCellMutationState>() )
        {
            cell_type = 1;
        }
        else
        {
            continue;
        }

        n_cells[cell_type]++;
        if ( height > max_height[cell_type] )
        {
            max_height[cell_type] = height;
        }
        mean_height[cell_type] += height;
    }

    // Calculate the mean
    mean_height[0] = mean_height[0]/n_cells[0];
    mean_height[1] = mean_height[1]/n_cells[1];

    // Output
    MAKE_PTR(WildTypeCellMutationState, wild_state);
    MAKE_PTR(TopOfTissueCellMutationState, tot_state);
    int wild_colour = (wild_state->GetColour());
    int tot_colour = (tot_state->GetColour());
    *this->mpOutStream <<  wild_colour;
    *this->mpOutStream << "\t" << mean_height[0] << "\t" << max_height[0] << "\t";
    *this->mpOutStream << tot_colour;
    *this->mpOutStream << "\t" << mean_height[1] << "\t" << max_height[1];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("Tissue height writer only valid for node based cell populations.");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("Tissue height writer only valid for node based cell populations.");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("Tissue height writer only valid for node based cell populations.");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TissueHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("Tissue height writer only valid for node based cell populations.");
}

// Explicit instantiation
template class TissueHeightWriter<1,1>;
template class TissueHeightWriter<1,2>;
template class TissueHeightWriter<2,2>;
template class TissueHeightWriter<1,3>;
template class TissueHeightWriter<2,3>;
template class TissueHeightWriter<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TissueHeightWriter)


/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "KLKOdeModelCellModifierWithAdaptiveHeight.hpp"

template<unsigned DIM>
KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::KLKOdeModelCellModifierWithAdaptiveHeight(double heightStartSC, unsigned vertical)
 : KLKOdeModelCellModifier<DIM>(heightStartSC,heightStartSC,vertical), 
    mSamplingHeightTimesteps(120), mLengthMovingAverage(120)
{
    // Note: expectedSCHeight parameter will be overwritten when setup solve is run
    // We just give it a dummy (non-zero) value for now

    // We start with an empty previous heights vector
    mPreviousHeights = std::vector<double>();
    // Calculate the length of the moving average from the timestep info
    mMovingAverageTimesteps = mLengthMovingAverage*mSamplingHeightTimesteps;
}

template<unsigned DIM>
KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::~KLKOdeModelCellModifierWithAdaptiveHeight()
{
}

template<unsigned DIM>
double KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::CalculateTissueHeight(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Iterate over cell population to calculate the sum of the heights of the top of tissue cells
    double sumHeights = 0.0;
    unsigned nTopCells = 0;
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter )
    {
        boost::shared_ptr<AbstractCellMutationState> cell_mut_state = (*cell_iter)->GetMutationState();
        if ( cell_mut_state->IsType<TopOfTissueCellMutationState>() )
        {
            double height = (rCellPopulation.GetLocationOfCellCentre(*cell_iter))[this->mVertical];
            sumHeights += height;
            nTopCells ++;
        }
    }
    // Calculate the mean of the tissue height and use to approximate the current SC height
    double meanHeight = sumHeights/( (double)nTopCells );
    // Check for case where there are no top cells, 
    // we set tissue height to zero so it will be dealt with in the UpdateCellData code
    if (nTopCells == 0)
    {
        meanHeight = 0.0;
    }
    return (meanHeight - (this->mHeightStartSC));
    //this->mExpectedSCHeight = meanHeight - (this->mHeightStartSC);
}

template<unsigned DIM>
void KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    unsigned currTimestep = SimulationTime::Instance() -> GetTimeStepsElapsed();

    // If we are in a sampling time step, calculate the current height
    if (currTimestep % mSamplingHeightTimesteps == 0)
    {
        // Calculate the current height
        double currHeight = CalculateTissueHeight(rCellPopulation);

        // Case where there is yet to be any cells in the SC, 
        // setting height=0 breaks the calculation of Z to give nan so we just set to 1.0 and don't add to the vector
        // Note, if the height is 0.0, then there will be no cells in the system so it shouldn't matter
        if (currHeight < 1.0e-6) 
        {
            (this->mExpectedSCHeight) = 1.0;
        }
        // Need to adjust for the first mLengthMovingAverage data points
        else if ( mPreviousHeights.size() < mLengthMovingAverage )
        {
            // Eg. formula: (x0 + x1 + x2)/3 = ((x0 + x1)/2)*2/3 + x2/3
            double currLength = (double) mPreviousHeights.size();
            // When vector is empty, should be (0 + currHeight)/1
            (this->mExpectedSCHeight) = ( (this->mExpectedSCHeight)*currLength + currHeight)/(currLength + 1.0);
            // Now add to the vector
            mPreviousHeights.push_back(currHeight);
        }
        // Otherwise, we swap out the values in the average and the vector
        else
        {
            // Determine the location in the vector
            unsigned vecIdx = (currTimestep % mMovingAverageTimesteps) / mSamplingHeightTimesteps;
            assert(vecIdx < mLengthMovingAverage);
            // Update the expected height by subtracting the previous value at this location and adding the new value
            (this->mExpectedSCHeight) = (this->mExpectedSCHeight) + (currHeight - mPreviousHeights[vecIdx])/( (double) mLengthMovingAverage );
            // Now replace the old value with the new height
            mPreviousHeights[vecIdx] = currHeight;
        }
    }

    // Run the parent class method to assign the z values to the cells
    KLKOdeModelCellModifier<DIM>::UpdateCellData(rCellPopulation);
}


template<unsigned DIM>
void KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::SetHeightSamplingDeltaTimesteps(unsigned samplingHeightTimesteps)
{
    // Check value is non-zero
    assert(samplingHeightTimesteps > 0);

    // Assign member variables
    mSamplingHeightTimesteps = samplingHeightTimesteps;
    mMovingAverageTimesteps = mLengthMovingAverage*samplingHeightTimesteps;

    // Do not need to change the length of the previous height vector
}

template<unsigned DIM>
unsigned KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::GetHeightSamplingDeltaTimesteps() const
{
    return mSamplingHeightTimesteps;
}

template<unsigned DIM>
void KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::SetNumberSamplesForMovingAverage(unsigned lengthMovingAverage)
{
    // Check value is non-zero
    assert(lengthMovingAverage > 0);

    // Assign member variables
    mLengthMovingAverage = lengthMovingAverage;
    mMovingAverageTimesteps = lengthMovingAverage*mSamplingHeightTimesteps;

    // Resize the previous heights vector if it is too long, using the current expected height
    if ( mPreviousHeights.size() > mLengthMovingAverage)
    {
        mPreviousHeights = std::vector<double>(mLengthMovingAverage,this->mExpectedSCHeight);
    }
}


template<unsigned DIM>
unsigned KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::GetNumberSamplesForMovingAverage() const
{
    return mLengthMovingAverage;
}


template<unsigned DIM>
void KLKOdeModelCellModifierWithAdaptiveHeight<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call parent class method
	KLKOdeModelCellModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class KLKOdeModelCellModifierWithAdaptiveHeight<1>;
template class KLKOdeModelCellModifierWithAdaptiveHeight<2>;
template class KLKOdeModelCellModifierWithAdaptiveHeight<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKOdeModelCellModifierWithAdaptiveHeight)

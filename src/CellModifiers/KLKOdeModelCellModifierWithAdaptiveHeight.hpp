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

#ifndef KLKODEMODELCELLMODIFIERWITHADAPTIVEHEIGHT
#define KLKODEMODELCELLMODIFIERWITHADAPTIVEHEIGHT

#include "ChasteSerialization.hpp"

#include "AbstractCellPopulation.hpp"
#include "KLKOdeModelCellModifier.hpp"
#include "TopOfTissueCellMutationState.hpp"

template<unsigned  DIM>
class KLKOdeModelCellModifierWithAdaptiveHeight : public KLKOdeModelCellModifier<DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
      archive & boost::serialization::base_object<KLKOdeModelCellModifier<DIM> >(*this);
      archive & mSamplingHeightTimesteps;
      archive & mMovingAverageTimesteps;
      archive & mLengthMovingAverage;
    }

    /** Parameters */
    unsigned mSamplingHeightTimesteps; // Number of time steps between each height sample
    unsigned mMovingAverageTimesteps; // Length of time to average the height over
    std::vector<double> mPreviousHeights; // Store the current heights included in the moving average
    unsigned mLengthMovingAverage; // The length of mPreviousHeights

public:

    /**
     * Default constructor.
     */
    KLKOdeModelCellModifierWithAdaptiveHeight(double heightStartSC,unsigned vertical=DIM-1);

    /**
     * Destructor.
     */
    virtual ~KLKOdeModelCellModifierWithAdaptiveHeight();

    /** Method to calculate the tissue height
     * @return the current height
     */
    double CalculateTissueHeight(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /** Override UpdateCellData */
    virtual void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /** Method to change the number of time steps between height sampling
     * @param samplingHeightTimesteps the new number of time steps
     */
    void SetHeightSamplingDeltaTimesteps(unsigned samplingHeightTimesteps);

    /** Get method for the sampling frequency
     * @return the value of mSamplingHeightTimesteps
     */
    unsigned GetHeightSamplingDeltaTimesteps() const;

    /** Method to change the number of samples in the moving average
     * @param 
     */
    void SetNumberSamplesForMovingAverage(unsigned lengthMovingAverage);

    /** Method to get the number of samples in the moving average
     * @return the value of mLengthMovingAverage
     */
    unsigned GetNumberSamplesForMovingAverage() const;
    /**
     * Output any simulation modifier parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKOdeModelCellModifierWithAdaptiveHeight)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const KLKOdeModelCellModifierWithAdaptiveHeight<DIM> * t, const unsigned int file_version)
		{
        double heightStartSC = t->GetHeightStartSC();
        unsigned vertical = t->GetVertical();
        unsigned sampling = t->GetHeightSamplingDeltaTimesteps();
        unsigned movavglength = t->GetNumberSamplesForMovingAverage();

        ar << heightStartSC;
        ar << vertical;
        ar << sampling;
        ar << movavglength;
    }

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, KLKOdeModelCellModifierWithAdaptiveHeight<DIM> * t, const unsigned int file_version)
		{
      double heightStartSC;
      unsigned vertical,sampling,movavglength;

      ar >> heightStartSC;
      ar >> vertical;
      ar >> sampling;
      ar >> movavglength;

			::new(t)KLKOdeModelCellModifierWithAdaptiveHeight<DIM>(heightStartSC, vertical);
      t->SetHeightSamplingDeltaTimesteps(sampling);
      t->SetNumberSamplesForMovingAverage(movavglength);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*KLKODEMODELCELLMODIFIERWITHADAPTIVEHEIGHT*/

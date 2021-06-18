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

#ifndef ROTATIONALDIVISIONFORCE_HPP_
#define ROTATIONALDIVISIONFORCE_HPP_

#include "AbstractForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SimulationTime.hpp"

/**
 * A class to apply a rotational force to two cells involved in the mitotic division phase from a single cell
 */
template<unsigned  DIM>
class RotationalDivisionForce : public AbstractForce<DIM, DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM,DIM> >(*this);
        archive & mTorsionCoefficient;
        archive & mNormalVector;
        archive & mGrowthDuration;
        archive & mPropSpringLength;
    }

    double DotProduct(c_vector<double,DIM> u, c_vector<double,DIM> v) const;

protected:

    // The scaling parameter - in radian^{-1}
    double mTorsionCoefficient;

    // The desired direction
    c_vector<double,DIM> mNormalVector;

    // The force duration
    double mGrowthDuration;

    // The vertical direction
    unsigned mVertical;

    // Whether the rotational spring force is proportional to spring length
    bool mPropSpringLength;

public:

    /**
     * Constructor.
     */
    RotationalDivisionForce(double torsionCoefficient, double growthDuration=1.0);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Set method for normal vector
     * @param the new value of normal vector
     */
    void SetNormalVector(c_vector<double, DIM> normalVec);

    /**
     * Set method for proportional to spring length
     * @param the new value of the boolean
     */
    void SetProportionalToSpringLength(bool propSpringLength);

    /**
     * Get method for torsion coefficient
     * 
     * @return the value of mTorsionCoefficient
     */
    virtual double GetTorsionCoefficient(double age = -1.0) const;

    /**
     * Get method for division vector
     * 
     * @return the value of mTorsionCoefficient
     */
    virtual c_vector<double,DIM> GetNormalVector(c_vector<double,DIM> p) const;

    /**
     * Get the distance from a point to the plane
     * 
     * @param the point of interest
     * @return the distance to the plane
     */
    virtual double GetDistanceToBasePlane(c_vector<double,DIM> p) const;

     /**
     * Get method for spring growth duration
     * 
     * @return the value of mGrowthDuration
     */
    double GetGrowthDuration() const;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Overridden WriteDataToVisualizerSetupFile() method.
     * Write any data necessary to a visualization setup file.
     * Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
     * 
     * @param pVizSetupFile a visualization setup file
     */
    virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RotationalDivisionForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const RotationalDivisionForce<DIM> * t, const unsigned int file_version)
		{
            double torsion_coefficient, growth_duration;
            c_vector<double,DIM> division_vector;
            torsion_coefficient = t->GetTorsionCoefficient();
            ar << torsion_coefficient;
            division_vector = t->GetNormalVector(zero_vector<double>(DIM));
            for (unsigned i =0; i < DIM; i++)
            {
                ar << division_vector[i];
            }
            growth_duration = t->GetGrowthDuration();
            ar << growth_duration;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, RotationalDivisionForce<DIM> * t, const unsigned int file_version)
		{
            double torsion_coefficient, growth_duration;
            c_vector<double,DIM> division_vector;
            ar >> torsion_coefficient;
            for (unsigned i =0; i < DIM; i++)
            {
                ar >> division_vector[i];
            }
            ar >> growth_duration;

            ::new(t)RotationalDivisionForce<DIM>(torsion_coefficient,growth_duration);
            t->SetNormalVector(division_vector);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*ROTATIONALDIVISIONFORCE_HPP_*/

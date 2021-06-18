
#ifndef PALSSONADHESIONFORCE_HPP_
#define PALSSONADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "FakePetscSetup.hpp"


template<unsigned DIM>
class PalssonAdhesionForce : public GeneralisedLinearSpringForce<DIM,DIM>
{
protected:
	// The spring strength of the force
	double mAlpha;

	// The calculation of the force
	double CalculatePalssonAdhesionForce(double bij);

private:
	// Friend the test class
	friend class TestPalssonAdhesionForce;

	/**
	 * Degrading factor computation
	 * @param hi the height of cell i
	 * @param hj the height of cell j
	 */
	//double CalculateDegradingFactor(double hi, double hj);

	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
		archive & mAlpha;
	}

public:
	/**
	 * Default constructor
	 * Default value taken from: Li et al. 2013, Skin Stem Cell Hypotheses and Long Term Clone Survival – Explored Using Agent-based Modelling
	 * @param strength the value of mStrength
	 * @param radius the value of mRadius
	 */
	PalssonAdhesionForce(double peakForce = 0.02);

	/**
	 * Overridden CalculateForceBetweenNodes method
	 *
	 * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
	 */

	virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation);

	/**
	 * Overridden AddForceContribution method
	 * @param rCellPopulation; a reference to the cell population
	 */
	void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

	/**
	 * Get methods for the member functions
	 * @return the value of mAlpha
	 */
	virtual double GetPeakForce() const;

	/**
	 * Overridden OutputForceParameters method
	 */
	void OutputForceParameters(out_stream& rParamsFile);

	/**
	 * Override the spring constants from generalised linear springa as well for safety
	 */
	void SetMeinekeSpringStiffness(double springStiffness);
	double GetMeinekeSpringStiffness() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PalssonAdhesionForce)

#endif /*PALSSONADHESIONFORCE_HPP_*/
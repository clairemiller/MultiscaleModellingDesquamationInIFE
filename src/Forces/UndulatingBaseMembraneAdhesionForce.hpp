#ifndef UNDULATINGBASEMEMBRANEADHESIONFORCE_HPP_
#define UNDULATINGBASEMEMBRANEADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AttachedCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PeriodicNdNodesOnlyMesh.hpp"

#include "FakePetscSetup.hpp"
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

template<unsigned DIM>
class UndulatingBaseMembraneAdhesionForce : public AbstractForce<DIM,DIM>
{
private:
	// Functions for getting vectors from CellData
	c_vector<double,DIM> GetVectorFromCellData(CellPtr p_cell);

	// Add archiving functions
	friend class boost::serialization::access;
	friend class TestUndulatingBaseMembraneForce;

	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractForce<DIM,DIM> >(*this);
		archive & mVert;
		for (unsigned i=0; i < (DIM-1); i++)
		{
			archive & mHorizVec[i];
		}
		archive & mAlphaProlifCell;
		archive & mAlphaDiffCell; 
		archive & mInteractionDistance;
	}
protected:
	// Member variables
	unsigned mVert;
	c_vector<unsigned,DIM-1> mHorizVec;
	double mInteractionDistance;
	double mAlphaProlifCell;
	double mAlphaDiffCell;
public:
	/**
	 * Default constructor
	 */
	UndulatingBaseMembraneAdhesionForce(unsigned verticalDirection = DIM-1);

	/**
	 * Default deconstructor
	*/
	virtual ~UndulatingBaseMembraneAdhesionForce(){}

	/**
     * Calculates the force on each node.
     * @param rCellPopulation reference to the cell population
     */
	virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

	
	/**
	 * Get function for member variable mAlphaProlifCell
	 * @return mAlphaProlifCell
	*/
	double GetAlphaProlifCell() const;

	/**
	 * Get function for member variable mAlphaDiffCell
	 * @return mAlphaDiffCell
	*/
	double GetAlphaDiffCell() const;

	/**
	 * Get the vertical direction
	 * @return mVert
	 */
	unsigned GetVerticalDirection() const;

	 /**
     * Outputs force parameters to file.
     * @param rParamsFile the file stream to which the parameters are output
     */
	 virtual void OutputForceParameters(out_stream& rParamsFile);
	 
	/**
	* Write any data necessary to a visualization setup file.
	* Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
	* 
	* @param pVizSetupFile a visualization setup file
	*/
	virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);

	/**
	 * Set the adhesive values for the force
	 * @param alpha the value of the force coefficient alpha
	 * @param prolif_multiplier a multiplier for proliferative cell types
	 * @param diff_multiplier a multiplier for differentiated cell types
	*/
	void SetAdhesionParameters(double alpha_prolif, double alpha_diff);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(UndulatingBaseMembraneAdhesionForce)

#endif /*UNDULATINGBASEMEMBRANEADHESIONFORCE_HPP_*/

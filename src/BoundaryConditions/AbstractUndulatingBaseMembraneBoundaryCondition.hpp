#ifndef ABSTRACTUNDULATINGBASEMEMBRABOUNDARYCONDITIONCE_HPP_
#define ABSTRACTUNDULATINGBASEMEMBRABOUNDARYCONDITIONCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "AttachedCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PeriodicNdNodesOnlyMesh.hpp"

#include "FakePetscSetup.hpp"
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

template<unsigned DIM>
class AbstractUndulatingBaseMembraneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM,DIM>
{
private:
	// Functions for setting/getting vectors in CellData
	void AddVectorToCellData(CellPtr p_cell, c_vector<double,DIM> vec);
	
	// Add archiving functions
	friend class boost::serialization::access;
	friend class TestAbstractUndulatingBaseMembraneBoundaryCondition;

	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		// The cell population will be archived in save_construct_data in daughter classes
		archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM,DIM> >(*this);
		archive & mVert;
		archive & mUseJiggledBottomCells;
		archive & mIsInitialised;
	}

protected:
	// Member variables
	unsigned mVert;
	bool mUseJiggledBottomCells;
	bool mIsInitialised;
public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    AbstractUndulatingBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection = DIM-1);

	void InitialisePopulation();
	
    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

	/** 
	 * A function for the membrane shape
	 * @param p a reference to a x,y(,z) location to determine the appropriate height given a point location
	*/
	virtual double BaseShapeFunction(c_vector<double,DIM> p) = 0;
	
	/**
	 * Calculates the derivatives vector to the membrane at a given point
	 * @param the point
	 * @return vector with the derivatives in each direction
	 */
	virtual c_vector<double,DIM> CalculateDerivativesAtPoint(c_vector<double,DIM> p) = 0;	

	/** 
	 * Determines closest point on the membrane to some point p1
	 * @param the point of interest
	 * @return the closest point on membrane
	 */
	c_vector<double,DIM> DetermineClosestMembraneLocation(c_vector<double,DIM> p1);
    
    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition(){return true;}

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);

	/**
	 * Default deconstructor
	*/
	virtual ~AbstractUndulatingBaseMembraneBoundaryCondition(){}

	/**
	 * Get the vertical direction
	 * @return mVert
	 */
	unsigned GetVerticalDirection() const;

	/**
	 * Set whether to use 'jiggled' bottom cells
	 * @param useJiggledBottomCells the new value of mUseJiggledBottomCells
	 */
	void SetUseJiggledBottomCells(bool useJiggledBottomCells);

	/**
	 * Get whether the boundary condition is using 'jiggled' bottom cells
	 * @return boolean whether using or not using 'jiggled' bottom cells
	 */
	bool GetUseJiggledBottomCells() const;

	/**
	 * Set method for the initialised flag
	 * @param isInitialised the new value of mIsInitialised
	 */
	void SetIsInitialised(bool isInitialised);

};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractUndulatingBaseMembraneBoundaryCondition)

#endif /*ABSTRACTUNDULATINGBASEMEMBRANEBOUNDARYCONDITION_HPP_*/

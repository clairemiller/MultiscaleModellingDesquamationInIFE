
#ifndef DETACHEDCELLKILLER_HPP_
#define DETACHEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AttachedCellMutationState.hpp"

/**
 * Cell Killer to remove cells that have no neighbours within a certain radius
 */

template< unsigned DIM >

class DetachedCellKiller : public AbstractCellKiller<DIM>
{
private:
	/* Member variables

	 * The threshold distance for what is defined as a neighbour
	 */
	double mNeighbourRadius; 

	/* 
	 * mIterationsToSkip: The number of iterations of the code to skip 
	 * mIterationsSkipped: The current iterations skipped
	*/
	unsigned mIterationsToSkip;
	unsigned mIterationsSkipped;

    friend class boost::serialization::access;
	friend class TestDetachedCellKiller;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
			archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
	}

public:
    /** 
     * Default constructor
     * @param pCellPopulation the cell population
     */
    DetachedCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius);

    /**
     * Default deconstructor
     */
    ~DetachedCellKiller(){}

    /**
     * Overridden virtual class
    */
    void CheckAndLabelCellsForApoptosisOrDeath();

	/**
	 * Calculate the detached cells
	 */
	std::vector<unsigned> GetDetachedCells(NodeBasedCellPopulation<DIM>* pNodePopulation);

	/**
	 * Get method for the distance threshold to define neighbours
	 * @return the value of mNeighbourDistance
	 */
	double GetNeighbourRadius() const;

	/**
	 * Set method for the distance threshold to define neighbours
	 * @param neighbourRadius the new value of mNeighbourRadius
	 */
	void SetNeighbourRadius(double neighbourRadius);

	/**
	 * Get method for the iterations to skip
	 * @return the value of mIterationsToSkip
	 */
	unsigned GetIterationsToSkip() const;

	/**
	 * Set method for the number of iterations to skip
	 * @param iterationsToSkip the new value of mIterationsToSkip
	 */
	void SetIterationsToSkip(unsigned iterationsToSkip);

    /**
     * Overridden virtual class
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DetachedCellKiller)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const DetachedCellKiller<DIM> * t, const unsigned int file_version)
		{
			const AbstractCellPopulation<DIM,DIM>* const p_cell_population = t->GetCellPopulation();
		    ar << p_cell_population;
			
			double neighbour_radius = t->GetNeighbourRadius();
			ar << neighbour_radius;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, DetachedCellKiller<DIM> * t, const unsigned int file_version)
		{
			AbstractCellPopulation<DIM,DIM>* p_cell_population;
			ar >> p_cell_population;
			double neighbour_radius;
			ar >> neighbour_radius;

		    ::new(t)DetachedCellKiller<DIM>(p_cell_population, neighbour_radius);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*DETACHEDCELLKILLER_HPP_*/


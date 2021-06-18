
#ifndef DETACHEDCELLKILLERWithWriterWITHWRITER_HPP_
#define DETACHEDCELLKILLERWithWriterWITHWRITER_HPP_

#include "AbstractCellKiller.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AttachedCellMutationState.hpp"

#include "DetachedCellKiller.hpp"
#include "CellAgeAtDeathWriter.hpp"

/**
 * Cell Killer to remove cells that have no neighbours within a certain radius
 */

template< unsigned DIM >

class DetachedCellKillerWithWriter : public DetachedCellKiller<DIM>
{
private:
    friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
			archive & boost::serialization::base_object<DetachedCellKiller<DIM> >(*this);
			archive & mpCellWriter;
	}

protected:
	// Member variables: pointer to the writer
	boost::shared_ptr<CellAgeAtDeathWriter<DIM,DIM> > mpCellWriter;

public:
    /** 
     * Default constructor
     * @param pCellPopulation the cell population
     */
    DetachedCellKillerWithWriter(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius, boost::shared_ptr<CellAgeAtDeathWriter<DIM,DIM> > rCellWriter = nullptr);

    /**
     * Default deconstructor
     */
    ~DetachedCellKillerWithWriter(){}

	/**
	 * Get function for the cell writer
	 */
	const boost::shared_ptr<CellAgeAtDeathWriter<DIM,DIM> > GetCellWriter() const;

    /**
     * Overridden virtual class
    */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden virtual class
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DetachedCellKillerWithWriter)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const DetachedCellKillerWithWriter<DIM> * t, const unsigned int file_version)
		{
			const AbstractCellPopulation<DIM,DIM>* const p_cell_population = t->GetCellPopulation();
		    ar << p_cell_population;
			
			double neighbour_radius = t->GetNeighbourRadius();
			ar << neighbour_radius;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, DetachedCellKillerWithWriter<DIM> * t, const unsigned int file_version)
		{
			AbstractCellPopulation<DIM,DIM>* p_cell_population;
			ar >> p_cell_population;

			double neighbour_radius;
			ar >> neighbour_radius;

		    ::new(t)DetachedCellKillerWithWriter<DIM>(p_cell_population, neighbour_radius);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*DETACHEDCELLKILLERWithWriterWITHWRITER_HPP_*/


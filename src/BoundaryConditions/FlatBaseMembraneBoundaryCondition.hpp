#ifndef FLATBASEMEMBRANEBOUNDARYCONDITION_HPP_
#define FLATBASEMEMBRANEBOUNDARYCONDITION_HPP_

#include "AbstractUndulatingBaseMembraneBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


template<unsigned DIM>
class FlatBaseMembraneBoundaryCondition : public AbstractUndulatingBaseMembraneBoundaryCondition<DIM>
{
private:
	// Add archiving functions
	friend class boost::serialization::access;
	friend class TestUndulatingBaseMembraneBoundaryCondition;
	
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractUndulatingBaseMembraneBoundaryCondition<DIM> >(*this);
	}
public:
	/**
	 * Constructor
	 & @param verticalDirection the 'up' direction
	 */
	FlatBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection = DIM-1);

	/**
	 * Default deconstructor
	*/
	virtual ~FlatBaseMembraneBoundaryCondition();

	/** 
	 * A function for the membrane shape
	 * @param p a x,y(,z) location to determine the appropriate height given a point location
	 * @return the height of the membrane (always 0)
	*/
	virtual double BaseShapeFunction(c_vector<double,DIM> p);

	/**
	 * Calculates the unit tangent vector to the membrane at a given point
	 * @param the point
	 * @return the tangent vector
	 */
	virtual c_vector<double,DIM> CalculateDerivativesAtPoint(c_vector<double,DIM> p);	
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FlatBaseMembraneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FlatBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const FlatBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a FlatBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, FlatBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)FlatBaseMembraneBoundaryCondition<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*FLATBASEMEMBRANEBOUNDARYCONDITION_HPP_*/

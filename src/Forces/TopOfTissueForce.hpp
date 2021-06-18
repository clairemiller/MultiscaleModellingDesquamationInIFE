
#ifndef TOPOFTISSUEFORCE_HPP_
#define TOPOFTISSUEFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "TopOfTissueCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MeshBasedCellPopulation.hpp"


#include "FakePetscSetup.hpp"


template<unsigned DIM>
class TopOfTissueForce : public AbstractForce<DIM,DIM>
{
private:
	// Member variables
	c_vector<double,DIM> mAppliedForce;
	double mMinHeight;
	unsigned mVertical;
	double mMinAge;

	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractForce<DIM,DIM> >(*this);
	}

public:
	/**
	 * Default constructor
	 * @param appliedForce the force vector to apply to the low neighbour cells
	 */
	TopOfTissueForce(c_vector<double,DIM> appliedForce, double minAge, double minHeight = 0.0, unsigned vertical = DIM-1);

	/**
	 * Default deconstructor
	*/
	~TopOfTissueForce(){}

	/**
     * Calculates the force on each node.
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Outputs force parameters to file.
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

	/**
	 * Set the neighbour threshold (if wanting to change from defaults)
	 * @param neighbourThreshold
	 */
	void SetNeighbourThreshold(unsigned neighbourThreshold);

	/**
	 * Get methods for the member variables
	 * @return the value of mNeighbourThreshold
	 */
	 unsigned GetNeighbourThreshold() const;

	/**
	 * Set method for the applied force
	 * @param the value of mAppliedForce
	 */
	void SetAppliedForce(c_vector<double, DIM> appliedForce);

	/**
	 * Get method for the applied force
	 * @return the value of mAppliedForce
	 */
	c_vector<double, DIM> GetAppliedForce() const;

	/**
	 * Get method for the minimum height variable
	 * @return the value of mMinHeight
	 */
	double GetMinimumHeight() const;

	/**
	 * Get method for the minimum age variable
	 */
	double GetMinimumAge() const;
	

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TopOfTissueForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const TopOfTissueForce<DIM> * t, const unsigned int file_version)
		{
			c_vector<double,DIM> applied_force = t->GetAppliedForce();
			for ( unsigned i = 0; i < DIM; i++ )
			{
				ar << applied_force[i];
			}
			double minimum_height = t->GetMinimumHeight();
			ar << minimum_height;
			double minimum_age = t->GetMinimumAge();
			ar << minimum_age;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, TopOfTissueForce<DIM> * t, const unsigned int file_version)
		{
			c_vector<double,DIM> applied_force;
			double minimum_height, minimum_age;

			for ( unsigned i = 0; i < DIM; i++ )
			{
				ar >> applied_force[i];
			}
			ar >> minimum_height;
			ar >> minimum_age;

		    ::new(t)TopOfTissueForce<DIM>(applied_force,minimum_age,minimum_height);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*TOPOFTISSUEFORCE_HPP_*/

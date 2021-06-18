
#ifndef KLKDRIVENPALSSONADHESIONFORCE_HPP_
#define KLKDRIVENPALSSONADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractGradientPalssonAdhesionForce.hpp"


#include "FakePetscSetup.hpp"


template<unsigned DIM>
class KLKDrivenPalssonAdhesionForce : public AbstractGradientPalssonAdhesionForce<DIM>
{
private:
	friend class TestKLKDrivenPalssonAdhesionForce;
	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractGradientPalssonAdhesionForce<DIM> >(*this);
	}

public:
	/**
	 * Default constructor
	 */
	KLKDrivenPalssonAdhesionForce(double peakForce);

    /**
     * The scaling function after decay start age
     * @param pCell a pointer to the cell to get scaling function value for
     */
    virtual double ScalingFunction(CellPtr pCell);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKDrivenPalssonAdhesionForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const KLKDrivenPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double peak_force = t->GetPeakForce();
			ar << peak_force;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, KLKDrivenPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double peak_force;
			ar >> peak_force;

			::new(t)KLKDrivenPalssonAdhesionForce<DIM>(peak_force);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*KLKDRIVENPALSSONADHESIONFORCE_HPP_*/
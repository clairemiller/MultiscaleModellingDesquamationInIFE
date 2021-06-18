
#ifndef HEIGHTDEPENDENTDIFFERENTIATIONCELLMODIFIER_HPP_
#define HEIGHTDEPENDENTDIFFERENTIATIONCELLMODIFIER_HPP_

#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractSimplePhaseBasedCellCycleModel.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

template<unsigned DIM>
class HeightDependentDifferentiationCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
		
		// Archive member variables
		archive & mMaxHeight;
	}

protected:
	double mMaxHeight;

public:
	// Default constructor
	HeightDependentDifferentiationCellModifier(double maxHeight=2.0);

	// Default destructor
	~HeightDependentDifferentiationCellModifier(){}

	// Override UpdatAtEndOfTimeStep method
	void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override UpdateCellData
	void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override output method
	void OutputSimulationModifierParameters(out_stream& rParamsFile);

	// Override virtual function
	virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

	/**
	 * Get mMaxHeight
	 * @return value of mMaxHeight
	 */
	double GetMaxHeight() const;

	/**
	 * Set mMaxHeight
	 * @param Height the value of mMaxHeight
	 */
	void SetMaxHeight(double maxHeight);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HeightDependentDifferentiationCellModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const HeightDependentDifferentiationCellModifier<DIM> * t, const unsigned int file_version)
		{
		     double height = t->GetMaxHeight();
		     ar << height;		    
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, HeightDependentDifferentiationCellModifier<DIM> * t, const unsigned int file_version)
		{
            double height;
			ar >> height;
			::new(t)HeightDependentDifferentiationCellModifier<DIM>(height);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*HEIGHTDEPENDENTDIFFERENTIATIONCELLMODIFIER_HPP_*/

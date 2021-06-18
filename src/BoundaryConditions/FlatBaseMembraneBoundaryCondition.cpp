#include "FlatBaseMembraneBoundaryCondition.hpp"


template<unsigned DIM>
FlatBaseMembraneBoundaryCondition<DIM>::FlatBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection)
    : AbstractUndulatingBaseMembraneBoundaryCondition<DIM>(pCellPopulation, verticalDirection)
{
}

template<unsigned DIM>
FlatBaseMembraneBoundaryCondition<DIM>::~FlatBaseMembraneBoundaryCondition()
{
}

template<unsigned DIM>
double FlatBaseMembraneBoundaryCondition<DIM>::BaseShapeFunction(c_vector<double,DIM> p)
{
    return(0.0);
}

template<unsigned DIM>
c_vector<double,DIM> FlatBaseMembraneBoundaryCondition<DIM>::CalculateDerivativesAtPoint(c_vector<double,DIM> p)
{
    c_vector<double, DIM> deriv = zero_vector<double>(DIM);
    deriv[this->mVert] = 1.0;
    return(deriv);
}


// Explicit instantiation
template class FlatBaseMembraneBoundaryCondition<1>;
template class FlatBaseMembraneBoundaryCondition<2>;
template class FlatBaseMembraneBoundaryCondition<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FlatBaseMembraneBoundaryCondition)
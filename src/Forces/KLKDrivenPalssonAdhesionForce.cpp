
#include "KLKDrivenPalssonAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
KLKDrivenPalssonAdhesionForce<DIM>::KLKDrivenPalssonAdhesionForce(double peakForce)
 : AbstractGradientPalssonAdhesionForce<DIM>(peakForce, 0.0)
{
    assert(peakForce > 0.0);
}

template<unsigned DIM>
double KLKDrivenPalssonAdhesionForce<DIM>::ScalingFunction(CellPtr cellA)
{
    double adhesive_proportion = cellA->GetCellData()->GetItem("AdhesiveProteinLevel");
    assert(adhesive_proportion > 0.0);
        
    return(adhesive_proportion);
}

// Explicit instantiation
template class KLKDrivenPalssonAdhesionForce<1>;
template class KLKDrivenPalssonAdhesionForce<2>;
template class KLKDrivenPalssonAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKDrivenPalssonAdhesionForce)

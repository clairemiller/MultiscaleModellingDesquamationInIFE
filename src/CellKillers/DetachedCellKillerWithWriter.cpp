
#include "DetachedCellKillerWithWriter.hpp"

template<unsigned DIM>
DetachedCellKillerWithWriter<DIM>::DetachedCellKillerWithWriter(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius, boost::shared_ptr<CellAgeAtDeathWriter<DIM,DIM> > rCellWriter)
: DetachedCellKiller<DIM>(pCellPopulation, neighbourRadius), mpCellWriter(rCellWriter)
{
}

template<unsigned DIM>
const boost::shared_ptr<CellAgeAtDeathWriter<DIM,DIM> > DetachedCellKillerWithWriter<DIM>::GetCellWriter() const
{
    return mpCellWriter;
}

template<unsigned DIM>
void DetachedCellKillerWithWriter<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
   // Cast to a node population
    NodeBasedCellPopulation<DIM>* pNodePopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    // Get the list of detached cells
    std::vector<unsigned> detached_cells = this->GetDetachedCells(pNodePopulation);

    // Store the cell ids and ages, then kill them
    for (std::vector<unsigned>::iterator id_it = detached_cells.begin(); id_it != detached_cells.end(); id_it++)
    {
        CellPtr p_deadcell = pNodePopulation->GetCellUsingLocationIndex(*id_it);
        if (!(p_deadcell->HasApoptosisBegun()) && !(p_deadcell->IsDead()))
        {
            double age = p_deadcell->GetAge();
            mpCellWriter->AddCellDeath((*id_it),age);
            p_deadcell->Kill();
        }
    }
}

template<unsigned DIM>
void DetachedCellKillerWithWriter<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{

    // Call parent method
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class DetachedCellKillerWithWriter<1>;
template class DetachedCellKillerWithWriter<2>;
template class DetachedCellKillerWithWriter<3>;

// Serialization
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DetachedCellKillerWithWriter)

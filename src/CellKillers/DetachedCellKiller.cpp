
#include "DetachedCellKiller.hpp"

template<unsigned DIM>
DetachedCellKiller<DIM>::DetachedCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius)
: AbstractCellKiller<DIM>(pCellPopulation), mNeighbourRadius(neighbourRadius), mIterationsToSkip(0)
{
    // NOTE: mIterationsToSkip is deprecated
    NodeBasedCellPopulation<DIM>* pNodePopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    if ( ! (bool) pNodePopulation )
    {
        EXCEPTION("No neighbour cell killer only implemented for node based cell populations.\n");
    }
    assert(pNodePopulation->GetMechanicsCutOffLength() >= neighbourRadius);
    mIterationsSkipped = mIterationsToSkip; 
}

template<unsigned DIM>
void DetachedCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    /*
     - Start with the set of attached cells
     - Get the set of their neighbours
     - Find the difference between the current set and the neighbours
     - Find their neighbours and iterate until the difference is an empty set
     - Any nodes not in this set are detached and removed
    */

   // Cast to a node population
    NodeBasedCellPopulation<DIM>* pNodePopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    // Get the vector of detached cells
    std::vector<unsigned> detached_cells = GetDetachedCells(pNodePopulation);

    // Kill anything detached
    for (std::vector<unsigned>::iterator cell_it = detached_cells.begin(); cell_it != detached_cells.end(); cell_it++)
    {
        CellPtr p_deadcell = pNodePopulation->GetCellUsingLocationIndex(*cell_it);
        if (!(p_deadcell->HasApoptosisBegun()) && !(p_deadcell->IsDead()))
        {
            p_deadcell->Kill();
        }
    }
}

template<unsigned DIM>
std::vector<unsigned> DetachedCellKiller<DIM>::GetDetachedCells(NodeBasedCellPopulation<DIM>* pNodePopulation)
{
       // Start with the set of attached cells
    std::vector<unsigned> tissue_set;
    for (typename NodeBasedCellPopulation<DIM>::Iterator cell_iter = pNodePopulation->Begin();
    cell_iter != pNodePopulation->End(); ++cell_iter)
    {
        if ((*cell_iter)->GetMutationState()->template IsType<AttachedCellMutationState>())
        {
            const unsigned index = pNodePopulation->GetLocationIndexUsingCell(*cell_iter);
            tissue_set.push_back(index);
        }
    }
    std::sort(tissue_set.begin(), tissue_set.end());

    // Iteratively find the attached sets neighbours
    std::vector<unsigned> set_neighbours = tissue_set; // The neighbours added during the previous loop
    unsigned iterator_check = 0;
    while (!set_neighbours.empty())
    {
        std::set<unsigned> new_neighbours;
        for ( std::vector<unsigned>::iterator it = set_neighbours.begin(); it != set_neighbours.end(); it++ )
        {
            std::set<unsigned> neighbours = pNodePopulation->GetNodesWithinNeighbourhoodRadius(*it,1.5);
            std::vector<unsigned> added_indices;
            // The set additions
            std::set_difference(neighbours.begin(), neighbours.end(), // Returns neighbours in this set...
                            tissue_set.begin(), tissue_set.end(), //... that are not in this set
                            std::back_inserter(added_indices));
            // Merge with the other new editions
            new_neighbours.insert(added_indices.begin(), added_indices.end());
        }
        // Merge the new neighbours into the full set
        tissue_set.insert(tissue_set.end(),new_neighbours.begin(), new_neighbours.end());
        std::sort(tissue_set.begin(), tissue_set.end());
        set_neighbours.resize(new_neighbours.size());
        std::copy(new_neighbours.begin(), new_neighbours.end(), set_neighbours.begin());
        iterator_check++;
        assert(iterator_check < pNodePopulation->GetNumNodes());
    }

    // Finally, determine those cells outside of the tissue
    std::vector<unsigned> all_indices = (pNodePopulation->rGetMesh()).GetAllNodeIndices();
    std::sort(all_indices.begin(), all_indices.end());
    std::vector<unsigned> detached_cells;
    std::set_difference(all_indices.begin(), all_indices.end(), 
        tissue_set.begin(), tissue_set.end(), 
        std::back_inserter(detached_cells));

    return(detached_cells);
}

template<unsigned DIM>
double DetachedCellKiller<DIM>::GetNeighbourRadius() const
{
    return(mNeighbourRadius);
}

template<unsigned DIM>
void DetachedCellKiller<DIM>::SetNeighbourRadius(double neighbourRadius)
{
    mNeighbourRadius = neighbourRadius;
}

template<unsigned DIM>
void DetachedCellKiller<DIM>::SetIterationsToSkip(unsigned iterationsToSkip)
{
    mIterationsToSkip = iterationsToSkip;
}

template<unsigned DIM>
unsigned DetachedCellKiller<DIM>::GetIterationsToSkip() const
{
    return(mIterationsToSkip);
}

template<unsigned DIM>
void DetachedCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NeighbourRadius>" << mNeighbourRadius << "</NeighbourRadius>\n";

    // Call parent method
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class DetachedCellKiller<1>;
template class DetachedCellKiller<2>;
template class DetachedCellKiller<3>;

// Serialization
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DetachedCellKiller)
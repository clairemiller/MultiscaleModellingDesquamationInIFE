
#include "PalssonAdhesionForce.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "Debug.hpp"

// Default constructor
template<unsigned DIM>
PalssonAdhesionForce<DIM>::PalssonAdhesionForce(double peakForce)
 : GeneralisedLinearSpringForce<DIM>(), mAlpha(peakForce)
{
	assert(peakForce > 0.0);
}

// Force calculation function
template<unsigned DIM>
double PalssonAdhesionForce<DIM>::CalculatePalssonAdhesionForce(double bij)
{
    // If they are calculate the force
    // Equation taken from: Li. et al. 2013. DOI: 10.1038/srep01904
	double lambda = 7.0;
	double c1 = sqrt(0.5/lambda);
	double c2 = c1*exp(-lambda*c1*c1);
    bij = 2.0*bij;
	double Fij = -1.0*mAlpha * ( (bij+c1)*exp(-lambda*pow(bij+c1,2.0)) - 
					c2*exp(-lambda*bij*bij) );
    return(Fij);
}

// Overridden AddForceContribution method
template<unsigned DIM>
void PalssonAdhesionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	// Check that it is a node based cell population
	NodeBasedCellPopulation<DIM>* pNodeBasedPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
	if ( pNodeBasedPopulation == NULL )
	{
		EXCEPTION("Currently only valid for node based population.\n");
	}
	
	// Call the base member class
	GeneralisedLinearSpringForce<DIM>::AddForceContribution(rCellPopulation);
}

// Overridden CalculateForceBetweenNodes method
// NOTE: Currently only implemented for node based cell populations
template<unsigned DIM>
c_vector<double, DIM> PalssonAdhesionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
	// We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM);
        }
    }

    // Calculate the spring length
    double rest_length = ( p_node_a->GetRadius() ) + ( p_node_b->GetRadius() );
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(p_node_a->GetIndex());
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(p_node_b->GetIndex());
    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();
    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));
    /*
    * If the cells are both newly divided, then the rest length of the spring
    * connecting them grows linearly with time, until 1 hour after division.
    */
    if (ageA < (this->mMeinekeSpringGrowthDuration) && ageB < (this->mMeinekeSpringGrowthDuration))
    {
        AbstractCentreBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);
        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a small value to the normal rest length over 1 hour
            double lambda = this->mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (rest_length - lambda) * ageA/(this->mMeinekeSpringGrowthDuration);
        }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= (this->mMeinekeSpringGrowthDuration) )
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    // Subtract the spring length and check it is not an overlap
    distance_between_nodes = distance_between_nodes - rest_length;
    if ( distance_between_nodes <= 0.0 )
    {
        return zero_vector<double>(DIM);
    }

    // Calculate force
    double Fij = CalculatePalssonAdhesionForce(distance_between_nodes/(rest_length));
    // Convert to the vector and return
    return ( Fij*unit_difference );
}

// Get methods for the peak force value
template<unsigned DIM>
double PalssonAdhesionForce<DIM>::GetPeakForce() const
{
	return mAlpha;
}

// Overridden OutputForceParameters method
template<unsigned DIM>
void PalssonAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<StrengthPeakForce>" << mAlpha << "</StrengthPeakForce>\n";
	GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned DIM>
void PalssonAdhesionForce<DIM>::SetMeinekeSpringStiffness(double springStiffness)
{
    mAlpha = springStiffness;
}

template<unsigned DIM>
double PalssonAdhesionForce<DIM>::GetMeinekeSpringStiffness() const
{
    return mAlpha;
}

// Explicit instantiation
template class PalssonAdhesionForce<1>;
template class PalssonAdhesionForce<2>;
template class PalssonAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PalssonAdhesionForce)
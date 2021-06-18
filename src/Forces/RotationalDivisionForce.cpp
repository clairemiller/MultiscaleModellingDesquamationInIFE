/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "RotationalDivisionForce.hpp"

#include "Debug.hpp"

template<unsigned DIM>
RotationalDivisionForce<DIM>::RotationalDivisionForce(double torsionCoefficient, double growthDuration)
   : AbstractForce<DIM>(), mTorsionCoefficient(torsionCoefficient), mGrowthDuration(growthDuration), mPropSpringLength(false)
{
    assert(torsionCoefficient > 0.0);
    assert(growthDuration > 0.0);
    if ( DIM == 1 )
    {
        EXCEPTION("Rotational division force non-sensical in 1D.\n");
    }
    // Create division vector
    mNormalVector = zero_vector<double>(DIM);
    mNormalVector[DIM-1] = 1.0;
    // Store a vertical direction
    mVertical = DIM-1;
}

template<unsigned DIM>
c_vector<double, DIM> RotationalDivisionForce<DIM>::GetNormalVector(c_vector<double,DIM> p) const
{
    return(mNormalVector);
}

template<unsigned DIM>
double RotationalDivisionForce<DIM>::GetDistanceToBasePlane(c_vector<double,DIM> p) const
{
    double d = DotProduct(mNormalVector,p);
    return( std::abs(d) );
}

template<unsigned DIM>
double RotationalDivisionForce<DIM>::DotProduct(c_vector<double,DIM> u, c_vector<double,DIM> v) const
{
    double dp = 0.0;
    for (unsigned i = 0; i < DIM; i++)
    {
        dp += u[i]*v[i];
    }
    return(dp);
}

template<unsigned DIM>
void RotationalDivisionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Currently only valid for node based cell populatoins
    if (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        AbstractCentreBasedCellPopulation<DIM>* p_population = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation);
        std::vector< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = p_population->rGetNodePairs();

        for (typename std::vector< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
            iter != r_node_pairs.end();
            iter++)
        {
            // Check age first (This will give better speeds than checking for a marked pair)
            unsigned node_a_index = (*iter).first->GetIndex();
            CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(node_a_index);
            double age = p_cell_A->GetAge();

            // Determine if the cell pair is still in division phase
            if ( age < mGrowthDuration )
            {
                unsigned node_b_index = (*iter).second->GetIndex();
                CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(node_b_index);
                std::pair<CellPtr,CellPtr> cell_pair = p_population->CreateCellPair(p_cell_A,p_cell_B);

                // Check for newly divided cell (marked)
                if ( p_population->IsMarkedSpring(cell_pair) )
                {
                    // Determine symmetric or asymmetric (we check for differentiated rather than proliferative, easier)
                    bool cell_a_prolif = !(p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>());
                    bool cell_b_prolif = !(p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>());
                    
                    // Required output from each case
                    Node<DIM>* force_node; // The node that experiences the force
                    c_vector<double, DIM> force_direction; // The direction of the force
                    double beta; // The angle of rotation from the desired angle
                    double spring_length; // The magnitude of the separation

                    // Symmetric proliferative division 
                    // Higher cell is rotated around the lower cell
                    //-----------------------------------------
                    if ( cell_a_prolif && cell_b_prolif )
                    {
                        // Get cell locations
                        const c_vector<double,DIM>& r_a_location = (*iter).first->rGetLocation();
                        const c_vector<double,DIM>& r_b_location = (*iter).second->rGetLocation();
                        
                        // Determine closest base locations
                        double dist_a = GetDistanceToBasePlane(r_a_location);
                        double dist_b = GetDistanceToBasePlane(r_b_location);

                        // Assign centre node and rotating node
                        Node<DIM>* centre_node = (dist_a <= dist_b) ? (*iter).first : (*iter).second;
                        force_node = (dist_a <= dist_b) ? (*iter).second : (*iter).first;

                        // Calculate the separation vector between the two cells
                        c_vector<double,DIM> separation_vec = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_node->rGetLocation(), force_node->rGetLocation());
                        spring_length = norm_2(separation_vec);
                        assert(spring_length > 0.0);

                        // Determine the projection of the separation vector onto the plane
                        c_vector<double,DIM> normal_vec = GetNormalVector(centre_node->rGetLocation());
                        c_vector<double, DIM> plane_vec = separation_vec*DotProduct(normal_vec,normal_vec)-normal_vec*DotProduct(normal_vec,separation_vec);
                        plane_vec = plane_vec / norm_2(plane_vec);

                        // Determine the angle between the separation and plane vector
                        double cosbeta = DotProduct(separation_vec,plane_vec)/(norm_2(separation_vec)*norm_2(plane_vec));
                        assert(std::abs(cosbeta) < (1.0+1e-6));
                        cosbeta = std::max(std::min(cosbeta,1.0),-1.0); // To prevent acos breaking with precision issues
                        beta = acos(cosbeta); // Return value in range [0,pi]

                        // Now determine the direction of the force:
                        // (v x (p x v)) where v is the separation vector and p is the plane vector
                        force_direction = DotProduct(separation_vec,separation_vec)*plane_vec - DotProduct(separation_vec,plane_vec)*separation_vec;
                    }
                    // Asymmetric division 
                    // Differentiated cell is rotated around the proliferative cell
                    //-----------------------------------------
                    else if ( cell_a_prolif || cell_b_prolif )
                    {  
                        Node<DIM>* stem_node = cell_a_prolif ? (*iter).first : (*iter).second;
                        Node<DIM>* diff_node = cell_a_prolif ? (*iter).second : (*iter).first;

                        // Get the difference vector
                        const c_vector<double,DIM>& r_stem_location = stem_node->rGetLocation();
                        const c_vector<double,DIM>& r_diff_location = diff_node->rGetLocation();
                        c_vector<double,DIM> location_vec = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_stem_location, r_diff_location);

                        // // Also determine the distance
                        spring_length = norm_2(location_vec);
                        assert(spring_length > 0.0);

                        // Determine the normal vector and the deformation angle, beta
                        c_vector<double,DIM> normal_vec = GetNormalVector(r_stem_location);
                        double denominator = norm_2(normal_vec) * norm_2(location_vec);
                        double dp_normal_location = DotProduct(normal_vec,location_vec);
                        double cosBeta = dp_normal_location/denominator;
                        assert(std::abs(cosBeta) < (1.0+1e-6));
                        cosBeta = std::max(std::min(cosBeta,1.0),-1.0); // To prevent acos breaking with precision issues
                        beta = acos(cosBeta); // Return value in range [0,pi]

                        // Calculate the triple product to determine direction: location_vec x (normal_vec x location_vec)
                        // This is equivalent to: nv(lv.lv) - lv(lv.nv)
                        force_direction = normal_vec*DotProduct(location_vec,location_vec) - location_vec*dp_normal_location;
                        
                        // Assign correct node
                        force_node = diff_node;
                    }
                    // Symmetric proliferative division 
                    // No force is applied
                    //-----------------------------------------
                    else
                    {
                        beta = 0.0;
                        force_direction = scalar_vector<double>(DIM,1.0);
                        force_node = (*iter).first;
                    }

                    // Calculate the force
                    double scalar_force = (this->GetTorsionCoefficient(age))*beta;
                    if ( mPropSpringLength ) 
                    {
                        scalar_force = scalar_force / spring_length;
                        assert(spring_length > 0.0);
                    }
                    // Convert to vector (first normalise force direction vector)
                    double force_dir_length = norm_2(force_direction);
                    if ( force_dir_length > 0.0 )
                    {
                        force_direction = force_direction / force_dir_length;
                    }
                    c_vector<double, DIM> force = scalar_force*force_direction;

                    // Add the force to the appropriate cell
                    force_node->AddAppliedForceContribution(force);

                    // Check for unmark spring next time step
                    if (age + SimulationTime::Instance()->GetTimeStep() >= mGrowthDuration)
                    {
                        p_population->UnmarkSpring(cell_pair);
                    }   
                }
            }
        }
    }
    else {
        EXCEPTION("The rotational division force is only implemented for node based cell populations.\n");
    }
}

template<unsigned DIM>
void RotationalDivisionForce<DIM>::SetNormalVector(c_vector<double,DIM> normalVec)
{
    assert(norm_2(normalVec)>0.0);
    mNormalVector=normalVec/norm_2(normalVec);
}

template<unsigned DIM>
void RotationalDivisionForce<DIM>::SetProportionalToSpringLength(bool propSpringLength)
{
    mPropSpringLength = propSpringLength;
}

template<unsigned DIM>
double RotationalDivisionForce<DIM>::GetTorsionCoefficient(double age) const
{
    return mTorsionCoefficient;
}

template<unsigned DIM>
double RotationalDivisionForce<DIM>::GetGrowthDuration() const
{
    return mGrowthDuration;
}

template<unsigned DIM>
void RotationalDivisionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TorsionCoefficient>" << mTorsionCoefficient << "</TorsionCoefficient>\n";
    for (unsigned i = 0; i < DIM; i++)
    {
        *rParamsFile << "\t\t\t<DirectionVector" << i << ">" << mNormalVector[i] << "</DirectionVector" << i << ">\n";
    }
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned DIM>
void RotationalDivisionForce<DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    *pVizSetupFile << "TorsionCoefficient\t" << mTorsionCoefficient << "\n";
    for (unsigned i = 0; i < DIM; i++)
    {
        *pVizSetupFile << "DirectionVector" << i << "\t" << mNormalVector[i] << "\n";
    }
    *pVizSetupFile << "GrowthDuration\t" << mGrowthDuration << "\n";
}

// Explicit instantiation
template class RotationalDivisionForce<1>;
template class RotationalDivisionForce<2>;
template class RotationalDivisionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RotationalDivisionForce)
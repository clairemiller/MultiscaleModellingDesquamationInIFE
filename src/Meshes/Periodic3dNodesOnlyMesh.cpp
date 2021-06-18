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

#include <map>
#include "Periodic3dNodesOnlyMesh.hpp"

template<unsigned SPACE_DIM>
Periodic3dNodesOnlyMesh<SPACE_DIM>::Periodic3dNodesOnlyMesh(double width)
        : NodesOnlyMesh<SPACE_DIM>(),
          mWidth(width)
{
    assert(width > 0.0);
    assert(SPACE_DIM==2);
}

template<unsigned SPACE_DIM>
void Periodic3dNodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(double cutOffLength, c_vector<double, SPACE_DIM*SPACE_DIM> domainSize, int numLocalRows, bool isPeriodic)
{
    // Ensure that the width is a multiple of cut-off length
    if( fmod( mWidth,cutOffLength ) > 1e-14 )
    {
        EXCEPTION("The periodic width must be a multiple of cut off length.");
    }
    else if ( mWidth/cutOffLength == 2.0 )
    {
        // A width of two boxes gives slightly different results. There is possibly some double counting happening somewhere.
        EXCEPTION( "The periodic domain width cannot be 2*CutOffLength." );
    }

    // We force the domain to the periodic width
    domainSize[0] = 0;
    domainSize[1] = mWidth;

    NodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(cutOffLength, domainSize, PETSC_DECIDE, true);    // Only difference is that this "true" makes the boxes periodic.
}

template<unsigned SPACE_DIM>
double Periodic3dNodesOnlyMesh<SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = MutableMesh<SPACE_DIM,SPACE_DIM>::GetWidth(rDimension);
    }
    return width;
}

template<unsigned SPACE_DIM>
void Periodic3dNodesOnlyMesh<SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point, bool concreteMove)
{
    // concreteMove should always be false for NodesOnlyMesh as no elements to check
    assert(!concreteMove);

    double x_coord = point.rGetLocation()[0];

    // Perform a periodic movement if necessary
    if (x_coord >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < 0.0)
    {
        double new_x_coord = x_coord + mWidth;
        double fudge_factor = 1e-14;
        // This is to ensure that the position is never equal to mWidth, which would be outside the box domain. 
        // This is due to the fact that mWidth-1e-16=mWidth
        if ( new_x_coord > mWidth-fudge_factor )
        {
            new_x_coord = mWidth-fudge_factor;
        }
        point.SetCoordinate(0, new_x_coord);
    }

    // Update the node's location
    this->GetNode(nodeIndex)->SetPoint(point);
}

template<unsigned SPACE_DIM>
unsigned Periodic3dNodesOnlyMesh<SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    // Call method on parent class
    unsigned node_index = NodesOnlyMesh<SPACE_DIM>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<SPACE_DIM> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;

}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> Periodic3dNodesOnlyMesh<SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2)
{
    c_vector<double, SPACE_DIM> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);

    /*
     * Handle the cylindrical condition here: if the points are more
     * than halfway around the cylinder apart, measure the other way.
     */
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }
    return vector;
}

// Refresh mesh -> Check and then move if outside box
template<unsigned SPACE_DIM>
void Periodic3dNodesOnlyMesh<SPACE_DIM>::RefreshMesh()
{

    // Check if the x values are in the domain, if not, get the fmod and relocate.
    unsigned num_nodes = mNodes.size();
    for (unsigned i=0; i<num_nodes; i++)
    {
        double& x_location = (mNodes[i]->rGetModifiableLocation())[0];
        if ( x_location < 0.0 )
        {
            x_location = fmod( x_location, mWidth) + mWidth;
        }
        else if ( x_location >= mWidth )
        {
            x_location = fmod( x_location, mWidth );
        }
    }

    // Now run the base class method
    NodesOnlyMesh<SPACE_DIM>::RefreshMesh();
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Periodic3dNodesOnlyMesh)

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

#ifndef CELLDEATHWRITER_HPP_
#define CELLDEATHWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * A modifier class which at each simulation time step calculates the volume of each cell
 * and stores it in in the CellData property as "volume". To be used in conjunction with
 * contact inhibition cell cycle models.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellDeathWriter : public AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

protected:

    unsigned mNumDeaths;

public:

    /**
     * Default constructor.
     */
    CellDeathWriter();

    /**
     * Destructor.
     */
    ~CellDeathWriter(){}

    
        /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation the population to write.
     */
    void VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation);
    
    /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
    * This line is appended to the output written by AbstractCellBasedWriter, which is a single
    * value [present simulation time], followed by a tab.
    *
    * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
    */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
    * This line is appended to the output written by AbstractCellBasedWriter, which is a single
    * value [present simulation time], followed by a tab.
    *
    * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
    */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the number of cells in the population that have each proliferative phase.
     *
    * This line is appended to the output written by AbstractCellBasedWriter, which is a single
    * value [present simulation time], followed by a tab.
    *
    * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
    */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellDeathWriter)

#endif /*CELLDEATHWRITER_HPP_*/

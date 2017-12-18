/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::linearMotionLayersFvMesh

Description
    A linear movement solver for a HPDC piston with layers

SourceFiles
    linearMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef linearMotionLayersFvMesh_H
#define linearMotionLayersFvMesh_H

#include "topoChangerFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class linearMotionLayersFvMesh Declaration
\*---------------------------------------------------------------------------*/

class linearMotionLayersFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearMotionLayersFvMesh(const linearMotionLayersFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const linearMotionLayersFvMesh&);

        //- Add linearMotionLayers zones and modifiers
        void addZonesAndModifiers();

        //- Make layering modifiers live
        void makeLayersLive();

        //- Return new point field
        tmp<pointField> newPoints() const;

public:

    //- Runtime type information
    TypeName("linearMotionLayersFvMesh");


    // Constructors

        //- Construct from database
        explicit linearMotionLayersFvMesh(const IOobject& io);

    //- Destructor
    virtual ~linearMotionLayersFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

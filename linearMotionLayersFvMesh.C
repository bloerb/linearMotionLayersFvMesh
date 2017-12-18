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

\*---------------------------------------------------------------------------*/

#include "linearMotionLayersFvMesh.H"
#include "Time.H"
#include "layerAdditionRemoval.H"
#include "pointFields.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

#include "IOdictionary.H"
#include "Function1.H"

/* Load balancing
#include "mapDistributePolyMesh.H"
#include "fvMeshDistribute.H"
*/
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearMotionLayersFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        linearMotionLayersFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearMotionLayersFvMesh::addZonesAndModifiers()
{
    // Check if point/face/cellZones or a meshModifier are present in polyMesh directory
    if
    (
        pointZones().size()
     || faceZones().size()
     || cellZones().size()
     || topoChanger_.size()
    )
    {
        Info<< "void linearMotionLayersFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    // Add zones and modifiers for motion action

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(0);
    List<faceZone*> fz(1);
    List<cellZone*> cz(0);


    // Add face zone for layer addition

    // read patch name from dynamicMeshDict
    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    // create Face Zone for patch
    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    labelList lpf(layerPatch.size());
    // read all patch faces
    forAll(lpf, i)
    {
        lpf[i] = layerPatch.start() + i;
    }

    fz[0] = new faceZone
    (
        "pistonFaceZone",
        lpf,
        boolList(layerPatch.size(), true),
        0,
        faceZones()
    );

    // adding them to polyMesh folder
    Info<< "Adding point and face zones" << endl;
    addZones(pz,fz,cz);

    // Add a topology modifier

    List<polyMeshModifier*> tm(1);

    tm[0] =
        new layerAdditionRemoval
        (
            "pistonModifier",
            0,
            topoChanger_,
            "pistonFaceZone",
            readScalar
            (
                motionDict_.subDict("layer").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("layer").lookup("maxThickness")
            )
        );

    // adding it to polyMesh folder
    Info<< "Adding topology modifiers" << endl;
    topoChanger_.addTopologyModifiers(tm);

    // Write mesh
    write();
}


void Foam::linearMotionLayersFvMesh::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll(topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void linearMotionLayersFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}

Foam::tmp<Foam::pointField> Foam::linearMotionLayersFvMesh::newPoints() const
{
    tmp<pointField> tnewPoints
    (
        new pointField(points())
    );
    
    pointField& np = tnewPoints.ref();

    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    const labelList& patchPoints = layerPatch.meshPoints();

    autoPtr<Function1<scalar>> vel
    (
        Function1<scalar>::New
        (
            "vel",
            motionDict_
        )
    );    
    
    const vectorField& nHat = layerPatch.pointNormals();
    
    forAll(patchPoints, ppI)
    {
        np[patchPoints[ppI]] += -nHat[ppI]*vel().value(time().value())*time().deltaTValue();
    }

    return tnewPoints;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearMotionLayersFvMesh::linearMotionLayersFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearMotionLayersFvMesh::~linearMotionLayersFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::linearMotionLayersFvMesh::update()
{
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true);

    makeLayersLive();

    movePoints(newPoints());

/* Load balancing 

    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        scalar nIdealCells =
            mesh_.globalData().nTotalCells()
          / Pstream::nProcs();

        scalar unbalance = returnReduce
        (
            mag(1.0-mesh_.nCells()/nIdealCells),
            maxOp<scalar>()
        );

        if (unbalance <= 0.5)
        {
            Info<< "Skipping balancing since max unbalance " << unbalance
                << " is less than allowable 0.5"
                << endl;
        }
        else
        {
            scalarField cellWeights(mesh_.nCells(), 1);

            distMap = balance
            (
                false,  //keepZoneFaces
                false,  //keepBaffles
                cellWeights,
                decomposer,
                distributor
            );
        }
    }
*/

    return true;
}


// ************************************************************************* //


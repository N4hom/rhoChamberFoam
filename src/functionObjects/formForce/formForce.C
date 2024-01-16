/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "formForce.H"
#include "fvcGrad.H"
//#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(formForce, 0);

    addToRunTimeSelectionTable(functionObject, formForce, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::functionObjects::formForce::formForce
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchNames_(wordList(dict.lookup("patches"))),
    totalForces_(wordList(dict.lookup("patches")).size()),
    forceDir_(Zero),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_(0),
    out_("force.dat")
{
    
    out_ << "time(s); ";

    forAll(patchNames_, patchi)
    {
        out_ << patchNames_[patchi] << ".x; " ;
        out_ << patchNames_[patchi] << ".y; " ;
        out_ << patchNames_[patchi] << ".z " ;


    }

    
    

    out_ << endl;

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::formForce::~formForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::formForce::read(const dictionary& dict)
{
    
    fvMeshFunctionObject::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    Info << "patch names " << pbm.names();
    //forceDir_ = dict.lookup("forceDirection");
    

    return true;
}


void Foam::functionObjects::formForce::calcForcesMoment()
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    //- SurfaceVectorField containing the area of each face at the boundary and its orientation
    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();


    
    label patchIndex = 0; 

    //- Calculate the force on each patch
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        Info << "patchIndex " << patchIndex << endl;
        const label patchi = iter.key();
        const vectorField force(Sfb[patchi]*(p.boundaryField()[patchi] - pRef_));
       
        //- Sum over all faces
        Vector<scalar> totalForce = gSum(force);

        //Vector<scalar> sumForce = sum(force);

        scalar patchArea = gSum(mag(Sfb[patchi]));
        
        Info << "totalForce "  << totalForce[0] << " " << totalForce[1] << "  " << totalForce[2]  << "N over " << patchArea << " m2" <<  endl;
//        Info << "sumForce "  << sumForce[0] << " " << sumForce[1] << "  " << sumForce[2]  << "N over " << patchArea << " m2" <<  endl;
        Info << "before totalForces_[patchIndex] " << endl;
        totalForces_[patchIndex] = totalForce;
        Info << "before patchIndex++ " << endl;

        patchIndex++;
    }

    
}





bool Foam::functionObjects::formForce::execute()
{
    calcForcesMoment();


    /*OFstream out("blastFoam.dat");

    out << "time(s);";
    out << endl << endl;*/
    

    return true;
}





bool Foam::functionObjects::formForce::write()
{

    out_.precision(5);
    out_ << mesh_.time().value() << "; ";

    forAll(totalForces_, patchi)
    {
        out_ << totalForces_[patchi][0] << ";  ";
        out_ << totalForces_[patchi][1] << ";  ";
        out_ << totalForces_[patchi][2] ;

        if (patchi < totalForces_.size() - 1)
        {
            out_ << "; ";
        }


    }
    

    out_ << endl;
    
   

   /* 
    logFiles::write();

    Log << type() << " " << name() << " write:" << nl;*/
    return true;
}


// ************************************************************************* //

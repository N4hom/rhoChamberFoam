/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "mywaveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mywaveTransmissiveFvPatchField<Type>::mywaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    speedOfSoundName_("thermo:speedOfSound"),
    gamma_(0.0)
{}


template<class Type>
Foam::mywaveTransmissiveFvPatchField<Type>::mywaveTransmissiveFvPatchField
(
    const mywaveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    speedOfSoundName_(ptf.speedOfSoundName_),
    gamma_(ptf.gamma_)
{}


template<class Type>
Foam::mywaveTransmissiveFvPatchField<Type>::mywaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    speedOfSoundName_(dict.getOrDefault<word>("speedOfSound", "thermo:speedOfSound")),
    gamma_(dict.get<scalar>("gamma"))
{}


template<class Type>
Foam::mywaveTransmissiveFvPatchField<Type>::mywaveTransmissiveFvPatchField
(
    const mywaveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    speedOfSoundName_(ptpsf.speedOfSoundName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
Foam::mywaveTransmissiveFvPatchField<Type>::mywaveTransmissiveFvPatchField
(
    const mywaveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    speedOfSoundName_(ptpsf.speedOfSoundName_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::mywaveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const auto& speedOfSoundp =
        this->patch().template lookupPatchField<volScalarField>(speedOfSoundName_);

    const auto& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    scalarField phip
    (
        this->patch().template
            lookupPatchField<surfaceScalarField>(this->phiName_)
    );

    if (phi.dimensions() == dimMass/dimTime)
    {
        const auto& rhop =
            this->patch().template
                lookupPatchField<volScalarField>(this->rhoName_);

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/speedOfSound)).
    return phip/this->patch().magSf() + speedOfSoundp;
}


template<class Type>
void Foam::mywaveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", this->rhoName_);
    os.writeEntryIfDifferent<word>("speedOfSound", "thermo:speedOfSound", speedOfSoundName_);

    os.writeEntry("gamma", gamma_);

    if (this->lInf_ > SMALL)
    {
        os.writeEntry("fieldInf", this->fieldInf_);
        os.writeEntry("lInf", this->lInf_);
    }

    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //

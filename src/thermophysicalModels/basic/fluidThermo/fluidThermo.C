/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName),
    f_
    (
        IOobject
        (
            "f",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    PV_
    (
        IOobject
        (
            "PV",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    PVs_
    (
        IOobject
        (
        	"PVs",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -3, -1, 0, 0)
    ),
    Sh_
    (
        IOobject
        (
        	"Sh",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    varf_
    (
        IOobject
        (
            "varf",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Chi_
    (
        IOobject
        (
            phasePropertyName("chi"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{}

Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicThermo(mesh, dict, phaseName),
    f_
    (
        IOobject
        (
            phasePropertyName("f"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    PVs_
    (
        IOobject
        (
        	phasePropertyName("PVs"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -3, -1, 0, 0)
    ),
    Sh_
    (
        IOobject
        (
        	"Sh",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    varf_
    (
        IOobject
        (
            phasePropertyName("varf"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Chi_
    (
        IOobject
        (
            phasePropertyName("chi"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}


Foam::volScalarField& Foam::fluidThermo::f()
{
    return f_;
}

const Foam::volScalarField& Foam::fluidThermo::f() const
{
    return f_;
}

// Progress Variable
Foam::volScalarField& Foam::fluidThermo::PV()
{
    return PV_;
}

const Foam::volScalarField& Foam::fluidThermo::PV() const
{
    return PV_;
}

// Progress Variable Source Term
Foam::volScalarField& Foam::fluidThermo::PVs()
{
    return PVs_;
}

const Foam::volScalarField& Foam::fluidThermo::PVs() const
{
    return PVs_;
}

/*
// mehdi 21032018
Foam::tmp<Foam::volScalarField> Foam::fluidThermo::PVs
(
	const volScalarField& f,
	const volScalarField& PV
) 
{

    const fvMesh& mesh = this->PVs_.mesh();

    tmp<volScalarField> tPVs
    (
        new volScalarField
        (
            IOobject
            (
                "PVs",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            PVs_.dimensions()
        )
    );
	tPVs() = this->getPVs(f,PV);
    return tPVs;
}
*/
// Enthalpy Source Term
Foam::volScalarField& Foam::fluidThermo::Sh()
{
    return Sh_;
}

const Foam::volScalarField& Foam::fluidThermo::Sh() const
{
    return Sh_;
}

Foam::volScalarField& Foam::fluidThermo::varf()
{
    return varf_;
}


const Foam::volScalarField& Foam::fluidThermo::varf() const
{
    return varf_;
}

Foam::volScalarField& Foam::fluidThermo::Chi()
{
    return Chi_;
}

const Foam::volScalarField& Foam::fluidThermo::Chi() const
{
    return Chi_;
}

// ************************************************************************* //

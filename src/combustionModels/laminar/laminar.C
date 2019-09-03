/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "laminar.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
laminar<Type>::laminar
(
    const word& modelType,
    const fvMesh& mesh
)
:
    Type(modelType, mesh),
    solver_(tableSolver(mesh, tables())),
    //he_(this->thermo().he()),
    f_(this->thermo().f()),
    mu_(this->thermo().mu()),
    alpha_(this->thermo().alpha()),
    ubIF_(mesh.cells().size()),
    ubP_(),
    posIF_(mesh.cells().size()),
    posP_(),
    hNorm_
    (
        IOobject
        (
            "hNorm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    integrateReactionRate_
    (
        this->coeffs().lookupOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }

    IOdictionary scalingParams
    (
       IOobject
       (
  	     "scalingParams",
          mesh.time().constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
       )
    );

    //hmaxTable_ = scalingParams.lookup("hmaxScalingTable");
    //hminTable_ = scalingParams.lookup("hminScalingTable");

	const polyBoundaryMesh& patches = mesh.boundaryMesh();
	int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }

    ubP_.setSize(patchSize);
    posP_.setSize(patchSize);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
laminar<Type>::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
hashedWordList laminar<Type>::tables()
{
	hashedWordList tableNames = this->thermo().composition().species();
	tableNames.append("T");
	tableNames.append("PVs");
    tableNames.append("mu");
    tableNames.append("alpha");

	return tableNames;
}

template<class Type>
tmp<Foam::volScalarField> laminar<Type>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class Type>
void laminar<Type>::correct()
{
	   const scalarField& fCells = f_.internalField();
	   //const scalarField& heCells = he_.internalField();
	   scalarField& muCells = mu_.internalField();
	   scalarField& alphaCells = alpha_.internalField();

       //scalarField& hNormCells = hNorm_.internalField();

       scalarList x(3, 0.0);
       //double hmin, hmax;

       // Interpolate for internal Field
          forAll(fCells, cellI)
          {

                 // Calculate normalized PV
                 x[0] = 1.;

                 //Calculate Zeta
                 x[1] = 0.;

        		 // f
        		 x[2] = fCells[cellI];

                 ubIF_[cellI] = solver_.upperBounds(x);
                 posIF_[cellI] = solver_.position(ubIF_[cellI], x);

                 //Lookup
                 muCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 2));
                 alphaCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));
          }

       // Interpolate for patches
       forAll(f_.boundaryField(), patchi)   //T direkt statt Umweg über h
       {
          const fvPatchScalarField& pf = f_.boundaryField()[patchi];
          //const fvPatchScalarField& phe = he_.boundaryField()[patchi];
          fvPatchScalarField& pmu = mu_.boundaryField()[patchi];
          fvPatchScalarField& palpha = alpha_.boundaryField()[patchi];

          //fvPatchScalarField& phNorm = hNorm_.boundaryField()[patchi];

              forAll(pf , facei)
              {

                     // Calculate normalized PV
                     x[0] = 1.;

                     //Calculate Zeta
                     x[1] = 0.;

            		 // f
            		 x[2] = pf[facei];

                     ubP_[facei] = solver_.upperBounds(x);
                     posP_[facei] = solver_.position(ubP_[facei], x);

                     //Lookup
                     pmu[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 2));
                     palpha[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));
              }
       }

        if (integrateReactionRate_)
        {
            word ddtScheme(this->mesh().ddtScheme("Yi"));

            if (ddtScheme == fv::localEulerDdtScheme<scalar>::typeName)
            {
                const scalarField& rDeltaT =
                    this->mesh().objectRegistry::
                    template lookupObject<volScalarField>
                    (
                        "rDeltaT"
                    );

                if (this->coeffs().found("maxIntegrationTime"))
                {
                    scalar maxIntegrationTime
                    (
                        readScalar(this->coeffs().lookup("maxIntegrationTime"))
                    );

                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxIntegrationTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
laminar<Type>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    fvScalarMatrix& Su = tSu();

    if (this->active())
    {
        const label specieI = this->thermo().composition().species()[Y.name()];

        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
laminar<Type>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        tdQ() = this->chemistryPtr_->dQ();
    }

    return tdQ;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
laminar<Type>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        tSh() = this->chemistryPtr_->Sh();
    }

    return tSh;
}

template<class Type>
bool laminar<Type>::correctDensity()
{
	return true;
}


template<class Type>
bool laminar<Type>::read()
{
    if (Type::read())
    {
        this->coeffs().lookup("integrateReactionRate")
            >> integrateReactionRate_;
        return true;
    }
    else
    {
        return false;
    }
}

}

}


// ************************************************************************* //

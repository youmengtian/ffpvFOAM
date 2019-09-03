/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    2018 Mehdi Jangi <m.jangi@bham.ac.uk> University of Birmingham

\*---------------------------------------------------------------------------*/

#include "rhoFPVModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
rhoFPVModel<CombThermoType, ThermoType>::rhoFPVModel
(
    const word& modelType, const fvMesh& mesh
)
:
    CombThermoType(modelType, mesh),
    solver_(tableSolver(mesh, tables())),
    Y_(this->thermo().composition().Y()),
    he_(this->thermo().he()),
    T_(this->thermo().T()),
    f_(this->thermo().f()),
    PV_(this->thermo().PV()),
    PVs_(this->thermo().PVs()),
    Sh_(this->thermo().Sh()),
    mu_(this->thermo().mu()),
    alpha_(this->thermo().alpha()),
    psi_(this->thermo().psi()),
    rho_(this->thermo().rho()),
    varf_(this->thermo().varf()),
    ubIF_(mesh.cells().size()),
    ubP_(),
    posIF_(mesh.cells().size()),
    posP_(),
	HRs_(this->thermo().Sh()),
    useMixtureFractionVariance_(this->coeffs().lookup("useMixtureFractionVariance")),
    Zeta_
    (
        IOobject
        (
            "Zeta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    PVNorm_
    (
        IOobject
        (
            "PVNorm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    hNorm_
    (
        IOobject
        (
            "hNorm",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    curY_(this->thermo().composition().Y().size()),
    curT_(this->thermo().T()),
    curPVs_(this->thermo().PVs()),
    curHRs_(this->thermo().Sh())   
{
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

    PVmaxTable_ = scalingParams.lookup("PVmaxScalingTable");


	const polyBoundaryMesh& patches = mesh.boundaryMesh();
	int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }

    ubP_.setSize(patchSize);
    posP_.setSize(patchSize);


    forAll(curY_, fieldi)
    {
        curY_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    "cur_"+Y_[fieldi].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Y_[fieldi]
            )
        );
    }

}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
rhoFPVModel<CombThermoType, ThermoType>::~rhoFPVModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
hashedWordList rhoFPVModel<CombThermoType, ThermoType>::tables()
{
	hashedWordList tableNames = this->thermo().composition().species();
	tableNames.append("T");
	tableNames.append("PVs");
	tableNames.append("mu");
	tableNames.append("alpha");
	tableNames.append("psi");
	tableNames.append("HeatRelease");

	return tableNames;
}

template<class CombThermoType, class ThermoType>
void rhoFPVModel<CombThermoType, ThermoType>::correct()
{

Info << "codeis here 00" <<endl;
    const scalarField& fCells = f_.internalField();
    const scalarField& varfCells = varf_.internalField();
    const scalarField& PVCells = PV_.internalField();

    scalarField& TCells = T_.internalField();
    scalarField& PVsCells = PVs_.internalField();
    scalarField& muCells = mu_.internalField();
    scalarField& alphaCells = alpha_.internalField();
    scalarField& psiCells = psi_.internalField();
    scalarField& HRsCells = HRs_.internalField();

    scalarField& ZetaCells = Zeta_.internalField();
    scalarField& PVNormCells = PVNorm_.internalField();
Info << "codeis here 01" <<endl;

    //- Update the species and temperature field
    if(this->active())
    {
       scalarList x(3, 0.0);
       double PVmax;

       // Interpolate for internal Field
       forAll(Y_, i)
       {
    	  scalarField& YCells = Y_[i].internalField();

          forAll(fCells, cellI)
          {
        	 if (i == 0)
        	 {
                 // Calculate normalized PV
                 PVmax = solver_.interp1D(fCells[cellI], PVmaxTable_);
                 PVNormCells[cellI] = max(min(PVCells[cellI]/max(PVmax,1e-10),1),0);
                 x[0] = PVNormCells[cellI];

                 //Calculate Zeta
                 ZetaCells[cellI] = min(sqrt(varfCells[cellI]/max(fCells[cellI]*(1 - fCells[cellI]), SMALL)),0.99);
                 if (useMixtureFractionVariance_) x[1] = ZetaCells[cellI];

        		 // f
        		 x[2] = fCells[cellI];

                 ubIF_[cellI] = solver_.upperBounds(x);
                 posIF_[cellI] = solver_.position(ubIF_[cellI], x);

                 //Lookup
            	 TCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 6));
            	 PVsCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 5));
            	 muCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 4));
            	 alphaCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 3));
            	 psiCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 2));
            	 HRsCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));

        	 }

        	 YCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], i);
          }
       }
Info << "codeis here 02" <<endl;

       // Interpolate for patches
       forAll(T_.boundaryField(), patchi)   //T direkt statt Umweg über h
       {
          const fvPatchScalarField& pvarf = varf_.boundaryField()[patchi];
          const fvPatchScalarField& pf = f_.boundaryField()[patchi];
          const fvPatchScalarField& pPV = PV_.boundaryField()[patchi];


          fvPatchScalarField& pT = T_.boundaryField()[patchi];
          fvPatchScalarField& pPVs = PVs_.boundaryField()[patchi];
          fvPatchScalarField& pmu = mu_.boundaryField()[patchi];
          fvPatchScalarField& palpha = alpha_.boundaryField()[patchi];
          fvPatchScalarField& ppsi = psi_.boundaryField()[patchi];
          fvPatchScalarField& pHRs = HRs_.boundaryField()[patchi];

          fvPatchScalarField& pZeta = Zeta_.boundaryField()[patchi];
          fvPatchScalarField& pPVNorm = PVNorm_.boundaryField()[patchi];

          forAll(Y_, i)
          {
        	  fvPatchScalarField& pY = Y_[i].boundaryField()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     // Calculate normalized PV
                     PVmax = solver_.interp1D(pf[facei], PVmaxTable_);
                     pPVNorm[facei] = max(min(pPV[facei]/max(PVmax,1e-10),1),0);
                     x[0] = pPVNorm[facei];

                     //Calculate Zeta
                     pZeta[facei] = min(sqrt(pvarf[facei]/max(pf[facei]*(1 - pf[facei]), SMALL)),0.99);
                     if (useMixtureFractionVariance_) x[1] = pZeta[facei];

            		 // f
            		 x[2] = pf[facei];

                     ubP_[facei] = solver_.upperBounds(x);
                     posP_[facei] = solver_.position(ubP_[facei], x);

                     //Lookup (T only if BC is not fixedValue)
                     pT[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 6));
                     pPVs[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 5));
                     pmu[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 4));
                     palpha[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 3));
                     ppsi[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 2));
                     pHRs[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));

             	 }

             	pY[facei] = solver_.interpolate(ubP_[facei], posP_[facei], i);
             }
          }

       }
    }

    rho_ = psi_*this->thermo().p();
}

//


// mehdi 21032018   Foam::tmp<Foam::volScalarField>
template<class CombThermoType, class ThermoType>
void rhoFPVModel<CombThermoType, ThermoType>::correctCurrentFields
(
	const volScalarField& f,
	const volScalarField& PV
)
{
    const scalarField& fCells = f.internalField();
    const scalarField& PVCells = PV.internalField();

    scalarField& TCells = curT_.internalField();
    scalarField& PVsCells = curPVs_.internalField();
    scalarField& HRsCells = curHRs_.internalField();

    scalarField PVNormCells = PVNorm_.internalField();

    //- Update the species and temperature field
    if(this->active())
    {
       scalarList x(3, 0.0);
       double PVmax;

       // Interpolate for internal Field
       forAll(Y_, i)
       {
    	  scalarField& YCells = curY_[i].internalField();

          forAll(fCells, cellI)
          {
        	 if (i == 0)
        	 {

                 // Calculate normalized PV
                 PVmax = solver_.interp1D(fCells[cellI], PVmaxTable_);
                 PVNormCells[cellI] = max(min(PVCells[cellI]/max(PVmax,1e-10),1),0);
                 x[0] = PVNormCells[cellI];

				 x[1] = 0;	

        		 // f
        		 x[2] = fCells[cellI];

//				 Info << "x= "<< x[0] <<" "<< x[1] <<" "<< x[2] << endl;
                 ubIF_[cellI] = solver_.upperBounds(x);

//				 Info << "ubIF_[cellI] "<< cellI <<" "<< ubIF_[cellI] << endl;

                 posIF_[cellI] = solver_.position(ubIF_[cellI], x);


//				 Info << "posIF_[cellI] "<< cellI <<" "<< posIF_[cellI] << endl;


//				 Info << "x= "<< x[0] <<" "<< x[1] <<" "<< x[2] << endl;

//				  Info <<cellI << " 1 useMixtureFractionVariance " <<" "<<  useMixtureFractionVariance_<< " " << x[1] << endl;	

                 //Lookup
                 TCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 6));
//				 Info << cellI <<" x[0] =" << x[0] <<" x[1] =" << x[1] <<" x[2] =" << x[2] << " in Tcells = " << TCells[cellI] << endl;

                 PVsCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 5));
//				 Info << cellI <<" x[0] =" << x[0] <<" x[1] =" << x[1] <<" x[2] =" << x[2] << " in PVs = " << PVsCells[cellI] << endl;

                 HRsCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));
//				 Info << cellI <<" x[0] =" << x[0] <<" x[1] =" << x[1] <<" x[2] =" << x[2] << " in HRs = " << HRsCells[cellI] << endl;


        	 }

        	 YCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], i);
//			 Info << i <<" "<< Y_[i].name() << " x[0] =" << x[0] <<" x[1] =" << x[1] <<" x[2] =" << x[2] << " in YCells = " << YCells[cellI] << endl;
			 		
          }
       }

       // Interpolate for patches
       forAll(curT_.boundaryField(), patchi)   //T direkt statt Umweg über h
       {
          const fvPatchScalarField& pf = f.boundaryField()[patchi];
          const fvPatchScalarField& pPV = PV.boundaryField()[patchi];

          fvPatchScalarField& pT = curT_.boundaryField()[patchi];
          fvPatchScalarField& pPVs = curPVs_.boundaryField()[patchi];
          fvPatchScalarField& pHRs = curHRs_.boundaryField()[patchi];
          fvPatchScalarField pPVNorm = PVNorm_.boundaryField()[patchi];

          forAll(Y_, i)
          {
        	  fvPatchScalarField& pY = curY_[i].boundaryField()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     // Calculate normalized PV
                     PVmax = solver_.interp1D(pf[facei], PVmaxTable_);
                     pPVNorm[facei] = max(min(pPV[facei]/max(PVmax,1e-10),1),0);
                     x[0] = pPVNorm[facei];

                     x[1] = 0;

            		 // f
            		 x[2] = pf[facei];

                     ubP_[facei] = solver_.upperBounds(x);
                     posP_[facei] = solver_.position(ubP_[facei], x);

                     pT[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 6));
                     pPVs[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 5));
                     pHRs[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));
             	 }

             	pY[facei] = solver_.interpolate(ubP_[facei], posP_[facei], i);
             }
          }
       }
       
	  }

}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::getHRs() const
{
    tmp<volScalarField> tHRs
    (
        new volScalarField
        (
            IOobject
            (
                "curHRs",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            curHRs_
        )
    );

    return tHRs;
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::getT() const
{
    tmp<volScalarField> tT
    (
        new volScalarField
        (
            IOobject
            (
                "curT",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            curT_
        )
    );

    return tT;
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::getPVs() const
{
    tmp<volScalarField> tPVs
    (
        new volScalarField
        (
            IOobject
            (
                "curPVs",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            curPVs_
        )
    );

    return tPVs;
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::getYi(const label& i) const
{
    tmp<volScalarField> tY
    (
        new volScalarField
        (
            IOobject
            (
                "tY",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            curY_[i]
        )
    );

    return tY;
}

// end mehdi 21032018

//

template<class CombThermoType, class ThermoType>
Switch rhoFPVModel<CombThermoType, ThermoType>::correctDensity()
{
	return true;
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::fvScalarMatrix>
rhoFPVModel<CombThermoType, ThermoType>::R
(
    volScalarField& Y              
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
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

    return tSh;
}

template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
rhoFPVModel< CombThermoType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
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

    return tdQ;
}

template<class CombThermoType, class ThermoType>
bool rhoFPVModel<CombThermoType, ThermoType>::read()
{
    if (CombThermoType::read())
    {
        this->coeffs().lookup("useScalarDissipation") >> useScalarDissipation_;
        this->coeffs().lookup("useMixtureFractionVariance") >> useMixtureFractionVariance_;
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

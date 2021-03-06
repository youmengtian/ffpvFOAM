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

Class
    Foam::combustionModels::modifiedFPVModel

Description
    Base class for modifiedFPVModel.

SourceFiles
    modifiedFPVModel.C

Contributors
    2015 Julian Zips <julian.zips@unibw.de> Universität der Bundeswehr München
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef modifiedFPVModel_H
#define modifiedFPVModel_H


#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "modifiedTableSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class modifiedFPVModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class modifiedFPVModel
:
    public CombThermoType
{

protected:

    // Protected data

        // tables
	    modifiedTableSolver solver_;

        //- Species mass fractions
        PtrList<volScalarField>& Y_;

        //- Enthalpy
        //const volScalarField& he_;

        //- Temperature
        volScalarField& T_;

        //- Mixture fraction
        const volScalarField& f_;

        //- Progress Variable
        const volScalarField& PV_;

        //- Progress Variable Source Term
        volScalarField& PVs_;

        //- Enthalpy Source Term
        volScalarField& Sh_;

        //- Viscosity
        volScalarField& mu_;

        //- Thermal conductivity
        volScalarField& alpha_;

        //- Mixture fraction variance
        const volScalarField& varf_;

        //- Upper Bounds for table interpolation
        List<List<int> > ubIF_, ubP_;

        //- Table Position
        List<scalarList> posIF_, posP_;

        //- Thermal conductivity
        volScalarField& HRs_;


private:

    // Private Member

        //- use scalar dissipation for interpolation switch
        Switch useScalarDissipation_;

        //- use mixture fraction variance for interpolation switch
        Switch useMixtureFractionVariance_;

        //- Disallow copy construct
	    modifiedFPVModel(const modifiedFPVModel&);

        //- Disallow default bitwise assignment
        void operator=(const modifiedFPVModel&);

        List<scalar> 	PVmaxTable_;
        //List<scalar> 	hmaxTable_;
        //List<scalar> 	hminTable_;

        //- Normalized Mixture fraction variance
        volScalarField Zeta_;

        //- Normalized Progress Variable
        volScalarField PVNorm_;

        //- Normalized Enthalpy
        //volScalarField hNorm_;
        
        
// mehdi 12042018
        PtrList<volScalarField> curY_;

        volScalarField curT_;
        
        volScalarField curPVs_;        

        volScalarField curHRs_;        

public:

    // Constructors

        //- Runtime type information
        TypeName("modifiedFPVModel");

        hashedWordList tables();

        //- Construct from components
        modifiedFPVModel(const word& modelType, const fvMesh& mesh);

        //- Destructor
        virtual ~modifiedFPVModel();

        //- Correct combustion rate
        virtual void correct();

		virtual tmp<volScalarField> getPVs() const;

		virtual tmp<volScalarField> getHRs() const;
		
		virtual tmp<volScalarField> getT() const;
		
		virtual tmp<volScalarField > getYi(label const& i) const;

		virtual void correctCurrentFields
        (
			const volScalarField& f,
            const volScalarField& PV
        );
        //- Correct density during PISO
        Switch correctDensity();

        //- Fuel consumption rate matrix  
        inline tmp<fvScalarMatrix> R(volScalarField& Y) const;         

        //- Enthalpy source term
        inline tmp<Foam::volScalarField> Sh() const;

        //- Heat release rate calculated from fuel consumption rate matrix
        inline tmp<volScalarField> dQ() const;

        // I-O
        //- Update properties
        virtual bool read();



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "modifiedFPVModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

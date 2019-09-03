/*---------------------------------------------------------------------------*\
  =========                 |

  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::fluidThermo

Description
    Fundamental fluid thermodynamic properties

SourceFiles
    fluidThermo.C

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    public basicThermo
{

protected:

    //- Mixture Fraction []
    volScalarField f_;

    //- Progress Variable [] //2DFGM
    volScalarField PV_;

    //- Progress Variable Source Term[kg/m3/s]
    volScalarField PVs_;

    //- Enthalpy Source Term[kg/m/s3]
    volScalarField Sh_;

    //- Mixture Fraction variance []
    volScalarField varf_;

     //- Scalar dissipation rate [1/s]
     volScalarField Chi_;
  
public:

    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~fluidThermo();


    // Member functions

        // Access to thermodynamic state variables

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual const scalarField& mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
	    	    
             //- Mixture Fraction []
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& f();
	    
	        //- Mixture Fraction []
	        virtual const volScalarField& f() const;

	        //- Progress Variable []
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& PV();

	        //- Progress Variable []
	        virtual const volScalarField& PV() const;

	        //- Progress Variable Source [kg/m^3/s]
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& PVs();

	        //- Progress Variable Source [kg/m^3/s] //2DFGM
	        virtual const volScalarField& PVs() const;



// mehdi 21032018

/*            //- Progress Variable Source [kg/m^3/s] //2DFGM
            //  for given mixture fraction and progress variable 
            virtual tmp<volScalarField> PVs
            (
                const volScalarField& f,
                const volScalarField& PV
            );


			virtual volScalarField& getPVs
            (
                const volScalarField& f,
                const volScalarField& PV
            );	
*/
//end mehdi 21032018



	        //- Enthalpy Source [kg/m/s^3]
	        //  Non-const access allowed for transport equations
	        virtual volScalarField& Sh();

	        //- Enthalpy Source [kg/m/s^3]
	        virtual const volScalarField& Sh() const;

	        //  Non-const access allowed for transport equations
	        virtual volScalarField& varf();

	        //- Mixture Fraction []
	        virtual const volScalarField& varf() const;

	        //  Non-const access allowed for transport equations
	        virtual volScalarField& Chi();

	        //- Scalar dissipation rate [1/s]
	        virtual const volScalarField& Chi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

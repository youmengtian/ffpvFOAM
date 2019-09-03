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

Application
    hfgmFoam

Description
    Solver for combustion with non-adiabatic FGM model.

@author Mehdi Jangi
@version 13.04.2018
@email m.jangi@bham.ac.uk

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoCombustionModel.H"

#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"

#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "IFstream.H"
#include "OFstream.H"
#include "Switch.H"
#include "Random.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createStochasticFields.H"
    #include "createSparkFields.H"        
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;


        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
    	{

		   if(sparkHotSpot)
		   {
			#include "sparkHotSpot.H"
		   }
		   else
		   {
			sparkingPVs *= 0.;
		   }

		   if(setReactingBox)
		   {
			#include "setReactingBox.H"
		   }

		   if(nonReacting)
		   {
			 meanPVs *= 0.;
		   }
		   else
		   {
			 meanPVs = thermo.PVs();
		   }
	        #include "UEqn.H"
		#include "PVEqn.H"
		#include "fEqn.H"
		turbulence->correctVarf();
		combustion->correct();
		meanf = thermo.f();
		meanPV = thermo.PV();
	
		Info << " T max/min: "<<  max(thermo.T()).value() <<" " << min(thermo.T()).value() << endl;

		// --- Pressure corrector loop

           	while (pimple.correct())
    	   	{
               		#include "pEqn.H"
           	}

           	turbulence->correct();

           	rho = thermo.rho();

        }

	if(runTime.write())
	{
		#include "computeRMS.H"
		forAll(Y, i)
		{
			Y[i].write();
		}

	}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	}

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

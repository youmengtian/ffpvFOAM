/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    translates cantera table Data to OpenFOAM

@author Gabriele Frank & Hagen Müller
@version 24.05.2013
@email hagen.mueller@unibw.de

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "fvCFD.H"
#include "rhoCombustionModel.H"
#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "matlabReader.H"
#include <stdio.h>
#include <stdlib.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"

   IOdictionary tableDict
   (
       IOobject
       (
          "tableProperties",
          runTime.constant(),
          runTime,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
       )
   );

   Info<< "Creating thermodynamics model\n" << endl;

   autoPtr<combustionModels::rhoCombustionModel> combustion
   (
      combustionModels::rhoCombustionModel::New
      (
         mesh
      )
   );

   rhoReactionThermo& thermo = combustion->thermo(); //2DFGM
   basicMultiComponentMixture& composition = thermo.composition();


   //create dummy tables
   hashedWordList dummytable(thermo.composition().species());
   List<scalar> dummy(0);

   for (int i=0; i<dummytable.size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
            dummytable[i]+"_table",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      word dummyName=dummytable[i]+"_table";
      dictionary.set(dummyName, dummy);
      OFstream output("constant/"+dummytable[i]+"_table");
      dictionary.writeHeader(output);
      dictionary.writeData(output);
   }

   //read the matlab output data
   matlabReader matlabRead(tableDict, thermo, composition, mesh);

   Info << nl << "Writing data..." << endl;

   for (int i=0; i<matlabRead.getNames().size(); i++)
   {
      IOdictionary dictionary
      (
         IOobject
         (
        	matlabRead.getNames()[i],
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         )
      );

      OFstream output("constant/"+matlabRead.getNames()[i]+"_table");
      matlabRead.write(i, dictionary, output);
   }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "End\n" << endl;

return 0;
}


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

Class
    Foam::matlabReader

Description
    Foam::matlabReader

SourceFiles
    matlabReader.C

\*---------------------------------------------------------------------------*/

#ifndef matlabReader_H
#define matlabReader_H

#include "fileName.H"
#include "speciesTable.H"
#include "IOdictionary.H"
#include "List.H"
#include "HashTable.H"
#include "scalar.H"
#include "rhoReactionThermo.H"
#include "OFstream.H"
#include <iostream>
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class chemkinReader Declaration
\*---------------------------------------------------------------------------*/

class matlabReader
{

private:

        void betaPDFIntegration(const scalar& varZ);

        //- Disallow default bitwise copy construct
        matlabReader(const matlabReader&);

        //- Disallow default bitwise assignment
        void operator=(const matlabReader&);

        int num_lines, num_chars, num_columns;

        basicMultiComponentMixture& composition;
        rhoReactionThermo& thermo;

        IOdictionary    tables_;
        hashedWordList  tableNames_;

        List<scalar> 	Zeta_param_;
        List<scalar>    PV_param_;
        List<scalar> 	f_param_;

        List<List<List<List<scalar> > > > sampledData_;
        List<List<List<scalar> > > readData_;

        List<List<scalar> > integratedData_;
        List<List<scalar> > unintegratedData_;



public:
        // Constructors
        matlabReader(const IOdictionary& matlabDict, rhoReactionThermo& thermo, basicMultiComponentMixture& composition, const fvMesh& mesh);

        //- Destructor
        virtual ~matlabReader();

        // Member functions
        int numberOfLines(){return num_lines;}

        int numberOfColumns(){return num_columns;}

        hashedWordList getNames();

        void write(const int& i,IOdictionary& dictionary,OFstream& output);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

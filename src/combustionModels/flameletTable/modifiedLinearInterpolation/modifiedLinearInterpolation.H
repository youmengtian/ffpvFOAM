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
    Foam::combustionModels::linearInter

Description
    Linear Table Interpolation

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "flameletTable.H"

namespace Foam
{

class fvMesh;

class modifiedLinearInterpolation
:
    public flameletTable
{

private:

	//- Disallow default bitwise copy construct
    modifiedLinearInterpolation(const modifiedLinearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const modifiedLinearInterpolation&);

    //- Table Values
    List<List<scalarList> > tableValues_ ; //Dimension 3DFGMFoam

public:

        static const List<List<scalarList> > defaultList; //Dimension 3DFGMFoam

        //- Runtime type information
        TypeName("modifiedLinearInterpolation");

        // Constructors

        //- Construct from components
        modifiedLinearInterpolation(const fvMesh& mesh, const word& tableName, const IOdictionary& tableDict);

	    //- Destructor
        virtual ~modifiedLinearInterpolation();

        //- Member functions

        //- Interpolate tableValues_ in a 3D-Table
        virtual inline scalar interpolate(const List<int>& ub, const scalarList& pos) const;
};

} // End Foam namespace

#endif /* LINEARINTERPOLATION_H_ */

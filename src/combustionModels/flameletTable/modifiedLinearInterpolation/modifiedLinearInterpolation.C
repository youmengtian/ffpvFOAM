/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/


#include "modifiedLinearInterpolation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<List<scalarList> > modifiedLinearInterpolation::defaultList(0.0); //Dimension 3DFGMFoam

defineTypeNameAndDebug(modifiedLinearInterpolation, 0);
addToRunTimeSelectionTable(flameletTable, modifiedLinearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

modifiedLinearInterpolation::modifiedLinearInterpolation(const fvMesh& mesh, const word& tableName, const IOdictionary& tableDict)
:
    flameletTable(mesh, tableName, tableDict),
	tableValues_(tableDict.lookupOrDefault<List<List<scalarList> > >(tableName, defaultList))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

modifiedLinearInterpolation::~modifiedLinearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

inline scalar modifiedLinearInterpolation::interpolate(const List<int>& ub, const scalarList& pos) const
{
   // Perform trilinear interpolation
	   // for PV

		label ub0m = ub[0] -1;
		if (ub0m < 0) ub0m =0;
		
		label ub1m = ub[1] -1;
		if (ub1m < 0) ub1m =0;

		label ub2m = ub[2] -1;
		if (ub2m < 0) ub2m =0;


	   scalar c00 = tableValues_[ub0m][ub1m][ub2m]*(1-pos[0]) + tableValues_[ub[0]][ub1m][ub2m]*pos[0];
	   scalar c10 = tableValues_[ub0m][ub[1]][ub2m]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub2m]*pos[0];
	   scalar c01 = tableValues_[ub0m][ub1m][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub1m][ub[2]]*pos[0];
	   scalar c11 = tableValues_[ub0m][ub[1]][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2]]*pos[0];
			
   
/*	   scalar c00 = tableValues_[ub[0] -1][ub[1] -1][ub[2] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2] -1]*pos[0];
	   scalar c10 = tableValues_[ub[0] -1][ub[1]][ub[2] -1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2] -1]*pos[0];
	   scalar c01 = tableValues_[ub[0] -1][ub[1] -1][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub[1] -1][ub[2]]*pos[0];
	   scalar c11 = tableValues_[ub[0] -1][ub[1]][ub[2]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]][ub[2]]*pos[0];*/
	   // for Zeta
	   scalar c0 = c00*(1-pos[1]) + c10*pos[1];
	   scalar c1 = c01*(1-pos[1]) + c11*pos[1];
	   // for Z
	   return c0*(1-pos[2]) + c1*pos[2];
}

} // End Foam namespace

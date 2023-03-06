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
    Foam::combustionModels::tableSolver

Description
    Abstract class for the interface between table look-up and combustion model.

SourceFiles
    tableSolver.C

\*---------------------------------------------------------------------------*/

#ifndef TABLESOLVER_H_
#define TABLESOLVER_H_

#include "Table.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "fvMesh.H"

namespace Foam
{
namespace combustionModels
{

class tableSolver
:
  public IOdictionary
{
public:
	
	enum searchMethod
	{
	    uniform,
	    bisect
	};

protected:

	//- List of the table used
	wordList tableNames_;

	//- List of the table parameters
	wordList paramNames_;

	//- Chemistry tables of all tabulated properties
	PtrList<Table> tables_;

	//- List of table parameter distributions
	List<scalarList> params_;

	label fp_;

	scalar PVMax_;

	searchMethod searchMethod_;


public:

    //- Runtime type information
    TypeName("tableSolver");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        tableSolver,
        dictionary,
        (
            const fvMesh &mesh
        ),
        (mesh)
    );

    //- Constructor
    tableSolver(const fvMesh& mesh, const word& tablePath, const wordList& tableNames);

    //- Destructor
    virtual ~tableSolver();

    //- Determine upper bounds for table interpolation
    List<int> upperBounds(const scalarList& x) const;

    int upperBounds(const scalar& x, const int& y) const;

    //- Determine table position
    scalarList position(const List<int>& ub, const scalarList& x) const;

    scalar position(const int& ub, const scalar& x, const int& y) const;

    //- Interpolate in lookup-table with progress variable, enthalpy, mixturefraction
    scalar interpolateTable(const List<int>& ub, const scalarList& pos, const label& i);

    //- Table name size
    int sizeTableNames() const;

    //- Return the search method as an enumeration
    searchMethod wordToSearchMethod(const word& searchMethod) const;

    scalar determineHt(const scalar& x) const;

    //- Determine the theta1 and theta2 values on the wall boundary based on CO2
    scalarList solveWallBoundary(const scalar& PV, const label& PVIndex);

};

} // End Foam namespace
} // End combustionModels namespace

#endif /* TABLESOLVER_H_ */

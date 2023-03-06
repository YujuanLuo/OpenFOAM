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
    Foam::combustionModels::tabulationSolver

Description
    Class that holds the table entries.

\*---------------------------------------------------------------------------*/

#ifndef TABLE2D_H_
#define TABLE2D_H_

#include "fvMesh.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "scalarList.H"

namespace Foam
{

class fvMesh;

class Table2D
:
 public IOdictionary

{
private:

//-the name of the table
    word tableName_;

public:

    //- Runtime type information
    TypeName("Table2D");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Table2D,
        dictionary,
        (
            const fvMesh &mesh,
	    const word& tablePath,
            const word& tableName
        ),
        (mesh, tablePath, tableName)
    );

    //- Constructor
    Table2D(const fvMesh& mesh, const word& tablePath, const word& tableName);

    //- Construct and return a clone
    autoPtr<Table2D> clone() const;

    //- Destructor
    virtual ~Table2D();

    // Selectors
    //- Return a reference to the selected tableSolver
    static autoPtr<Table2D> New(const fvMesh& mesh, const word& tablePath, const word& tableName);

    //- Return the table values
    virtual const List<scalarList> tableValues2D() const = 0;

    //- Return the interpolated value
    virtual scalar interpolate2D(const List<int>& ub, const scalarList& pos) const = 0;

};

} // End Foam namespace

#endif /* TABLE2D_H_ */

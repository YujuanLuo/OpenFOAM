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

\*---------------------------------------------------------------------------*/

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include "Table.H"
#include "typeInfo.H"
//#include "Interpolation.H"

namespace Foam
{


class linearInterpolation
:
    public Table
{

private:

        //- Disallow default bitwise copy construct
        linearInterpolation(const linearInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const linearInterpolation&);

	//- Table values
	List<scalarList> tableValues_;

public:

	static const List<scalarList> defaultList;

        //- Runtime type information
        TypeName("linearInterpolation");

        // Constructors

        //- Construct from components
        linearInterpolation(const fvMesh& mesh,const word& tablePath, const word& tableName);

	//- Destructor
        virtual ~linearInterpolation();

        //- Member functions
	
	//- tableValues of a 2D-Table
	const List<scalarList> tableValues() const;

	//- Interpolate tableValues_ in a 2D-Table
	scalar interpolate2D(const List<int>& ub, const scalarList& pos) const;

};

} // End Foam namespace

#endif /* LINEARINTERPOLATION_H_ */

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
    Foam::combustionModels::Interpolation

Description
    Table Interpolation

\*---------------------------------------------------------------------------*/

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "fvMesh.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "scalarList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{


class Interpolation
{

private:

        //- Disallow default bitwise copy construct
//        Interpolation(const Interpolation&);

        //- Disallow default bitwise assignment
        void operator=(const Interpolation&);

        //- Table Values
//        List<List<scalarList> > tableValues_;

public:


        //- Runtime type information
        TypeName("Interpolation");

	//- Declare run-time constructor selection table
	declareRunTimeSelectionTable
	(
	    autoPtr,
	    Interpolation,
	    interpo,
	    (const fvMesh& mesh),
	    (mesh)
	);

	//- Selectors
	
	//- Return a reference to the selected interpolation scheme
	static autoPtr<Interpolation> New(const fvMesh& mesh);

        // Constructors

        //- Construct from components
        Interpolation(const fvMesh& mesh);

	//- Destructor
        virtual ~Interpolation();

	virtual scalar interpolateValues(List<scalarList> tableValue2D, const List<int>& ub, const scalarList& pos) const =0;

	virtual scalar interpolateValues(List<List<scalarList> > tableValue3D, const List<int>& ub, const scalarList& pos) const =0;

//	virtual scalar interpolateValues(const List<List<List<scalarList> > >& tableValue4D,const scalarList& area, const List<int>& ub, const scalarList& pos, const scalarList& x) const =0;

};

} // End Foam namespace

#endif /* INTERPOLATION_H_ */

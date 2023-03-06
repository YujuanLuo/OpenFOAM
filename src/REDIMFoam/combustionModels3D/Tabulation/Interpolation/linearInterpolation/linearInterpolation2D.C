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

\*---------------------------------------------------------------------------*/

#include "linearInterpolation2D.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<scalarList> linearInterpolation2D::defaultList2D(0.0);

defineTypeNameAndDebug(linearInterpolation2D, 0);
addToRunTimeSelectionTable(Table2D, linearInterpolation2D, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolation2D::linearInterpolation2D(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
Table2D(mesh, tablePath, tableName),
tableValues2D_(this->lookupOrDefault<List<scalarList> > (tableName, defaultList2D))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearInterpolation2D::~linearInterpolation2D()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //
const List<scalarList> linearInterpolation2D::tableValues2D() const
{
    return tableValues2D_;
}

scalar linearInterpolation2D::interpolate2D(const List<int>& ub, const scalarList& pos) const
{
      scalar c0 = tableValues2D_[ub[0] -1][ub[1] -1]*(1-pos[0]) + tableValues2D_[ub[0]][ub[1] -1]*pos[0];
      scalar c1 = tableValues2D_[ub[0] -1][ub[1]]*(1-pos[0]) + tableValues2D_[ub[0]][ub[1]]*pos[0];

      return c0*(1-pos[1]) + c1*pos[1];
}

} // End Foam namespace

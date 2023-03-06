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

#include "linearInterpolation3D.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<List<scalarList> > linearInterpolation3D::defaultList3D(0.0);

defineTypeNameAndDebug(linearInterpolation3D, 0);
addToRunTimeSelectionTable(Table3D, linearInterpolation3D, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolation3D::linearInterpolation3D(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
Table3D(mesh, tablePath, tableName),
tableValues3D_(this->lookupOrDefault<List<List<scalarList> > >(tableName, defaultList3D))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearInterpolation3D::~linearInterpolation3D()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //
const List<List<scalarList> > linearInterpolation3D::tableValues3D() const
{
    return tableValues3D_;
}

scalar linearInterpolation3D::interpolate3D(const List<int>& ub, const scalarList& pos) const
{
     scalar c00 = tableValues3D_[ub[0] -1][ub[1] -1][ub[2] -1]*(1-pos[0]) + tableValues3D_[ub[0]][ub[1] -1][ub[2] -1]*pos[0];
     scalar c10 = tableValues3D_[ub[0] -1][ub[1]][ub[2] -1]*(1-pos[0]) + tableValues3D_[ub[0]][ub[1]][ub[2] -1]*pos[0];
     scalar c01 = tableValues3D_[ub[0] -1][ub[1] -1][ub[2]]*(1-pos[0]) + tableValues3D_[ub[0]][ub[1] -1][ub[2]]*pos[0];
     scalar c11 = tableValues3D_[ub[0] -1][ub[1]][ub[2]]*(1-pos[0]) + tableValues3D_[ub[0]][ub[1]][ub[2]]*pos[0];

     scalar c0 = c00*(1-pos[1]) + c10*pos[1];
     scalar c1 = c01*(1-pos[1]) + c11*pos[1];

     return c0*(1-pos[2]) + c1*pos[2];
}

} // End Foam namespace

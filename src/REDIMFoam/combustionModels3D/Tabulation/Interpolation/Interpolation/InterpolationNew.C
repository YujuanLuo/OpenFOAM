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

\*---------------------------------------------------------------------------*/

#include "Interpolation.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::Interpolation> Foam::Interpolation::New(const fvMesh& mesh)
{
    const word interpolationType
    (
        IOdictionary
	(
	   IOobject
	   (
	      "tableProperties",
	      mesh.time().constant(),
	      mesh,
	      IOobject::MUST_READ,
	      IOobject::NO_WRITE,
	      false
	   )
	).lookup("InterpolationScheme")
    );

    interpoConstructorTable::iterator cstrIter =
    interpoConstructorTablePtr_->find(interpolationType);

    if (cstrIter == interpoConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Interpolation::New(const word&, "
            "const fvMesh& mesh)"
        )   << "Unknown interpolation type " << interpolationType
            << "Valid interpolation types : " << endl
            << interpoConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}


// ************************************************************************* //

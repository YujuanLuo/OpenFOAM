/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "combustionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combustionModel, 0);
}

const Foam::word Foam::combustionModel::combustionPropertiesName
(
    "combustionProperties"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::combustionModel::createIOobject
(
    basicThermo& thermo,
    const word& combustionProperties
) const
{
    IOobject io
    (
        thermo.phasePropertyName(combustionProperties),
        thermo.db().time().constant(),
        thermo.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const word& modelType,
    basicThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    IOdictionary(createIOobject(thermo, combustionProperties)),
    mesh_(thermo.p().mesh()),
    turb_(turb),
    active_(getOrDefault<Switch>("active", true)),
    coeffs_(optionalSubDict(modelType + "Coeffs")),
    modelType_(modelType),
    nDims_(readInt(this->lookup("nDims"))),
    tableForm_(this->lookupOrDefault<word>("tableForm","flip")),
    nSpec_(readInt(this->lookup("nSpec"))),
    theta_(),
    stheta_()
{
    Info << "Dimensions: " << nDims_ << endl;

    theta_.setSize(nDims_);

    for (int m = 0; m < nDims_; m++)
    {
	word fileName=Foam::name(m+1);
        theta_.set
        (
          m,
          new volScalarField
          (
            IOobject
            (
              "theta"+fileName,
              mesh_.time().timeName(),
              mesh_,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
            ),
            mesh_
          )
        );
    }

    stheta_.setSize(nDims_);

    for (int m = 0; m < nDims_; m++)
    {
	word fileName=Foam::name(m);
        stheta_.set
        (
          m,
          new surfaceScalarField
          (
            IOobject
            (
              "stheta"+fileName,
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            fvc::interpolate(theta_[m])()
          )
        );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModel::~combustionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::combustionModel::read()
{
    if (regIOobject::read())
    {
        this->readEntry("active", active_);
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }

    return false;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "basicThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, fvMesh);
    defineRunTimeSelectionTable(basicThermo, fvMeshDictPhase);
}

const Foam::word Foam::basicThermo::dictName("thermophysicalProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf = this->T_.boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::Boundary& tbf = this->T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& name,
    bool& isOwner
)
{
    volScalarField* ptr =
        mesh.objectRegistry::getObjectPtr<volScalarField>(name);

    isOwner = !ptr;

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        // Transfer ownership of this object to the objectRegistry
        ptr->store();
    }

    return *ptr;
}


Foam::basicThermo::basicThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phaseName_(phaseName),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),

//    T_(lookupOrConstruct(mesh, phasePropertyName("T"), TOwner_)),
    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T",dimensionSet(0,0,0,1,0,0,0),300.0)
    ),

    TOwner_(getOrDefault<Switch>("updateT", TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),

    dpdt_(getOrDefault<Switch>("dpdt", true))
{}


Foam::basicThermo::basicThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),

    phaseName_(phaseName),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),

//    T_(lookupOrConstruct(mesh, phasePropertyName("T"), TOwner_)),
    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T",dimensionSet(0,0,0,1,0,0,0),300.0)
    ),

    TOwner_(getOrDefault<Switch>("updateT", TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),

    dpdt_(getOrDefault<Switch>("dpdt", true))
{}


Foam::basicThermo::basicThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    IOdictionary
    (
        IOobject
        (
            dictionaryName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phaseName_(phaseName),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),

//    T_(lookupOrConstruct(mesh, "T", TOwner_)),
    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T",dimensionSet(0,0,0,1,0,0,0),300.0)
    ),

    TOwner_(getOrDefault<Switch>("updateT", TOwner_)),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), Zero)
    ),

    dpdt_(getOrDefault<Switch>("dpdt", true))
{
    if (debug)
    {
        Pout<< "Constructed shared thermo : mesh:" << mesh.name()
            << " phase:" << phaseName
            << " dictionary:" << dictionaryName
            << " T:" << T_.name()
            << " updateT:" << TOwner_
            << " alphaName:" << alpha_.name()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<basicThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{
    db().checkOut("p");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::basicThermo& Foam::basicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    const basicThermo* thermo = pf.db().findObject<basicThermo>(dictName);

    if (thermo)
    {
        return *thermo;
    }

    HashTable<const basicThermo*> thermos =
        pf.db().lookupClass<basicThermo>();

    forAllConstIters(thermos, iter)
    {
        if
        (
            &(iter()->he().internalField())
         == &(pf.internalField())
        )
        {
            return *iter();
        }
    }

    return pf.db().lookupObject<basicThermo>(dictName);
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
         || he().name() == phasePropertyName(d)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList::null();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList::null();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


Foam::volScalarField& Foam::basicThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::basicThermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


bool Foam::basicThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //

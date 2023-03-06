#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

//pressureControl pressureControl(p, rho, pimple.dict(), false);

//mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> combustion
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);

const word gridRes(combustion().get<word>("gridResolution"));
const word tableForm(combustion().get<word>("tableForm"));
int  nSpec(combustion().get<int>("nSpec"));
//*******************************************************************
volScalarField& theta1 = combustion->theta(0);
volScalarField& theta2 = combustion->theta(1);

surfaceScalarField& stheta1 = combustion->stheta(0);
surfaceScalarField& stheta2 = combustion->stheta(1);

combustion->correct();

theta1.correctBoundaryConditions();
theta2.correctBoundaryConditions();

stheta1 = fvc::interpolate(theta1)();
stheta2 = fvc::interpolate(theta2)();

combustion->correctBoundary();

thermo.correct();

rho = thermo.rho(); 

phi = linearInterpolate(rho*U) & mesh.Sf();

volScalarField divTerm
(
    IOobject
    (
        "divTerm",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0) , 0.0)
);

//********************************************************************
label nproc = Pstream::nProcs();

//*******************************************************************

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"



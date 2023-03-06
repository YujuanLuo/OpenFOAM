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
//*******************************************************************
volScalarField& theta1 = combustion->theta(0);
volScalarField& theta2 = combustion->theta(1);
volScalarField& theta3 = combustion->theta(2);

surfaceScalarField& stheta1 = combustion->stheta(0);
surfaceScalarField& stheta2 = combustion->stheta(1);
surfaceScalarField& stheta3 = combustion->stheta(2);

labelField& nAtomsC_ = combustion->nAtomsC();
labelField& nAtomsH_ = combustion->nAtomsH();
labelField& nAtomsO_ = combustion->nAtomsO();

scalar& beta1 = combustion->beta1();
scalar& beta0 = combustion->beta0();

//*******************************************************************
speciesTable species;
dictionary fuel(thermo.subDict("fuel"));
dictionary oxidizer(thermo.subDict("oxidizer"));
const word fractionBasis(thermo.lookup("fractionBasis"));
if ((fractionBasis != "mass") && (fractionBasis != "mole"))
{
    FatalError << "Unknown fractionBasis type"
               << token::SPACE
               << fractionBasis
               << nl
               << "Valid types are: mass or mole."
               << abort(FatalError);
}

fileName chemkinFile(thermo.lookup("CHEMKINFile"));
fileName transFile(thermo.lookup("CHEMKINTransportFile"));
fileName thermoFile(thermo.lookup("CHEMKINThermoFile"));
chemkinReader cr(species, chemkinFile.expand(), transFile.expand(), thermoFile.expand(), false);
const HashPtrTable<gasHThermoPhysics>& speciesThermo = cr.speciesThermo();
const HashTable<List<specieElement> >& specieComposition = cr.specieComposition();

scalarField yFuel(species.size(), 0.0);
scalarField yOx(species.size(), 0.0);
scalarField xFuel(species.size(), 0.0);
scalarField xOx(species.size(), 0.0);
forAll(species, specieI)
{
    if (fractionBasis == "mass")
    {
            if (fuel.found(species[specieI]))
            {
                    yFuel[specieI] = readScalar(fuel.lookup(species[specieI]));
            }
            if (oxidizer.found(species[specieI]))
            {
                    yOx[specieI] = readScalar(oxidizer.lookup(species[specieI]));
            }
    }
    else
    {
            if (fuel.found(species[specieI]))
            {
                    xFuel[specieI] = readScalar(fuel.lookup(species[specieI]));
            }
            if (oxidizer.found(species[specieI]))
            {
                    xOx[specieI] = readScalar(oxidizer.lookup(species[specieI]));
            }
    }
}
scalar mwf = 0.0;
scalar mwo = 0.0;

if (fractionBasis == "mole")
{
    const scalar mTotf = sum(xFuel);
    const scalar mToto = sum(xOx);

    forAll(species, i)
    {
            xFuel[i] /= mTotf;
            xOx[i] /= mToto;

            mwf += speciesThermo[species[i]]->W()*xFuel[i];
            mwo += speciesThermo[species[i]]->W()*xOx[i];
    }

    forAll(species, i)
    {
            yFuel[i] = xFuel[i]*speciesThermo[species[i]]->W()/mwf;
            yOx[i] = xOx[i]*speciesThermo[species[i]]->W()/mwo;
    }
}

scalar C=0;
scalar H=0;
scalar O=0;
scalar C0=0;
scalar H0=0;
scalar O0=0;
forAll(species,specieI)
{
    forAll(specieComposition[species[specieI]], elemI)
    {
            word name = specieComposition[species[specieI]][elemI].name();
            label n = specieComposition[species[specieI]][elemI].nAtoms();

            if (name == "C")
            {
                    C = C + n*atomicWeights["C"]*yFuel[specieI]/speciesThermo[species[specieI]]->W();
                    C0 = C0 + n*atomicWeights["C"]*yOx[specieI]/speciesThermo[species[specieI]]->W();
                    nAtomsC_[specieI] = n;
            }
            if (name == "H")
            {
                    H = H + n*atomicWeights["H"]*yFuel[specieI]/speciesThermo[species[specieI]]->W();
                    H0 = H0 + n*atomicWeights["H"]*yOx[specieI]/speciesThermo[species[specieI]]->W();
                    nAtomsH_[specieI] = n;
            }
            if (name == "O")
            {
                    O = O + n*atomicWeights["O"]*yFuel[specieI]/speciesThermo[species[specieI]]->W();
                    O0 = O0 + n*atomicWeights["O"]*yOx[specieI]/speciesThermo[species[specieI]]->W();
                    nAtomsO_[specieI] = n;
            }

    }
}

beta1 = 2.0/atomicWeights["C"]*C+1.0/(2.0*atomicWeights["H"])*H-1.0/atomicWeights["O"]*O;
beta0 = 2.0/atomicWeights["C"]*C0+1.0/(2.0*atomicWeights["H"])*H0-1.0/atomicWeights["O"]*O0;

//*******************************************************************
//
combustion->correct();

theta1.correctBoundaryConditions();
theta2.correctBoundaryConditions();
theta3.correctBoundaryConditions();

stheta1 = fvc::interpolate(theta1)();
stheta2 = fvc::interpolate(theta2)();
stheta3 = fvc::interpolate(theta3)();

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



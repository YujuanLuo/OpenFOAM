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
    Foam::combustionModels::REDIMModel

Description
    Base class for REDIMModel.

SourceFiles
    REDIMModel.C

\*---------------------------------------------------------------------------*/

#ifndef REDIMModel_H
#define REDIMModel_H

#include "ThermoCombustion.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "tableSolver2D.H"
#include "tableSolver3D.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class REDIMModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class REDIMModel
:
    public ThermoCombustion<ReactionThermo>
{
protected:

    // Protected data

	const fvMesh& mesh_;

	//- Dimension of the REDIM
	int nDim_;

        // tables
        tableSolver3D solver_;

        tableSolver2D BoundarySolver_;

        //- Species mass fractions
        PtrList<volScalarField>& Y_;

        //- Species mass fractions
        PtrList<volScalarField> PM_;

        //- DPSI
        PtrList<volScalarField> DPSI_;

        PtrList<surfaceScalarField> DPSIS_;

        //- DPSE
        PtrList<volScalarField> DPSE_;

        PtrList<surfaceScalarField> DPSES_;

        //- DPSP
        volScalarField DPSP_;

        surfaceScalarField DPSPS_;

	//- rETA
	PtrList<volScalarField> rETA_;

        //- Enthalpy
        volScalarField& he_;

        //- Temperature
        volScalarField& T_;

        //- Progress variable production rate
        tmp<volScalarField> tProdRate_;

	scalarList wallTheta2_;

        scalarList wallTheta3_;

	//- Max and min of theta
	List<scalarList> thetalim_;

        //- species are saved for paraview
	wordList list_of_species_;

	label CO2Index, H2OIndex;

        word gridRes_;

        labelField nAtomsC_;

        labelField nAtomsH_;

        labelField nAtomsO_;

        scalar beta1_;

        scalar beta0_;

        scalar iterativeTolerance_;



private:

    // Private Member

        //- Disallow copy construct
        REDIMModel(const REDIMModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const REDIMModel&) = delete;

public:

    // Constructors

        //- Runtime type information
        TypeName("REDIMModel");

        hashedWordList tables();

	//- Flamelet tables selection
	hashedWordList selectedTables();

	hashedWordList selectedBoundaryTables();

        //- Construct from components
        REDIMModel
	(
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Destructor
        virtual ~REDIMModel();

        //- Correct combustion rate
        virtual void correct();

	virtual void correctBoundary();
        //virtual void correctPV();

        //- Correct boundary values based on normalized enthalpy and normalized CO2
        tmp<scalarField> correctWallTheta1(const label patchI);

        tmp<scalarField> correctWallTheta2(const label patchI);
        
	tmp<scalarField> correctWallTheta3(const label patchI);

        //- Correct density during PISO
        //Switch correctDensity();

        //- Fuel consumption rate matrix  
        tmp<fvScalarMatrix> R(volScalarField& Y) const;         

        //- Enthalpy source term
        //tmp<Foam::volScalarField> Sh() const;

        //- Heat release rate calculated from fuel consumption rate matrix
        tmp<volScalarField> Qdot() const;

        //- Return rETA
        const volScalarField& rETA(const label i);

        //- Return PM
        const volScalarField& PM(const label i);

        //- Return DPSI
        const volScalarField& DPS(const label i);

	const surfaceScalarField& DPSS(const label i);

	List<scalarList> getDPS(const scalar theta1R, const scalar theta1L, const scalar theta2R, const scalar theta2L, const scalar theta3R, const scalar theta3L);

        // I-O
        //- Update properties
        virtual bool read();

        //- Progress variable production rate
        const volScalarField& prodRate();

        labelField& nAtomsC();

        labelField& nAtomsH();

        labelField& nAtomsO();

        scalar& beta1();

        scalar& beta0();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "REDIMModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

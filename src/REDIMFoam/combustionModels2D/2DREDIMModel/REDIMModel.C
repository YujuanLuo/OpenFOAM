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

\*---------------------------------------------------------------------------*/

#include "REDIMModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class ReactionThermo>
Foam::combustionModels::REDIMModel<ReactionThermo>::REDIMModel
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    mesh_(this->mesh()),
    nDim_(this->nDims()),
    solver_(tableSolver(this->mesh(), "InternalTable", selectedTables())),
    Outletsolver_(tableSolver(this->mesh(), "OutletTable", selectedOutletTables())),
    Y_(thermo.composition().Y()),
    PM_(),
    DPSI_(),
    DPSIS_(),
    DPSE_(),
    DPSES_(),
    DPSP_
    (
	IOobject
	(
	   "DPSP",
	   this->mesh().time().timeName(),
	   this->mesh(),
	   IOobject::NO_READ,
	   IOobject::NO_WRITE
	),
	this->mesh(),
	dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0.0)
    ),
    DPSPS_(fvc::interpolate(DPSP_)),
    rETA_(),
    he_(thermo.he()),
    T_(thermo.T()),
    tProdRate_
    (
       new volScalarField
       (
           IOobject
           (
              "tProdRate",
               this->mesh().time().timeName(),
               this->mesh(),
               IOobject::NO_READ,
               IOobject::NO_WRITE
            ),
           this->mesh(),
           dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
        )
    ),
    wallTheta2_(),
    thetalim_(),
    gridRes_(),
    tableForm_(this->tableForm()),
    mode_(),
    Tw_()
{
    const polyBoundaryMesh& patches = this->mesh().boundaryMesh();
    int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }

    wallTheta2_.setSize(patchSize);

    IOdictionary tableProperties_
    (
       IOobject
       (
          "tableProperties",
          this->mesh().time().constant(),
          this->mesh(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
       )
    );


    list_of_species_.setSize(0);

    list_of_species_.append(wordList(tableProperties_.lookup("species")));

    PM_.setSize(nDim_*(list_of_species_.size()+2));

    for (int m = 0; m < (nDim_*(list_of_species_.size()+2)); m++)
    {
	PM_.set
	(
	  m,
	  new volScalarField
	  (
	    IOobject
	    (
	      "PM_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0),0.0)
          )
	);
    }

    DPSI_.setSize(nDim_*list_of_species_.size());

    for (int m = 0; m < (nDim_*list_of_species_.size()); m++)
    {
	DPSI_.set
	(
	  m,
	  new volScalarField
	  (
	    IOobject
	    (
	      "DPSI_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0.0)
	  )
	);
    }

    DPSIS_.setSize(nDim_*list_of_species_.size());

    for (int m = 0; m < (nDim_*list_of_species_.size()); m++)
    {
	DPSIS_.set
	(
	  m,
	  new surfaceScalarField
	  (
	    IOobject
	    (
	      "DPSIS_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0.0)
	  )
	);
    }

    DPSE_.setSize(nDim_);

    for (int m = 0; m < nDim_; m++)
    {
	DPSE_.set
	(
	  m,
	  new volScalarField
	  (
	    IOobject
	    (
	      "DPSE_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0.0)
	  )
	);
    }

    DPSES_.setSize(nDim_);

    for (int m = 0; m < nDim_; m++)
    {
	DPSES_.set
	(
	  m,
	  new surfaceScalarField
	  (
	    IOobject
	    (
	      "DPSES_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(1,-1,-1,0,0,0,0),0.0)
	  )
	);
    }

    rETA_.setSize(nDim_);

    for (int m = 0; m < nDim_; m++)
    {
	rETA_.set
	(
	  m,
	  new volScalarField
	  (
	    IOobject
	    (
	      "rETA_"+name(m),
	      this->mesh().time().timeName(),
	      this->mesh(),
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	  this->mesh(),
	  dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0)
          )
	);
    }

    thetalim_.setSize(0);
    thetalim_.append(scalarList(tableProperties_.lookup("thetaMax")));
    thetalim_.append(scalarList(tableProperties_.lookup("thetaMin")));

 
    forAll(list_of_species_,i)
    {
	label k = this->thermo().composition().species()[list_of_species_[i]];

	if (Y_[k].name()=="CO2")

	{
	    CO2Index = k;
	    break;
	}

    }

    gridRes_ = tableProperties_.lookupOrDefault<word>("gridResolution","coarse");
    //tableForm_ = tableProperties_.lookupOrDefault<word>("tableForm","flip");
    mode_ = tableProperties_.lookupOrDefault<word>("mode","mass");
    Tw_ = tableProperties_.lookupOrDefault<scalar>("Tw",300);
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::REDIMModel<ReactionThermo>::~REDIMModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class ReactionThermo>
Foam::hashedWordList Foam::combustionModels::REDIMModel<ReactionThermo>::selectedTables()
{
   IOdictionary tableProperties_
   (
      IOobject
      (
         "tableProperties",
         this->mesh().time().constant(),
         this->mesh(),
         IOobject::MUST_READ,
         IOobject::NO_WRITE
      )
   );
   hashedWordList tableNamesSelected(tableProperties_.lookup("species"));

   int list_of_species = tableNamesSelected.size();

   tableNamesSelected.append("PV");
   tableNamesSelected.append("T");
   tableNamesSelected.append("He");

   //-source terms
   for (int m = 1; m <= nDim_; m++)
   {
	tableNamesSelected.append("rETA"+name(m));
   }

   //-transport coefficient

//-PM1,3,5...109 for theta1, PM2,4,6...110 for theta2
   for (int m = 1; m <= (nDim_*(list_of_species+2)); m++)
   {
	tableNamesSelected.append("PM"+name(m));
   }
//-DPSE1 for theta1, DPSE2 for theta2
   for (int m = 1; m <= nDim_; m++)
   {
	tableNamesSelected.append("DPSE"+name(m));
   }
//-DPSI1,2,...53 for theta1, DPSI54,55,...106 for theta2
   for (int m = 1; m <= (nDim_*list_of_species); m++)
   {
	tableNamesSelected.append("DPSI"+name(m));
   }

   Info << "tableNamesSelected: " << tableNamesSelected << endl;

   return tableNamesSelected;
}

template<class ReactionThermo>
Foam::hashedWordList Foam::combustionModels::REDIMModel<ReactionThermo>::selectedOutletTables()
{

   hashedWordList tableNamesSelected;

   for (int m = 1; m <= nDim_; m++)
   {
	tableNamesSelected.append("theta"+name(m));
   }

   return tableNamesSelected;
}

template<class ReactionThermo>
void Foam::combustionModels::REDIMModel<ReactionThermo>::correct()
{
    scalarField& TCells = T_.primitiveFieldRef();
    scalarField& heCells = he_.primitiveFieldRef();

    //- Update the species and enthalpy field
    if(this->active())
    {

       scalarList cv(nDim_, 0.0);
       List<int> ubIF_(nDim_, 0);
       scalarList posIF_(nDim_, 0.0);
       List<int> ubIS_(nDim_, 0);
       scalarList posIS_(nDim_, 0.0);


       forAll(TCells, cellI)
       {
	    int m=0;

	    for (int ii = 0; ii < nDim_; ii++)
	    {
		cv[ii] = this->theta_[ii].primitiveField()[cellI];
		if (cv[ii] < thetalim_[1][ii]) cv[ii] = thetalim_[1][ii];
		if (cv[ii] > thetalim_[0][ii]) cv[ii] = thetalim_[0][ii];
	    }

            ubIF_ = solver_.upperBounds(cv);
	    posIF_ = solver_.position(ubIF_,cv);

            forAll(list_of_species_, i)
            {
		if (list_of_species_[i] != "ersatz1")
		{
                    label k = this->thermo().composition().species()[list_of_species_[i]];
                    Y_[k].primitiveFieldRef()[cellI] = max(solver_.interpolateTable(ubIF_, posIF_, i), 0.);
		}
            }

            TCells[cellI]= max(solver_.interpolateTable(ubIF_, posIF_, list_of_species_.size()+1), Tw_);

            heCells[cellI] = solver_.interpolateTable(ubIF_, posIF_, list_of_species_.size()+2);

	    int j = list_of_species_.size() + 3;

	    for (m=0; m < nDim_; m++)
	    {
		 rETA_[m].primitiveFieldRef()[cellI] = solver_.interpolateTable(ubIF_, posIF_, j+m);
	    }

	    j = j+m;
	    //-PM1-110
	    for (m=0; m < (nDim_*(list_of_species_.size()+2)); m++)
	    {
		 PM_[m].primitiveFieldRef()[cellI] = solver_.interpolateTable(ubIF_, posIF_, j+m);
	    }

	    j = j+m;
	    for (m=0; m < nDim_; m++)
	    {
		 DPSE_[m].primitiveFieldRef()[cellI] = -1.0 * solver_.interpolateTable(ubIF_, posIF_, j+m);
	    }

	    //-DPSI1-106
	    j = j+m;
	    for (m=0; m < (nDim_*list_of_species_.size()); m++)
	    {
		 DPSI_[m].primitiveFieldRef()[cellI] = -1.0 * solver_.interpolateTable(ubIF_, posIF_, j+m);
	    }

       }

       if (gridRes_ == "fine")
       {

            forAll(DPSPS_.internalField(), facei)
            {
               for (int ii = 0; ii < nDim_; ii++)
               {
                     cv[ii] = this->stheta_[ii].primitiveField()[facei];
                     if (cv[ii] < thetalim_[1][ii]) cv[ii] = thetalim_[1][ii];
                     if (cv[ii] > thetalim_[0][ii]) cv[ii] = thetalim_[0][ii];
               }


               ubIS_ = solver_.upperBounds(cv);
	       posIS_ = solver_.position(ubIS_,cv);

	       int j = list_of_species_.size() + 3 + nDim_ + nDim_ * (list_of_species_.size() + 2);
	       
	       int m=0;

	       for (m=0; m < nDim_; m++)
	       {
		     DPSES_[m].primitiveFieldRef()[facei] = -1.0 * solver_.interpolateTable(ubIS_, posIS_, j+m);
	       }

	       j = j+m;
	       for (m=0; m < (nDim_*list_of_species_.size()); m++)
	       {
		     DPSIS_[m].primitiveFieldRef()[facei] = -1.0 * solver_.interpolateTable(ubIS_, posIS_, j+m);
	       }

            }
       }
    }
}


template<class ReactionThermo>
void Foam::combustionModels::REDIMModel<ReactionThermo>::correctBoundary()
{
    if(this->active())
    {
       scalarList cv(nDim_, 0.0);
       List<int> ubP_(nDim_, 0);
       scalarList posP_(nDim_, 0.0);


       forAll(this->theta_[0].boundaryField(),patchi)
       {
	   fvPatchScalarField& pT = T_.boundaryFieldRef()[patchi];
	   fvPatchScalarField& pHe = he_.boundaryFieldRef()[patchi];

           forAll(pT, facei)
           {
		int m=0;

            	for (int ii = 0; ii < nDim_; ii++)
            	{
                	cv[ii] = this->theta_[ii].boundaryField()[patchi][facei];
                	if (cv[ii] < thetalim_[1][ii]) cv[ii] = thetalim_[1][ii];
                	if (cv[ii] > thetalim_[0][ii]) cv[ii] = thetalim_[0][ii];
            	}

                ubP_ = solver_.upperBounds(cv);
		posP_ = solver_.position(ubP_,cv);

                forAll(list_of_species_, i)
                {
		     if (list_of_species_[i] != "ersatz1")
		     {
                         label k = this->thermo().composition().species()[list_of_species_[i]];
                         Y_[k].boundaryFieldRef()[patchi][facei] = max(solver_.interpolateTable(ubP_, posP_, i), 0.);
		     }
                }
	
                pT[facei]= max(solver_.interpolateTable(ubP_, posP_, list_of_species_.size()+1), Tw_);

                pHe[facei] = solver_.interpolateTable(ubP_, posP_, list_of_species_.size()+2);

		int j = list_of_species_.size() + 3;

                for (m=0; m < nDim_; m++)
                {
                     rETA_[m].boundaryFieldRef()[patchi][facei] = solver_.interpolateTable(ubP_, posP_, j+m);
                }

		j = j+m;
		for (m=0; m < (nDim_*(list_of_species_.size()+2)); m++)
		{
		     PM_[m].boundaryFieldRef()[patchi][facei] = solver_.interpolateTable(ubP_, posP_, j+m);
		}

		j=j+m;
		for (m=0; m < nDim_; m++)
		{
		     DPSE_[m].boundaryFieldRef()[patchi][facei] = -1.0 * solver_.interpolateTable(ubP_, posP_, j+m);
		     if (gridRes_ == "fine") DPSES_[m].boundaryFieldRef()[patchi][facei] = -1.0 * solver_.interpolateTable(ubP_, posP_, j+m);
		}

		j=j+m;
		for (m=0; m < (nDim_*list_of_species_.size()); m++)
		{
		     DPSI_[m].boundaryFieldRef()[patchi][facei] = -1.0 * solver_.interpolateTable(ubP_, posP_, j+m);
		     if (gridRes_ == "fine") DPSIS_[m].boundaryFieldRef()[patchi][facei] = -1.0 * solver_.interpolateTable(ubP_, posP_, j+m);
		}
		
	   }

       }

    }
}


template<class ReactionThermo>
Foam::tmp<Foam::scalarField> Foam::combustionModels::REDIMModel<ReactionThermo>::correctWallTheta1(const label patchI)
{

   const fvPatchScalarField& pppv = Y_[CO2Index].boundaryField()[patchI];

   scalarList theta(2,0.0);

   tmp<scalarField> ttheta1(new scalarField(pppv.size()));
   scalarField& theta1 = ttheta1.ref();

   //-zeroGradient for CO2, only update before theta1 is solved
   scalar coeff_ = 1;
   if (mode_ == "mole") coeff_ = 0.0440098;

   const tmp<scalarField>& ppv = Y_[CO2Index].boundaryField()[patchI].patchInternalField()/coeff_;

   forAll(pppv,facei)
   {
	theta = solver_.solveWallBoundary(ppv()[facei],list_of_species_.size());

	theta1[facei] = theta[0];

        wallTheta2_[facei] = theta[1];
   }

   return ttheta1;
}



template<class ReactionThermo>
Foam::tmp<Foam::scalarField> Foam::combustionModels::REDIMModel<ReactionThermo>::correctWallTheta2(const label patchI)
{
    const fvPatchScalarField& pppv = Y_[CO2Index].boundaryField()[patchI];

    tmp<scalarField> ttheta2(new scalarField(pppv.size()));
    scalarField& theta2 = ttheta2.ref();

    forAll(pppv,facei)
    {
        theta2[facei] = wallTheta2_[facei];
    }

    return ttheta2;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::REDIMModel<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
}


template<class ReactionThermo>
Foam::List<Foam::scalarList>
Foam::combustionModels::REDIMModel<ReactionThermo>::getDPS(const scalar theta1R,const scalar theta1L,const scalar theta2R,const scalar theta2L)
{
    //-need to be changed for different tables

    List<scalarList> DPS;
    DPS.setSize(list_of_species_.size()+2);
    scalarList cv(2,0.0);
    List<int> ubR(2,0), ubL(2,0);
    scalarList posR(2,0.0), posL(2,0.0);
    scalar DPS1R, DPS1L, DPS2R, DPS2L;

    cv[0] = theta1R;
    cv[1] = theta2R;
    for (int ii = 0; ii < 2; ii++)
    {
          if (cv[ii] < thetalim_[1][ii]) cv[ii] = thetalim_[1][ii];
          if (cv[ii] > thetalim_[0][ii]) cv[ii] = thetalim_[0][ii];
    }
    ubR = solver_.upperBounds(cv);
    posR = solver_.position(ubR,cv);

    cv[0] = theta1L;
    cv[1] = theta2L;
    for (int ii = 0; ii < 2; ii++)
    {
          if (cv[ii] < thetalim_[1][ii]) cv[ii] = thetalim_[1][ii];
          if (cv[ii] > thetalim_[0][ii]) cv[ii] = thetalim_[0][ii];
    }
    ubL = solver_.upperBounds(cv);
    posL = solver_.position(ubL,cv);

    for(int m=0; m<list_of_species_.size()+2; m++)
    {
	//-DPS from interpolation

	if (m == 0)
	{
                if (tableForm_ == "flip")
                {
                        DPS1R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+1);
                        DPS2R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3);
                        DPS1L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+1);
                        DPS2L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3);
                }
                else
                {
                        DPS1R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3);
                        DPS2R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+1);
                        DPS1L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3);
                        DPS2L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+1);
                }
	}
	else if (m == 1)
	{
	   DPS1R = 0.0;
	   DPS2R = 0.0;
	   DPS1L = 0.0;
	   DPS2L = 0.0;
	}
	else
	{
		if (tableForm_ == "flip")
		{
	   		DPS1R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+m+list_of_species_.size());
	   		DPS2R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+m);
	   		DPS1L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+m+list_of_species_.size());
	   		DPS2L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+m);
		}
		else
		{
                        DPS1R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+m);
                        DPS2R = -1.0 * solver_.interpolateTable(ubR,posR,(nDim_+1)*(list_of_species_.size()+2)+3+m+list_of_species_.size());
                        DPS1L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+m);
                        DPS2L = -1.0 * solver_.interpolateTable(ubL,posL,(nDim_+1)*(list_of_species_.size()+2)+3+m+list_of_species_.size());
		}

	}

	DPS[m].setSize(4);

	DPS[m][0] = DPS1L;
	DPS[m][1] = DPS1R;
	DPS[m][2] = DPS2L;
	DPS[m][3] = DPS2R;

    }

    return DPS;
}


template<class ReactionThermo>
const Foam::volScalarField&
Foam::combustionModels::REDIMModel<ReactionThermo>::rETA(const label i)
{
    return rETA_[i];	
}


template<class ReactionThermo>
const Foam::volScalarField&
Foam::combustionModels::REDIMModel<ReactionThermo>::PM(const label i)
{
    return PM_[i];
}


template<class ReactionThermo>
const Foam::volScalarField&
Foam::combustionModels::REDIMModel<ReactionThermo>::DPS(const label i)
{
    if (i == 0)
    {
	return DPSE_[0];
    }
    else if (i==1)
    {
	return DPSP_;
    }
    else if (i<(list_of_species_.size()+2))
    {
	return DPSI_[i-2];
    }
    else if (i<(list_of_species_.size()+3))
    {
	return DPSE_[1];
    }
    else if (i<(list_of_species_.size()+4))
    {
	return DPSP_;
    }
    else
    {
	return DPSI_[i-4];
    }
}


template<class ReactionThermo>
const Foam::surfaceScalarField&
Foam::combustionModels::REDIMModel<ReactionThermo>::DPSS(const label i)
{
    if (i == 0)
    {
	return DPSES_[0];
    }
    else if (i==1)
    {
	return DPSPS_;
    }
    else if (i<(list_of_species_.size()+2))
    {
	return DPSIS_[i-2];
    }
    else if (i<(list_of_species_.size()+3))
    {
	return DPSES_[1];
    }
    else if (i<(list_of_species_.size()+4))
    {
	return DPSPS_;
    }
    else
    {
	return DPSIS_[i-4];
    }
}

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::REDIMModel<ReactionThermo>::R(volScalarField& VC) const
{
   const dimensionedScalar VCSMALL("CSMALL", dimless, SMALL);

   return fvm::SuSp(tProdRate_()/(VC+VCSMALL), VC);

}

template<class ReactionThermo>
bool Foam::combustionModels::REDIMModel<ReactionThermo>::read()
{
       return true;
}

template<class ReactionThermo>
const Foam::volScalarField&
Foam::combustionModels::REDIMModel<ReactionThermo>::prodRate()
{
    return tProdRate_();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

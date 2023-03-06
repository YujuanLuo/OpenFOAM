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
#include "tableSolver3D.H"
#include "OFstream.H"

namespace Foam
{
namespace combustionModels
{

  defineTypeNameAndDebug(tableSolver3D, 0);
  defineRunTimeSelectionTable(tableSolver3D, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tableSolver3D::tableSolver3D(const fvMesh& mesh, const word& tablePath, const wordList& tableNames)
:
   IOdictionary
   (
      IOobject
      (
         "tableProperties",
         mesh.time().constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
      )
   ),
   tableNames_(tableNames),
   paramNames_(3),
   tables_(tableNames_.size()),
   fp_(readLabel(this->lookup("fp"))),
   searchMethod_
   (
	wordToSearchMethod
	(
	    this->lookupOrDefault<word>("searchMethod","uniform")
	)
   )
{

    forAll(tableNames_,i)
    {
    	tableNames_[i] = tableNames_[i] + "_table";

    	tables_.set(i,Table3D::New(mesh,tablePath,tableNames_[i]));
    }

    paramNames_[0] = "theta1_param";
    paramNames_[1] = "theta2_param";
    paramNames_[2] = "theta3_param";
    	
    forAll(paramNames_, i)
    {
       	params_.append(scalarList(this->lookup(paramNames_[i])));
    }
  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tableSolver3D::~tableSolver3D()
{}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

//- Locate the position of flame tables
List<int> tableSolver3D::upperBounds(const scalarList& x) const
{
    // Upper, lower bound for table Interpolation and temp Value for bisection method
    List<int> ub(paramNames_.size(), 0);

    // Determine upper bounds and interpolation weights for table interpolation
    for (register label j=0; j<paramNames_.size(); j++)
    {
	label jj = j;

	ub[j] = upperBounds(x[j],jj);
    }

    // Return upper bounds
    return ub;
}

scalarList tableSolver3D::position(const List<int>& ub, const scalarList& x) const
{
    scalarList pos(paramNames_.size(), 0.0);

    for (register label j=0; j<paramNames_.size(); j++)
    {
	label jj = j;

	pos[j] = position(ub[j], x[j], jj);
    }

    return pos;
}

int tableSolver3D::upperBounds(const scalar& x, const int& y) const
{
    // Upper, lower bound for table Interpolation and temp Value for bisection method
    int ub(0);
    int tlb, tub;
    int newVal = 0;

    // Determine upper bounds and interpolation weights for table interpolation
    switch(searchMethod_)
    {
	case tableSolver3D::uniform:
	{
	    scalar last = params_[y].size()-1;
	    if (x == params_[y][0]) 
	    {
		ub = 1;
	    }
	    else
	    {
		scalar psi = (x - params_[y][0])/(params_[y][last]-params_[y][0]);
		psi = std::max(std::min(psi,1.0),0.0);
		ub = ceil(psi*last);
	    }
	    break;
	}
	case tableSolver3D::bisect:
	{
       	    tub = params_[y].size() - 1;
            tlb = 0;
            while (tub-1 != tlb)
            {
          	newVal = (tub + tlb)/2;
          	if (x < params_[y][newVal])
              		tub = newVal;
          	else
              		tlb = newVal;
       	    }
       	    ub = tub;
	    break;
	}
    } 

    // Return upper bounds
    return ub;
}

scalar tableSolver3D::position(const int& ub, const scalar& x, const int& y) const
{
    scalar pos(0.0);

    pos = (x - params_[y][ub-1]) / (params_[y][ub] - params_[y][ub-1]);

    return pos;
}

scalar tableSolver3D::interpolateTable3D(const List<int>& ub, const scalarList& pos, const label& i)
{
    return tables_[i].interpolate3D(ub, pos);
}


int tableSolver3D::sizeTableNames() const
{
    return tableNames_.size();
}


typename tableSolver3D::searchMethod
tableSolver3D::wordToSearchMethod
(
    const word& searchMethod
) const
{

    if (searchMethod == "uniform")
    {
	return tableSolver3D::uniform;
    } 
    else if (searchMethod == "bisect")
    {
	return tableSolver3D::bisect;
    }
    else
    {
	Info
	    << "bad searchMethod specifier "
	    << searchMethod << "using 'uniform' " << endl;
	
	return tableSolver3D::uniform;
    }
}

} // End combustionModels namespace
} // End Foam namespace

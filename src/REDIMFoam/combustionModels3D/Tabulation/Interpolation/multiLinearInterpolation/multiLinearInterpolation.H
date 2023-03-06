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
    Foam::combustionModels::linearInter

Description
    Linear Table Interpolation

\*---------------------------------------------------------------------------*/

#ifndef MULTILINEARINTERPOLATION_H_
#define MULTILINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "Interpolation.H"
#include "OpenSMOKE_ChiDistribution.hpp"
#include "OpenSMOKE_External_Functions.hpp"

namespace Foam
{

class fvMesh;

class multiLinearInterpolation
:
    public Interpolation
{

private:

    //- Disallow default bitwise copy construct
    multiLinearInterpolation(const multiLinearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const multiLinearInterpolation&);

    //- chi distribution kind
    //chiPDF_kinds CHIPDF_;
    string CHIPDF_;

    //- Log Normal chi distribution: variance of PDF sigma
    scalar chi_logNormal_sigma_;

    //- Log Normal chi distribution: number of integration intervals
    scalar chi_logNormal_number_of_points_;

    //- chi_st List
    scalarList chi_st_List_;

    //- scalar dissipation rate probability distribution function
    OpenSMOKE_ChiDistribution  chi_pdf;

    //- cell values: A, chi, Z, Z
    scalarList cellValues_;

public:

        static const List<List<List<scalarList> > > defaultList;

        //- Runtime type information
        TypeName("multiLinearInterpolation");

   // Constructors

        //- Construct from components
        multiLinearInterpolation();

	//- Destructor
        virtual ~multiLinearInterpolation();

   //- Member functions

        //- Interpolate tableValues_ in a 4D-Table
        virtual scalar interpolateValues(const scalarList& area, const List<int>& ub, const scalarList& pos, const scalarList& x);

	//- Interpolate tableValues_ in a 3D-Table: Dirac distribution
	void interpolateDirac(List<List<scalarList> >& table3D, const scalarList& area1, const List<int>& ub1, const scalarList& pos1, scalar& extracted);

        //- Interpolate tableValues_ in a 3D-Table: Log-Normal distribution
        void interpolateLogNormal(List<List<scalarList> >& table3D, const scalarList& area1, const List<int>& ub1, const scalarList& pos1, scalar& extracted);

	//- Interpolate tableValues_ in a 2D-Table
	void interpolatePDF(List<scalarList>& table2D, const scalarList& area2, const List<int>& ub2, const scalarList& pos2, scalar& extractedValue);

};
} // End Foam namespace

#endif /* MULTILINEARINTERPOLATION_H_ */

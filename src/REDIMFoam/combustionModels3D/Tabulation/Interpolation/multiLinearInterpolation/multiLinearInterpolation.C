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


#include "multiLinearInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "OpenSMOKE_ChiDistribution.hpp"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<List<List<scalarList> > > multiLinearInterpolation::defaultList(0.0);

defineTypeNameAndDebug(multiLinearInterpolation, 0);
addToRunTimeSelectionTable(Interpolation, multiLinearInterpolation, interpo);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiLinearInterpolation::multiLinearInterpolation()
:
    Interpolation(),
    chi_pdf(),
{
       cellValues_.resize(4);

	IOdictionary tableProperties_
	(
		IOobject
		(
		    "tableProperties",
		    mesh.time().constant(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);

       CHIPDF_ = tableProperties_.lookup("CHIPDF");

     if (CHIPDF_ == "logNormal")
     {
        chi_logNormal_sigma_ = readScalar(tableProperties_.lookup("sigma"));
        chi_logNormal_number_of_points_ = readScalar(tableProperties_.lookup("points"));
        chi_st_List_ = scalarList(tableProperties_.lookup("chi_param"));
        chi_pdf.SetSigma(chi_logNormal_sigma_);
        chi_pdf.SetNumberOfPoints(chi_logNormal_number_of_points_);
        chi_pdf.BuildGrid();
        std::vector<double> chi_stVector(chi_st_List_.size());
        forAll(chi_st_List_,i)
        {
          chi_stVector[i] = chi_st_List_[i];
        }
        chi_pdf.AssignScalarDissipationRates(chi_stVector);
     }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

multiLinearInterpolation::~multiLinearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

scalar multiLinearInterpolation::interpolateValues(const List<List<List<scalarList> > > & table4D, const scalarList& area, const List<int>& ub, const scalarList& pos, const scalarList& x) const
{
      cellValues_ = x;
   
      scalar extracted1(0.);
      scalar extracted2(0.);
      
      List<List<scalarList> > paraTable3DDW_(table4D[ub[0]-1]);
      List<List<scalarList> > paraTable3DUP_(table4D[ub[0]]);

      if (CHIPDF_ == "dirac")
      {
      interpolateDirac(paraTable3DDW_, area, ub, pos, extracted1);
      interpolateDirac(paraTable3DUP_, area, ub, pos, extracted2);
      }
      else
      {
      interpolateLogNormal(paraTable3DDW_, area, ub, pos, extracted1);
      interpolateLogNormal(paraTable3DUP_, area, ub, pos, extracted2);
      }

      return extracted1 + pos[0]*extracted2;
}

void multiLinearInterpolation::interpolateDirac(List<List<scalarList> >& table3D, const scalarList& area1, const List<int>& ub1, const scalarList& pos1, scalar& extracted)
{
     scalar extracted3(0.);
     scalar extracted4(0.);
     
     List<scalarList> paraTable2DDW_(table3D[ub1[1]-1]);
     List<scalarList> paraTable2DUP_(table3D[ub1[1]]);

     interpolatePDF(paraTable2DDW_, area1, ub1, pos1, extracted3);
     interpolatePDF(paraTable2DUP_, area1, ub1, pos1, extracted4);
     
     extracted = extracted3 + pos1[1]*extracted4;

}

void multiLinearInterpolation::interpolateLogNormal(List<List<scalarList> >& table3D, const scalarList& area1, const List<int>& ub1, const scalarList& pos1, scalar& extracted)
{ 
        double chi_st(cellValues_[1]);

	chi_pdf.AssignMeanScalarDissipationRate(chi_st);

        std::vector<double> extractedValues(2);

        for(label j=0; j<2;j++)
        {
               interpolatePDF(table3D[j], area1, ub1, pos1, extractedValues[j]);
        }

        extracted = chi_pdf.ExtractMeanValue(extractedValues);
        
}

void multiLinearInterpolation::interpolatePDF(List<scalarList>& table2D, const scalarList& area2, const List<int>& ub2, const scalarList& pos2, scalar& extractedValue)
{

      extractedValue = (table2D[ub2[2]-1][ub2[3]-1]*area2[1] + table2D[ub2[2]-1][ub2[3]]*area2[3] +table2D[ub2[2]][ub2[3]-1]*area2[2] + table2D[ub2[2]][ub2[3]]*area2[4])/area2[0];

}

} // End Foam namespace

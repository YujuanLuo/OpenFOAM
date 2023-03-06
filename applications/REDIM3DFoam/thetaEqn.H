{
     volScalarField diffusion1(0.0*divTerm);
     volScalarField diffusion2(0.0*divTerm);
     volScalarField diffusion3(0.0*divTerm);

     for(int m=0; m<22; m++)
     {

	divTerm = fvc::laplacian(combustion->DPS(m+22), theta1) + fvc::laplacian(combustion->DPS(m), theta2) + fvc::laplacian(combustion->DPS(m+44), theta3);

        diffusion1 = diffusion1 + combustion->PM(3*m+1) * divTerm;

        diffusion2 = diffusion2 + combustion->PM(3*m) * divTerm;
        	
	diffusion3 = diffusion3 + combustion->PM(3*m+2) * divTerm;

     }

//-Special treatment to the first cell layer near the wall

     if (gridRes == "fine")
     {
     	forAll (theta1.boundaryField(),patchI)
     	{

        	if (mesh.boundary()[patchI].type() == "wall")
        	{

                	const polyPatch& cPatch = mesh.boundaryMesh()[patchI];
                	const labelUList& faceCells = cPatch.faceCells();

                	if (faceCells.size()>0)
                	{
                        	scalar deltax(mesh.deltaCoeffs().boundaryField()[patchI][0]);
                        	scalar deltay(1.0);
                        	scalar theta1FS, theta2FS, theta3FS;
                        	scalar theta1R, theta1L, theta2R, theta2L, theta3R, theta3L;
                        	scalar theta1U, theta1B, theta2U, theta2B, theta3U, theta3B;
                        	scalar theta1UU, theta1BB, theta2UU, theta2BB, theta3UU, theta3BB;
                        	scalar gradxtheta1R, gradxtheta1L, gradxtheta2R, gradxtheta2L, gradxtheta3R, gradxtheta3L;
                        	scalar gradytheta1U, gradytheta1B, gradytheta2U, gradytheta2B, gradytheta3U, gradytheta3B;
                        	List<scalarList> DPSx, DPSy;
                        	label  patchBottom(7), patchTop(8); //-have to predefine the ID of bottom and top patch
                        	label  patchOutlet(4), patchInlet(0);
                        	label  patchTopBottom(7), patchBottomTop(7);

                        	forAll(faceCells,i)
                        	{
                                	label cellID = faceCells[i];
                                	const labelUList& faces = mesh.cells()[cellID];
                                	label facexID(-1), faceyuID(-1), faceybID(-1);
                                	bool bottom(true), top(true);

                                	scalar xFace = 2.0 / mesh.deltaCoeffs().boundaryField()[patchI][i];
                                	scalar halfxFace = 1.0 / mesh.deltaCoeffs().boundaryField()[patchI][i];
                                	scalar yCell = mesh.C().primitiveField()[cellID].component(vector::Y);

                                	forAll(faces,j)
                                	{
                                        	if (faces[j] < mesh.Cf().size())
                                        	{
                                                	scalar xCoor = mesh.Cf().primitiveField()[faces[j]].component(vector::X);
                                                	scalar yCoor = mesh.Cf().primitiveField()[faces[j]].component(vector::Y);

                                                	if (mag(xCoor - xFace) < 1e-16) facexID = faces[j];
                                                	if ((mag(xCoor - halfxFace) < 1e-16) && (yCoor < yCell) )
                                                	{
                                                        	faceybID = faces[j];
                                                        	bottom = false;
                                                	}
                                                	if ((mag(xCoor - halfxFace) < 1e-16) && (yCoor > yCell) )
                                                	{
                                                        	faceyuID = faces[j];
                                                        	top = false;
                                                	}
                                        	}
                                	}

                                	theta1FS = stheta1.primitiveField()[facexID];
                                	theta2FS = stheta2.primitiveField()[facexID];
                                	theta3FS = stheta3.primitiveField()[facexID];

                                	gradxtheta1R = (theta1FS - theta1.primitiveField()[cellID]) * deltax;
                                	gradxtheta1L = (theta1.primitiveField()[cellID] - theta1.boundaryField()[patchI][i]) * deltax;
                                	gradxtheta2R = (theta2FS - theta2.primitiveField()[cellID]) * deltax;
                                	gradxtheta2L = (theta2.primitiveField()[cellID] - theta2.boundaryField()[patchI][i]) * deltax;
                                	gradxtheta3R = (theta3FS - theta3.primitiveField()[cellID]) * deltax;
                                	gradxtheta3L = (theta3.primitiveField()[cellID] - theta3.boundaryField()[patchI][i]) * deltax;

                                	theta1R = (theta1.primitiveField()[cellID] + theta1FS) * 0.5;
                                	theta1L = (theta1.boundaryField()[patchI][i] + theta1.primitiveField()[cellID]) * 0.5;
                                	theta2R = (theta2.primitiveField()[cellID] + theta2FS) * 0.5;
                                	theta2L = (theta2.boundaryField()[patchI][i] + theta2.primitiveField()[cellID]) * 0.5;
                                	theta3R = (theta3.primitiveField()[cellID] + theta3FS) * 0.5;
                                	theta3L = (theta3.boundaryField()[patchI][i] + theta3.primitiveField()[cellID]) * 0.5;

                                	DPSx = combustion->getDPS(theta1R,theta1L,theta2R,theta2L,theta3R,theta3L);

                                	if ((bottom == false) && (top == false)) //-internal cell without boundary faces
                                	{
                                        	theta1U = stheta1.primitiveField()[faceyuID];
                                        	theta1B = stheta1.primitiveField()[faceybID];
                                        	theta2U = stheta2.primitiveField()[faceyuID];
                                        	theta2B = stheta2.primitiveField()[faceybID];
                                        	theta3U = stheta3.primitiveField()[faceyuID];
                                        	theta3B = stheta3.primitiveField()[faceybID];

                                        	deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));

                                	}
                                	else if ((bottom == true) && (top == false)) //boundary at the Bottom
                                	{
                                        	theta1U = stheta1.primitiveField()[faceyuID];
                                        	theta2U = stheta2.primitiveField()[faceyuID];
                                        	theta3U = stheta3.primitiveField()[faceyuID];

                                        	if (Pstream::myProcNo() == 0)
                                        	{
                                                	theta1B = stheta1.boundaryField()[patchInlet][0];
                                                	theta2B = stheta2.boundaryField()[patchInlet][0];
                                                	theta3B = stheta3.boundaryField()[patchInlet][0];

                                                	deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchInlet][0].component(vector::Y));
                                        	}
                                        	else if (Pstream::myProcNo() == (nproc - 1))
                                        	{
                                                	theta1B = stheta1.boundaryField()[patchTopBottom][0];
                                                	theta2B = stheta2.boundaryField()[patchTopBottom][0];
                                                	theta3B = stheta3.boundaryField()[patchTopBottom][0];

                                                	deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchTopBottom][0].component(vector::Y));
                                        	}
                                        	else
                                        	{
                                                	theta1B = stheta1.boundaryField()[patchBottom][0];
                                                	theta2B = stheta2.boundaryField()[patchBottom][0];
                                                	theta3B = stheta3.boundaryField()[patchBottom][0];

                                                	deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchBottom][0].component(vector::Y));
                                        	}

                                	}
                                	else if ((bottom == false) && (top == true))//boundary at the Top
                                	{
                                        	theta1B = stheta1.primitiveField()[faceybID];
                                        	theta2B = stheta2.primitiveField()[faceybID];
                                        	theta3B = stheta3.primitiveField()[faceybID];

                                        	if (Pstream::myProcNo() == (nproc -1))
                                        	{
                                                	theta1U = stheta1.boundaryField()[patchOutlet][0];
                                                	theta2U = stheta2.boundaryField()[patchOutlet][0];
                                                	theta3U = stheta3.boundaryField()[patchOutlet][0];

                                                	deltay = 2.0 / (mesh.Cf().boundaryField()[patchOutlet][0].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));
                                        	}
                                        	else if (Pstream::myProcNo() == 0)
                                        	{
                                                	theta1U = stheta1.boundaryField()[patchBottomTop][0];
                                                	theta2U = stheta2.boundaryField()[patchBottomTop][0];
                                                	theta3U = stheta3.boundaryField()[patchBottomTop][0];

                                                	deltay = 2.0 / (mesh.Cf().boundaryField()[patchBottomTop][0].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));
                                        	}
                                        	else
                                        	{
                                                	theta1U = stheta1.boundaryField()[patchTop][0];
                                                	theta2U = stheta2.boundaryField()[patchTop][0];
                                                	theta3U = stheta3.boundaryField()[patchTop][0];

                                                	deltay = 2.0 / (mesh.Cf().boundaryField()[patchTop][0].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));
                                        	}

                                	}
                                	else
                                	{
                                        	FatalErrorIn
                                        	(
                                           	"REDIMFDFoam solver"
                                        	)
                                        	<< "no such a cell"
                                        	<< abort(FatalError);
                                	}


                                	gradytheta1U = (theta1U - theta1.primitiveField()[cellID]) * deltay;
                               		gradytheta1B = (theta1.primitiveField()[cellID] - theta1B) * deltay;
                                	gradytheta2U = (theta2U - theta2.primitiveField()[cellID]) * deltay;
                                	gradytheta2B = (theta2.primitiveField()[cellID] - theta2B) * deltay;
                                	gradytheta3U = (theta3U - theta3.primitiveField()[cellID]) * deltay;
                                	gradytheta3B = (theta3.primitiveField()[cellID] - theta3B) * deltay;


                                	theta1UU = (theta1U + theta1.primitiveField()[cellID]) * 0.5;
                                	theta1BB = (theta1B + theta1.primitiveField()[cellID]) * 0.5;
                                	theta2UU = (theta2U + theta2.primitiveField()[cellID]) * 0.5;
                                	theta2BB = (theta2B + theta2.primitiveField()[cellID]) * 0.5;
                                	theta3UU = (theta3U + theta3.primitiveField()[cellID]) * 0.5;
                                	theta3BB = (theta3B + theta3.primitiveField()[cellID]) * 0.5;

                                	DPSy = combustion->getDPS(theta1UU,theta1BB,theta2UU,theta2BB,theta3UU,theta3BB);

                                	scalar diffusion1FG(0.0);
                                	scalar diffusion2FG(0.0);
                                	scalar diffusion3FG(0.0);

                                	for(int m=0; m<22; m++)
                                	{
                                        	scalar divterm =  (gradxtheta1R * DPSx[m][1] - gradxtheta1L * DPSx[m][0]) * deltax
                                                        	+ (gradxtheta2R * DPSx[m][3] - gradxtheta2L * DPSx[m][2]) * deltax
                                                        	+ (gradxtheta3R * DPSx[m][5] - gradxtheta3L * DPSx[m][4]) * deltax
                                                        	+ (gradytheta1U * DPSy[m][1] - gradytheta1B * DPSy[m][0]) * deltay
                                                        	+ (gradytheta2U * DPSy[m][3] - gradytheta2B * DPSy[m][2]) * deltay
                                                        	+ (gradytheta3U * DPSy[m][5] - gradytheta3B * DPSy[m][4]) * deltay;

                                        	diffusion1FG = diffusion1FG + combustion->PM(3*m+1).primitiveField()[cellID] * divterm;
                                        	diffusion2FG = diffusion2FG + combustion->PM(3*m).primitiveField()[cellID] * divterm;
                                        	diffusion3FG = diffusion3FG + combustion->PM(3*m+2).primitiveField()[cellID] * divterm;

                                	}

                                	diffusion1.primitiveFieldRef()[cellID] = diffusion1FG;
                                	diffusion2.primitiveFieldRef()[cellID] = diffusion2FG;
                                	diffusion3.primitiveFieldRef()[cellID] = diffusion3FG;

                        	}

                	}
        	}
     	}
     }
//-----------------------------------------------------------------------------------------------------------
     fvScalarMatrix theta1Eqn
     (
         fvm::ddt(rho, theta1)
       + fvm::div(phi, theta1)
       ==
         rho * combustion->rETA(1)
       + diffusion1
     );

     theta1Eqn.relax();

     theta1Eqn.solve(mesh.solver("theta1"));

     Info << "theta1 max/min from equation:" << max(theta1).value() << ", " << min(theta1).value() << endl;
//-----------------------------------------------------------------------------------------------------------

     fvScalarMatrix theta2Eqn
     (
         fvm::ddt(rho, theta2)
       + fvm::div(phi, theta2)
       ==
         rho * combustion->rETA(0)
       + diffusion2
     );

     theta2Eqn.relax();

     theta2Eqn.solve();
     //theta2Eqn.solve(mesh.solver("theta2"));

     Info << "theta2 max/min from equation:" << max(theta2).value() << ", " << min(theta2).value() << endl;

//-----------------------------------------------------------------------------------------------------------

     fvScalarMatrix theta3Eqn
     (
         fvm::ddt(rho, theta3)
       + fvm::div(phi, theta3)
       ==
         rho * combustion->rETA(2)
       + diffusion3
     );

     theta3Eqn.relax();

     theta3Eqn.solve();
     //theta3Eqn.solve(mesh.solver("theta3"));

     Info << "theta3 max/min from equation:" << max(theta3).value() << ", " << min(theta3).value() << endl;

//-----------------------------------------------------------------------------------------------------------

     combustion->correct();

     theta1.correctBoundaryConditions();
     theta2.correctBoundaryConditions();
     theta3.correctBoundaryConditions();


//     p.correctBoundaryConditions();

     if (gridRes == "fine")
     {
         stheta1 = fvc::interpolate(theta1)();
         stheta2 = fvc::interpolate(theta2)();
         stheta3 = fvc::interpolate(theta3)();
     }

     combustion->correctBoundary();

     thermo.correct();

     Info<< "T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
}

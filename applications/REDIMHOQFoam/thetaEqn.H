{
     volScalarField diffusion1(0.0*divTerm);
     volScalarField diffusion2(0.0*divTerm);

     for(int m=0; m<nSpec+2; m++)
     {

//-Attention: D is negative in original REDIM Table, be careful about the sign. Theta1 and theta2 are exchanged, be careful about the coefficient order.
	if (tableForm == "flip")
        {
		divTerm = fvc::laplacian(combustion->DPS(m+nSpec+2), theta1) + fvc::laplacian(combustion->DPS(m), theta2);

        	diffusion1 = diffusion1 + combustion->PM(2*m+1) * divTerm;

        	diffusion2 = diffusion2 + combustion->PM(2*m) * divTerm;
	}
	else
	{
        	divTerm = fvc::laplacian(combustion->DPS(m), theta1) + fvc::laplacian(combustion->DPS(m+nSpec+2), theta2);

        	diffusion1 = diffusion1 + combustion->PM(2*m) * divTerm;

        	diffusion2 = diffusion2 + combustion->PM(2*m+1) * divTerm;
	}

     }

//-Special treatment to the first cell layer near the wall
     forAll (theta1.boundaryField(),patchI)
     {

        if (mesh.boundary()[patchI].type() == "wall")
        {

                const polyPatch& cPatch = mesh.boundaryMesh()[patchI];
                const labelUList& faceCells = cPatch.faceCells();

                if (faceCells.size()>0)
                {
                        scalar deltax(mesh.deltaCoeffs().boundaryField()[patchI][0]);
                        scalar theta1FS, theta2FS;
			scalar theta1R, theta1L, theta2R, theta2L;
                        scalar gradxtheta1R, gradxtheta1L, gradxtheta2R, gradxtheta2L; 
			List<scalarList> DPS;

			label cellID = faceCells[0];
			const labelUList& faces = mesh.cells()[cellID];
			label facexID=faces[0];

			theta1FS = stheta1.primitiveField()[facexID];
			theta2FS = stheta2.primitiveField()[facexID];
			gradxtheta1R = (theta1FS - theta1.primitiveField()[cellID])*deltax;
			gradxtheta1L = (theta1.primitiveField()[cellID] - theta1.boundaryField()[patchI][0])*deltax;
			gradxtheta2R = (theta2FS - theta2.primitiveField()[cellID])*deltax;
			gradxtheta2L = (theta2.primitiveField()[cellID] - theta2.boundaryField()[patchI][0])*deltax;

			theta1R = (theta1.primitiveField()[cellID] + theta1FS) * 0.5;
			theta1L = (theta1.boundaryField()[patchI][0] + theta1.primitiveField()[cellID]) * 0.5;
			theta2R = (theta2.primitiveField()[cellID] + theta2FS) * 0.5;
			theta2L = (theta2.boundaryField()[patchI][0] + theta2.primitiveField()[cellID]) * 0.5;

			DPS = combustion->getDPS(theta1R,theta1L,theta2R,theta2L);

			scalar diffusion1FG(0.0);
			scalar diffusion2FG(0.0);

			for(int m=0; m<nSpec+2; m++)
			{
				scalar divterm =  (gradxtheta1R * DPS[m][1] - gradxtheta1L * DPS[m][0])*deltax
						+ (gradxtheta2R * DPS[m][3] - gradxtheta2L * DPS[m][2])*deltax;

				diffusion1FG = diffusion1FG + combustion->PM(2*m+1).internalField()[cellID] * divterm;
				diffusion2FG = diffusion2FG + combustion->PM(2*m).internalField()[cellID] * divterm;
			}

			diffusion1.primitiveFieldRef()[cellID] = diffusion1FG;
			diffusion2.primitiveFieldRef()[cellID] = diffusion2FG;

               }
        }
     }
//-----------------------------------------------------------------------------------------------------------
     int mm = 0;
     int nn = 0;

     if (tableForm == "flip")
     {
	mm = 1;
     }
     else
     {
	nn = 1;
     }

     fvScalarMatrix theta1Eqn
     (
         fvm::ddt(rho, theta1)
       + fvm::div(phi, theta1)
       ==
         rho * combustion->rETA(mm)
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
         rho * combustion->rETA(nn)
       + diffusion2
     );

     theta2Eqn.relax();

     theta2Eqn.solve();

     Info << "theta2 max/min from equation:" << max(theta2).value() << ", " << min(theta2).value() << endl;

//-----------------------------------------------------------------------------------------------------------

     combustion->correct();

     theta1.correctBoundaryConditions();
     theta2.correctBoundaryConditions();


//     p.correctBoundaryConditions();

     if (gridRes == "fine")
     {
         stheta1 = fvc::interpolate(theta1)();
         stheta2 = fvc::interpolate(theta2)();
     }

     combustion->correctBoundary();

     thermo.correct();

     Info<< "T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
}

{
     volScalarField diffusion1(0.0*divTerm);
     volScalarField diffusion2(0.0*divTerm);

     for(int m=0; m<29; m++)
     {

//-Attention: D is negative in original REDIM Table, be careful about the sign. Theta1 and theta2 are exchanged, be careful about the coefficient order.

	divTerm = fvc::laplacian(combustion->DPSS(m+29), theta1) + fvc::laplacian(combustion->DPSS(m), theta2);

        diffusion1 = diffusion1 + combustion->PM(2*m+1) * divTerm;

        diffusion2 = diffusion2 + combustion->PM(2*m) * divTerm;

     }

//-Special treatment to the first cell layer near the wall

//     surfaceScalarField stheta1 = fvc::interpolate(theta1)();
//     surfaceScalarField stheta2 = fvc::interpolate(theta2)();

     //const labelUList& neighbour = mesh.neighbour();

     forAll(theta1.boundaryField(),patchI)
     {

     	if (mesh.boundary()[patchI].type() == "wall")
     	{
     		const polyPatch& cPatch = mesh.boundaryMesh()[patchI];
		const labelUList& faceCells = cPatch.faceCells();

		//-to make sure that wall is included in the processor
		if (faceCells.size()>0)
		{
     			scalar deltax(mesh.deltaCoeffs().boundaryField()[patchI][0]);
			scalar deltay(1.0);
     			scalar theta1FS, theta2FS;
     			scalar theta1R, theta1L, theta2R, theta2L;
			scalar theta1U, theta1B, theta2U, theta2B;
			scalar theta1UU, theta1BB, theta2UU, theta2BB;
     			scalar gradxtheta1R, gradxtheta1L, gradxtheta2R, gradxtheta2L; 
     			scalar gradytheta1U, gradytheta1B, gradytheta2U, gradytheta2B;
			List<scalarList> DPSx, DPSy;
			label  patchBottom(7), patchTop(8); //-have to predefine the ID of bottom and top patch
			label  patchOutlet(4), patchInlet(0);
			label  patchTopBottom(7), patchBottomTop(7);

	   		//-loop over cells near the wall
	   		forAll(faceCells,i)
	   		{
				//-determine cell ID
				label cellID = faceCells[i];
				//-determine all face IDs of the cell
	      			const labelUList& faces = mesh.cells()[cellID];
			
				//- to determine whether the cell includes face that on the boundary 
			      	//-Find the face near the wall (in X direction), and the face at the top and bottom of the cell (in Y direction)
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

				//label neighbourID = neighbour[faces[0]];

				//-Calculate gradient of theta1 and theta2 in y direction
				
				//-Calculation in x direction
				theta1FS = stheta1.primitiveField()[facexID];//theta1FS = (theta1.internalField()[cellID] + theta1.internalField()[neighbourID]) * 0.5;
				theta2FS = stheta2.primitiveField()[facexID];//theta1FS = (theta2.internalField()[cellID] + theta2.internalField()[neighbourID]) * 0.5;

				//-Calculate gradient of theta1 and theta2 in x direction
				gradxtheta1R = (theta1FS - theta1.primitiveField()[cellID]) * deltax;
				gradxtheta1L = (theta1.primitiveField()[cellID] - theta1.boundaryField()[patchI][i]) * deltax;
				gradxtheta2R = (theta2FS - theta2.primitiveField()[cellID]) * deltax;
				gradxtheta2L = (theta2.primitiveField()[cellID] - theta2.boundaryField()[patchI][i]) * deltax;

				//-theta1 and theta2 value on the right and left side of the first grid point
				theta1R = (theta1.primitiveField()[cellID] + theta1FS) * 0.5;
				theta1L = (theta1.boundaryField()[patchI][i] + theta1.primitiveField()[cellID]) * 0.5;
				theta2R = (theta2.primitiveField()[cellID] + theta2FS) * 0.5;
				theta2L = (theta2.boundaryField()[patchI][i] + theta2.primitiveField()[cellID]) * 0.5;

				DPSx = combustion->getDPS(theta1R,theta1L,theta2R,theta2L);

				//-Calculation in y direction 
				//- 1: the cell is totally within the inner domain; 2. the bottom of the cell is on the boundary (process or physical); 3. the top of the cell is on the boundary
				if ((bottom == false) && (top == false)) //-internal cell without boundary faces
				{	
					theta1U = stheta1.primitiveField()[faceyuID]; 
					theta1B = stheta1.primitiveField()[faceybID];
					theta2U = stheta2.primitiveField()[faceyuID];
					theta2B = stheta2.primitiveField()[faceybID];

					deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));

				}
				else if ((bottom == true) && (top == false)) //boundary at the Bottom
				{		
					theta1U = stheta1.primitiveField()[faceyuID]; 
					theta2U = stheta2.primitiveField()[faceyuID];
						
					if (Pstream::myProcNo() == 0)
					{	
						theta1B = stheta1.boundaryField()[patchInlet][0];
						theta2B = stheta2.boundaryField()[patchInlet][0];

						deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchInlet][0].component(vector::Y));
					}
					else if (Pstream::myProcNo() == (nproc - 1))
					{
						theta1B = stheta1.boundaryField()[patchTopBottom][0];
						theta2B = stheta2.boundaryField()[patchTopBottom][0];

						deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchTopBottom][0].component(vector::Y));	
					}
					else
					{
						theta1B = stheta1.boundaryField()[patchBottom][0];
						theta2B = stheta2.boundaryField()[patchBottom][0];

						deltay = 2.0 / (mesh.Cf().primitiveField()[faceyuID].component(vector::Y) - mesh.Cf().boundaryField()[patchBottom][0].component(vector::Y));
					}
					
				}
				else if ((bottom == false) && (top == true))//boundary at the Top
				{
					theta1B = stheta1.primitiveField()[faceybID];
					theta2B = stheta2.primitiveField()[faceybID];

					if (Pstream::myProcNo() == (nproc -1))
					{
						theta1U = stheta1.boundaryField()[patchOutlet][0]; 
						theta2U = stheta2.boundaryField()[patchOutlet][0];

						deltay = 2.0 / (mesh.Cf().boundaryField()[patchOutlet][0].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));
					}
					else if (Pstream::myProcNo() == 0)
					{
						theta1U = stheta1.boundaryField()[patchBottomTop][0]; 
						theta2U = stheta2.boundaryField()[patchBottomTop][0];

						deltay = 2.0 / (mesh.Cf().boundaryField()[patchBottomTop][0].component(vector::Y) - mesh.Cf().primitiveField()[faceybID].component(vector::Y));
					}
					else
					{
						theta1U = stheta1.boundaryField()[patchTop][0]; 
						theta2U = stheta2.boundaryField()[patchTop][0];

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


				theta1UU = (theta1U + theta1.primitiveField()[cellID]) * 0.5;
				theta1BB = (theta1B + theta1.primitiveField()[cellID]) * 0.5;
				theta2UU = (theta2U + theta2.primitiveField()[cellID]) * 0.5;
                                theta2BB = (theta2B + theta2.primitiveField()[cellID]) * 0.5;

				DPSy = combustion->getDPS(theta1UU,theta1BB,theta2UU,theta2BB);
				
				scalar diffusion1FG(0.0);
				scalar diffusion2FG(0.0);

				for(int m=0; m<29; m++)
				{
					scalar divterm =  (gradxtheta1R * DPSx[m][1] - gradxtheta1L * DPSx[m][0]) * deltax
						        + (gradxtheta2R * DPSx[m][3] - gradxtheta2L * DPSx[m][2]) * deltax
							+ (gradytheta1U * DPSy[m][1] - gradytheta1B * DPSy[m][0]) * deltay
							+ (gradytheta2U * DPSy[m][3] - gradytheta2B * DPSy[m][2]) * deltay;
					
					diffusion1FG = diffusion1FG + combustion->PM(2*m+1).primitiveField()[cellID] * divterm;
					diffusion2FG = diffusion2FG + combustion->PM(2*m).primitiveField()[cellID] * divterm;
				}

				diffusion1.primitiveFieldRef()[cellID] = diffusion1FG;
				diffusion2.primitiveFieldRef()[cellID] = diffusion2FG;

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
         rho * combustion->rETA(2)
       + diffusion1
     );

     theta1Eqn.relax();

     theta1Eqn.solve(mesh.solver("theta1"));

     Info << "theta1 max/min from equation:" << max(theta1).value() << ", " << min(theta1).value() << endl;

     theta1.max(1.0);
     theta1.min(101.0);

     Info << "theta1 max/min after correction:" << max(theta1).value() << ", " << min(theta1).value() << endl;
//-----------------------------------------------------------------------------------------------------------

     fvScalarMatrix theta2Eqn
     (
         fvm::ddt(rho, theta2)
       + fvm::div(phi, theta2)
       ==
         rho * combustion->rETA(1)
       + diffusion2
     );

     theta2Eqn.relax();

     theta2Eqn.solve();

     Info << "theta2 max/min from equation:" << max(theta2).value() << ", " << min(theta2).value() << endl;

     theta2.max(1.0);
     theta2.min(97.0);

     Info << "theta2 max/min after correction:" << max(theta2).value() << ", " << min(theta2).value() << endl;
//-----------------------------------------------------------------------------------------------------------

     combustion->correctPV();

     theta1.correctBoundaryConditions();
     theta2.correctBoundaryConditions();

     p.correctBoundaryConditions();

     stheta1 = fvc::interpolate(theta1)();
     stheta2 = fvc::interpolate(theta2)();

     combustion->correct();

     Info<< "T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
}

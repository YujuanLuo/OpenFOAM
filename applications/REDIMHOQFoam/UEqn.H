// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
    rho*g
);

/*
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
    rho * g
);
fvVectorMatrix& UEqn = tUEqn.ref();
*/

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}

Info << "Ux max/min: " << max(U.component(vector::X)).value() << " " << min(U.component(vector::X)).value() << endl;
Info << "Uy max/min: " << max(U.component(vector::Y)).value() << " " << min(U.component(vector::Y)).value() << endl;

fvScalarMatrix PVEqn
(


	(
	  fvm::ddt(rho, PV)
	+ fvm::div(phi, PV)
 	- fvm::laplacian(turbulence->DfEff(), PV)
	==
	  //thermo.PVs()
	  meanPVs
	+ sparkPVs       	  
	  
	)
);

PVEqn.relax();
PVEqn.solve(mesh.solver("PV"));


forAll(PV,icell)
{
	PV[icell] = min(PV[icell],maxProgVar.value());
	PV[icell] = max(PV[icell],0.0);

}

forAll(PV.boundaryField(), patchi)
{
	forAll(PV.boundaryField()[patchi],facei)
	{
		PV.boundaryField()[patchi][facei] = min(PV.boundaryField()[patchi][facei],maxProgVar.value());
		PV.boundaryField()[patchi][facei] = max(PV.boundaryField()[patchi][facei],0.0);
	}
}



Info <<"max/min meanf " << max(f).value() << " " << min(f).value() <<endl; 
Info <<"max/min meanPV " << max(PV).value() << " " << min(PV).value() <<endl; 

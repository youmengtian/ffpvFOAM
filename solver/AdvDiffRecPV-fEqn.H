tmp<fv::convectionScheme<scalar> > mvConvectionPVSF
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,PVSF)")
    )
);
tmp<fv::convectionScheme<scalar> > mvConvectionMFSF
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,fSF)")
    )
);

/*
if(turbulenceModelType == "LESModel")
{

    volScalarField deltaLES
    (
        IOobject
        (
            "deltaLES",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("deltaLES", dimLength, SMALL)
    );

    if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
    {
        const volScalarField& delta =
            mesh.objectRegistry::lookupObject<volScalarField>("delta");

        deltaLES.internalField() = delta.internalField();
    }

	const volScalarField gammaEff = turbulence->muEff();	

	forAll(Omega_t,icell)
	{
		Omega_t[icell] = gammaEff[icell]/(rho[icell]*deltaLES[icell]*deltaLES[icell]);
	}

							
	const fvPatchList& patches = mesh.boundary();
	forAll(patches, patchi)
	{
		const fvPatch& currPatch = patches[patchi];
		forAll(currPatch,facei)
		{
			label faceCelli = currPatch.faceCells()[facei];
			Omega_t.boundaryField()[patchi][facei] = Omega_t[faceCelli];
		}
	}

Info << " Omega_t (nuEff/delatLES^2) max/min: "<<  max(Omega_t).value() <<" " << min(Omega_t).value() << endl;

	
}
else
{*/

	dimensionedScalar smallK("smallK",dimensionSet(0, 2, -2, 0, 0, 0, 0), 1e-6);

	const volScalarField e_turb = turbulence->epsilon();

	Info << " e_turb max/min: "<<  max(e_turb).value() <<" " << min(e_turb).value() << endl;

	const volScalarField k_turb = turbulence->k();

	Info << " k_turb max/min: "<<  max(k_turb).value() <<" " << min(k_turb).value() << endl;

	forAll(Omega_t, icell)
	{ 
		Omega_t[icell] = max(e_turb[icell],0.0)/max(k_turb[icell],SMALL);
	}

	forAll(Omega_t.boundaryField(), patchi)
	{
		forAll(Omega_t.boundaryField()[patchi],facei)
		{
			Omega_t.boundaryField()[patchi][facei] = max(e_turb.boundaryField()[patchi][facei],0)/max(k_turb.boundaryField()[patchi][facei],SMALL);
		}
	}

Info << " Omega_t (k/epsilon) max/min: "<<  max(Omega_t).value() <<" " << min(Omega_t).value() << endl;

//}




label tmpDebug = solverPerformance::debug; 
solverPerformance::debug = 0;
forAll(sfPV,nsf)
{	
//  for PV
    volScalarField& PVSF = sfPV[nsf];
    volScalarField& fSF = sff[nsf];

	volScalarField curPVs
	(
		IOobject
		(
			"curPVs",
			 mesh.time().timeName(),
			 mesh,
			 IOobject::NO_READ,
			 IOobject::NO_WRITE
		),
        mesh,
        dimensionedScalar("zero", meanPVs.dimensions(), 0.0)
	);

	if(!nonReacting)
	{
		curPVs = sfPVs[nsf];
	}
    if(setReactingBox)
    {
	   curPVs *= inReactingBox;
	}


	fvScalarMatrix PVEqn
	(

		(
		  fvm::ddt(rho, PVSF)
		  + mvConvectionPVSF->fvmDiv(phi,PVSF)
		//+ fvm::div(phi, PVSF)
		- fvm::laplacian(turbulence->DfEff(), PVSF)
		==
		  curPVs
		- fvm::Sp(0.5*Cphi*rho*Omega_t,PVSF)
		+ 0.5*Cphi*rho*Omega_t*meanPV
		+ sparkPVs       
		)
	);
	PVEqn.relax();
	PVEqn.solve(mesh.solver("PVSF"));


//  for f
	fvScalarMatrix fEqn
	(
		(
		  fvm::ddt(rho, fSF)
	    + mvConvectionMFSF->fvmDiv(phi,fSF)
//		+ fvm::div(phi, fSF)
        - fvm::laplacian(turbulence->DfEff(), fSF)
		==
		- fvm::Sp(0.5*Cphi*rho*Omega_t,fSF)
		+ 0.5*Cphi*rho*Omega_t*meanf 
        )
    );

	fEqn.relax();
    fEqn.solve(mesh.solver("fSF"));

	
}

solverPerformance::debug = tmpDebug; 

Info <<"For " << nStochasticFields << " stochastic fields: transport of f and PV done " <<endl; 
Info << " meanPVs max/min: "<<  max(meanPVs).value() <<" " << min(meanPVs).value() << endl;


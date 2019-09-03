forAll(PV, icell)
{ 
	meanPVs[icell]  = 0;
}

forAll(PV.boundaryField(), patchi)
{
	forAll(PV.boundaryField()[patchi],facei)
	{
		meanPVs.boundaryField()[patchi][facei] = 0;
	}
}


dimensionedScalar smallK("smallK",dimensionSet(0, 2, -2, 0, 0, 0, 0), 1e-6);
volScalarField Omega_t
(
	IOobject
	(
		"Omega_t",
		 mesh.time().timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
	),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)
);

const volScalarField& e_turb = turbulence->epsilon();
const volScalarField& k_turb = turbulence->k();

forAll(Omega_t, icell)
{ 
	Omega_t[icell]  = max(e_turb[icell],0.0)/max(k_turb[icell],SMALL);
}

forAll(Omega_t.boundaryField(), patchi)
{
	forAll(Omega_t.boundaryField()[patchi],facei)
	{
		Omega_t.boundaryField()[patchi][facei] = max(e_turb.boundaryField()[patchi][facei],0)/max(k_turb.boundaryField()[patchi][facei],SMALL);
	}
}


//
if(dichotomicVectorsSym)
{
#include "DichotomicVectorsSym.H"
}
else
{
#include "DichotomicVectors.H"
}

label nLastField = -1;
scalar sLastField = -1;
//
if (Pstream::myProcNo() == 0)
{
    //scalar slastField = rndGen.scalar01();
    sLastField =rndGen.scalar01()*nStochasticFields;
    for(label nsf = 0; nsf < nStochasticFields; nsf++)
    {
		if(sLastField >= nsf && sLastField < nsf+1)
		{
			nLastField = nsf;
		}
	}
    
   for(label iproc=1; iproc<Pstream::nProcs(); iproc++)
   {   
    OPstream sendStream(Pstream::blocking,iproc);
    sendStream   << nLastField;
   }
 
}
else
{
       IPstream reciveStream(Pstream::blocking,0);
       reciveStream >> nLastField;	
}

vector vecSumm (0,0,0);
for(label nsf = 0; nsf < nStochasticFields-1; nsf++)
{
	if(nsf != nLastField)
	{
		vecSumm += eta_sf[nsf]; 
	}
}
eta_sf[nLastField] = -vecSumm;
Info << "field " <<  nLastField <<" sLastField " << sLastField << " is selected to sum and balance the wiener term "<<" eta_sf[nLastField] = " << eta_sf[nLastField] << endl;
//

label tmpDebug = solverPerformance::debug; 
solverPerformance::debug = 0;
forAll(sfPV,nsf)
{	
    volScalarField& PVSF = sfPV[nsf];
    volScalarField& fSF = sff[nsf];

	volScalarField curPV
	(
		IOobject
		(
			"curPV",
			 mesh.time().timeName(),
			 mesh,
			 IOobject::NO_READ,
			 IOobject::NO_WRITE
		),
        PVSF
	);

	volScalarField curf
	(
		IOobject
		(
			"curf",
			 mesh.time().timeName(),
			 mesh,
			 IOobject::NO_READ,
			 IOobject::NO_WRITE
		),
        fSF
	);

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

////
	const dimensionedVector dW_sf = eta_sf[nsf]*sqrt(runTime.deltaT()); 
//	Info << "code is here 0" << endl;

	volScalarField   wienerSourcePV
	(      
		IOobject
		(
			"wienerSourcePV_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", PV.dimensions()*dimDensity/dimTime, 0.0),
		zeroGradientFvPatchScalarField::typeName
	);
	if(writeWienerSources)
	{
		wienerSourcePV.writeOpt() = IOobject::AUTO_WRITE;
	}

	volScalarField   wienerSourceMixFrc
	(      
		IOobject
		(
			"wienerSourceMixFrc_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", PV.dimensions()*dimDensity/dimTime, 0.0),
		zeroGradientFvPatchScalarField::typeName
	);
	if(writeWienerSources)
	{
		wienerSourceMixFrc.writeOpt() = IOobject::AUTO_WRITE;
	}

    volVectorField gradPV
	(      
		IOobject
		(
			"gradPV_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero",  PV.dimensions()/dimLength, vector (0, 0, 0) ),
		zeroGradientFvPatchScalarField::typeName
	);		
	gradPV = fvc::grad(curPV);
	gradPV.correctBoundaryConditions();
    
    volVectorField gradf
	(      
		IOobject
		(
			"gradf_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero",  f.dimensions()/dimLength, vector (0, 0, 0) ),
		zeroGradientFvPatchScalarField::typeName
	);		
	gradf = fvc::grad(curf);
	gradf.correctBoundaryConditions();
    
	const volScalarField gamma = turbulence->DfEff();

//	Info << "code is here 1" << endl;


	forAll(wienerSourcePV, icell)
	{
		const scalar gamTerm = max(double(2)*gamma[icell],SMALL);
		const scalar gradTermPV = (gradPV[icell] & dW_sf).value();
		const scalar gradTermf  = (gradf[icell] & dW_sf).value();

		wienerSourcePV[icell] = gradTermPV * rho[icell]*Foam::sqrt(gamTerm);
		wienerSourceMixFrc[icell] = gradTermf * rho[icell]*Foam::sqrt(gamTerm);

	}
//	Info << "code is here 02" << endl;
 	 forAll(wienerSourcePV.boundaryField(), patchi)
	 {
		forAll(wienerSourcePV.boundaryField()[patchi],facei)
		{
			scalar gamTerm = max(double(2)*gamma.boundaryField()[patchi][facei],SMALL);
			const vector vecPVi = gradPV.boundaryField()[patchi][facei];
			const vector vecfi = gradf.boundaryField()[patchi][facei];
			const scalar gradTermPV = (vecPVi & dW_sf).value();
			const scalar gradTermf  = (vecfi & dW_sf).value();

			wienerSourcePV.boundaryField()[patchi][facei] = gradTermPV * rho.boundaryField()[patchi][facei]*Foam::sqrt(gamTerm);
			wienerSourceMixFrc.boundaryField()[patchi][facei] = gradTermf * rho.boundaryField()[patchi][facei]*Foam::sqrt(gamTerm);

		}
			
	}

//	Info << "code is here 2" << endl;

//  for PV


	
	if(!nonReacting)
	{
		curPVs = combustion->getPVs(fSF,PVSF)();
	}
    if(setReactingBox)
    {
	   curPVs *= inReactingBox;
	}



// for PV
//Info << "before trans curPV max/min: "<< nsf <<" "<<  max(curPV).value() <<" " << min(curPV).value() << endl;
	
	fvScalarMatrix PVEqn
	(

		(
		  fvm::ddt(rho, curPV)
		+ fvm::div(phi, curPV)
		- fvm::laplacian(turbulence->DfEff(), curPV)
		==
		  curPVs
		- fvm::Sp(0.5*Cphi*rho*Omega_t,curPV)
		+ wienerSourcePV
		+ 0.5*Cphi*rho*Omega_t*meanPV
		+ sparkPVs       
		)
	);
	PVEqn.relax();
	PVEqn.solve(mesh.solver("curPV"));
//Info << "after trans curPV max/min: "<< nsf <<" "<<  max(curPV).value() <<" " << min(curPV).value() << endl;
	
//  for f
	fvScalarMatrix fEqn
	(
		(
		  fvm::ddt(rho, curf)
		+ fvm::div(phi, curf)
        	- fvm::laplacian(turbulence->DfEff(), curf)
		==
		- fvm::Sp(0.5*Cphi*rho*Omega_t,curf)
		+ wienerSourceMixFrc
		+ 0.5*Cphi*rho*Omega_t*meanf 
        )
    );

	fEqn.relax();
    fEqn.solve(mesh.solver("curf"));

//Info << "after trans curf max/min: "<< nsf <<" "<<  max(curf).value() <<" " << min(curf).value() << endl;

	forAll(PV, icell)
	{ 
		meanPVs[icell]  += curPVs[icell]/nStochasticFields;
		PVSF[icell]     = curPV[icell];
		fSF[icell] 		= curf[icell];
		
	}

	forAll(PV.boundaryField(), patchi)
	{
		forAll(PV.boundaryField()[patchi],facei)
		{
			meanPVs.boundaryField()[patchi][facei] += curPVs.boundaryField()[patchi][facei]/nStochasticFields;
			PVSF.boundaryField()[patchi][facei] = curPV.boundaryField()[patchi][facei];
			fSF.boundaryField()[patchi][facei] = curf.boundaryField()[patchi][facei];
		}
	}
	

}

solverPerformance::debug = tmpDebug; 

Info <<"For " << nStochasticFields << " stochastic fields: transport of f and PV done " <<endl; 
Info << " meanPVs max/min: "<<  max(meanPVs).value() <<" " << min(meanPVs).value() << endl;


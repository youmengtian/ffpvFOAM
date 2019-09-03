volScalarField wFactorPV
(
	IOobject
	(
		"wFactorPV",
		 mesh.time().timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
	),
    mesh,
    dimensionedScalar("wFactorPV", dimless, 1.0)
);

volScalarField wFactorf
(
	IOobject
	(
		"wFactorf",
		 mesh.time().timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
	),
    mesh,
    dimensionedScalar("wFactorf", dimless, 1.0)
);

volScalarField Tstar
(
	IOobject
    (
		"T_start",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     0*T
);
PtrList <volScalarField> Ystar(Y.size());
forAll(Ystar,i)
{
		Ystar.set
        (
			i,
            new volScalarField
            (
				IOobject
                (
					Y[i].name()+"_start",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*Y[i]
             )
        );
}

forAll(PV, icell)
{
	scalar& fCelli = f[icell];
	scalar& meanfCelli = meanf[icell];
	scalar& PVCelli = PV[icell];
	scalar& meanPVCelli = meanPV[icell];
	
	wFactorPV[icell] = Foam::tanh(C_wFactor.value()*mag((meanPVCelli-PVCelli)/(meanPVCelli+SMALL) ) ); 
	wFactorf[icell] = Foam::tanh(C_wFactor.value()*mag((meanfCelli-fCelli)/(meanfCelli+SMALL) ) ); 

	
	fCelli = (1. - wFactorf[icell])*meanfCelli + wFactorf[icell]*fCelli;
	PVCelli = (1. - wFactorPV[icell])*meanPVCelli + wFactorPV[icell]*PVCelli;

}

forAll(PV.boundaryField(), patchi)
{
	forAll(PV.boundaryField()[patchi],facei)
	{

		scalar& ffacei = f.boundaryField()[patchi][facei];
		scalar& meanffacei = meanf.boundaryField()[patchi][facei];
		scalar& PVfacei = PV.boundaryField()[patchi][facei];
		scalar& meanPVfacei = meanPV.boundaryField()[patchi][facei];

		wFactorPV.boundaryField()[patchi][facei] = 
		Foam::tanh(C_wFactor.value()*mag((meanPVfacei-PVfacei)/(meanPVfacei+SMALL) ) ); 

		wFactorf.boundaryField()[patchi][facei] = 
		Foam::tanh(C_wFactor.value()*mag( (meanffacei-ffacei)/(meanffacei+SMALL) ) ); 

		ffacei  = (1. - wFactorf.boundaryField()[patchi][facei])*meanffacei + wFactorf.boundaryField()[patchi][facei]*ffacei;
		PVfacei = (1. - wFactorPV.boundaryField()[patchi][facei])*meanPVfacei + wFactorPV.boundaryField()[patchi][facei]*PVfacei;
	}
			
}

combustion->correctCurrentFields(f,PV);
Tstar = combustion->getT();
forAll(Y,i)
{
	Ystar[i] = combustion->getYi(i);
}
forAll(T, icell)
{
	scalar wFactor = Foam::tanh(C_wFactor.value()*mag((T[icell]-Tstar[icell])/(T[icell]+SMALL) ) ); 
	T[icell] = (1. - wFactor)*T[icell] + wFactor*Tstar[icell];
	forAll(Y,i)
	{
		scalar wFactor = Foam::tanh(C_wFactor.value()*mag((Y[i][icell]-Ystar[i][icell])/(Y[i][icell]+SMALL) ) ); 
		Y[i][icell] = (1. - wFactor)*Y[i][icell] + wFactor*Ystar[i][icell];
	}

}

forAll(T.boundaryField(), patchi)
{
	forAll(T.boundaryField()[patchi],facei)
	{
		scalar& Tfacei = T.boundaryField()[patchi][facei];
		scalar& sTfacei = Tstar.boundaryField()[patchi][facei];
		
		scalar wFactor = Foam::tanh(C_wFactor.value()*mag((Tfacei-sTfacei)/(Tfacei+SMALL) ) ); 
		Tfacei  = (1. - wFactor)*Tfacei + wFactor*sTfacei;

		forAll(Y,i)
		{		

			scalar& Yifacei = Y[i].boundaryField()[patchi][facei];
			scalar& sYifacei = Ystar[i].boundaryField()[patchi][facei];
			scalar wFactor = Foam::tanh(C_wFactor.value()*mag((Yifacei-sYifacei)/(Yifacei+SMALL) ) ); 
			Yifacei = (1. - wFactor)*Yifacei + wFactor*sYifacei;
		}
	}
			
}

//thermo.correct();

if(runTime.write())
{
	wFactorf.write();
	wFactorPV.write();

}

Info <<"For " << nStochasticFields << " stochastic fields: remove stochastic noise done " <<endl; 
Info <<"max/min wFactorPV " << max(wFactorPV).value() << " " << min(wFactorPV).value() <<endl; 
Info <<"max/min wFactorf " << max(wFactorf).value() << " " << min(wFactorf).value() <<endl; 


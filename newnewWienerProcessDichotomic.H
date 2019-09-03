#include "DichotomicVectors.H"

forAll(WienerTermPV_t, icell)
{ 
	WienerTermPV_t[icell] = 0;
	WienerTermMF_t[icell] = 0;
//	WienerTermH_t[icell]  = 0; //deleted by jianhong
}

forAll(WienerTermPV_t.boundaryField(), patchi)
{
	forAll(WienerTermPV_t.boundaryField()[patchi],facei)
	{
		WienerTermPV_t.boundaryField()[patchi][facei] = 0;
		WienerTermMF_t.boundaryField()[patchi][facei] = 0;
	}
}

// for PV		
forAll(sfPV,nsf)
{
	const dimensionedVector dW_sf = eta_sf[nsf]*sqrt(runTime.deltaT()); 

    volScalarField   WienerTerm
    (      
      IOobject
      (
       "WienerTerm_PV_"+Foam::name(nsf),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
	  mesh,
	  dimensionedScalar("zero", PV.dimensions(), 0.0),
	  zeroGradientFvPatchScalarField::typeName
    );
    if(writeWienerSources)
    {
		WienerTerm.writeOpt() = IOobject::AUTO_WRITE;
    }

    volScalarField & PVSF = sfPV[nsf];
    volScalarField   varWiener = PVSF;
    volScalarField   gamma = turbulence->muEff()/rho;
    
	volVectorField gradTerm
	(      
		IOobject
		(
			"gradTerm_PV_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero",  PV.dimensions()/dimLength, vector (0, 0, 0) ),
		zeroGradientFvPatchScalarField::typeName
	);		
	gradTerm = fvc::grad(varWiener);
	gradTerm.correctBoundaryConditions();
		
	volScalarField  varTerm  = gradTerm & dW_sf; 

     forAll(varTerm, icell)
     {
		varTerm[icell]  = gradTerm[icell] & dW_sf.value();

		dimensionedScalar igamma("igamma",dimless,2*gamma[icell]);
		WienerTerm[icell] = sqrt(igamma).value()*varTerm[icell];

     }

 	 forAll(varTerm.boundaryField(), patchi)
	 {
		forAll(varTerm.boundaryField()[patchi],facei)
		{
			vector vecTerm = gradTerm.boundaryField()[patchi][facei];
			scalar sclTerm = gamma.boundaryField()[patchi][facei];

			varTerm.boundaryField()[patchi][facei] = vecTerm & dW_sf.value();

			dimensionedScalar igamma("igamma",dimless,2*sclTerm);
			WienerTerm.boundaryField()[patchi][facei] = sqrt(igamma).value()*varTerm.boundaryField()[patchi][facei];
		}
		
	 }

//make sure the Wiener Term for fixedValue boundaries must be zero
		const fvPatchList& patches = mesh.boundary();
		const volScalarField::GeometricBoundaryField& PVSFbf = PVSF.boundaryField();
		wordList PVSFBoundaryTypes(PVSF.boundaryField().size());
		PVSFBoundaryTypes = PVSFbf.types();
		forAll(patches, patchi)
		{

			if(PVSFBoundaryTypes[patchi] == "fixedValue")
			{
				const fvPatch& currPatch = patches[patchi];			
//				Info << currPatch.name() <<" " << PVSFBoundaryTypes[patchi] << endl;

				forAll(currPatch,facei)
				{
					WienerTerm.boundaryField()[patchi][facei] = 0;
					label faceCelli = currPatch.faceCells()[facei];
					WienerTerm[faceCelli] = 0;
					
				}
			}
		}
//


	forAll(PVSF,icell)
	{
		
	    if(setReactingBox)
		{
			WienerTerm[icell] *= inReactingBox[icell];
		}		
		scalar PVSFicell = PVSF[icell];
		
		PVSF[icell] += WienerTerm[icell];
		if(!balanceWienerTerm)
		{
			PVSF[icell] = min(PVSF[icell],maxProgVar.value());
			PVSF[icell] = max(PVSF[icell],0.0);
			WienerTerm[icell] = PVSF[icell] - PVSFicell;
		}
		WienerTermPV_t[icell] += WienerTerm[icell];


	}

	forAll(PVSF.boundaryField(), patchi)
	{
		forAll(PVSF.boundaryField()[patchi],facei)
		{

			if(setReactingBox)
			{
				WienerTerm.boundaryField()[patchi][facei] *= inReactingBox.boundaryField()[patchi][facei];
			}	
			
			scalar PVSFiface = PVSF.boundaryField()[patchi][facei];
			
			PVSF.boundaryField()[patchi][facei] += WienerTerm.boundaryField()[patchi][facei];
			if(!balanceWienerTerm)
			{
				PVSF.boundaryField()[patchi][facei] = min(PVSF.boundaryField()[patchi][facei],maxProgVar.value());
				PVSF.boundaryField()[patchi][facei] = max(PVSF.boundaryField()[patchi][facei],0.0);
				WienerTerm.boundaryField()[patchi][facei] = PVSF.boundaryField()[patchi][facei] - PVSFiface;
			}

			WienerTermPV_t.boundaryField()[patchi][facei] += WienerTerm.boundaryField()[patchi][facei];


		}
	}

	if(runTime.write() && writeWienerSources)
	{
		WienerTerm.write();
		gradTerm.write();
	}
	
}

//for f
forAll(sff,nsf)
{
	const dimensionedVector dW_sf = eta_sf[nsf]*sqrt(runTime.deltaT()); 

    volScalarField   WienerTerm
    (      
      IOobject
      (
       "WienerTerm_f_"+Foam::name(nsf),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
	  mesh,
	  dimensionedScalar("zero", f.dimensions(), 0.0),
	  zeroGradientFvPatchScalarField::typeName
    );
    if(writeWienerSources)
    {
		WienerTerm.writeOpt() = IOobject::AUTO_WRITE;
    }

    volScalarField & fSF = sff[nsf];
    volScalarField   varWiener = fSF;
    volScalarField   gamma = turbulence->muEff()/rho;
    
	volVectorField gradTerm
	(      
		IOobject
		(
			"gradTerm_f_"+Foam::name(nsf),
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("zero",  f.dimensions()/dimLength, vector (0, 0, 0) ),
		zeroGradientFvPatchScalarField::typeName
	);		
	gradTerm = fvc::grad(varWiener);
	gradTerm.correctBoundaryConditions();
		
	volScalarField  varTerm  = gradTerm & dW_sf; 
	//varTerm.correctBoundaryConditions();

     forAll(varTerm, icell)
     {
		varTerm[icell]  = gradTerm[icell] & dW_sf.value();

		dimensionedScalar igamma("igamma",dimless,2*gamma[icell]);
		WienerTerm[icell] = sqrt(igamma).value()*varTerm[icell];
     }

 	 forAll(varTerm.boundaryField(), patchi)
	 {
		forAll(varTerm.boundaryField()[patchi],facei)
		{
			vector vecTerm = gradTerm.boundaryField()[patchi][facei];
			scalar sclTerm = gamma.boundaryField()[patchi][facei];

			varTerm.boundaryField()[patchi][facei] = vecTerm & dW_sf.value();

			dimensionedScalar igamma("igamma",dimless,2*sclTerm);
			WienerTerm.boundaryField()[patchi][facei] = sqrt(igamma).value()*varTerm.boundaryField()[patchi][facei];
		}
		
	 }

//make sure the Wiener Term for fixedValue boundaries must be zero
		const fvPatchList& patches = mesh.boundary();
		const volScalarField::GeometricBoundaryField& fSFbf = fSF.boundaryField();
		wordList fSFBoundaryTypes(fSF.boundaryField().size());
		fSFBoundaryTypes = fSFbf.types();
		forAll(patches, patchi)
		{

			if(fSFBoundaryTypes[patchi] == "fixedValue")
			{
				const fvPatch& currPatch = patches[patchi];			
//				Info << currPatch.name() <<" " << fSFBoundaryTypes[patchi] << endl;

				forAll(currPatch,facei)
				{
					WienerTerm.boundaryField()[patchi][facei] = 0;
					label faceCelli = currPatch.faceCells()[facei];
					WienerTerm[faceCelli] = 0;
				}
			}
		}
//
	forAll(fSF,icell)
	{
	    if(setReactingBox)
		{
			WienerTerm[icell] *= inReactingBox[icell];
		}	

		scalar fSFicell = fSF[icell];
		fSF[icell] += WienerTerm[icell];

		if(!balanceWienerTerm)
		{
			fSF[icell] = min(fSF[icell],maxMixFrc.value());
			fSF[icell] = max(fSF[icell],0.0);
			WienerTerm[icell] = fSF[icell] - fSFicell;
		}	

		WienerTermMF_t[icell] += WienerTerm[icell];
			
	}

	forAll(fSF.boundaryField(), patchi)
	{
		forAll(fSF.boundaryField()[patchi],facei)
		{

			if(setReactingBox)
			{
				WienerTerm.boundaryField()[patchi][facei] *= inReactingBox.boundaryField()[patchi][facei];
			}	
			scalar fSFiface = fSF.boundaryField()[patchi][facei];
			
			fSF.boundaryField()[patchi][facei] += WienerTerm.boundaryField()[patchi][facei];
			if(!balanceWienerTerm)
			{
				fSF.boundaryField()[patchi][facei] = min(fSF.boundaryField()[patchi][facei],maxMixFrc.value());
				fSF.boundaryField()[patchi][facei] = max(fSF.boundaryField()[patchi][facei],0.0);
				WienerTerm.boundaryField()[patchi][facei] = fSF.boundaryField()[patchi][facei]  - fSFiface;
			}
			WienerTermMF_t.boundaryField()[patchi][facei] += WienerTerm.boundaryField()[patchi][facei];
		}
	}

	if(runTime.write() && writeWienerSources)
	{
		WienerTerm.write();
		gradTerm.write();
	}

}

//

if(balanceWienerTerm)
{

	for(label nsf = 0; nsf < nStochasticFields; nsf++)
	{
		volScalarField & PVSF = sfPV[nsf];
		volScalarField & fSF  = sff[nsf];
		

		forAll(WienerTermPV_t,icell)
		{
			fSF[icell]  -= WienerTermMF_t[icell]/nStochasticFields;			
			PVSF[icell] -= WienerTermPV_t[icell]/nStochasticFields;	
			
			fSF[icell] = min(fSF[icell],maxMixFrc.value());
			fSF[icell] = max(fSF[icell],0.0);

			PVSF[icell] = min(PVSF[icell],maxProgVar.value());
			PVSF[icell] = max(PVSF[icell],0.0);				
		}

		forAll(WienerTermPV_t.boundaryField(), patchi)
		{
			forAll(WienerTermPV_t.boundaryField()[patchi],facei)
			{

				fSF.boundaryField()[patchi][facei]  -= WienerTermMF_t.boundaryField()[patchi][facei]/nStochasticFields;
				PVSF.boundaryField()[patchi][facei] -= WienerTermPV_t.boundaryField()[patchi][facei]/nStochasticFields;

				fSF.boundaryField()[patchi][facei] = min(fSF.boundaryField()[patchi][facei],maxMixFrc.value());
				fSF.boundaryField()[patchi][facei] = max(fSF.boundaryField()[patchi][facei],0.0);

				PVSF.boundaryField()[patchi][facei] = min(PVSF.boundaryField()[patchi][facei],maxProgVar.value());
				PVSF.boundaryField()[patchi][facei] = max(PVSF.boundaryField()[patchi][facei],0.0);

			}
		}

	}
}

Info <<"For " << nStochasticFields << " stochastic fields: Wiener process done " << endl;
Info << " WienerTermPV_t max/min: "<<  max(WienerTermPV_t).value() <<" " << min(WienerTermPV_t).value() << endl;
Info << " WienerTermMF_t max/min: "<<  max(WienerTermMF_t).value() <<" " << min(WienerTermMF_t).value() << endl;
//Info << " WienerTermH_t max/min: "<<  max(WienerTermH_t).value() <<" " << min(WienerTermH_t).value() << endl;//deleted by jianhong

// step 0: create thermophysic of the current stochastic field
IOdictionary stochasticFieldsProperties
(
	IOobject
    (
		"stochasticFieldsProperties",
		runTime.constant(),
		runTime,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

fileName XYpath = "XY"+args.globalCaseName();
mkDir(XYpath);

Info<< nl << "reading stochastic fields properties" << endl;
label nStochasticFields = 1;
stochasticFieldsProperties.lookup("nStochasticFields") >> nStochasticFields;
Switch writeWienerSources(stochasticFieldsProperties.lookup("writeWienerSources"));
Switch readStochasticFields(stochasticFieldsProperties.lookup("readStochasticFields"));
Switch writeStochasticFields(stochasticFieldsProperties.lookup("writeStochasticFields"));
Switch transportPVf(stochasticFieldsProperties.lookup("transportPVf"));
Switch nonReacting(stochasticFieldsProperties.lookup("nonReacting"));
Switch removeStochasticNoise(stochasticFieldsProperties.lookup("removeStochasticNoise"));
Switch balanceWienerTerm(stochasticFieldsProperties.lookup("balanceWienerTerm"));

Switch writeSpeciesStochasticFields(stochasticFieldsProperties.lookup("writeSpeciesStochasticFields"));

Switch setReactingBox(stochasticFieldsProperties.lookup("setReactingBox"));


dimensionedVector reactingBoxCenter = 0*mesh.C()[0];
dimensionedScalar reactingBoxL("reactingBoxL",dimLength, 0.0);
dimensionedScalar reactingBoxW("reactingBoxW",dimLength, 0.0);
dimensionedScalar reactingBoxH("reactingBoxH",dimLength, 0.0);

volScalarField inReactingBox
(
	IOobject
	(
		"inReactingBox",
		 mesh.time().timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::NO_WRITE
	),
    mesh,
   dimensionedScalar("zero", dimless, 1.0)
);

if(setReactingBox)
{
stochasticFieldsProperties.lookup("reactingBoxCenter")   >> reactingBoxCenter;
stochasticFieldsProperties.lookup("reactingBoxL")        >> reactingBoxL;
stochasticFieldsProperties.lookup("reactingBoxW")        >> reactingBoxW;
stochasticFieldsProperties.lookup("reactingBoxH")        >> reactingBoxH;
inReactingBox.writeOpt() = IOobject::AUTO_WRITE;
}


dimensionedScalar C_wFactor
(
stochasticFieldsProperties.lookup("C_wFactor")
);


dimensionedScalar Cphi
(
stochasticFieldsProperties.lookup("Cphi")
);


dimensionedScalar maxMixFrc
(
stochasticFieldsProperties.lookup("maxMixFrc")
);


dimensionedScalar maxProgVar
(
stochasticFieldsProperties.lookup("maxProgVar")
);

const  label nSFields = nStochasticFields;
Info<< "*** creating " << nSFields << " f and PV stochastic fields " << endl;

Random rndGen(0);

List < vector > eta_sf(nSFields);
forAll(eta_sf,nsf) eta_sf[nsf] = vector (0,0,0);

vector sum_eta (0,0,0);
List < vector > sum_eta_sf(nSFields);
forAll(sum_eta_sf,nsf) sum_eta_sf[nsf] = vector (0,0,0);

scalar tot_eta = 0;
scalar mean_eta =0;

PtrList< volScalarField > sff(nSFields);
PtrList< volScalarField > sfPV(nSFields);
PtrList< volScalarField > sfPVs(nSFields);

PtrList < PtrList<volScalarField>  > Ysf(nSFields);
PtrList< volScalarField > Tsf(nSFields);

OFstream WienerWalk_out(XYpath+"/WienerWalkAllFields_"+runTime.timeName());

List<OFstream*> WienerWalkEachField_Out(nSFields);

forAll(WienerWalkEachField_Out, nsf)
{
    WienerWalkEachField_Out[nsf] = new OFstream(XYpath+"/WienerWalk_"+Foam::name(nsf)+"_"+runTime.timeName());
}


WienerWalk_out << "# time " <<tab<< " x-direction "<<tab<< "y-direction"<<tab<<"z-direction"<< tab<<" mean " <<endl; 



if(nSFields < 2)
{
readStochasticFields  = false;
writeStochasticFields = false;
writeWienerSources    = false;
}


    volScalarField WienerTermPV_t
    (
      IOobject
      (
       "WienerTermPV_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
          mesh,
          dimensionedScalar("zero", PV.dimensions(), 0.0),
          zeroGradientFvPatchScalarField::typeName
    );

    volScalarField   WienerTermMF_t
    (
      IOobject
      (
       "WienerTermMF_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
          mesh,
          dimensionedScalar("zero", f.dimensions(), 0.0),
          zeroGradientFvPatchScalarField::typeName
    );

forAll(sff,nsf)
{
	IOobject header
	(
		f.name() +"_" + Foam::name(nsf)+"_SF",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);
        // Check if field exists and can be read
    if (header.headerOk())
    {
		sff.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					f.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
             )
        );
     }
     else
     {

        sff.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					f.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                f
            )
        );
     }
	 if(!writeStochasticFields)
	 {
		sff[nsf].writeOpt() = IOobject::NO_WRITE;
	 }
}

forAll(sfPV,nsf)
{
	IOobject header
	(
		PV.name() +"_" + Foam::name(nsf)+"_SF",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);
        // Check if field exists and can be read
    if (header.headerOk())
    {
		sfPV.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					PV.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
             )
        );
     }
     else
     {
        sfPV.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					PV.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                PV
            )
        );
     }
	 if(!writeStochasticFields)
	 {
		sfPV[nsf].writeOpt() = IOobject::NO_WRITE;
	 }
}

forAll(sfPVs,nsf)
{
	sfPVs.set
    (
		nsf,
        new volScalarField
        (
			IOobject
            (
				"PVs_" + Foam::name(nsf)+"_SF",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            meanPVs
        )
     );
	if(!writeStochasticFields)
	{
		sfPVs[nsf].writeOpt() = IOobject::NO_WRITE;
	}

}

//
const label  nSpecies = composition.Y().size();
forAll(Ysf,nsf)
{
    Ysf.set
    (
        nsf,
        new PtrList<volScalarField>(nSpecies)
    );


    forAll(Ysf[nsf], i)
    {
        IOobject header
        (
            Y[i].name() +"_" + Foam::name(nsf)+"_SF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // Check if field exists and can be read
        if (header.headerOk())
        {
            Ysf[nsf].set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        Y[i].name() +"_" + Foam::name(nsf)+"_SF",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {

				Ysf[nsf].set
				(
					i,
					new volScalarField
					(
						IOobject
						(
							Y[i].name() +"_" + Foam::name(nsf)+"_SF",
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
						),
						Y[i]
					)
				);

        }
		if(!writeSpeciesStochasticFields)
		{
			Ysf[nsf][i].writeOpt() = IOobject::NO_WRITE;
		}		

    }
}
 
forAll(Ysf,nsf)
{
	volScalarField Yt(0.0*Y[0]);
	label inertIndex = -1;
	
	forAll(Ysf[nsf],i)
	{
		if (Y[i].name() != inertSpecie )
		{
			Yt += Ysf[nsf][i]; 
		}
		else
		{
			inertIndex = i;
		}
		

	}

	Ysf[nsf][inertIndex] = scalar(1) - Yt;
	Ysf[nsf][inertIndex].max(0.0);
}

forAll(Tsf,nsf)
{
	IOobject header
	(
		T.name() +"_" + Foam::name(nsf)+"_SF",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);
        // Check if field exists and can be read
    if (header.headerOk())
    {
		Tsf.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					T.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
             )
        );
     }
     else
     {
		volScalarField Tii
        (
			IOobject
            (
				T.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        Tsf.set
        (
			nsf,
            new volScalarField
            (
				IOobject
                (
					T.name() +"_" + Foam::name(nsf)+"_SF",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Tii
            )
        );
     }
	 if(!writeStochasticFields)
	 {
		Tsf[nsf].writeOpt() = IOobject::NO_WRITE;
	 }
}



volScalarField Omega_t
(
	IOobject
	(
		"Omega_t",
		 mesh.time().timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::AUTO_WRITE
	),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)
);


const word turbulenceModelType
(
	IOdictionary
    (
		IOobject
        (
			"turbulenceProperties",
			runTime.constant(),
			runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ).lookup("simulationType")
);




//
/*	volScalarField curPV
	(
		IOobject
		(
			"curPV",
			 mesh.time().timeName(),
			 mesh,
			 IOobject::NO_READ,
			 IOobject::NO_WRITE
		),
        PV
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
        f
	);
*/
/*fields.add(curf);
fields.add(curPV);*/

forAll(sff,nsf)
{
fields.add(sff[nsf]);
fields.add(sfPV[nsf]);

forAll(Y,ns) fields.add(Ysf[nsf][ns]);
fields.add(Tsf[nsf]);
}

// RMS Fields only for output

volScalarField           rmsf
(
    IOobject
    (
        "rmsf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
      mesh,
      dimensionedScalar("rmsf", f.dimensions(), 0),
      zeroGradientFvPatchScalarField::typeName
);
if(nSFields < 2)
{
	rmsf.writeOpt() = IOobject::NO_WRITE;
}

volScalarField           rmsPV
(
    IOobject
    (
        "rmsPV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
      mesh,
      dimensionedScalar("rmsPV", PV.dimensions(), 0),
      zeroGradientFvPatchScalarField::typeName
);
     if(nSFields < 2)
     {
      rmsPV.writeOpt() = IOobject::NO_WRITE;
     }


volScalarField           rmsT
(
    IOobject
    (
        "rmsT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
      mesh,
      dimensionedScalar("rmsT", T.dimensions(), 0),
      zeroGradientFvPatchScalarField::typeName
);
     if(nSFields < 2)
     {
      rmsT.writeOpt() = IOobject::NO_WRITE;
     }

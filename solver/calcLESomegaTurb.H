volScalarField LESdelta
(
	IOobject
    (
		"LESdelta",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("LESdelta", dimLength, SMALL)
);



if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
{
	const volScalarField& delta =
    mesh.objectRegistry::lookupObject<volScalarField>("delta");

    LESdelta.internalField() = delta.internalField();
}

Info << " LESdelta max/min: "<<  max(LESdelta).value() <<" " << min(LESdelta).value() << endl;

const volScalarField muSgsCell = turbulence->DfEff();
forAll(Omega_t,icell)
{
	Omega_t[icell] = (lesC_e/lesC_k)*muSgsCell[icell]/(rho[icell]*LESdelta[icell]*LESdelta[icell]);
}

Info << " Omega_t max/min: "<<  max(Omega_t).value() <<" " << min(Omega_t).value() << endl;

Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::rhoCombustionModel> combustion
(
    combustionModels::rhoCombustionModel::New(mesh)
);
combustion->correct();

rhoReactionThermo& thermo = combustion->thermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

Info <<" Y.size() "<<Y.size()<< endl;
    

word inertSpecie(thermo.lookup("inertSpecie"));

PtrList<gasHsThermoPhysics> specieData(Y.size());
forAll(specieData, i)
{
	specieData.set
    (
		i,
        new gasHsThermoPhysics
        (
			dynamic_cast<const reactingMixture<gasHsThermoPhysics>&>
            (thermo).speciesData()[i]
        )
    );
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

// Mixture Fraction
volScalarField& f = thermo.f();

// Progress Variable
volScalarField& PV = thermo.PV();

// Enthalpy
//volScalarField& he = thermo.he();

volScalarField& T = thermo.T();

// Compressibility
const volScalarField& psi = thermo.psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField meanPVs
(
	IOobject
    (
     	"meanPVs",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    thermo.PVs()
);


volScalarField meanPV
(
	IOobject
    (
     	"meanPV",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    thermo.PV()
);


volScalarField meanf
(
	IOobject
    (
     	"meanf",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    thermo.f()
);


#include "compressibleCreatePhi.H"

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "readGravitationalAcceleration.H"
Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());

surfaceScalarField ghf("ghf", g & mesh.Cf());


//#include "readhRef.H"
//#include "gh.H"


volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

/*
surfaceScalarField heatFlux
(
       fvc::interpolate(turbulence->alphaEff()*thermo.Cp())*fvc::snGrad(thermo.T())
);

volScalarField wallHeatFlux
(
    IOobject
    (
        "wallHeatFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
);

#include "heatFlux.H"
*/

heatFlux = fvc::interpolate(turbulence->alphaEff()*thermo.Cp())*fvc::snGrad(thermo.T());

const surfaceScalarField::GeometricBoundaryField& patchHeatFlux = heatFlux.boundaryField();

forAll(wallHeatFlux.boundaryField(), patchi)
{
      wallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
}

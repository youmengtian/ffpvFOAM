Info<< nl << "compute RMS from stochastic fields" << endl;

volScalarField rms2PV = 0*PV*PV;
volScalarField rms2f = 0*f*f;
volScalarField rms2T = 0*T*T;

forAll(sfPV,nsf)
{
	rms2PV +=(sfPV[nsf]-PV)*(sfPV[nsf]-PV)/nStochasticFields;
	rms2f +=(sff[nsf]-f)*(sff[nsf]-f)/nStochasticFields;
	rms2T +=(Tsf[nsf]-T)*(Tsf[nsf]-T)/nStochasticFields;
}
rmsPV =sqrt(rms2PV);
rmsf =sqrt(rms2f);
rmsT =sqrt(rms2T);


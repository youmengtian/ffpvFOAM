bool isSparking = runTime.value() > sparkBeg.value() && runTime.value()<= sparkBeg.value() + sparkDuration.value();

label nSparkingCell = 0;
scalarField sparkingList(rho.size(),-1);

if(isSparking)
{
	volVectorField cellCenter = mesh.C();

	if(spherSpark)
	{

		forAll(PV,icell)
		{
			if(mag(cellCenter[icell] - sparkCenter.value() ) <= sparkRadi.value() &&
				f[icell] > minSparkingF.value() && f[icell] < maxSparkingF.value() &&
				T[icell] < sparkingT.value() )
			{
				sparkingList[nSparkingCell] = icell;
				nSparkingCell ++;
			}
		}

	}
	if(boxSpark)
	{

		forAll(PV,icell)
		{
			bool cellInBox = false;
			if( mag(cellCenter[icell].x() - sparkCenter.value().x() ) <= 0.5*sparkL.value()     &&
				mag(cellCenter[icell].y() - sparkCenter.value().y() ) <= 0.5*sparkW.value() &&
				mag(cellCenter[icell].z() - sparkCenter.value().z() ) <= 0.5*sparkH.value() && 
                                f[icell] > minSparkingF.value() && f[icell] < maxSparkingF.value()          &&
				T[icell] < sparkingT.value()	
			)
			{
				cellInBox = true;
			}
			if(cellInBox)
			{
				sparkingList[nSparkingCell] = icell;
				nSparkingCell ++;
			}
		}


	}

}

label totSparkingCell = nSparkingCell;
reduce(totSparkingCell, sumOp<label>());


sparkPVs *= 0;

for(label ispark = 0; ispark < nSparkingCell; ispark ++)
{
	label icell     = sparkingList[ispark];
	sparkPVs[icell] = sparkingPVs.value();
}
Info << "set sparking cells: " << " added PVs " << sparkingPVs.value() << " for detected igniting cells of: " << totSparkingCell<< endl;


/*
if(nStochasticFields > 1 )
{
forAll(sfPV,nsf)
{
forAll(sparkingList,spCell)
{
	label icell = sparkingList[spCell];
	sfPV[nsf][icell]   = sparkingPV.value();
}
Info << "set sparking cells: " << nsf << " number of igniting cells detected: " << nSparkingCell << endl;
}
}
else
{

forAll(sparkingList,spCell)
{
	label icell = sparkingList[spCell];
	PV[icell]   = sparkingPV.value();
}
Info << "set sparking cells: " << " PV " << "number of igniting cells detected: " << nSparkingCell << endl;
}*/

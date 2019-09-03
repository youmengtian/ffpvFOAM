label nInReactingBox = 0;

volVectorField cellCenter = mesh.C();

forAll(inReactingBox,icell)
{
		if( mag(cellCenter[icell].x() - reactingBoxCenter.value().x() ) <= 0.5*reactingBoxL.value() &&
			mag(cellCenter[icell].y() - reactingBoxCenter.value().y() ) <= 0.5*reactingBoxW.value() && 
			mag(cellCenter[icell].z() - reactingBoxCenter.value().z() ) <= 0.5*reactingBoxH.value() 	
		  )
			{
				inReactingBox[icell] = 1;
				nInReactingBox ++;
			}
			else
			{
				inReactingBox[icell] = 0;
			}
}

forAll(inReactingBox.boundaryField(), patchi)
{
	forAll(inReactingBox.boundaryField()[patchi],facei)
	{
		if( mag(cellCenter.boundaryField()[patchi][facei].x() - reactingBoxCenter.value().x() ) <= 0.5*reactingBoxL.value() &&
			mag(cellCenter.boundaryField()[patchi][facei].y() - reactingBoxCenter.value().y() ) <= 0.5*reactingBoxW.value() && 
			mag(cellCenter.boundaryField()[patchi][facei].z() - reactingBoxCenter.value().z() ) <= 0.5*reactingBoxH.value() 	
		  )
		  {
			inReactingBox.boundaryField()[patchi][facei] = 1;
			nInReactingBox ++;
		  }
		  else
		  {
			inReactingBox.boundaryField()[patchi][facei] = 0;		
		  }

	}
}

reduce(nInReactingBox, sumOp<label>());
Info << " reacting box detected in " <<  nInReactingBox <<" number of cells " << endl;
Info <<	" reactingBoxCenter "  << reactingBoxCenter.value() <<
"\t" <<" reactingBoxL "    << reactingBoxL.value() <<
"\t"<< " reactingBoxW "    << reactingBoxW.value() <<
"\t" <<" reactingBoxH "    << reactingBoxH.value() << endl;



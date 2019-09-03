{
	List<label> patches(3,0);
	scalar flux = 0;

	label patchI = mesh.boundaryMesh().findPatchID("inlet_o2");
	label patchII = mesh.boundaryMesh().findPatchID("inlet_ch4");
	label patchIII = mesh.boundaryMesh().findPatchID("outlet");

    flux = sum(phi.boundaryField()[patchI]);
        reduce(flux, sumOp<scalar>());
    Info<< "Mass flux at inlet_o2: " << flux << endl;

    flux = sum(phi.boundaryField()[patchII]);
        reduce(flux, sumOp<scalar>());
    Info<< "Mass flux at inlet_ch4: " << flux << endl;

    flux = sum(phi.boundaryField()[patchIII]);
        reduce(flux, sumOp<scalar>());
    Info<< "Mass flux at outlet: " << flux << endl;

}

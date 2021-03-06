label cur_eta = -1;
forAll(eta_sf,nsf) eta_sf[nsf] = vector (0,0,0);

for(label nsf =0; nsf <nStochasticFields; nsf +=2)
{
// step : Wiener process
  
if (Pstream::myProcNo() == 0)
{
    //scalar kisi_0 = rndGen.GaussNormal();
    cur_eta ++;
    tot_eta ++; 
    scalar kisi_0=rndGen.GaussNormal();
    if(kisi_0 < 0)
    {
	 kisi_0=-1;
	}
	else
	{
	 kisi_0=1; 
	}
	cur_eta ++;
    tot_eta ++;
	scalar kisi_1=rndGen.GaussNormal();
    if(kisi_1 < 0)
    {
	 kisi_1=-1;
	}
	else
	{
	 kisi_1=1; 
	}
	cur_eta ++;
    tot_eta ++; 
	scalar kisi_2=rndGen.GaussNormal();
    if(kisi_2 < 0)
    {
	 kisi_2=-1;
	}
	else
	{
	 kisi_2=1; 
	}
    
    eta_sf[nsf]   = vector (kisi_0,kisi_1,kisi_2);
    eta_sf[nsf+1] = -eta_sf[nsf];

    sum_eta             += eta_sf[nsf]+eta_sf[nsf+1];  
    sum_eta_sf[nsf]     +=eta_sf[nsf];
    sum_eta_sf[nsf+1]   +=eta_sf[nsf+1];
    mean_eta             = mag(sum_eta)/tot_eta;
    
   for(label iproc=1; iproc<Pstream::nProcs(); iproc++)
   {   
    OPstream sendStream(Pstream::blocking,iproc);
    sendStream   << eta_sf[nsf] << eta_sf[nsf+1];
   }
   WienerWalk_out << runTime.timeName() << tab << sum_eta.x() << tab << sum_eta.y() << tab << sum_eta.z() << tab << mean_eta << endl; 
   *WienerWalkEachField_Out[nsf]<< runTime.timeName() << tab << sum_eta_sf[nsf].x() << tab << sum_eta_sf[nsf].y() << tab << sum_eta_sf[nsf].z() << endl; 
   *WienerWalkEachField_Out[nsf+1]<< runTime.timeName() << tab << sum_eta_sf[nsf+1].x() << tab << sum_eta_sf[nsf+1].y() << tab << sum_eta_sf[nsf+1].z() << endl; 

}
else
{
       IPstream reciveStream(Pstream::blocking,0);
       reciveStream >> eta_sf[nsf] >> eta_sf[nsf+1];	
} 
 

}

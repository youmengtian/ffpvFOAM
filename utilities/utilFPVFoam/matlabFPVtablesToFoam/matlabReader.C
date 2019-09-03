/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
#include "matlabReader.H"
#include "IFstream.H"
#include "fileName.H"
#include "OFstream.H"
#include "IOList.H"
#include "Gamma.h"
#include "SortableList.H"
#include <sstream>
#include <cmath>

void Foam::matlabReader::betaPDFIntegration(const scalar& Zeta)
{
if (Zeta != 0)
{
   scalar varf, pdfAlpha, pdfBeta;
   List<scalar> PDF(f_param_.size(), 0.0);

   for (int i=0;i<f_param_.size();i++)
   {
      varf = sqr(Zeta) * (f_param_[i]*(1.0 - f_param_[i]));

      if (varf > 1e-7)
      {
          pdfAlpha = f_param_[i] * ((f_param_[i]*(1.0-f_param_[i]))/varf-1);
          pdfBeta = (1.0 - f_param_[i]) * ((f_param_[i]*(1.0-f_param_[i]))/varf-1);

          // Limit alpha and beta but keep their ratio
          if (pdfAlpha > 500)
          {
             pdfBeta = pdfBeta/pdfAlpha * 500;
             pdfAlpha = 500;
          }
          if (pdfBeta > 500)
          {
             pdfAlpha = pdfAlpha/pdfBeta * 500;
             pdfBeta = 500;
          }

          int    gridPoints = 250;
          List<long double> hf(gridPoints, 0.0);
          List<long double> helpZ(gridPoints, 0.0);
          List<long double> delta(gridPoints, 0.0);

          if ((pdfAlpha > 1) && (pdfBeta > 1))
          {
             // Allocation of Z for PDF integration
             scalar Zmax = 0;
             int    n1 = 0;
             int    n2 = 0;
             PDF.clear();
             PDF.resize(f_param_.size(), 0.0);
             
             for (int j=0;j<f_param_.size();j++)
             {
                PDF[j] = std::pow(f_param_[j],(pdfAlpha-1.0)) * std::pow((1.0 - f_param_[j]),(pdfBeta-1.0)) * Gamma(pdfAlpha + pdfBeta) / (min(Gamma(pdfAlpha),1e17)*min(Gamma(pdfBeta),1e17));
                if (PDF[j] > PDF[max(j-1, 0)]) Zmax = f_param_[j];
             }

             if(pdfAlpha/pdfBeta <= 0.5)
             {
                n1 = 0.2*gridPoints;
                n2 = 0.8*gridPoints+1;
             }
             else if (pdfAlpha/pdfBeta >= 2)
             {
                n1 = 0.8*gridPoints;
                n2 = 0.2*gridPoints+1;
             }
             else
             {
                n1 = 0.5*gridPoints;
                n2 = 0.5*gridPoints+1;
             }

             //  Allocate Z for 0 < Z < Zmax
             scalar ex1 = 0.9;
             delta[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             for (int j=0; j<n1; j++)
             {
                delta[j] = pow(ex1,j) * delta[0];
                hf[j+1] = hf[j] + delta[j];
             }
             for (int j=1; j<n1; j++)
             {
                hf[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 1.1;
             delta[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
               delta[j] = pow(ex2,j) * delta[0];
               helpZ[j+1] = helpZ[j] + delta[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hf[n1+j] = hf[n1-1]+helpZ[j];
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hf[j] /= hf[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hf.size(), 0.0);

             for (int j=0;j<hf.size();j++)
             {
                PDF[j] = (std::pow(hf[j],(pdfAlpha-1e0))) * std::pow((1e0 - hf[j]),(pdfBeta-1e0)) * Gamma(pdfAlpha + pdfBeta)/(min(Gamma(pdfAlpha),1e17)*min(Gamma(pdfBeta),1e17));
             }
          }

          else if ((pdfAlpha <= 1) && (pdfBeta > 1))
          {
             // PDF Singularity at Z = 0
             // Allocation of Z for PDF integration
             int    n1 = 0;
             int    n2 = 0;

             if (pdfAlpha/pdfBeta > 0.5)
             {
                scalar Zmax = 0.5;
                scalar ex1 = 1.1;
                n1 = 0.7*gridPoints;
                delta[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

                // Allocate Z for 0 < Z < Zmax
                for (int j=0; j<n1; j++)
                {
                   delta[j] = pow(ex1,j) * delta[0];
                   hf[j+1] = hf[j] + delta[j];
                }
                for (int j=1; j<n1; j++)
                {
                   hf[j] *= Zmax;
                }

                // Allocate Z for Zmax < Z < 1
                scalar ex2 = 1.1;
                n2 = 0.3*gridPoints+1;
                delta[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

                for (int j=0; j<n2-1; j++)
                {
                   delta[j] = pow(ex2,j) * delta[0];
                   helpZ[j+1] = helpZ[j] + delta[j];
                }
                for (int j=0;j<n2;j++)
                {
                   helpZ[j] *= (1.0 - Zmax);
                }
                for (int j=0;j<n2-1;j++)
                {
                   hf[n1+j] = hf[n1-1]+helpZ[j];
                }
             }

             else
	         {
                scalar ex2 = 1.05;
                delta[0] = (1.0 - ex2)/(1.0 - pow(ex2,(gridPoints-1)));
                for (int j=0; j<gridPoints-1; j++)
                {
                   delta[j] = pow(ex2,j) * delta[0];
                   hf[j+1] = hf[j] + delta[j];
                }
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hf[j] /= hf[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hf.size(), 0.0);
             for (int j=1;j<hf.size();j++)
             {
                PDF[j] = std::pow(hf[j],(pdfAlpha-1e0)) * std::pow((1e0 - hf[j]),(pdfBeta-1.0)) * min(Gamma(pdfAlpha + pdfBeta),1e17)/(min(Gamma(pdfAlpha),1e17)*min(Gamma(pdfBeta),1e17));
             }
             PDF[0] = 1.5 * PDF[1] / pdfAlpha;

          }

          else if ((pdfAlpha > 1) && (pdfBeta <= 1))
          {
          // PDF Singularity at Z = 1
          // Allocation of Z for PDF integration

          int    n1 = 0;
          int    n2 = 0;

          if (pdfAlpha/pdfBeta < 2)
          {
             scalar Zmax = 0.5;
             scalar ex1 = 1.1;
             n1 = 0.3*gridPoints;
             delta[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             // Allocate Z for 0 < Z < Zmax
             for (int j=0; j<n1; j++)
             {
                delta[j] = pow(ex1,j) * delta[0];
                hf[j+1] = hf[j] + delta[j];
             }
             for (int j=1; j<n1; j++)
             {
                hf[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 0.9;
             n2 = 0.7*gridPoints+1;
             delta[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
                delta[j] = pow(ex2,j) * delta[0];
                helpZ[j+1] = helpZ[j] + delta[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hf[n1+j] = hf[n1-1]+helpZ[j];
             }
          }

          else
          {
             scalar ex1 = 0.95;
             delta[0] = (1.0 - ex1)/(1.0 - pow(ex1,(gridPoints-1)));

             for (int j=0; j<gridPoints-1; j++)
             {
                 delta[j] = pow(ex1,j) * delta[0];
                 hf[j+1] = hf[j] + delta[j];
             }
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hf[j] /= hf[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hf.size(), 0.0);

          for (int j=0;j<hf.size()-1;j++)
          {
              PDF[j] = std::pow((hf[j]),(pdfAlpha- 1e0)) * std::pow((1e0 - hf[j]),(pdfBeta-1.0)) * min(Gamma(pdfAlpha + pdfBeta),1e17)/(min(Gamma(pdfAlpha),1e17)*min(Gamma(pdfBeta),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta;

       }

       else if ((pdfAlpha <= 1) && (pdfBeta <= 1))
       {
          // PDF Singularity at Z = 1 and Z = 0
          // Allocation of Z for PDF integration
          int    n1 = 0;
          int    n2 = 0;
          scalar Zmax = 0.5;
          scalar ex1 = 1.1;
          n1 = 0.5*gridPoints;
          delta[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

          // Allocate Z for 0 < Z < Zmax
          for (int j=0; j<n1; j++)
          {
              delta[j] = pow(ex1,j) * delta[0];
              hf[j+1] = hf[j] + delta[j];
          }
          for (int j=1; j<n1; j++)
          {
             hf[j] *= Zmax;
          }

          // Allocate Z for Zmax < Z < 1
          scalar ex2 = 0.9;
          n2 = 0.5*gridPoints+1;
          delta[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

          for (int j=0; j<n2-1; j++)
          {
              delta[j] = pow(ex2,j) * delta[0];
              helpZ[j+1] = helpZ[j] + delta[j];
          }
          for (int j=0;j<n2;j++)
          {
              helpZ[j] *= (1.0 - Zmax);
          }
          for (int j=0;j<n2-1;j++)
          {
             hf[n1+j] = hf[n1-1]+helpZ[j];
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hf[j] /= hf[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hf.size(), 0.0);

          for (int j=1;j<hf.size()-1;j++)
          {
             PDF[j] = std::pow(hf[j],(pdfAlpha- 1e0)) * std::pow((1e0 - hf[j]),(pdfBeta-1.0)) * min(Gamma(pdfAlpha + pdfBeta),1e17)/(min(Gamma(pdfAlpha),1e17)*min(Gamma(pdfBeta),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta;
          PDF[0] = 1.5 * PDF[1] / pdfAlpha;
       }

       // Calculate the area of the PDF for scaling
       scalar intPDF = 0;
       for (int j=1;j<hf.size();j++)
       {
    	  intPDF += (hf[j-1] - hf[j]) * (PDF[j-1] + PDF[j])/2; //area calculation for trapezoid, denoted by jianhong
       }

       // Interpolate singleData entries to the new mixture fraction space
       List<scalar> hY(gridPoints, 0.0);
       scalar intY = 0;

       for (int j=0;j<integratedData_.size();j++)
       {
          hY = 0.0;
          hY[0] = unintegratedData_[j][0];
          hY[hY.size()-1] = unintegratedData_[j][unintegratedData_[j].size()-1];
          intY = 0;

          for (int k=1;k<hf.size()-1;k++)
          {
             int ubZ = 0;
             for (int l=0;l<f_param_.size();l++)
             {
                ubZ = l;
                if (hf[k] < f_param_[l])
                break;
             }
             int lbZ = ubZ -1;

             // Interpolation to hf space
             hY[k] = (unintegratedData_[j][ubZ] - unintegratedData_[j][lbZ])/max(f_param_[ubZ] - f_param_[lbZ], SMALL) * (hf[k] - f_param_[lbZ]) + unintegratedData_[j][lbZ];
             // PDF Integration using the trapezoidal rule
             intY += (hf[k-1] - hf[k]) * (hY[k-1]*PDF[k-1] + hY[k]*PDF[k])/(2.0 * intPDF);
          }

          // Special treatment for the boundaries
          intY += (hf[hf.size()-2] - hf[hf.size()-1]) * (hY[hf.size()-2]*PDF[hf.size()-2] + hY[hf.size()-1]*PDF[hf.size()-1])/(2.0 * intPDF);
          if (i != 0 && i != f_param_.size()-1)
             integratedData_[j][i] = intY;
       }
     }
   }
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matlabReader::matlabReader(const IOdictionary& matlabDict, rhoReactionThermo& thermo, basicMultiComponentMixture& composition, const fvMesh& mesh)
:  composition(composition),
   thermo(thermo),
   tables_
   (
       IOdictionary
       (
          IOobject
          (
    	     "tables",
             mesh.time().constant(),
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE,
             false
          )
       )
   ),
   tableNames_(thermo.composition().species()),
   Zeta_param_(matlabDict.lookup("Zeta_param")),
   PV_param_(matlabDict.lookup("PV_param")),
   f_param_(matlabDict.lookup("f_param"))
{
   tableNames_.append("T");
   tableNames_.append("PVs");
   tableNames_.append("mu");
   tableNames_.append("alpha");
   tableNames_.append("psi");
   tableNames_.append("HeatRelease");
   tableNames_.append("rho"); //added by jianhong, to print the rho_table"

   readData_.resize(tableNames_.size());
   sampledData_.resize(tableNames_.size());
   integratedData_.resize(tableNames_.size());
   unintegratedData_.resize(tableNames_.size());

   // Initialize readData_
   for (int i=0; i<readData_.size(); i++)
       {
   	   readData_[i].resize(PV_param_.size());
       for (int k=0; k<readData_[i].size();k++)
           {
           		readData_[i][k].resize(f_param_.size());
           }

      }

   // Initialize sampledData_
   for (int i=0; i<sampledData_.size(); i++)
       {
            sampledData_[i].resize(PV_param_.size());
            for (int k=0; k<sampledData_[i].size();k++)
               {
                   sampledData_[i][k].resize(Zeta_param_.size());
                   for (int l=0; l<sampledData_[i][k].size();l++)
                   {
                 	  sampledData_[i][k][l].resize(f_param_.size());
                   }
               }
      }

   for (int i=0; i<tableNames_.size(); i++)
   {
	   // Correct size of integratedData_ and unintegratedData_
	   integratedData_[i].resize(f_param_.size());
	   unintegratedData_[i].resize(f_param_.size());

	   // Read all tables
	   Info << "Reading " << tableNames_[i] << "-table" << endl;
	   List<List<scalar> > tmpRead(tables_.lookup(tableNames_[i]));
	   readData_[i] = tmpRead;
   }

   Info << nl << "Start PDF integration:" << nl << endl;

   // Integration


       // for all PV
       for (int numPV=0; numPV<PV_param_.size();numPV++)
       {
    	   Info << "Progress: " << 100*numPV/(PV_param_.size()-1) << " %" <<  endl;
          // for all species
    	  for (int numSpec=0; numSpec<tableNames_.size();numSpec++)
    	  {
    		  unintegratedData_[numSpec] = readData_[numSpec][numPV];
    	  }
    	  // for all Zeta
    	  for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
    	  {
    		  // Write Integrated Data in sampledData
    		  integratedData_ = unintegratedData_;
    		  betaPDFIntegration(Zeta_param_[numZeta]);

    		  for (int k=0; k<integratedData_.size();k++)
    		  {
    			  sampledData_[k][numPV][numZeta] = integratedData_[k];
    		  }
    	  }
      }

   Info << nl << "PDF integration done!" << endl;

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::matlabReader::~matlabReader()
{}

// * * * * * * * * * * * * * * * * MemberFunctions  * * * * * * * * * * * * * * * //

Foam::hashedWordList Foam::matlabReader::getNames()
{
	return tableNames_;
}

void	Foam::matlabReader::write(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_table";
	List<List<List<scalar> > >lists=sampledData_[i];
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}

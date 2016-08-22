#ifndef _eocascadeInit2_h
#define _eocascadeInit2_h
#include <eoInit.h>
#include <vector>
#include "Systeme.h"

template <class GenotypeT>
class eocascadeInit2: public eoInit<GenotypeT> {
public:
  eocascadeInit2(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
  {
	  V=_V;
	int i,j,h;
	systeme=_systeme;
	nbReservoirs=systeme->getNbReservoirs();
	nbHeures= _nbHeures;
	for(i=0;i<nbReservoirs;i++)
	  {
		Reservoir* R=systeme->getReservoir(i);
		qTot.push_back(0);
		for(h=0;h<nbHeures;h++)
		{
			qTot[i]=qTot[i]+R->getApport(h);
		}
	        for(j=0;j<R->getNbParents();j++)
		{
		
			 int parent=R->getParents()[j];
			 qTot[i]=qTot[i]+qTot[parent];
		}
		
	  }
	
  }


  void operator()(GenotypeT & _genotype)
  {

	int i,j,h,t;
	vector<double> quantite;
	_genotype.adEtat(qTot);
	//Autres états :
	eoUniformGenerator<double> cb1(0,40);
	double b1=cb1();
	eoUniformGenerator<double> cb2(70,95);
	double b2=cb2();
	eoUniformGenerator<double> cb3(b1,b2);
	double b3=cb3();
	for(h=nbHeures-2;h>-1;h--)
	{
		  quantite.clear();
		  for(i=0;i<nbReservoirs;i++)
		  {
				Reservoir* R=systeme->getReservoir(i);
				quantite.push_back(0);
				double Vi=V[h][i];
			        int nbP=R->getNbParents();
				for(j=0;j<nbP;j++)
				{
					  
					int parent=R->getParents()[j];
					Vi=Vi+quantite[parent]*3600;
				
				}
				  //calcul de qmin:
				
				double qminC=R->getQmin();
				if(qminC<0) qminC=0;
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				double qmin=(Vi-Vmax)/3600;
				double qmax=_genotype.getQuantite(nbHeures-h-2,i)-qminC;
				if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
				double qmaxD=R->getQmax();
				if(R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					int Int=systeme->getTurbine(turbine)->getNbInt()-1;
					qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
				}
				if(qmin<_genotype.getQuantite(nbHeures-h-2,i)-qmaxD&&-qmaxD+_genotype.getQuantite(nbHeures-h-2,i)<=qmax&&qmaxD>0) qmin=-qmaxD+_genotype.getQuantite(nbHeures-h-2,i);
				if(qmin<qminC*(1+h)&&qminC*(1+h)<=qmax)qmin=qminC*(1+h);;
				if(nbP==0 &&qminC>0)
				{
					for(t=0;t<h;t++)
					{
						if(V[t][i]/3600+(h-t)*qminC-Vmax/3600>qmin && V[t][i]/3600+(h-t)*qminC-Vmax/3600<qmax)qmin=V[t][i]/3600+(h-t)*qminC-Vmax/3600;
					}
				}
				if(qmax>qTot[i]) qmax=qTot[i];
				if(qmax<qmin)qmax=qmin;
				//choix aléatoire
				eoUniformGenerator<double> random(qmin,qmax);
			        eoUniformGenerator<double> rim(0,100);
				double p=rim();
				quantite[i]=random();
				if(p<b1) quantite[i]=qmax;
			        if(p>b2)quantite[i]=qmin;
				
			  }
			  _genotype.adEtat(quantite);
	  	}
		//Inversion:
		for(h=0;h<nbHeures/2;h++)
		{
			for(j=0;j<nbReservoirs;j++)
			{
				double temp=_genotype.getQuantite(h,j);
				_genotype.setQuantite(h,j,_genotype.getQuantite(nbHeures-1-h,j));
				_genotype.setQuantite(nbHeures-1-h,j,temp);
			}
		}
	
       		_genotype.invalidate();	   
  }

private:

vector< vector<double> >V;
Systeme* systeme;
int nbHeures;
int nbReservoirs;
 vector<double> qTot;

};

#endif

#ifndef _eocascadeInit5_h
#define _eocascadeInit5_h
#include <eoInit.h>
#include <vector>
#include "Systeme.h"

template <class GenotypeT>
class eocascadeInit5: public eoInit<GenotypeT> {
public:

  eocascadeInit5(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures,int _hFin,vector<double> _qTot)
{
	V=_V;
	int i,j,h;
	systeme=_systeme;
	nbReservoirs=systeme->getNbReservoirs();
	nbHeures= _nbHeures;
	hFin=_hFin;
	qTot=_qTot;
  }


  void operator()(GenotypeT & _genotype)
  {
	int i,j,h,t;
	_genotype.setNbReservoirs(nbReservoirs);
	
	eoUniformGenerator<double> cb1(0,40);
	double b1=cb1();
	eoUniformGenerator<double> cb2(70,95);
	double b2=cb2();
	eoUniformGenerator<double> cb3(b1,b2);
	double b3=cb3();
	for(h=hFin-1;h>-1;h--)
	{
		for(i=0;i<nbReservoirs;i++)
		{
			Reservoir* R=systeme->getReservoir(i);
			double Vi=V[h][i];
			int nbP=R->getNbParents();
			for(j=0;j<nbP;j++)
			{
					  
				int parent=R->getParents()[j];
				Vi=Vi+_genotype.getQuantite(h,parent)*3600;
			}
			double qminC=R->getQmin();
			if(qminC<0) qminC=0;
			double Vmin=R->getVmin(h);
			double Vmax=R->getVmax();
			double qmin=(Vi-Vmax)/3600;
			double qmax=_genotype.getQuantite(h+1,i)-qminC;
			if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
			//calcul de qmax
			double qmaxD=R->getQmax();
			//cas avec turbine
			if(R->getNbTurbines()>0)
			{
				int turbine=R->getTurbine(0);
				int Int=systeme->getTurbine(turbine)->getNbInt()-1;
				qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
			}
			if(qmin<_genotype.getQuantite(h+1,i)-qmaxD&&-qmaxD+_genotype.getQuantite(h+1,i)<=qmax&&qmaxD>0) qmin=-qmaxD+_genotype.getQuantite(h+1,i);
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
			eoUniformGenerator<double> random(qmin,qmax);
			eoUniformGenerator<double> rim(0,100);
			double p=rim();
			double qte=random();
			if(p<b1) qte=qmax;
			if(p>b2)qte=qmin;
			_genotype.setQuantite(h,i,qte);
		}
			
	  }
	
	//tests croissance:
	for(h=1;h<hFin;h++)
	{
		for(j=0;j<nbReservoirs;j++)
		{
			if(_genotype.getQuantite(h-1,j)>_genotype.getQuantite(h,j)) cout<<"erreur ini5 à l'h "<<h<<" "<<j<<endl;
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
int hFin;

};

#endif

#ifndef _eocascadeInit4_h
#define _eocascadeInit4_h
#include <eoInit.h>
#include <vector>
#include "Systeme.h"

template <class GenotypeT>
class eocascadeInit4: public eoInit<GenotypeT> {
public:
  eocascadeInit4(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures,int _hDeb,vector<double> _qDeb)
  {

	V=_V;
	int i,j,h;
	systeme=_systeme;
	nbReservoirs=systeme->getNbReservoirs();
	hDeb=_hDeb;
	qDeb=_qDeb;
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
	_genotype.setNbReservoirs(nbReservoirs);
	double qte;
	//premier état:
	for(i=0;i<nbReservoirs;i++)
	{
		_genotype.setQuantite(hDeb,i,qDeb[i]);
	}
	//Autres
	eoUniformGenerator<double> cb1(0,50);
	double b1=cb1();
	eoUniformGenerator<double> cb2(70,95);
	double b2=cb2();
	eoUniformGenerator<double> cb3(b1,b2);
	double b3=cb3();
	for(h=hDeb+1;h<nbHeures-1;h++)
	{
		for(i=0;i<nbReservoirs;i++)
		{
			Reservoir* R=systeme->getReservoir(i);
			double Vinit=V[h-1][i]-_genotype.getQuantite(h-1,i)*3600;
			double Vi=V[h][i];
			int nbP=R->getNbParents();
			for(j=0;j<nbP;j++)
			{
				int parent=R->getParents()[j];
				Vi=Vi+_genotype.getQuantite(h,parent)*3600;
				Vinit=Vinit+_genotype.getQuantite(h-1,parent)*3600;
			}
				  //calcul de qmin:
				
			double qminC=R->getQmin();
			if(qminC<0) qminC=0;
			double qmin=_genotype.getQuantite(h-1,i)+qminC;
			double Vmin=R->getVmin(h);
			double Vmax=R->getVmax();
			if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
			double qmax=(Vi-Vmin)/3600;
			if(i==1)
			{
				double Vh=Vi-(h+1)*qminC*3600;
				for(j=h+1;j<nbHeures;j++)
				{
					Vh=Vh+R->getApport(j)*3600-3600*qminC;
					double Vminh=R->getVmin(j);
					if((Vh-Vminh)/3600<qmax)qmax=(Vh-Vminh)/3600;
				}

			}
			//calcul de qmax
			double qmaxD=R->getQmax();
			//cas avec turbine
			if(R->getNbTurbines()>0)
			{
				int turbine=R->getTurbine(0);
				int Int=systeme->getTurbine(turbine)->getIntervalle(Vinit);
				qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
			}
			if(qmax>qmaxD+_genotype.getQuantite(h-1,i)&&qmaxD+_genotype.getQuantite(h-1,i)>qmin&&qmaxD>0) qmax=qmaxD+_genotype.getQuantite(h-1,i);
			if(qmax>qTot[i]-qminC*(nbHeures-1-h)&&qTot[i]-qminC*(nbHeures-1-h)>=qmin)qmax=qTot[i]-qminC*(nbHeures-1-h);
			if(qminC>0&&nbP==0)
			{
				for(t=h+1;t<nbHeures-1;t++)
				{
					double	Vmint=R->getVmin(t)/3600;
					if(V[t][i]/3600-(t-h)*qminC-Vmint<qmax &&V[t][i]/3600-(t-h)*qminC-Vmint>qmin)qmax=V[t][i]/3600-(t-h)*qminC-Vmint;
				}
			}
			if(qmax>qTot[i]) qmax=qTot[i];
			if(qmax<qmin)
			{
				//Correction:
				double Vh=Vinit;
				if(V[h-1][i]-(qmax-qminC)*3600<Vmax &&_genotype.getQuantite(h-1,i)-qmax+qminC<_genotype.getQuantite(h-1,i)-_genotype.getQuantite(h-2,i))
				{
					_genotype.setQuantite(h-1,i,qmax-qminC);
				}
					//else cout<<_genotype.getQuantite(h-1,i)-qmax+qminC<<" "<<_genotype.getQuantite(h-1,i)-_genotype.getQuantite(h-2,i)<<endl;
			}
			if(qmax<qmin)qmin=qmax;
			//choix aléatoire
			eoUniformGenerator<double> random(qmin,qmax);
			eoUniformGenerator<double> rim(0,100);
			double p=rim();
			qte=random();
			if(p<b1) qte=qmin;
			if(p>b2)qte=qmax;
	                _genotype.setQuantite(h,i,qte);
		}
			
	}
	//DernierEtat:
	_genotype.adEtat(qTot);

	for(i=0;i<nbReservoirs;i++)
	{
		for(h=hDeb+1;h<nbHeures;h++)
		{
			if(_genotype.getQuantite(h,i)<_genotype.getQuantite(h-1,i))cout<<"pb grave init4 "<<h<<" "<<i<<endl;
		}
	}
	
	_genotype.invalidate();	
  }

private:
// START Private data of an eocascadeInit object
vector< vector<double> >V;
Systeme* systeme;
int nbHeures;
int hDeb;
vector<double> qDeb;
int nbReservoirs;
 vector<double> qTot;
// END   Private data of an eocascadeInit object
};

#endif

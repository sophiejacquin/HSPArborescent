#ifndef _eocascadeInit_h
#define _eocascadeInit_h
#include <eoInit.h>
#include <vector>
#include "Systeme.h"

template <class GenotypeT>
class eocascadeInit: public eoInit<GenotypeT> {
public:
	eocascadeInit(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
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
		_genotype.setNbReservoirs(nbReservoirs);
	  	vector<double> quantite;
		//premier état:
	  	for(i=0;i<nbReservoirs;i++)
	  	{
		 	Reservoir* R=systeme->getReservoir(i);
		 	quantite.push_back(0);
		  	//calcul de l'interval dans lequel se trouve Vinit:
		  	double Vinit=R->getVinit();
			
		  	double Vi=V[0][i];
			int nbP=R->getNbParents();
			for(j=0;j<nbP;j++)
		  	{
				int parent=R->getParents()[j];
				Vi=Vi+quantite[parent]*3600;
		  	}
		  	//calcul de qmin:
		  	double qmin=0;
		  	double qminC=R->getQmin();
			
		  	if(qminC<0) qminC=0;
		  	if(qmin<qminC)qmin=qminC;
		  	double Vmax=R->getVmax();
		  	if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
		  	//calcul de qmax
		  	double Vmin=R->getVmin(0);
		  	double qmax=(Vi-Vmin)/3600;
		  	double qmaxD=R->getQmax();
		  
		  	//cas avec turbine
			if(R->getNbTurbines()>0)
			{
				int turbine=R->getTurbine(0);
				int Int=systeme->getTurbine(turbine)->getIntervalle(Vinit);
				
				if(qmaxD>0)  qmaxD=qmaxD+systeme->getTurbine(turbine)->getQMax(Int);
				else qmaxD=systeme->getTurbine(turbine)->getQMax(Int);
				
			}
		 	if(qmax>qmaxD&&qmaxD>qmin) qmax=qmaxD;
			if(qmax>qTot[i])qmax=qTot[i];
			if(qmax>qTot[i]-qminC*(nbHeures-1))qmax=qTot[i]-qminC*(nbHeures-1);
			if(qminC>0&&nbP==0)
			{
				for(t=1;t<nbHeures-1;t++)
				{
					double Vmint=R->getVmin(t)/3600;
					if(V[t][i]/3600-t*qminC-Vmint<qmax)
						qmax=V[t][i]/3600-t*qminC-Vmint;
				}
			}
			eoUniformGenerator<double> random(qmin,qmax);
			quantite[i]=random();
	 	};
		_genotype.adEtat(quantite);
	 	//Autres états :
		eoUniformGenerator<double> cb1(0,45);
		double b1=60;
		eoUniformGenerator<double> cb2(70,95);
		double b2=95;
		eoUniformGenerator<double> cb3(b1,b2);
		double b3=cb3();
		
		for(h=1;h<nbHeures-1;h++)
		{
			
			quantite.clear();
			for(i=0;i<nbReservoirs;i++)
			{
				
				Reservoir* R=systeme->getReservoir(i);
				quantite.push_back(0);
				//calcul de l'interval dans lequel se trouve Vinit:
				eoUniformGenerator<int> tMarche(0,1); 
				int marche=tMarche();
				double Vinit=V[h-1][i]-_genotype.getQuantite(h-1,i)*3600;
				double Vi=V[h][i];
			        int nbP=R->getNbParents();
				for(j=0;j<nbP;j++)
				{
					int parent=R->getParents()[j];
					
					Vi=Vi+quantite[parent]*3600;
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
				
				if(nbP==0 &&qminC>0)
				{
					double Vh=Vi-(h+1)*qminC*3600;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600-3600*qminC;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qmax)qmax=(Vh-Vminh)/3600;
					}

				}
				
				double qmaxD=R->getQmax();
				if(R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					int Int=systeme->getTurbine(turbine)->getIntervalle(Vinit);
					qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
				}
				
				if(qmax>qmaxD+_genotype.getQuantite(h-1,i)&&qmaxD+_genotype.getQuantite(h-1,i)>=qmin&&qmaxD>0) qmax=qmaxD+_genotype.getQuantite(h-1,i);
				
				if(qmax>qTot[i]-qminC*(nbHeures-1-h)&&qTot[i]-qminC*(nbHeures-1-h)>=qmin)qmax=qTot[i]-qminC*(nbHeures-1-h);
				/*TODO*/
				else{
					
					if(qTot[i]-qminC*(nbHeures-1-h)<qmin &&h<8000){
						qmax=qmin;
					}
				}
				
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
					
					//cout<<"correction "<<i<<" h "<<h<<" qmax="<<qmax<<" qTot="<<qTot[i]<<" (Vi-Vmin)/3600="<<(Vi-Vmin)/3600<<" (Vi-Vmax)/3600="<<(Vi-Vmax)/3600<<" qmin="<<qmin<<endl;//" "<<qminC<<"   "<<V[h-1][i]-(qmax-qminC)*3600<<"   "<<Vmax<<"   "<<_genotype.getQuantite(h-2,i)<<"      "<<qmax-qminC<<endl;
				//Correction:
					double Vh=Vinit;
					//if(V[h-1][i]-(qmax-qminC)*3600<Vmax && qmax-qminC>=_genotype.getQuantite(h-2,i))
					int htemp=h;
					while(htemp>0 && _genotype.getQuantite(htemp-1,i)>qmax)
					{
						//cout<<"correction "<<_genotype.getQuantite(htemp-1,i)<<" "<<qmax<<" htemp-1 "<<htemp-1<<" i="<<i<<endl;
						_genotype.setQuantite(htemp-1,i,qmax);
						htemp--;
					}
					
				}
				//cout<<"init vmax borne "<<(Vi-Vmax)/3600<<endl;
				
				
				  if(qmax<qmin)qmin=qmax;
				  if(qmax<(Vi-Vmax)/3600|| qmin<(Vi-Vmax)/3600) cout<<"*****************PROBLEME****"<<endl;
				  //choix aléatoire
				  eoUniformGenerator<double> random(qmin,qmax);
			          eoUniformGenerator<double> rim(0,100);
				  double p=rim();
				  quantite[i]=random();
				  if(p<b1) quantite[i]=qmin;
			          if(p>b2)quantite[i]=qmax;
				  if((Vi-Vmax)/3600>quantite[i]) cout<<"mystere"<<endl;
			  }
			  _genotype.adEtat(quantite);
	  	}
	  	//DernierEtat:
		
	  	_genotype.adEtat(qTot);

		for(i=0;i<nbReservoirs;i++)
		{
			for(h=1;h<nbHeures;h++)
			{
				if(_genotype.getQuantite(h,i)<_genotype.getQuantite(h-1,i))cout<<"*************************pb grave "<<h<<" "<<i<<"   "<<_genotype.getQuantite(h,i)-_genotype.getQuantite(h-1,i)<<endl;
			}
		}
		_genotype.invalidate();
			
		_genotype.check_Vmax(V,systeme);  
		
	}

private:

	vector< vector<double> >V;
	Systeme* systeme;
	int nbHeures;
	int nbReservoirs;
	vector<double> qTot;

};

#endif

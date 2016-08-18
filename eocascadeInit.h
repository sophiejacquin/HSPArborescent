/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for EO objects initialization in EO
============================================
*/

#ifndef _eocascadeInit_h
#define _eocascadeInit_h

// include the base definition of eoInit
#include <eoInit.h>
#include <vector>
#include "Systeme.h"
/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen
 *
 * There is NO ASSUMPTION on the class GenoypeT.
 * In particular, it does not need to derive from EO (e.g. to initialize
 *    atoms of an eoVector you will need an eoInit<AtomType>)
 */
template <class GenotypeT>
class eocascadeInit: public eoInit<GenotypeT> {
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeInit(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
  //  eocascadeInit( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeInit object
	  V=_V;
	int i,j,h;
	//cout<<"ici"<<endl;
	  systeme=_systeme;
	nbReservoirs=systeme->getNbReservoirs();

	  nbHeures= _nbHeures;
	 for(i=0;i<nbReservoirs;i++)
	  {
		//cout<<i<<" "<<nbReservoirs<<endl;
		Reservoir* R=systeme->getReservoir(i);
		  qTot.push_back(0);
		  for(h=0;h<nbHeures;h++)
		  {
			//cout<<h<<endl;
			  qTot[i]=qTot[i]+R->getApport(h);
		  }
		 for(j=0;j<R->getNbParents();j++)
		  {
			//cout<<j<<endl;
			  int parent=R->getParents()[j];
			 qTot[i]=qTot[i]+qTot[parent];
		  }
		  //
	  }
	//cout<<"ici fin"<<endl;
    // END   Code of Ctor of an eocascadeInit object
  }


  /** initialize a genotype
   *
   * @param _genotype  generally a genotype that has been default-constructed
   *                   whatever it contains will be lost
   */
  void operator()(GenotypeT & _genotype)
  {
	//cout<<"debut ini"<<endl;
    // START Code of random initialization of an eocascade object
	  int i,j,h,t;

	  _genotype.setNbReservoirs(nbReservoirs);
	  vector<double> quantite;
	 
	  //Calcul qTot;manque les qte parentes !
	 
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
		 // if(i==1)cout<<" qminC1 "<<qminC<<endl;
		  if(qmin<qminC)qmin=qminC;
		  double Vmax=R->getVmax();
		  if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
		  //calcul de qmax
		  double Vmin=R->getVmin(0);
		  double qmax=(Vi-Vmin)/3600;
		  double qmaxD=R->getQmax();
		  
		  //cas avec turbine
		 // if(qmaxD>=0)
		  //{
			  if(R->getNbTurbines()>0)
			  {
				  	 
				  int turbine=R->getTurbine(0);
				  int Int=systeme->getTurbine(turbine)->getIntervalle(Vinit);
				if(qmaxD>0)  qmaxD=qmaxD+systeme->getTurbine(turbine)->getQMax(Int);
				else qmaxD=systeme->getTurbine(turbine)->getQMax(Int);
			  }
		 	if(qmax>qmaxD&&qmaxD>qmin) qmax=qmaxD;
		//  }
		  if(qmax>qTot[i])qmax=qTot[i];
		if(qmax>qTot[i]-qminC*(nbHeures-1))qmax=qTot[i]-qminC*(nbHeures-1);
		//TRUC REvolutionaire? non :(
		if(qminC>0&&nbP==0)
		{
			for(t=1;t<nbHeures-1;t++)
			{
					double Vmint=R->getVmin(t)/3600;
				if(V[t][i]/3600-t*qminC-Vmint<qmax)qmax=V[t][i]/3600-t*qminC-Vmint;
			}
		}
		//	 if(qmax<qmin)qmin=qmax;
		  //choix aléatoire
//			cout<<"random : "<<qmin<<" "<<qminC<<" "<<((Vi-Vmax)/3600)<<" "<<qmax<<" "<<qmaxD<<" "<<(Vi-Vmin)/3600<<" "<<qTot[i]<<" "<<endl;
		  eoUniformGenerator<double> random(qmin,qmax);
//		cout<<"apres random"<<endl;
		  quantite[i]=random();
	  }
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
				/*if(marche==1 && R->getNbTurbines()>0)
				{
						  	 
						  int turbine=R->getTurbine(0);
						
						  qminC=qminC+systeme->getTurbine(turbine)->getQmin(Vinit);
					  }*/
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
				/*if(qmax<qmin)
				{
					qmin=qmax;
					_genotype.setQuantite(h-1,i,qmin-qminC);
					int k=h-1;
					while(_genotype.getQuantite(k-1,i)>_genotype.getQuantite(k,i))
					{

						_genotype.setQuantite(k-1,i,_genotype.getQuantite(k,i)-qminC); k--;
					}	
				}*/
				 

				  //calcul de qmax
			
				  double qmaxD=R->getQmax();
				
				  //cas avec turbine
				  //if(qmaxD>=0)
				  //{
					  if(R->getNbTurbines()>0)
					  {
						  	 
						  int turbine=R->getTurbine(0);
						  int Int=systeme->getTurbine(turbine)->getIntervalle(Vinit);
						 // if(qmaxD>0)qmaxD=qmaxD+systeme->getTurbine(turbine)->getQMax(Int);
						//else 
						qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
					  }
				 	if(qmax>qmaxD+_genotype.getQuantite(h-1,i)&&qmaxD+_genotype.getQuantite(h-1,i)>=qmin&&qmaxD>0) qmax=qmaxD+_genotype.getQuantite(h-1,i);
				/*  else {
					if(R->getQmax()>0 && qmaxD+_genotype.getQuantite(h-1,i)<qmin)cout<<"P2 EN PERSPECTIVE "<<i<<" "<<qmaxD+_genotype.getQuantite(h-1,i)<<" "<<qmin<<" "<<(Vi-Vmax)/3600<<" "<<qmaxD<<" Heure "<<h<<endl;}*/
					if(qmax>qTot[i]-qminC*(nbHeures-1-h)&&qTot[i]-qminC*(nbHeures-1-h)>=qmin)qmax=qTot[i]-qminC*(nbHeures-1-h);
				else{
					
					if(qTot[i]-qminC*(nbHeures-1-h)<qmin &&h<8000){//cout<<"ini1 qmin "<<qmin<<" qtot-qminC "<<qTot[i]-qminC*(nbHeures-1-h)<<" i "<<i<<" qminC "<<qminC<<" Vi-Vmax "<<(Vi-Vmax)/3600<<" qte+qminC "<<_genotype.getQuantite(h-1,i)+qminC<<" h "<<h<<endl;
					qmax=qmin;}
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
				//Correction:
					double Vh=Vinit;
					if(V[h-1][i]-(qmax-qminC)*3600<Vmax &&_genotype.getQuantite(h-1,i)-qmax+qminC<_genotype.getQuantite(h-1,i)-_genotype.getQuantite(h-2,i))
					{
						//cout<<i<<" "<<h<<endl;
						
						_genotype.setQuantite(h-1,i,qmax-qminC);
						//qmax=qmin;
					}
					
				}
				
				  if(qmax<qmin)qmin=qmax;
				//if(qmin<_genotype.getQuantite(h-1,i)+qminC) cout<<"Reservoir "<<i<<" heure "<<h<<" "<<_genotype.getQuantite(h-1,i)+qminC<<" "<<qmax<<"apport "<<R->getApport(h)<<endl;
				  //choix aléatoire
				  eoUniformGenerator<double> random(qmin,qmax);
			          eoUniformGenerator<double> rim(0,100);
				  double p=rim();
				  quantite[i]=random();
				  if(p<b1) quantite[i]=qmin;
			          if(p>b2)quantite[i]=qmax;
			//	if(p>b1&&p<b3)quantite[i]=(qmin+qmax)/2;
			  }
			  _genotype.adEtat(quantite);
	  }
	  //DernierEtat:
	  _genotype.adEtat(qTot);

	for(i=0;i<nbReservoirs;i++)
	{
		for(h=1;h<nbHeures;h++)
		{
			if(_genotype.getQuantite(h,i)<_genotype.getQuantite(h-1,i))cout<<"pb grave "<<h<<" "<<i<<endl;
		}
	}
		
    // END   Code of random initialization of an eocascade object
    _genotype.invalidate();	   // IMPORTANT in case the _genotype is old
//cout<<"FIN INI"<<endl;
  }

private:
// START Private data of an eocascadeInit object
vector< vector<double> >V;
Systeme* systeme;
int nbHeures;
int nbReservoirs;
 vector<double> qTot;
// END   Private data of an eocascadeInit object
};

#endif
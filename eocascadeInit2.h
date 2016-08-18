/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for EO objects initialization in EO
============================================
*/

#ifndef _eocascadeInit2_h
#define _eocascadeInit2_h

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
class eocascadeInit2: public eoInit<GenotypeT> {
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeInit2(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
  {
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
	  //cout<<"debut ini2"<<endl;
    
	  int i,j,h,t;

	  _genotype.setNbReservoirs(nbReservoirs);
	  vector<double> quantite;
	 
	  //DernierEtat:
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
				/*if(marche==1 && R->getNbTurbines()>0)
				{
						  	 
						  int turbine=R->getTurbine(0);
						
						  qminC=qminC+systeme->getTurbine(turbine)->getQmin(Vinit);
					  }*/
				   // double qmin=_genotype.getQuantite(h-1,i)+qminC;
				double Vmin=R->getVmin(h);
				 double Vmax=R->getVmax();
				//  if(qmin<(Vi-Vmax)/3600)
				double qmin=(Vi-Vmax)/3600;

				double qmax=_genotype.getQuantite(nbHeures-h-2,i)-qminC;
				if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
				/*if(i==1)
				{
					double Vh=Vi-(h+1)*qminC*3600;
					for(j=h-1;j>-1;j--)
					{
						Vh=Vh-R->getApport(j)*3600+3600*qminC;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600>qmin)qmin=(Vh-Vminh)/3600;
					}

				}No need??*/
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
						  int Int=systeme->getTurbine(turbine)->getNbInt()-1;
						 // if(qmaxD>0)qmaxD=qmaxD+systeme->getTurbine(turbine)->getQMax(Int);
						//else 
						qmaxD=systeme->getTurbine(turbine)->getQMax(Int)+qminC;
					  }
				 	if(qmin<_genotype.getQuantite(nbHeures-h-2,i)-qmaxD&&-qmaxD+_genotype.getQuantite(nbHeures-h-2,i)<=qmax&&qmaxD>0) qmin=-qmaxD+_genotype.getQuantite(nbHeures-h-2,i);
				//  else cout<<"P2 EN PERSPECTIVE "<<qmaxD+_genotype.getQuantite(h-1,i)<<" "<<qmin<<" "<<(Vi-Vmax)/3600<<" "<<qmaxD<<endl;//}
					if(qmin<qminC*(1+h)&&qminC*(1+h)<=qmax)qmin=qminC*(1+h);;
				/*else{
					if(i==1&&qTot[i]-qminC*(nbHeures-1-h)<qmin)cout<<" qmin "<<qmin<<" qtot-qminC "<<qTot[i]-qminC*(nbHeures-1-h)<<" i "<<i<<" qminC "<<qminC<<" Vi-Vmax "<<(Vi-Vmax)/3600<<" qte+qminC "<<_genotype.getQuantite(h-1,i)+qminC<<" h "<<h<<endl;
				}*/
					if(nbP==0 &&qminC>0)
									{
										for(t=0;t<h;t++)
										{
											if(V[t][i]/3600+(h-t)*qminC-Vmax/3600>qmin && V[t][i]/3600+(h-t)*qminC-Vmax/3600<qmax)qmin=V[t][i]/3600+(h-t)*qminC-Vmax/3600;
										}
									}
				if(qmax>qTot[i]) qmax=qTot[i];
				/*if(qmax<qmin)
				{
				//Correction:
					double Vh=Vinit;
					if(V[h-1][i]-(qmax-qminC)*3600<Vmax &&_genotype.getQuantite(h-1,i)-qmax+qminC<_genotype.getQuantite(h-1,i)-_genotype.getQuantite(h-2,i))
					{
						//cout<<i<<" "<<h<<endl;
						
						_genotype.setQuantite(h-1,i,qmax-qminC);
						//qmax=qmin;
					}
					//else cout<<_genotype.getQuantite(h-1,i)-qmax+qminC<<" "<<_genotype.getQuantite(h-1,i)-_genotype.getQuantite(h-2,i)<<endl;
				}*/
					
				  if(qmax<qmin)qmax=qmin;
				  //choix aléatoire
				  eoUniformGenerator<double> random(qmin,qmax);
			          eoUniformGenerator<double> rim(0,100);
				  double p=rim();
				  quantite[i]=random();
				  if(p<b1) quantite[i]=qmax;
			          if(p>b2)quantite[i]=qmin;
				  //if(p>b1&&p<b3)quantite[i]=(qmin+qmax)/2;
			  }
			  _genotype.adEtat(quantite);
	  }
	  
	//	cout<<"FIN INI"<<endl;
    // END   Code of random initialization of an eocascade object
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
	//tests croissance:
	/*for(h=1;h<nbHeures;h++)
	{
		for(j=0;j<nbReservoirs;j++)
		{
			if(_genotype.getQuantite(h-1,j)>_genotype.getQuantite(h,j)) 					cout<<"erreur "<<h<<" "<<j<<endl;
		}
	}*/
        _genotype.invalidate();	   // IMPORTANT in case the _genotype is old
	//cout<<"fin init2"<<endl;
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

/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for EO objects initialization in EO
============================================
*/

#ifndef _eocascadeInit4_h
#define _eocascadeInit4_h

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
class eocascadeInit4: public eoInit<GenotypeT> {
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeInit4(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures,int _hDeb,vector<double> _qDeb)
  //  eocascadeInit( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeInit object
	  V=_V;
	int i,j,h;
	//cout<<"ici"<<endl;
	  systeme=_systeme;
	nbReservoirs=systeme->getNbReservoirs();
	hDeb=_hDeb;
	qDeb=_qDeb;
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
	//cout<<"debut ini4"<<endl;
    // START Code of random initialization of an eocascade object
	  int i,j,h,t;

	  _genotype.setNbReservoirs(nbReservoirs);
	 double qte;
	 
	  //Calcul qTot;manque les qte parentes !
	 
	  //premier état:
	  for(i=0;i<nbReservoirs;i++)
	  {
		
		_genotype.setQuantite(hDeb,i,qDeb[i]);

	  }
	 // _genotype.adEtat(quantite);
	  //Autres états :
	eoUniformGenerator<double> cb1(0,50);
	double b1=cb1();
	eoUniformGenerator<double> cb2(70,95);
	double b2=cb2();
	eoUniformGenerator<double> cb3(b1,b2);
	double b3=cb3();
	  for(h=hDeb+1;h<nbHeures-1;h++)
	  {
		  //quantite.clear();
		  for(i=0;i<nbReservoirs;i++)
			  {
				Reservoir* R=systeme->getReservoir(i);
				//  quantite.push_back(0);
				  //calcul de l'interval dans lequel se trouve Vinit:
				// eoUniformGenerator<int> tMarche(0,1); 
				//	int marche=tMarche();
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
				 	if(qmax>qmaxD+_genotype.getQuantite(h-1,i)&&qmaxD+_genotype.getQuantite(h-1,i)>qmin&&qmaxD>0) qmax=qmaxD+_genotype.getQuantite(h-1,i);
				//  else cout<<"P2 EN PERSPECTIVE "<<qmaxD+_genotype.getQuantite(h-1,i)<<" "<<qmin<<" "<<(Vi-Vmax)/3600<<" "<<qmaxD<<endl;//}
					if(qmax>qTot[i]-qminC*(nbHeures-1-h)&&qTot[i]-qminC*(nbHeures-1-h)>=qmin)qmax=qTot[i]-qminC*(nbHeures-1-h);
				//else{
					//if(qTot[i]-qminC*(nbHeures-1-h)<qmin)cout<<"inoit4 qmin "<<qmin<<" qtot-qminC "<<qTot[i]-qminC*(nbHeures-1-h)<<" i "<<i<<" qminC "<<qminC<<" Vi-Vmax "<<(Vi-Vmax)/3600<<" qte+qminC "<<_genotype.getQuantite(h-1,i)+qminC<<" h "<<h<<endl;
				//}
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

				//if(p>b1&&p<b3)quantite[i]=(qmin+qmax)/2;
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
	//	cout<<"FIN INI"<<endl;
    // END   Code of random initialization of an eocascade object
    _genotype.invalidate();	   // IMPORTANT in case the _genotype is old
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

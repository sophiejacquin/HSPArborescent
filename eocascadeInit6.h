/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for EO objects initialization in EO
============================================
*/

#ifndef _eocascadeInit6_h
#define _eocascadeInit6_h

// include the base definition of eoInit
#include <eoInit.h>
#include <vector>
#include "Systeme.h"
#include "eocascadeInit4.h"
#include "eocascadeInit5.h"
#include"eocascade.h"
/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen
 *
 * There is NO ASSUMPTION on the class GenoypeT.
 * In particular, it does not need to derive from EO (e.g. to initialize
 *    atoms of an eoVector you will need an eoInit<AtomType>)
 */
template <class GenotypeT>
class eocascadeInit6: public eoInit<GenotypeT> {
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeInit6(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
  {
	int i,j,h;
    // START Code of Ctor of an eocascadeInit object
	  V=_V;
	  systeme=_systeme;
	  nbHeures= _nbHeures;
	for(i=0;i<systeme->getNbReservoirs();i++)
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
    // END   Code of Ctor of an eocascadeInit object
  }


  /** initialize a genotype
   *
   * @param _genotype  generally a genotype that has been default-constructed
   *                   whatever it contains will be lost
   */
  void operator()(GenotypeT & _genotype)
  {
	//cout<<"deb ini6"<<endl;
	int i,j;
	vector<double> qNul;
	//creation vecteur nul :
	for(i=0;i<systeme->getNbReservoirs();i++)
	{
		qNul.push_back(0);
	}
	//creation genotype :
	for(i=0;i<nbHeures-1;i++)
	{
		_genotype.adEtat(qNul);
	}
	//choix heure:
	eoUniformGenerator<int> choix(nbHeures/3,2*nbHeures/3);

	int h=choix();
	//qte :
	vector<double> qte;
	for(i=0;i<systeme->getNbReservoirs();i++)
	{
		double qmin= (V[h][i]-systeme->getReservoir(i)->getVmax())/3600;
		for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
		{
			int p=systeme->getReservoir(i)->getParents()[j];
			qmin=qmin+qte[p];
		}
		//cout<<qmin<<" "<<systeme->getReservoir(i)->getVmax()<<" "<<V[i][h]<<" "<<qTot[i]*3600<<endl;
		double qminC=systeme->getReservoir(i)->getQmin();
		if(qminC<0) qminC=0;
		if(qmin<qminC*h)qmin=h*qminC;
		double qmax= (V[h][i]-systeme->getReservoir(i)->getVmin(h))/3600;
		for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
		{
			int p=systeme->getReservoir(i)->getParents()[j];
			qmax=qmax+qte[p];
		}
		if(qmax>qTot[i]-qminC*(nbHeures-h))qmax=qTot[i]-qminC*(nbHeures-h);
		//cout<<qmin<<" "<<qmax<<" "<<h<<" "<<qTot[i]-qminC*(nbHeures-h)<<" "<<h*qminC<<endl;
		eoUniformGenerator<double> g(qmin,qmax);
		double q=g();
		qte.push_back(q);
	}
	//creation des inis :
	eocascadeInit4<eocascade<double> > init4(V,systeme,nbHeures,h,qte);
	 eocascadeInit5<eocascade<double> > init5 (V, systeme,nbHeures,h,qte);
	
	//ini :
	init4(_genotype);
	init5(_genotype);
        _genotype.invalidate();	
	//cout<<"fin ini6"<<endl;
	
  }

private:
// START Private data of an eocascadeInit object
vector< vector<double> >V;
Systeme* systeme;
int nbHeures;
vector<double> qTot;
// END   Private data of an eocascadeInit object
};

#endif

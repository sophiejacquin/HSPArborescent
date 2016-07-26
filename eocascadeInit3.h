/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for EO objects initialization in EO
============================================
*/

#ifndef _eocascadeInit3_h
#define _eocascadeInit3_h

// include the base definition of eoInit
#include <eoInit.h>
#include <vector>
#include "Systeme.h"
#include "eocascadeInit.h"
#include "eocascadeInit2.h"
#include "eocascadeInit6.h"
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
class eocascadeInit3: public eoInit<GenotypeT> {
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeInit3(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
  //  eocascadeInit( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeInit object
	  V=_V;
	  systeme=_systeme;
	  nbHeures= _nbHeures;
    // END   Code of Ctor of an eocascadeInit object
  }


  /** initialize a genotype
   *
   * @param _genotype  generally a genotype that has been default-constructed
   *                   whatever it contains will be lost
   */
  void operator()(GenotypeT & _genotype)
  {
	//cout<<"deb ini"<<endl;
	 eocascadeInit2<eocascade<double> > init2(V, systeme,nbHeures);
	eocascadeInit6<eocascade<double> > init6(V, systeme,nbHeures);
	 eocascadeInit<eocascade<double> > init1(V, systeme,nbHeures);

	eoUniformGenerator<double> choix(0,300);
	double c=choix();
	/*vector<double> qTot;
	for(i=0;nbR;i++)
	{
		double qmin=(V[c][i]-Vmax)/3600;
		if(qmin<c*qminC)qmin=c*qminC;
		qmax=(V[c][i]-Vmin)/3600;
		if(qmax>Q[i]-(nbHeures-c-1)*qminC) qmax=Q[i]-(nbHeures-c-1)*qminC;
	}*/
	if(c>250)init6(_genotype);
	else 
	{
	if(c<150) init2(_genotype);
		else 
		init1(_genotype);
	}

    _genotype.invalidate();	   // IMPORTANT in case the _genotype is old
	
  }

private:
// START Private data of an eocascadeInit object
vector< vector<double> >V;
Systeme* systeme;
int nbHeures;
// END   Private data of an eocascadeInit object
};

#endif

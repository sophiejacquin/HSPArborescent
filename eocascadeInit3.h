#ifndef _eocascadeInit3_h
#define _eocascadeInit3_h
#include <eoInit.h>
#include <vector>
#include "Systeme.h"
#include "eocascadeInit.h"
#include "eocascadeInit2.h"
#include "eocascadeInit6.h"
#include"eocascade.h"

template <class GenotypeT>
class eocascadeInit3: public eoInit<GenotypeT> {
public:

eocascadeInit3(vector< vector<double> > _V, Systeme* _systeme,int _nbHeures)
{

	V=_V;
	systeme=_systeme;
	nbHeures= _nbHeures;
}



void operator()(GenotypeT & _genotype)
{

	eocascadeInit2<eocascade<double> > init2(V, systeme,nbHeures);
	eocascadeInit6<eocascade<double> > init6(V, systeme,nbHeures);
	eocascadeInit<eocascade<double> > init1(V, systeme,nbHeures);
	eoUniformGenerator<double> choix(0,300);
	double c=choix();
	if(c>250){
		init6(_genotype);	
	}
	else 
	{
		if(c<150) init2(_genotype);
		else 
			init1(_genotype);
	}

        _genotype.invalidate();	   
  }

private:

vector< vector<double> >V;
Systeme* systeme;
int nbHeures;

};

#endif

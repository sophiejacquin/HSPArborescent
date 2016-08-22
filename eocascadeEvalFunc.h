#ifndef _eocascadeEvalFunc_h
#define _eocascadeEvalFunc_h
#define P1 1000000
#define P2 1000000
#include <stdexcept>
#include <fstream>
#include "eoEvalFunc.h"


template <class EOT>
class eocascadeEvalFunc : public eoEvalFunc<EOT>
{
public:
	eocascadeEvalFunc(Systeme* _systeme, vector< vector<double> > _V)
	{
		systeme=_systeme;
		V= _V;
	}
	void operator()(EOT & _eo)
	{
		if (_eo.invalid())
		{
			double fit=0;//_eo.getLast_fitness();	
			int h,i,j;
			for(h=0;h<_eo.getNbEtats();h++)
			{
				if(_eo.getModif(h))
				{
					double profit=0;
					for(i=0;i<_eo.getNbReservoirs();i++)
					{
					//calcul quantite :
						double qte;
						if(h==0)
						{
							qte=_eo.getQuantite(0,i);
						}
						else qte=_eo.getQuantite(h,i)-_eo.getQuantite(h-1,i);
					//qminC:
						double qminC=systeme->getReservoir(i)->getQmin();
						if(qminC<0)qminC=0;
						if(qte<qminC){
						if(qminC-qte>0.00001)
						{
			 	 			profit= profit- P1*(qminC-qte);
					
						}
					}
				
					else
					{
						double qteT=0;
						double Vi;
						int turbine=systeme->getReservoir(i)->getTurbine(0);
						if(systeme->getReservoir(i)->getNbTurbines()>0)
						{
						
							if(h==0) Vi=systeme->getReservoir(i)->getVinit();
							else
							{
								Vi=V[h-1][i]-_eo.getQuantite(h-1,i)*3600;
								for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
								{
									int parent=systeme->getReservoir(i)->getParents()[j];
									Vi=Vi+_eo.getQuantite(h-1,parent)*3600;
								}
									
							}
							int Int=systeme->getTurbine(turbine)->getIntervalle(Vi);
							double qmaxT=systeme->getTurbine(turbine)->getQMax(Int);
							double qminT=systeme->getTurbine(turbine)->getQmin(Vi);
							//calcul qteT :
							qteT=qte-qminC;
							if(qteT>qmaxT)qteT=qmaxT;
							if(qteT<qminT-0.000001)qteT=0;
							//calcl profit:
							profit=profit+systeme->getTurbine(turbine)->getBenefice(Vi,qteT,h);
						
						}
						double qmaxC=systeme->getReservoir(i)->getQmax();
						if(qmaxC>0 && qte-qteT>qmaxC+0.000001){
							 profit=profit-P2*(qte-qteT-qmaxC);//
						}
					}
				}
				_eo.setEval(h,profit);
				_eo.setModif(h,false);
			}
			fit=fit+_eo.getEval(h);
		}

		_eo.fitness(fit);
     	}
}

private:

  Systeme* systeme;
  vector< vector<double> > V;
};


#endif

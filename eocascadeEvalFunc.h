/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for evaluator in EO, a functor that computes the fitness of an EO
==========================================================================
*/

#ifndef _eocascadeEvalFunc_h
#define _eocascadeEvalFunc_h
#define P1 1000000
#define P2 1000000
// include whatever general include you need
#include <stdexcept>
#include <fstream>

// include the base definition of eoEvalFunc
#include "eoEvalFunc.h"

/**
  Always write a comment in this format before class definition
  if you want the class to be documented by Doxygen
*/
template <class EOT>
class eocascadeEvalFunc : public eoEvalFunc<EOT>
{
public:
	/// Ctor - no requirement
// START eventually add or modify the anyVariable argument
  eocascadeEvalFunc(Systeme* _systeme, vector< vector<double> > _V)
  //  eocascadeEvalFunc( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeEvalFunc object
	  systeme=_systeme;
	  V= _V;
    // END   Code of Ctor of an eocascadeEvalFunc object
  }

  /** Actually compute the fitness
   *
   * @param EOT & _eo the EO object to evaluate
   *                  it should stay templatized to be usable
   *                  with any fitness type
   */
  void operator()(EOT & _eo)
  {
	//cout<<"debut eval"<<endl;
    // test for invalid to avoid recomputing fitness of unmodified individuals
    if (_eo.invalid())
      {
	//cout<<_eo.getLast_fitness();
	double fit=0;//_eo.getLast_fitness();	
	int h,i,j;
	//bool realisable=true;
	// to hold fitness value
    // START Code of computation of fitness of the eocascade object
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
					//cout<<"h=0 qte "<<i<<" "<<_eo.getQuantite(0,i)<<endl;
					qte=_eo.getQuantite(0,i);
				}
				else qte=_eo.getQuantite(h,i)-_eo.getQuantite(h-1,i);
				//qminC:
				double qminC=systeme->getReservoir(i)->getQmin();
				if(qminC<0)qminC=0;
				//if(h>0)cout<<i<<" "<<qte<<" "<<qminC<<" "<<_eo.getQuantite(h,i)<<" "<<_eo.getQuantite(h-1,i)<<endl;
				if(qte<qminC){
					if(qminC-qte>0.00001)
					{
				/*	if(h<_eo.getNbEtats()-1&&_eo.getQuantite(h,i)+qminC-qte<=_eo.getQuantite(h+1,i)&&V[h][i]-(_eo.getQuantite(h,i)+qminC-qte)*3600>=systeme->getReservoir(i)->getVmin(h))
					{
						_eo.setQuantite(h,i,_eo.getQuantite(h,i)+qminC-qte);
						_eo.setModif(h+1,true);
					}//ajout avec deconcentration
					else*/  profit= profit- P1*(qminC-qte);
					//realisable=false;
			 //cout<<"P1 "<<h<<" "<<qte<<" "<<qminC<<" "<<qminC-qte<<endl;
					}
				}
				
				else
				{
					//calcul qteT:
					//cout<<"dans else "<<i<<endl;
					double qteT=0;
					double Vi;
					int turbine=systeme->getReservoir(i)->getTurbine(0);
					if(systeme->getReservoir(i)->getNbTurbines()>0)
					{
						//cout<<"il y a une turbine"<<endl;
						
					       // cout<<i<<": "<<turbine<<endl;
						//calcul qmaxT et qminT:
						//calcul de Vi:
						
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
						//cout<<"Vi "<<Vi<<endl;
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
					/*	if(h>0 &&h+1<_eo.getNbEtats()&& _eo.getQuantite(i,h)-qte+qteT+qmaxC>=_eo.getQuantite(i,h-1)&&V[i][h]-3600*(_eo.getQuantite(i,h)-qte+qteT+qmaxC)<=systeme->getReservoir(i)->getVmax())
						{
							_eo.setQuantite(i,h,_eo.getQuantite(i,h)-qte+qteT+qmaxC);
							_eo.setModif(h+1,true);
						}

						else*/ profit=profit-P2*(qte-qteT-qmaxC);//cout<<"P2"<<qte-qteT<<" "<<qmaxC<<" "<<qte-qteT-qmaxC<<" "<<h<<" "<<" "<<qteT<<" "<<systeme->getTurbine(turbine)->getQmin(Vi)<<endl;
						//realisable=false;
					}
				}
			}
			//fit=fit+profit-_eo.getEval(h);
			_eo.setEval(h,profit);
			_eo.setModif(h,false);
			
		}
		fit=fit+_eo.getEval(h);
	}
    // END   Code of computation of fitness of the eocascade object
	//_eo.setLast_fitness(fit);
	//if(realisable) cout<<"SOLUTION réalisable!!!! :D youpiii"<<endl; 
	_eo.fitness(fit);
      }
	//cout<<"fin eval"<<endl;
  }

private:
// START Private data of an eocascadeEvalFunc object
  //  varType anyVariable;		   // for example ...
  Systeme* systeme;
  vector< vector<double> > V;
// END   Private data of an eocascadeEvalFunc object
};


#endif

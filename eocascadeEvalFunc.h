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
void details(EOT & _eo)
	{
		
			cout<<"contraintes violées:"<<endl;
			int h,i,j;
			vector<vector <double> > VolumesReservoirs;
			vector<vector<double> > q;
			vector< vector<double> > reserve;
			for(h=0;h<_eo.getNbEtats()+1;h++)
			{
				
				vector<double> rh;
				vector<double> Vh;
				for(i=0;i<_eo.getNbReservoirs();i++)
				{
		
					double ri;
					double Vi;
					if(h==0) {
						Vi=systeme->getReservoir(i)->getVinit();
						ri=_eo.getQuantite(h,i)*3600;
					}
					
					else{
						if (h<_eo.getNbEtats()) ri=(_eo.getQuantite(h,i)-_eo.getQuantite(h-1,i))*3600;
						Vi=V[h-1][i]-_eo.getQuantite(h-1,i)*3600;
						for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
						{
							int parent=systeme->getReservoir(i)->getParents()[j];
							Vi=Vi+_eo.getQuantite(h-1,parent)*3600;
						}

									
					} 
					if(Vi<systeme->getReservoir(i)->getVmin(h)-0.01)
						cout<<"Contrainte de volume min non respectée pour le reservoir "<<i<<" à l'heure "<<h<<" écart ="<<Vi-systeme->getReservoir(i)->getVmin(h)<<endl;
					if(Vi>systeme->getReservoir(i)->getVmax()+0.01)
						cout<<"Contrainte de volume max non respectée pour le reservoir "<<i<<" à l'heure "<<h<<" écart ="<<Vi-systeme->getReservoir(i)->getVmax()<<endl;
					Vh.push_back(Vi);
					rh.push_back(ri);
				}
				reserve.push_back(rh);
				VolumesReservoirs.push_back(Vh);
				if (h<_eo.getNbEtats())
				{
					vector<double> qh;
					for(i=0; i<systeme->getNbTurbines(); i++)
					{
						double qi;
						double qte;
						int r=systeme->getTurbine(i)->getReservoirParent();
						if(h==0)
						{
							qte=_eo.getQuantite(0,r);
						}
						else qte=_eo.getQuantite(h,r)-_eo.getQuantite(h-1,r);
						int Int=systeme->getTurbine(i)->getIntervalle(VolumesReservoirs[h][r]);
						double qmaxT=systeme->getTurbine(i)->getQMax(Int);
						double qminT=systeme->getTurbine(i)->getQmin(VolumesReservoirs[h][r]);
						double qminC=systeme->getReservoir(r)->getQmin();
						if(qminC<0)qminC=0;
						qi=qte-qminC;
						if(qi>qmaxT)qi=qmaxT;
						if(qi<qminT)qi=0;
						reserve[h][r]-=qi*3600;
						qh.push_back(qi*3600);
					
					
					}
					q.push_back(qh);
				}
				for(i=0;i<_eo.getNbReservoirs();i++)
				{
					if(reserve[h][i]<systeme->getReservoir(i)->getQmin()-0.01)
						cout<<"quantite min reserve non respectee pour le reservoir "<<i<<" à l'heure "<<h<<" ecart "<<reserve[h][i]-systeme->getReservoir(i)->getQmin()<<endl;
					if(reserve[h][i]>systeme->getReservoir(i)->getQmax()+0.01)
						cout<<"quantite max reserve non respectee pour le reservoir "<<i<<" à l'heure "<<h<<" ecart "<<reserve[h][i]-systeme->getReservoir(i)->getQmax()<<endl;
				}
			}
			
			
							


			_eo.setV(VolumesReservoirs);
			_eo.setReserve(reserve);
			_eo.setQ(q);
			
				
     	}

  void operator()(EOT & _eo)
  {
	//cout<<"debut eval"<<endl;
    // test for invalid to avoid recomputing fitness of unmodified individuals
    if (_eo.invalid())
      {
	//cout<<_eo.getLast_fitness();
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
						  profit= profit- P1*(qminC-qte);
					
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

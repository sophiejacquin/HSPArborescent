/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is useful in Emacs-like editors
 */

/*
Template for simple mutation operators
======================================
*/

#ifndef eocascadeMutation3_H
#define eocascadeMutation3_H
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include<ilconcert/ilomodel.h>
#include <iostream>
#include <eoOp.h>
#include "Sommet.h"
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen
 *
 * THere is NO ASSUMPTION on the class GenoypeT.
 * In particular, it does not need to derive from EO
 */
template<class GenotypeT>
class eocascadeMutation3: public eoMonOp<GenotypeT>
{
public:
  /**
   * Ctor - no requirement
   */
// START eventually add or modify the anyVariable argument
  eocascadeMutation3(vector< vector<double> > _V, Systeme* _systeme,int _nbH,vector<double> _delta, vector<double> _nbDeltas)
  //  eocascadeMutation( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeEvalFunc object
	  V=_V;
	  systeme=_systeme;
	  nbH=_nbH;
	  delta=_delta;
	 nbDeltas=_nbDeltas;

    // END   Code of Ctor of an eocascadeEvalFunc object
  }

  /// The class name. Used to display statistics
  string className() const { return "eocascadeMutation4"; }

  /**
   * modifies the parent
   * @param _genotype The parent genotype (will be modified)
   */
  bool operator()(GenotypeT & _genotype)
  {
	int i,h;
	int nbReservoirs=systeme->getNbReservoirs();
	cout<<"debut mutation4"<<endl;
      bool isModified(true);//a tester
    // START code for mutation of the _genotype object
	//Declarations :
	vector<vector<Sommet2*> > graphe;
      //choix de hDeb :
	eoUniformGenerator<int> choixHdeb(0,_genotype.getNbEtats()-nbH);
	int hDeb=choixHdeb();
      //etat ini :
	vector<double> null;
	if(hDeb>0)
	{
		for(i=0;i<nbReservoirs;i++)
		{
			null.push_back(_genotype.getQuantite(hDeb-1,i));
		}
	}
	else
	{
		for(i=0;i<nbReservoirs;i++)
		{
			null.push_back(0);
		}
		
	}
	Sommet2* racine=new Sommet2(null,0);
	vector<Sommet2*> tab0;
	tab0.push_back(racine);
	graphe.push_back(tab0);
	//Remplissage des tables :
	cout<<"remplissage des tables :"<<endl;
	for(h=1;h<nbH;h++)
	{
		
		Sommet2* deb=new Sommet2();
			
		vector<Sommet2*> tab;
		tab.push_back(deb);
		//Calcul de V :l;
		remplir(tab,graphe[h-1],0,h+hDeb,hDeb+nbH-1,_genotype);
		cout<<"heure "<<h<<" "<<tab.size()<<endl;
		graphe.push_back(tab);
	}
	//restitution résultats:
	Sommet2 * s=graphe[nbH-1][0];
	for(i=hDeb+nbH;i<hDeb;i--)
		{
			
			for(int n=0;n<nbReservoirs;n++)
			{
				_genotype.setQuantite(i,n,s->contenu[n]);
				
			}
			_genotype.setEval(i,s->valeur);
			s=s->pred;
		}
	cout<<"fin mutation4"<<endl;
    return isModified;

    // END code for mutation of the _genotype object
  }
void remplir(vector<Sommet2*> & t, vector<Sommet2*> & tpred,int n, int h,int hFin,GenotypeT & _genotype)
{
	//Declarations :
	int i,j,k,l,m;
	bool trouve;
	double Vmin=systeme->getReservoir(n)->getVmin(h);
	double qmaxConduite=systeme->getReservoir(n)->getQmax();
	Reservoir* R=systeme->getReservoir(n);
	int nbParents=R->getNbParents();
	double Vmax=systeme->getReservoir(n)->getVmax();
	double qminC=systeme->getReservoir(n)->getQmin();
	double qTot= _genotype.getQuantite(hFin,n);
	double q;
	int nbReservoirs=systeme->getNbReservoirs();
	//Recopie de t dans tc :
	vector<Sommet2*> tc=t;
	int tcSize=tc.size();
	//cout<<"reservoir "<<n<<" "<<tcSize<<endl;
	//vider t :
	t.clear();
	//Reconstruction : 
	for(i=0;i<tcSize;i++)
	{
			//cout<<"i "<<i<<endl;
		int tcTaille=tc[i]->contenu.size();
		//calcul de vR:
		double vR=V[h][n];
		for(j=0;j<nbParents;j++)
		{
			int parent = R->getParents()[j];
			vR=vR+tc[i]->contenu[parent]*3600;
		}
		//calcul qmin qmax et pas :
		double pas=delta[n]/nbDeltas[n];
		if (pas==0) pas=0.01;
		
		double qmin=0;
			
		if(qmin<(vR-Vmax)/3600)
		{
			qmin=(vR-Vmax)/3600;
		}
		trouve=false;
		double cont=delta[n];
		//cout<<"debut while1"<<endl;
		while(trouve==false && cont>-delta[n]-pas)
		{
			//if(n==4&&i==3)cout<<qmin<<" "<<_genotype.getQuantite(h,n)-cont<<" pas "<<pas<<" j "<<cont<<endl;
			if(qmin<=_genotype.getQuantite(h,n)-cont)
			{
				//cout<<"ici"<<endl;
				trouve=true;
				qmin=_genotype.getQuantite(h,n)-cont;

			}
			else cont=cont-pas;
		}
		double qmax=(vR-Vmin)/3600;
		if(qmax> qTot) qmax=qTot;
		if(trouve==true)
		{
		trouve=false;
		cont=delta[n];
		//cout<<"debut while"<<endl;
		while(trouve==false && cont>-delta[n]-pas)
		{
			//cout<<qmax<< " "<<_genotype.getQuantite(h,n)+j<<" pas "<<pas<<" j "<<j<<endl;
			if(qmax>=_genotype.getQuantite(h,n)+cont)
			{
				trouve=true;
				qmax=_genotype.getQuantite(h,n)+cont;
			}
			else cont=cont-pas;
		}
		}
		if(h==hFin){qmin=qTot;qmax=qTot;}
		//Remplissage :
		//cout<<"reservoir "<<n<<" qmin "<<qmin<<" qmax "<<qmax<<" pas "<<pas <<" tcTaile "<<tcTaille<<endl;
		if(trouve==true)
		{
		for(q=qmin;q<=qmax;q=q+pas)
		{
			//cout<<"q "<<q<<endl;
			Sommet2* s=new Sommet2();
			for(j=0;j<tcTaille;j++)
			{
				s->contenu.push_back(tc[i]->contenu[j]);
			}
			s->contenu.push_back(q);
			bool trouve=true;
			//cout<<"là"<<endl;
			if(n==nbReservoirs-1)
			{
				trouve=false;
				//then add a links:
				for(k=0;k<tpred.size();k++)
				{
					//regarder si lien est possible:
					bool pos=true;
					l=0;
					while(pos==true && l<nbReservoirs)
						{
							if(tpred[k]->contenu[l]>s->contenu[l]){ pos=false;}
							l++;
						}
					//cout<<"apres while"<<endl;
					//si lien possible alors calculer son prix :	
					if(pos==true)
					{
						double eval=0;
						trouve=true;
						for(l=0;l<nbReservoirs;l++)
						{
							Reservoir* Rl=systeme->getReservoir(l);
							double conT=tpred[k]->contenu[l];
							//calcul qteT et qteC
							double qte=s->contenu[l]-conT;
							double qminCl=Rl->getQmin();
							if(qminCl<0) qminCl=0;
							double qteT= qte-qminCl;
							if(qteT<0)
							{
								eval=eval+qteT*P1;
								qteT=0;
							}
							else
							{
								//calcul vini:
								
								//calcul qminT et qmaxT:
								if(Rl->getNbTurbines()>0){
									double vini=V[h][l]-conT*3600;
									for(m=0;m<Rl->getNbParents();m++)
									{
										int parent=Rl->getParents()[m];
										vini=vini+tpred[k]->contenu[parent]*3600;
									}
									int tourbine =Rl->getTurbine(0);
									Turbine* Tl=systeme->getTurbine(tourbine);
									double qminTurbine=Tl->getQmin(vini);
									if(qminTurbine<0)qminTurbine=0;
									int inter=Tl->getIntervalle(vini);
									double qmaxTurbine=Tl->getQMax(inter);
									if(qteT<qminTurbine) qteT=0;
									else
									{
										if(qmaxTurbine>0&&qteT>qmaxTurbine) qteT=qmaxTurbine;
									}
									//calcul profit:
									eval=eval+Tl->getBenefice(vini,q,h);
								}
								else qteT=0;
								//respect borne max:
								double qmaxCl=Rl->getQmax();
								double qteC=qte-qteT;
								if(qmaxCl>0&& qteC>qmaxCl) eval=eval-P2*(qteC-qmaxCl);
								
							}
							
						
						}
						if(s->pred==NULL ||s->valeur<eval)
						{
							s->valeur=eval;
							s->pred=tpred[k];
						}
						
					}
					
				}
				//cout<<"dans cette boucle?"<<endl;
			}
			//cout<<"ici"<<endl;
			if(trouve==true) t.push_back(s); //si sommet raccordé au graphe on l'ajoute ! ou si sommet incomplet
			else delete s; //sinon on l'élimine
			//cout<<"apres delete"<<endl;

			
		}
		}
			
	}
	//iteration:
			if(n<nbReservoirs-1)
			{
				//cout<<"on va iterer"<<endl;
				remplir(t,tpred,n+1,h,hFin,_genotype);
			}

}

private:
// START Private data of an eocascadeMutation object
 	vector< vector<double> > V;
 	Systeme* systeme;
 	int nbH;
	vector< double> delta;
	vector<double> nbDeltas;

		   // for example ...
// END   Private data of an eocascadeMutation object
};

#endif

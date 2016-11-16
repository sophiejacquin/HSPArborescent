#ifndef eocascadeMutation_H
#define eocascadeMutation_H


#include <eoOp.h>
template<class GenotypeT>
class eocascadeMutation: public eoMonOp<GenotypeT>
{
public:
  
  eocascadeMutation(vector< vector<double> > _V, Systeme* _systeme)
  {
	  V=_V;
	  systeme=_systeme;
  }

  string className() const { return "eocascadeMutation"; }


  bool operator()(GenotypeT & _genotype)
  {
cout<<"deb mut"<<endl;
      bool isModified(true);
      int i,j;
      int nbHeures=_genotype.getNbEtats();
      int nbReservoirs=_genotype.getNbReservoirs();  
      eoUniformGenerator<int> choixEtat(0,nbHeures-2);
      eoUniformGenerator<double> exacte(0,1);
      int etat=choixEtat();
      for(i=0;i<nbReservoirs;i++)
      {
	  cout<<"mut i="<<i<<endl;
	  double ex=exacte();
	  Reservoir* R = systeme->getReservoir(i);
    	  double Vmin= R->getVmin(etat);
	  cout<<"mut Vmin="<<Vmin<<endl;
    	  double Vmax= R->getVmax();
    	  double qmin,qmax;
    	  //calcul qmin:
    	  if(etat>0) qmin=_genotype.getQuantite(etat-1,i);
    	  else qmin=0;
    	  double qminC= R->getQmin();
    	  double Vi=V[etat][i];
    	  for(j=0;j<R->getNbParents();j++)
    	  {
    	  	int parent=R->getParents()[j];
    	  	Vi=Vi+_genotype.getQuantite(etat,parent)*3600;
    	  }
	  if(ex<0.99 &&qminC>0&&(Vi-Vmin)/3600>=qmin+qminC)qmin=qmin+qminC;//chercher à ne pas faire à tout les coups?
    	  if((Vi-Vmax)/3600>qmin)qmin=(Vi-Vmax)/3600;
	  int Rdeversement=R->getDeversement();
	  cout<<"DEVERSEMENT "<<Rdeversement<<endl;
	  if(Rdeversement>0)
	  {
	  	double Vsuc= V[etat][Rdeversement] -_genotype.getQuantite(etat,Rdeversement)*3600;
	  	if(qmin<systeme->getReservoir(Rdeversement)->getVmin(etat)-Vsuc)
			qmin=systeme->getReservoir(Rdeversement)->getVmin(etat)-Vsuc;
	  }	
    	  //calcul de qmax :
    	  qmax=_genotype.getQuantite(etat+1,i);
    	  double qmaxC=R->getQmax();
    	  if(qmaxC>0)
    	  {
    		  if(R->getNbTurbines()>0)
    		  {
    			  int turbine=R->getTurbine(0);
    			  double Vini;
    			  if(etat==0) Vini=R->getVinit();
    			  else
    			  {
    				  Vini=V[etat-1][i]-_genotype.getQuantite(etat-1,i);
    				  for(j=0;j<R->getNbParents();j++)
    				  {
    					 int parent=R->getParents()[j];
    					  Vini=Vini+_genotype.getQuantite(etat-1,parent)*3600;
    				  } 
    			  }
    			  int Int=systeme->getTurbine(turbine)->getIntervalle(Vini);
    			  qmaxC= qmaxC+systeme->getTurbine(turbine)->getQMax(Int);
    		  }
    		  if(ex<0.9&&etat>0&&qmax>qmaxC+_genotype.getQuantite(etat-1,i))qmax=qmaxC+_genotype.getQuantite(etat-1,i);
		if(etat==0 &&qmax>qmaxC) qmax=qmaxC;
    	  }
	  cout<<"mut ap if"<<endl;
    	  if((Vi-Vmin)/3600<qmax) qmax=(Vi-Vmin)/3600;
	  if(ex<0.8&&qmax>_genotype.getQuantite(etat+1,i)-qminC&&_genotype.getQuantite(etat+1,i)-qminC>=qmin)qmax=_genotype.getQuantite(etat+1,i)-qminC;//chercher à ne pas faire à tout les coups?
	if(Rdeversement>0)
	  {
	  	double Vsuc= V[etat][Rdeversement] -_genotype.getQuantite(etat,Rdeversement)*3600;
	 	if(qmax>systeme->getReservoir(Rdeversement)->getVmax()-Vsuc)
			qmax=systeme->getReservoir(Rdeversement)->getVmax()-Vsuc;
	}
    	  if(qmin>qmax){
			
			qmin=qmax;
	  }
    	  //Choix aléatoire d'une nouvelle quantite:
    	   eoUniformGenerator<double> choix(qmin,qmax);
    	   double qte=choix();
    	   _genotype.setQuantite(etat,i,qte);
	   _genotype.setModif(etat,true);
	   _genotype.setModif(etat+1,true);
    	
      }
cout<<"fin mut"<<endl;
      return isModified;
  
  }

private:

 	vector< vector<double> > V;
 	Systeme* systeme;

};

#endif

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

      bool isModified(true);
      int i,j;
      int nbHeures=_genotype.getNbEtats();
      int nbReservoirs=_genotype.getNbReservoirs();  
      eoUniformGenerator<int> choixEtat(0,nbHeures-2);
      eoUniformGenerator<double> exacte(0,1);
      int etat=choixEtat();
      for(i=0;i<nbReservoirs;i++)
      {
	  double ex=exacte();
	  Reservoir* R = systeme->getReservoir(i);
    	  double Vmin= R->getVmin(etat);
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
    	  if((Vi-Vmin)/3600<qmax) qmax=(Vi-Vmin)/3600;
	  if(ex<0.8&&qmax>_genotype.getQuantite(etat+1,i)-qminC&&_genotype.getQuantite(etat+1,i)-qminC>=qmin)qmax=_genotype.getQuantite(etat+1,i)-qminC;//chercher à ne pas faire à tout les coups?
	
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
      return isModified;
    // END code for mutation of the _genotype object
  }

private:

 	vector< vector<double> > V;
 	Systeme* systeme;

};

#endif

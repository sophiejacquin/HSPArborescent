/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is useful in Emacs-like editors
 */

/*
Template for simple mutation operators
======================================
*/

#ifndef eocascadeMutation_H
#define eocascadeMutation_H


#include <eoOp.h>

/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen
 *
 * THere is NO ASSUMPTION on the class GenoypeT.
 * In particular, it does not need to derive from EO
 */
template<class GenotypeT>
class eocascadeMutation: public eoMonOp<GenotypeT>
{
public:
  /**
   * Ctor - no requirement
   */
// START eventually add or modify the anyVariable argument
  eocascadeMutation(vector< vector<double> > _V, Systeme* _systeme)
  //  eocascadeMutation( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeEvalFunc object
	  V=_V;
	  systeme=_systeme;
    // END   Code of Ctor of an eocascadeEvalFunc object
  }

  /// The class name. Used to display statistics
  string className() const { return "eocascadeMutation"; }

  /**
   * modifies the parent
   * @param _genotype The parent genotype (will be modified)
   */
  bool operator()(GenotypeT & _genotype)
  {
	//cout<<"debut mutation"<<endl;
      bool isModified(true);
    // START code for mutation of the _genotype object
      //choix aléatoire de l'état à modifier:
      int i,j;
	   int nbHeures=_genotype.getNbEtats();
      int nbReservoirs=_genotype.getNbReservoirs();
	/*for(i=0;i<nbReservoirs;i++){
		for(int h=1;h<nbHeures;h++)
		{
			if(_genotype.getQuantite(h-1,i)>_genotype.getQuantite(h,i))cout<<"avant mut i : "<<i<<" h "<<h<<" qte "<<_genotype.getQuantite(h,i)<<" qte pred "<<_genotype.getQuantite(h-1,i)<<" "<<endl;

		}}*/
   
      eoUniformGenerator<int> choixEtat(0,nbHeures-2);
	eoUniformGenerator<double> exacte(0,1);
      int etat=choixEtat();
      //eoUniformGenerator<int> choixSup(etat+1,nbHeures-2);
      //int eSup=choixSup();
     // for(int e=etat;e<etat+500 && e<nbHeures-2;e++)
    	 // if(_genotype.getEval(e)<0 && _genotype.getEval(e)<_genotype.getEval(etat))etat=e;
	//cout<<"etat "<<etat<<endl;
      for(i=0;i<nbReservoirs;i++)
      {
	double ex=exacte();
	Reservoir* R = systeme->getReservoir(i);
    	  double Vmin= R->getVmin(etat);
    	  double Vmax= R->getVmax();
    	  double qmin,qmax;
    	  //calcul qmin:
    	  if(etat>0)qmin=_genotype.getQuantite(etat-1,i);
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
	//if((Vi-Vmin)/3600<qmin) cout<<"pb respect vmin imp"<<endl;
    	  if(qmin>qmax){
			//if(etat>0)cout<<"mut qmin>qmax "<<i<<" "<<qmin<<" "<<qmax<<" "<<(Vi-Vmax)/3600<<" "<<_genotype.getQuantite(etat-1,i)<<" "<<(Vi-Vmin)/3600<<" "<<_genotype.getQuantite(etat+1,i)<<endl;
			qmin=qmax;
			//isModified=false;
		}
    	  //Choix aléatoire d'une nouvelle quantite:
    	   eoUniformGenerator<double> choix(qmin,qmax);
    	   double qte=choix();
    	   _genotype.setQuantite(etat,i,qte);
	   _genotype.setModif(etat,true);
	   _genotype.setModif(etat+1,true);
    	
      }
       /** Requirement
	* if (_genotype has been modified)
	*     isModified = true;
	* else
	*     isModified = false;
	*/
	/*for(i=0;i<nbReservoirs;i++){
		for(int h=1;h<nbHeures;h++)
		{
			if(_genotype.getQuantite(h-1,i)>_genotype.getQuantite(h,i))cout<<"mut i : "<<i<<" h "<<h<<" qte "<<_genotype.getQuantite(h,i)<<" qte pred "<<_genotype.getQuantite(h-1,i)<<" "<<etat<<endl;

		}}*/
	//
//cout<<"fin mutation"<<endl;
    return isModified;
    // END code for mutation of the _genotype object
  }

private:
// START Private data of an eocascadeMutation object
 	vector< vector<double> > V;
 	Systeme* systeme;
		   // for example ...
// END   Private data of an eocascadeMutation object
};

#endif

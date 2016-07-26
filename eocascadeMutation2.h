/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is useful in Emacs-like editors
 */

/*
Template for simple mutation operators
======================================
*/

#ifndef eocascadeMutation2_H
#define eocascadeMutation2_H


#include <eoOp.h>
typedef struct Sommet{
	double valeur;
	double quantite;
	int pred;
}Sommet;
/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen
 *
 * THere is NO ASSUMPTION on the class GenoypeT.
 * In particular, it does not need to derive from EO
 */
//#define P1 100000
//#define P2 100000
template<class GenotypeT>
class eocascadeMutation2: public eoMonOp<GenotypeT>
{
public:
  /**
   * Ctor - no requirement
   */
// START eventually add or modify the anyVariable argument
  eocascadeMutation2(vector< vector<double> > _V, Systeme* _systeme,int _longueur,int _nbPas)
  //  eocascadeMutation( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
    // START Code of Ctor of an eocascadeEvalFunc object
	  V=_V;
	  systeme=_systeme;
	  longueur=_longueur;
	  nbPas=_nbPas;
    // END   Code of Ctor of an eocascadeEvalFunc object
  }

  /// The class name. Used to display statistics
  string className() const { return "eocascadeMutation2"; }

  /**
   * modifies the parent
   * @param _genotype The parent genotype (will be modified)
   */
  bool operator()(GenotypeT & _genotype)
  {
	//cout<<"debut mutation2 "<<longueur<<endl;
      bool isModified(true);
    // START code for mutation of the _genotype object
      //déclarations:
      vector< vector<Sommet> > graphe;
      int i,j,h,hi,hf,r,rSuc;
      int nbReservoirs=_genotype.getNbReservoirs();
      int nbHeures=_genotype.getNbEtats();

     // cout<<"choix du reservoir:"<<endl;
      eoUniformGenerator<int> choixReservoir(1,nbReservoirs-1);
      r=choixReservoir();
	if(r==2)r=5;
	//cout<<nbReservoirs-1<<" R :"<<r<<endl;
      Reservoir* R=systeme->getReservoir(r);
	int bT=R->getNbTurbines();
	Turbine* turbine;
	 if(R->getNbTurbines()>0)
    	  {
    		 int t=R->getTurbine(0);
    		turbine=systeme->getTurbine(t);
    		
    	  }
      rSuc=R->getDeversement();
      Reservoir* Rsuc=systeme->getReservoir(rSuc);
      //choix hi hf
      eoUniformGenerator<int> choixHi(0,nbHeures-longueur);
      hi=choixHi();
     if(r==0) hf=hi+24;
	if(r==3&&hi+800<nbHeures)hf=hi+800;
	else hf=hi+longueur;
//cout<<"hi "<<hi<<" hf: "<<hf<<endl;
      //calcul pas:
      double Vmax=R->getVmax();
      double VmaxS=Rsuc->getVmax();
      double pas =(Vmax/3600)/nbPas;
      if(pas==0)pas=0.1;
	//cout<<hf<<endl;
      double qTot=_genotype.getQuantite(hf-1,r);
	//if(hi>0)cout<<"qte "<<qTot-_genotype.getQuantite(hi-1,r)<<endl;
      double qmaxC=R->getQmax();
      double qminC=R->getQmin();
	  double qminCs=Rsuc->getQmin();
      double qmaxCs=Rsuc->getQmax();
      if(qminC<0) qminC=0;
      double qmaxT=0;
     


    	  if(bT>0)
    	  {
    		  
    		  int derIntervalle=turbine->getNbInt()-1;
    		  qmaxT=turbine->getQMax(derIntervalle);
    	  }
    	  if(qmaxC>0&&(qmaxT+qmaxC)/longueur<pas)pas=(qmaxT+qmaxC)/longueur;
    	  int turbineS=-1;
	Turbine* Ts;
    	  if(Rsuc->getNbTurbines()>0)
    	  {
    		  turbineS=Rsuc->getTurbine(0);
			Ts=systeme->getTurbine(turbineS);

    	  }

//	cout<<"debut DP"<<endl;
      //Algorithme dynamique:
      double qmin=0;
      double qmax=0;
      if(hi>0){qmin=_genotype.getQuantite(hi-1,r);
      //qmax=_genotype.getQuantite(hi-1,r);

      }
	//cout<<"heures intermediaires"<<endl;
      for(h=hi;h<hf-1;h++)
      {
		//cout<<"heure "<<h<<endl;
    	  vector<Sommet> vec;
    	  //calcul qmin:
    	   //qmin=qmin;
	//	cout<<qmin<<endl;
    	   //respect vmax:
    	   double Vr=V[h][r];
    	   for(i=0;i<R->getNbParents();i++)
    	   {
    		   int p=R->getParents()[i];
    		   Vr=Vr+_genotype.getQuantite(h,p)*3600;
    	   }
    	   if((Vr-Vmax)/3600>qmin)qmin=(Vr-Vmax)/3600;
    	   //respect Vmin succ:
    	   double VminS=Rsuc->getVmin(h);
    	   double VrSuc=V[h][rSuc]-_genotype.getQuantite(h,rSuc)*3600;
    	   for(i=0;i<Rsuc->getNbParents();i++)
    	   {
    		   int p=Rsuc->getParents()[i];
    		   if(p!=r)VrSuc=VrSuc+_genotype.getQuantite(h,p)*3600;
    	   }
    	   if((VminS-VrSuc)/3600>qmin)qmin=(VminS-VrSuc)/3600;

    	  //calcul qmax:
    	 /* if(qmaxC>0)
    	  {
    		 qmax=qmax+qmaxC+qmaxT;
    		 if(qmax>qTot)qmax=qTot;
    	  }*/
    	  qmax=qTot;
    	  double Vmin=R->getVmin(h);
    	  if((Vr-Vmin)/3600<qmax) qmax=(Vr-Vmin)/3600;
    	  if((VmaxS-VrSuc)/3600<qmax)qmax=(VmaxS-VrSuc)/3600;
    	  if(qmin>qmax)//respect qminC impossible
    	  {
		//	cout<<"cas limite "<<qmin<<" "<<qmax<<" "<<(Vr-Vmin)/3600<<" "<<(VmaxS-VrSuc)/3600<<" "<<qTot<<endl;
    		  qmin=qmax;
		
    	  }
	//cout<<"qmin "<<qmin<<" qmax "<<qmax<<endl;
	//cout<<"? "<<qmin<<" "<<qmax<<" "<<(Vr-Vmin)/3600<<" "<<(VmaxS-VrSuc)/3600<<" "<<_genotype.getQuantite(h,r)<<" "<<(Vr-Vmax)/3600<<" "<<(VminS-VrSuc)/3600<<" "<<qTot<<endl;
    	  double qte;
    	  for(qte=qmin;qte<=qmax;qte=qte+pas)
    	  {
    		  //creer sommet :
    		  Sommet s;
    		  s.quantite=qte;
    		  s.pred=-1;
    		  s.valeur=-100000;
		if(qte>qTot) cout<<"wtf?"<<endl;
    		  //recherche pred :
    		  if(h==hi)
    		  {
    			  //qte deversée:
    			  double qtePred=0;
    			  if(hi>0)qtePred=_genotype.getQuantite(hi-1,r);
    			  double qteD= qte-qtePred;
			//cout<<"h "<<h<<" qteD "<<qteD<<" qte "<<qte<<" pas "<<pas<<" qmaxT "<<qmaxT<<" turbine :"<<turbine<<" qmaxC "<<qmaxC<<" qmax "<<qmax<<" Vmax "<<Vmax<<endl;
    			  if(qteD>=0)
    			  {
			//	cout<<"respect croissance"<<endl;
    				    //vini:
    				  double Vini;
    				  if(hi==0) Vini=R->getVinit();
    				  else{
    					  Vini=V[hi-1][r]-qtePred*3600;
    					  for(i=0;i<R->getNbParents();i++)
    					 {
    					     	int p=R->getParents()[i];
    					     	Vini=Vini+_genotype.getQuantite(h-1,p)*3600;
    					  }

    				  }
    			       //réalisabilité:
		//		cout<<qteD<<" "<<qminC<<endl;
    				  if(qteD>=qminC-0.00001)
    				  {
				//	cout<<"respect qminC"<<endl;
    					  double val=0;
    					  double qteT=0;
    					  if(bT>0)
    					  {

    				  //calcul qteT :
    					      qteT=qteD-qminC;
    				  //calcul qmaxT:
    						  int Int=turbine->getIntervalle(Vini);
    						  double qmaxTi=turbine->getQMax(Int);
    						  double qminTi=turbine->getQmin(Vini);
    						  if(qteT>qmaxTi)qteT=qmaxTi;
    						  if(qteT<qminTi)qteT=0;
    						  val=turbine->getBenefice(Vini,qteT,h);
    					  }
    			      //la le profit fait avec le successeur ne change pas donc pas besoin de le calculer
    			     //ajout de s?
    					  if(qmaxC<0||qteD-qteT<=qmaxC)
    					  {
    						  if(val>s.valeur)
    						  {
    							  s.valeur=val;
    							  vec.push_back(s);
							//cout<<"ici 1 tout respecté "<<val<<endl;
    						  }
    					  }
    				  }
    			  }

    		  }
    		  else
    		  {
    			  bool plusGrand=false;
    			  j=0;
		//cout<<"taille vecteur pred else"<<graphe[h-hi-1].size()<<endl;
    			  while(j<graphe[h-hi-1].size())
    			  {
				 double val=0;
    				  double qtePred=graphe[h-hi-1][j].quantite;
    				  double qteD=qte-qtePred;
				//cout<<"qteD "<<qteD<<" "<<qte<<" "<<qtePred<<endl;
    				  if(qteD>=0)
    				  {
				//	cout<<"respect croissance"<<endl;
					//cout<<qteD<<" "<<qminC<<endl;
    					  if(qteD>=qminC+0.000001)
    					  {
						//cout<<"respectqminC 2"<<endl;
    						 
    						  double qteT=0;
    						  if(bT>0)
    						  {
    							  	qteT=qteD-qminC;
    							  	//calcul Vini :
    							  	double Vini;

    							    Vini=V[h-1][r]-qtePred*3600;
    							  	for(i=0;i<R->getNbParents();i++)
    							  	 {
    							  	    	int p=R->getParents()[i];
    							  	    	Vini=Vini+_genotype.getQuantite(h-1,p)*3600;
    							  	  }
    							  	 int Int=turbine->getIntervalle(Vini);
    							  	  double qmaxTi=turbine->getQMax(Int);
    							  	  double qminTi=turbine->getQmin(Vini);
    							  	   if(qteT>qmaxTi)qteT=qmaxTi;
    							  	    if(qteT<qminTi)qteT=0;
    							  	    val=turbine->getBenefice(Vini,qteT,h);

    						  }
    						  //calcul profit avec rSuc:
    						  double qteDs=_genotype.getQuantite(h,rSuc)-_genotype.getQuantite(h-1,rSuc);
    						  double qteTs=0;
    						  if(turbineS>-1)
    						  {
    							 qteTs=qteDs-qminCs;
    							  double VrSucPred=V[h-1][rSuc]-(_genotype.getQuantite(h-1,rSuc)-qtePred)*3600 ;
    							  for(i=0;i<Rsuc->getNbParents();i++)
    							  {
    							      int p=Rsuc->getParents()[i];
    							       if(p!=r)VrSucPred=VrSuc+_genotype.getQuantite(h-1,p)*3600;
    							   }
    							  int Int=Ts->getIntervalle(VrSucPred);
    							  double qmaxTi=Ts->getQMax(Int);
    							  double qminTi=Ts->getQmin(VrSucPred);
    							  if(qteTs>qmaxTi)qteTs=qmaxTi;
    						     if(qteTs<qminTi)qteTs=0;
    						     val=val+Ts->getBenefice(VrSucPred,qteTs,h);
    						    // else val=val-10000000;

    						  }
							 if(qmaxC>0&&qteD-qteT>qmaxC)
							{
								val=val-P2*(qteD-qteT-qmaxC);    			
							}
							if(qmaxCs>0&&qteDs-qteTs>qmaxCs)
							{
								val=val-P2*(qteDs-qteTs-qmaxCs);
							}
			 /* if((qmaxC<0||qteD-qteT<=qmaxC)&&(qmaxCs<0||qteDs-qteTs<=qmaxCs))
    						   {
							//cout<<"respect qmaxC"<<endl;

    						      	if(val+graphe[h-hi-1][j].valeur>s.valeur)
    						      	{
								//cout<<"ici"<<endl;
    						      		s.valeur=val+graphe[h-hi-1][j].valeur;
    						      		s.pred=j;
							//	cout<<" j "<<j<<" "<<s.valeur<<endl;
    						      		//vec.push_back(s);
    						      	}
    						   }*/
    					  }
						else val=val-P1*(qteD-qminC);
					if(val+graphe[h-hi-1][j].valeur>s.valeur)//MODIF
    						      	{
								//cout<<"ici"<<endl;
    						      		s.valeur=val+graphe[h-hi-1][j].valeur;
    						      		s.pred=j;
							//	cout<<" j "<<j<<" "<<s.valeur<<endl;
    						      		//vec.push_back(s);
    						      	}
    				  }
    				 // else plusGrand=true;
    				  j++;
    			  }
			if(s.pred>-1)vec.push_back(s);

    		  }


    	  }
    	  //gestion sommet max

    	  //creation:
    	  Sommet s;
    	  s.quantite=_genotype.getQuantite(h,r);
    	  s.valeur=-100000000000000000;
    	  s.pred=-1;
	qte=_genotype.getQuantite(h,r);
	if(qte>qTot) cout<<"euh...?"<<endl;
    	  bool plusGrand=false;
	if(h==hi)
    		  {
    			  //qte deversée:
					double val=0;
    			  double qtePred=0;
    			  if(hi>0)qtePred=_genotype.getQuantite(hi-1,r);
    			  double qteD= qte-qtePred;
			//cout<<"h "<<h<<" qteD "<<qteD<<" qte "<<qte<<" pas "<<pas<<" qmax "<<qmax<<endl;
    			  if(qteD>=0)
    			  {
		//		cout<<"respect croissance qmax"<<endl;
    				    //vini:
    				  double Vini;
    				  if(hi==0) Vini=R->getVinit();
    				  else{
    					  Vini=V[hi-1][r]-qtePred*3600;
    					  for(i=0;i<R->getNbParents();i++)
    					 {
    					     	int p=R->getParents()[i];
    					     	Vini=Vini+_genotype.getQuantite(h-1,p)*3600;
    					  }

    				  }
    			       //réalisabilité:
    			//	  cout<<qteD<<" "<<qminC<<endl;
    				  if(qteD>=qminC-0.000001)
    				  {
		//			cout<<"respect qminC qmax"<<endl;

    					  double qteT=0;
    					  if(bT>0)
    					  {
						

    				  //calcul qteT :
    					      qteT=qteD-qminC;
    				  //calcul qmaxT:
    						  int Int=turbine->getIntervalle(Vini);
    						  double qmaxTi=turbine->getQMax(Int);
    						  double qminTi=turbine->getQmin(Vini);
    						  if(qteT>qmaxTi)qteT=qmaxTi;
    						  if(qteT<qminTi)qteT=0;
    						  val=turbine->getBenefice(Vini,qteT,h);
    					  }
    			      //la le profit fait avec le successeur ne change pas donc pas besoin de le calculer
    			     //ajout de s?
    					  if(qmaxC>0 && qteD-qteT>qmaxC) val=val-100000*( qteD-qteT-qmaxC)*( qteD-qteT-qmaxC);
    					  
    						  if(val>s.valeur)
    						  {
    							  s.valeur=val;
    							  vec.push_back(s);
						//	cout<<"ici 1"<<" "<<s.valeur<<endl;
    						  }
    					  
    				  }
    				  else
    				  {
    			//		  cout<<"respect qmin impossible"<<endl;
    					  val=val-P1*(qminC-qteD)*(qminC-qteD);
    					  if(val>s.valeur || s.pred<0)
    					  {
    						  s.valeur=val;
    						  vec.push_back(s);
    			//			  cout<<"ici"<<endl;
    					  }
    				  }
    			  }

    		  }
			else{
    	     			  j=0;
    	     			  while(j<graphe[h-hi-1].size())
    	     			  {
    	     				  double val=0;
    	     				  double qtePred=graphe[h-hi-1][j].quantite;
    	     				  double qteD=qte-qtePred;
					//	cout<<qteD<<" qtePred "<<qtePred<<" qte "<<qte<<endl;
    	     				  if(qteD>=0)
    	     				  {
					//	cout<<"respect positivité qtemax "<<qteD<<" "<<qminC<<endl;
    	     					  if(qteD>=qminC-0.000001)
    	     					  {
						//	cout<<"respect qminC qmax"<<endl;

    	     						  double qteT=0;
    	     						  if(bT>0)
    	     						  {
									
    	     							  	qteT=qteD-qminC;
    	     							  	//calcul Vini :
    	     							  	double Vini;

    	     							    Vini=V[h-1][r]-qtePred*3600;
    	     							  	for(i=0;i<R->getNbParents();i++)
    	     							  	 {
    	     							  	    	int p=R->getParents()[i];
    	     							  	    	Vini=Vini+_genotype.getQuantite(h-1,p)*3600;
    	     							  	  }
    	     							  	 int Int=turbine->getIntervalle(Vini);
    	     							  	  double qmaxTi=turbine->getQMax(Int);
    	     							  	  double qminTi=turbine->getQmin(Vini);
    	     							  	   if(qteT>qmaxTi)qteT=qmaxTi;
    	     							  	    if(qteT<qminTi)qteT=0;
    	     							  	    val=turbine->getBenefice(Vini,qteT,h);

    	     						  }
    	     						  //calcul profit avec rSuc:
    	     						  double qteDs=_genotype.getQuantite(h,rSuc)-_genotype.getQuantite(h-1,rSuc);
    	     						  double qteTs=0;
    	     						  if(turbineS>-1)
    	     						  {
								
    	     							 qteTs=qteDs-qminCs;
    	     							  double VrSucPred=V[h-1][rSuc]-(_genotype.getQuantite(h-1,rSuc)-qtePred)*3600 ;
    	     							  for(i=0;i<Rsuc->getNbParents();i++)
    	     							  {
    	     							      int p=Rsuc->getParents()[i];
    	     							       if(p!=r)VrSucPred=VrSuc+_genotype.getQuantite(h-1,p)*3600;
    	     							   }
    	     							  int Int=Ts->getIntervalle(VrSucPred);
    	     							  double qmaxTi=Ts->getQMax(Int);
    	     							  double qminTi=Ts->getQmin(VrSucPred);
    	     							  double val=0;	  if(qteTs>qmaxTi)qteTs=qmaxTi;
    	     						     if(qteTs<qminTi)qteTs=0;
    	     						     val=val+Ts->getBenefice(VrSucPred,qteTs,h);
    	     						    // else val=val-10000000;

    	     						  }
    	     						  if(qmaxC>0&&qteD-qteT>qmaxC)val=val-P2*(qteD-qteT-qmaxC);
    	     								  if(qmaxCs>0&&qteDs-qteTs>qmaxCs)val=val-P2*(qteDs-qteTs-qmaxCs);
    	     							//	  	cout<<val+graphe[h-hi-1][j].valeur<<" sval "<<s.valeur<<endl;
    	     						      	if(val+graphe[h-hi-1][j].valeur>s.valeur||s.pred<0)
    	     						      	{
    	     						      		s.valeur=val+graphe[h-hi-1][j].valeur;
    	     						      		s.pred=j;
    	     						      		
								//	cout<<"ici der "<<s.valeur<<endl;
    	     						      	}


    	     					  }
    	     					 else
    	     					 {
    	     					    val=val-P1*(qminC-qteD);
    	     					    if(val+graphe[h-hi-1][j].valeur>s.valeur||s.pred<0)
    	     					    {
    	     					    	   s.valeur=val+graphe[h-hi-1][j].valeur;
    	     					    	    s.pred=j;

    	     					    //	    cout<<"ici"<<endl;
    	     					    }
    	     					}
    	     				  }

    	     				  j++;
    	     			  }
				if(s.pred>-1)vec.push_back(s);
				else cout<<"pb de merde, was passiert?m2 "<<h<<" "<<hi<<" "<<hf<<" "<<_genotype.getQuantite(h,r)<<" "<<_genotype.getQuantite(h-1,r)<<" "<<graphe[h-hi-1][graphe[h-hi-1].size()-1].quantite<<endl;
			}

    	  //recherche pred avec possibilité de sol non realisable si vec vide:
	//cout<<"taille vecteur "<<vec.size()<<endl;
    	  graphe.push_back(vec);
      }
      //derniere heure:
	//cout<<"derniere heure"<<endl;
      Sommet s;
        	  s.quantite=qTot;
        	  s.valeur=-1000000000000000000;
        	  s.pred=-1;
        	  h=hf-1;
        	  //bool plusGrand=false;
		double qte=qTot;
        	     			  j=0;
        	     			  while(j<graphe[h-hi-1].size())
        	     			  {
        	     				  double val=0;
        	     				  double qtePred=graphe[h-hi-1][j].quantite;
        	     				  double qteD=qte-qtePred;
        	     				  if(qteD>=0)
        	     				  {
        	     					  if(qteD>=qminC)
        	     					  {

        	     						  double qteT=0;
        	     						  if(bT>0)
        	     						  {
        	     							  	qteT=qteD-qminC;
        	     							  	//calcul Vini :
        	     							  	double Vini;

        	     							    Vini=V[h-1][r]-qtePred*3600;
        	     							  	for(i=0;i<R->getNbParents();i++)
        	     							  	 {
        	     							  	    	int p=R->getParents()[i];
        	     							  	    	Vini=Vini+_genotype.getQuantite(h-1,p)*3600;
        	     							  	  }
        	     							  	 int Int=turbine->getIntervalle(Vini);
        	     							  	  double qmaxTi=turbine->getQMax(Int);
        	     							  	  double qminTi=turbine->getQmin(Vini);
        	     							  	   if(qteT>qmaxTi)qteT=qmaxTi;
        	     							  	    if(qteT<qminTi)qteT=0;
        	     							  	    val=turbine->getBenefice(Vini,qteT,h);

        	     						  }
        	     						  //calcul profit avec rSuc:
        	     						  double qteDs=_genotype.getQuantite(h,rSuc)-_genotype.getQuantite(h-1,rSuc);
        	     						  double qteTs=0;
        	     						  if(turbineS>-1)
        	     						  {
        	     							 qteTs=qteDs-qminCs;
        	     							  double VrSucPred=V[h-1][rSuc]-(_genotype.getQuantite(h-1,rSuc)-qtePred)*3600 ;
        	     							  for(i=0;i<Rsuc->getNbParents();i++)
        	     							  {
        	     							      int p=Rsuc->getParents()[i];
        	     							       if(p!=r)VrSucPred=VrSucPred+_genotype.getQuantite(h-1,p)*3600;
        	     							   }
        	     							  int Int=Ts->getIntervalle(VrSucPred);
        	     							  double qmaxTi=Ts->getQMax(Int);
        	     							  double qminTi=Ts->getQmin(VrSucPred);
        	     							  if(qteTs>qmaxTi)qteTs=qmaxTi;
        	     						     if(qteTs<qminTi)qteTs=0;
        	     						     val=val+Ts->getBenefice(VrSucPred,qteTs,h);
        	     						    // else val=val-10000000;

        	     						  }
        	     						  if(qmaxC>0&&qteD-qteT>qmaxC)val=val-100000*(qteD-qteT-qmaxC)*(qteD-qteT-qmaxC);
        	     								  if(qmaxCs>0&&qteDs-qteTs>qmaxCs)val=val-100000*(qteDs-qteTs-qmaxCs)*(qteDs-qteTs-qmaxCs);

        	     						      	if(val+graphe[h-hi-1][j].valeur>s.valeur)
        	     						      	{
        	     						      		s.valeur=val+graphe[h-hi-1][j].valeur;
        	     						      		s.pred=j;
        	     						      		
        	     						      	}


        	     					  }
        	     					  else{
        	     						 val=val-100000*(qminC-qteD)*(qminC-qteD);
        	     						 if(val+graphe[h-hi-1][j].valeur>s.valeur)
        	     						 {
        	     							 s.valeur=val+graphe[h-hi-1][j].valeur;
        	     						          s.pred=j;

        	     						    	     					    //	    cout<<"ici"<<endl;
        	     						    	     					    }
        	     					  }
        	     				  }
        	     				 /// else plusGrand=true;
        	     				  j++;
        	     			  }
      //modification genotype:
//	cout<<"modification genotype "<<longueur<<endl;
	  _genotype.setModif(hf-1,true);
	//cout<<"val "<<s.valeur<<endl;
        	     					  for(h=hf-2;h>=hi;h--)
        	     					  {
								//cout<<s.pred<<" "<<h<<" "<<hi<<endl;
								//cout<<graphe[h-hi].size()<<endl;
        	     					//	  cout<<_genotype.getQuantite(h,r)<<endl;
        	     						  _genotype.setQuantite(h,r,graphe[h-hi][s.pred].quantite);
        	     						  _genotype.setModif(h,true);
        	     						  if(_genotype.getQuantite(h,r)>_genotype.getQuantite(h+1,r)){
        	     							  cout<<"PROBLEME "<<_genotype.getQuantite(h,r)<<" "<<_genotype.getQuantite(h+1,r)<<" "<<h<<" "<<hi<<" hf "<<hf<<" qTot "<<qTot<<" "<<" "<<s.pred<<graphe[h-hi].size()<<endl;

        	     						  }
        	     						  s=graphe[h-hi][s.pred];
        	     					  }
	//cout<<"fin mutation2 verif :"<<endl;

/*for(i=0;i<nbReservoirs;i++){
	for(int h=1;h<nbHeures;h++)
	{
		if(_genotype.getQuantite(h-1,i)>_genotype.getQuantite(h,i))cout<<"mut 2 i : "<<i<<" h "<<h<<" qte "<<_genotype.getQuantite(h,i)<<" qte pred "<<_genotype.getQuantite(h-1,i)<<endl;


	}}*/
    return isModified;
    // END code for mutation of the _genotype object
  }

private:
// START Private data of an eocascadeMutation object
 	vector< vector<double> > V;
 	Systeme* systeme;
 	int longueur;
 	int nbPas;
		   // for example ...
// END   Private data of an eocascadeMutation object
};

#endif

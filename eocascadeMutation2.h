#ifndef eocascadeMutation2_H
#define eocascadeMutation2_H


#include <eoOp.h>
typedef struct Sommet{
	double valeur;
	double quantite;
	int pred;
}Sommet;
template<class GenotypeT>
class eocascadeMutation2: public eoMonOp<GenotypeT>
{
public:
  eocascadeMutation2(vector< vector<double> > _V, Systeme* _systeme,int _longueur,int _nbPas)
  {
	  V=_V;
	  systeme=_systeme;
	  longueur=_longueur;
	  nbPas=_nbPas;
  }

  string className() const { return "eocascadeMutation2"; }

  /**
  à revoir, beaucoup de choix sont instance dépendant...
   */
  bool operator()(GenotypeT & _genotype)
  {
      bool isModified(true);
      vector< vector<Sommet> > graphe;
      int i,j,h,hi,hf,r,rSuc;
      int nbReservoirs=_genotype.getNbReservoirs();
      int nbHeures=_genotype.getNbEtats();

      eoUniformGenerator<int> choixReservoir(1,nbReservoirs-1);
      r=choixReservoir();
      if(r==2)r=5; //revoir si Vmin=Vmax...
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
      //if(r==0) hf=hi+24;
      //if(r==3&&hi+800<nbHeures)hf=hi+800;
      //else 
      hf=hi+longueur;
      //calcul pas:
      double Vmax=R->getVmax();
      double VmaxS=Rsuc->getVmax();
      double pas =(Vmax/3600)/nbPas;
      if(pas==0)pas=0.1;
      double qTot=_genotype.getQuantite(hf-1,r);
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
      if(qmaxC>0&&(qmaxT+qmaxC)/longueur<pas) pas=(qmaxT+qmaxC)/longueur;
      int turbineS=-1;
      Turbine* Ts;
      if(Rsuc->getNbTurbines()>0)
      {
	turbineS=Rsuc->getTurbine(0);
	Ts=systeme->getTurbine(turbineS);
      }

//	cout<<"debut DP"<<endl;

      double qmin=0;
      double qmax=0;
      if(hi>0){
	qmin=_genotype.getQuantite(hi-1,r);
      }
	
      for(h=hi;h<hf-1;h++)
      {

    	  vector<Sommet> vec;
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
    	  qmax=qTot;
    	  double Vmin=R->getVmin(h);
    	  if((Vr-Vmin)/3600<qmax) qmax=(Vr-Vmin)/3600;
    	  if((VmaxS-VrSuc)/3600<qmax)qmax=(VmaxS-VrSuc)/3600;
    	  if(qmin>qmax)//respect qminC impossible
    	  {

    		  qmin=qmax;
    	  }
    	  double qte;
    	  for(qte=qmin;qte<=qmax;qte=qte+pas)
    	  {
    		  //creer sommet :
		Sommet s;
		s.quantite=qte;
		s.pred=-1;
		s.valeur=-100000;
		//if(qte>qTot) cout<<"wtf?"<<endl;
		//recherche pred :
		if(h==hi)
		{
			//qte deversée:
			double qtePred=0;
			if(hi>0)qtePred=_genotype.getQuantite(hi-1,r);
			double qteD= qte-qtePred;
			if(qteD>=0)
    			{
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

    				  if(qteD>=qminC-0.00001)
    				  {
    					  double val=0;
    					  double qteT=0;
    					  if(bT>0)
    					  {
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
					if(qteD>=qminC+0.000001)
    					{
						double qteT=0;
    						if(bT>0)
    						{
							qteT=qteD-qminC;
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
    						  }
					          if(qmaxC>0&&qteD-qteT>qmaxC)
						  {
							val=val-P2*(qteD-qteT-qmaxC);    			
						  }
						  if(qmaxCs>0&&qteDs-qteTs>qmaxCs)
						  {
							val=val-P2*(qteDs-qteTs-qmaxCs);
						  }

    					}
					else val=val-P1*(qteD-qminC);
					if(val+graphe[h-hi-1][j].valeur>s.valeur)//MODIF
    					{
						s.valeur=val+graphe[h-hi-1][j].valeur;
						s.pred=j;

    				        }
    				  }
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
    		if(qteD>=0)
    		{
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
    			if(qteD>=qminC-0.000001)
    			{
				double qteT=0;
    				if(bT>0)
    				{
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
    			        if(qmaxC>0 && qteD-qteT>qmaxC) val=val-100000*( qteD-qteT-qmaxC)*( qteD-qteT-qmaxC);
    				if(val>s.valeur)
    				{
					s.valeur=val;
    					vec.push_back(s);
				}
    					  
    	                }
    			else
    			{
				val=val-P1*(qminC-qteD)*(qminC-qteD);
    				if(val>s.valeur || s.pred<0)
    				{
					 s.valeur=val;
    					 vec.push_back(s);
    			
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
			if(qteD>=0)
    	     		{
				if(qteD>=qminC-0.000001)
    	     			{
					 double qteT=0;
    	     				 if(bT>0)
    	     				 {
					 	qteT=qteD-qminC;
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
    	     					double val=0;
						if(qteTs>qmaxTi)qteTs=qmaxTi;
    	     					if(qteTs<qminTi)qteTs=0;
    	     					val=val+Ts->getBenefice(VrSucPred,qteTs,h);
    	     				}
    	     				if(qmaxC>0&&qteD-qteT>qmaxC)val=val-P2*(qteD-qteT-qmaxC);
    	     				if(qmaxCs>0&&qteDs-qteTs>qmaxCs)val=val-P2*(qteDs-qteTs-qmaxCs);
    	     				if(val+graphe[h-hi-1][j].valeur>s.valeur||s.pred<0)
    	     				{
    	     					s.valeur=val+graphe[h-hi-1][j].valeur;
    	     					s.pred=j;
    	     						      		
					}
				}
    	     			else
    	     			{
					val=val-P1*(qminC-qteD);
    	     				if(val+graphe[h-hi-1][j].valeur>s.valeur||s.pred<0)
    	     				{
    	     					s.valeur=val+graphe[h-hi-1][j].valeur;
    	     					s.pred=j;

    	     				}
    	     			}
    	     		}
		        j++;
    	     	}
		if(s.pred>-1)vec.push_back(s);
		//else cout<<"pb de merde, was passiert?m2 "<<h<<" "<<hi<<" "<<hf<<" "<<_genotype.getQuantite(h,r)<<" "<<_genotype.getQuantite(h-1,r)<<" "<<graphe[h-hi-1][graphe[h-hi-1].size()-1].quantite<<endl;
	}
	//recherche pred avec possibilité de sol non realisable si vec vide:
	graphe.push_back(vec);
      }
      //derniere heure:
	;
      Sommet s;
      s.quantite=qTot;
      s.valeur=-1000000000000000000;
      s.pred=-1;
      h=hf-1;
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
        	     	}
        	}
        }
        j++;
  }
  _genotype.setModif(hf-1,true);
  for(h=hf-2;h>=hi;h--)
  {
	_genotype.setQuantite(h,r,graphe[h-hi][s.pred].quantite);
        _genotype.setModif(h,true);
        if(_genotype.getQuantite(h,r)>_genotype.getQuantite(h+1,r)){
        	cout<<"PROBLEME "<<_genotype.getQuantite(h,r)<<" "<<_genotype.getQuantite(h+1,r)<<" "<<h<<" "<<hi<<" hf "<<hf<<" qTot "<<qTot<<" "<<" "<<s.pred<<graphe[h-hi].size()<<endl;
	}
        s=graphe[h-hi][s.pred];
 }
 return isModified;
    
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

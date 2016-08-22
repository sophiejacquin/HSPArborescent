#ifndef eocascadeMutation3_H
#define eocascadeMutation3_H
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include<ilconcert/ilomodel.h>
#include <iostream>
#include <eoOp.h>

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
template<class GenotypeT>
class eocascadeMutation3: public eoMonOp<GenotypeT>
{
public:

  eocascadeMutation3(vector< vector<double> > _V, Systeme* _systeme,int _nbR,int _nbH)

  {
	  V=_V;
	  systeme=_systeme;
	  nbH=_nbH;
	  nbR=_nbR;
  }

 
  string className() const { return "eocascadeMutation3"; }

  bool operator()(GenotypeT & _genotype)
  {
	cout<<"debut mutation3"<<endl;
      bool isModified(true);
      int i,j,k,l,m,t,h,s;
      IloEnv env;
      IloModel model(env);
      IloInt nbHeures=nbH;
      //choix des réservoirs :
      int nbReservoirs=_genotype.getNbReservoirs();
      int reservoirs[nbReservoirs];
      for(i=0;i<nbReservoirs;i++)
      {
    	  reservoirs[i]=-1;
      }
      eoUniformGenerator<int> choixR(0,nbReservoirs-1);
      int cont=0;
      IloInt nbT=0;
      vector<int> turbines;
      while(cont<nbR)
      {
    	  i=choixR();
    	  if(reservoirs[i]==-1)
    	  {

    		  reservoirs[i]=cont;
    		  cont++;
    		  if(systeme->getReservoir(i)->getNbTurbines()>0)
    		  {
    			int turbine=systeme->getReservoir(i)->getTurbine(0);
    			turbines.push_back(turbine);
			nbT++;
    		  }
    	  }
      }
      //choix heure début :
      int nbHe=_genotype.getNbEtats();
      eoUniformGenerator<int> choixHdeb(1,nbHe-nbH);
      int hDeb=choixHdeb();
      cout<<hDeb<<endl;
      //calcul nbB:
      IloInt nbB=0;
      IloInt nbS=0;
      vector<int> suc;
      for(i=0;i<nbReservoirs;i++)
      {
    	  int dev=systeme->getReservoir(i)->getDeversement();
    	  if(reservoirs[i]>-1 &&dev>-1&& reservoirs[dev]==-1)
    	  {
    		  //recherche si dev deja dans la liste:
    		  bool trouve=false;
    		  for(j=0;j<nbS;j++)
    		  {
    			if(suc[j]==dev)trouve=true;
    		  }
    		  if(trouve==false){
    			  suc.push_back(dev);
    			  if(systeme->getReservoir(dev)->getNbTurbines()>0)
    			  {
    				  int turbine=systeme->getReservoir(dev)->getTurbine(0);
    				  turbines.push_back(turbine);
    				  nbB++;
    			  }
    			  nbS++;
    		  }

    	  }
      }
      IloInt sum=nbT+nbB;
            IloIntArray nbIntervals(env,sum);
            IloIntArray nbPieces(env,sum);
            for(i=0;i<sum;i++)
            {
		//cout<<"turbine "<<turbines[i]<<endl;
            	IloInt inter=systeme->getTurbine(turbines[i])->getNbInt();
            	nbIntervals[i]=inter;
            }

            for(i=0;i<sum;i++)
            {

            	IloInt piece=systeme->getTurbine(turbines[i])->getNbPieces();
            	nbPieces[i]=piece;

            }
	//cout<<"debut déclarations"<<endl;
	IloInt nbRe=nbR;
  	IloNumVarArray4 x(env,nbT);
      	IloNumVarArray4 d(env,nbT);
      	IloNumVarArray2 r(env,nbRe);
      	IloNumVarArray2 R1(env,nbRe);
      	IloNumVarArray2 R2(env,nbRe);
      	IloNumVarArray2 tot(env,nbRe);
      	IloNumVarArray3 b(env,nbB);
      //Declaration des x :
      for(i=0;i<nbT;i++)
      {
    	  x[i]=IloNumVarArray3(env,nbIntervals[i]);
    	  for(j=0;j<nbIntervals[i];j++)
    	  {
    		  x[i][j]=IloNumVarArray2(env,nbHeures);
    		  for(k=0;k<nbH;k++)
    		  {
    			  x[i][j][k]=IloNumVarArray(env,nbPieces[i],0,IloInfinity);
    		  }
    	  }
      }
      //pour d:
      	for(i=0;i<nbT;i++)
          	{

              	d[i]=IloNumVarArray3(env,nbIntervals[i]);
              	for(j=0;j<nbIntervals[i];j++)
              	{
              		d[i][j]=IloNumVarArray2(env,nbHeures);
              		for(k=0;k<nbHeures;k++)
              		{
              			d[i][j][k]=IloNumVarArray(env,nbPieces[i],0,1,ILOBOOL);
              		}
              	}
          	}
      	//pour r
      	for(i=0;i<nbR;i++)
      	    	{

      	            		r[i]=IloNumVarArray(env,nbHeures,0,IloInfinity);
      	    	}
    	//pour tot
    	for(i=0;i<nbR;i++)
        {

                	tot[i]=IloNumVarArray(env,nbHeures,0,IloInfinity);

        }
	//cout<<"tot ini"<<endl;
    	//pour b
    	   for(i=nbT;i<sum;i++)
    	   {
				//cout<<sum<<endl;
			//	cout<<turbines[i]<<endl;
			 //  cout<<nbIntervals[i]<<endl;
			
    		   	   b[i-nbT]=IloNumVarArray2(env,nbIntervals[i]);
			//cout<<"là"<<endl;
    		   	   for(j=0;j<nbIntervals[i];j++)
    		   	   {
    		   		   b[i-nbT][j]=IloNumVarArray(env,nbHeures,0,1,ILOBOOL);

    		   	   }
			//cout<<"ici"<<endl;
    	   }
	//	cout<<"bini"<<endl;
    	   //pour R1 :
    	   for(i=0;i<nbR;i++)
    	   {
    		   R1[i]=IloNumVarArray(env,nbHeures,0,1,ILOBOOL);
    	   }
	//cout<<"R1 ini"<<endl;
    	   //pour R2 :
    	      	   for(i=0;i<nbR;i++)
    	      	   {
    	      		   R2[i]=IloNumVarArray(env,nbHeures,0,1,ILOBOOL);
    	      	   }
	//cout<<"fin des déclarations"<<endl;
    	   //fonction objective :
    	   IloExpr objective(env);
    	   for(h=0;h<nbH;h++)
    	   {

    		   for(s=0;s<nbR;s++)
    		   {
    			   objective+= -10000*R1[s][h];
			  objective+= -10000*R2[s][h];
    		   }
		//	cout<<"fonction obj Part one"<<endl;
    		   for(i=0;i<nbT;i++)
    		   {
    			   IloNum prixSpot= systeme->getTurbine(turbines[i])->getPrix(h+hDeb) ;
    			   for(j=0;j<nbIntervals[i];j++)
    			   {
    				   for(k=0;k<nbPieces[i];k++)
    				   {
    					   IloNum w = systeme->getTurbine(turbines[i])->getProduction(k,j)*prixSpot ;
    					   objective += w * x[i][j][h][k];

    				   }
    			   }
    		   }
		//cout<<"fonction obj Part two"<<endl;
    		   for(i=nbT;i<nbT+nbB;i++)
    		   {
    			   IloNum prixSpot= systeme->getTurbine(turbines[i])->getPrix(h+hDeb) ;
    			   for(j=0;j<nbIntervals[i];j++)
    			   {

    					   //calcul qteT:
    					   //num reservoir associé à la turbine:
    					   int res=systeme->getTurbine(turbines[i])->getReservoirParent();
    					   double Qtot;
    					   Qtot=_genotype.getQuantite(res,h+hDeb)-_genotype.getQuantite(res,h+hDeb-1);
    					 //  else Qtot=_genotype.getQuantite(res,0);
    					   //calcul qminC
    					   double qminC=systeme->getReservoir(res)->getQmin();
    					   if(qminC<0)qminC=0;
    					   //calcul qmaxC
    					   double qmaxC=systeme->getReservoir(res)->getQmax();
    					   //calcul qmin
    					   double qminT=systeme->getTurbine(turbines[i])->getQminInt(j);
    					   //calcul qmax:
    					   double qmaxT=systeme->getTurbine(turbines[i])->getQMax(j);
    					   //calcul qte :
    					   double qte=Qtot-qminC;
    					   if(qte<0)
    					   {
    						   objective+=-10000000*b[i-nbT][j][h];//*qte*qte;
    					   }
    					   else
    					   {
    						   if(qte<qminT)qte=0;
    						   if(qte>qmaxT)qte=qmaxT;
    						   if(Qtot-qte>qmaxC)
    						   {
    							   objective +=-10000000*b[i-nbT][j][h];//*(Qtot-qte-qmaxC)*(Qtot-qte-qmaxC)*b[i-nbT][j][h];
    						   }
    						   else  objective +=systeme->getTurbine(turbines[i])->getBeneficeInt(j,qte,h+hDeb) *b[i-nbT][j][h];
    					   }
    			   }
    		   }
		//cout<<"fonction obj Part 3"<<endl;

    	   }
    	   model.add(IloMaximize(env,objective));
	//cout<<"fonction objective créée"<<endl;
    	   //LES CONTRAINTES :
    	   	//Contraintes 1: débits maximaux :
    	   	for(i=0;i<nbT;i++)
    	   	{
    	   		for(t=0;t<nbHeures;t++)
    	   		{
    	   			for(j=0;j<nbIntervals[i];j++)
    	   			{
    	   				IloExpr contrainte1(env);
    	   				IloNum qmax=systeme->getTurbine(turbines[i])->getQMax(j);

    	   				//cout<<i<<" "<<j<<" "<<qmax<<endl;
    	   				contrainte1 += -qmax*d[i][j][t][0];
    	   				for(k=0;k<nbPieces[i];k++)
    	   				{
    	   					contrainte1+=x[i][j][t][k];
    	   				}
    	   				model.add(contrainte1<=0);
    	   				contrainte1.end();
    	   			}
    	   		}
    	   	}
	//	cout<<"contraintes1 créées"<<endl;
    		//Contraintes 2: débits minimaux :
    	   	for(i=0;i<nbT;i++)
    		{
    			IloNum prodMin=systeme->getTurbine(turbines[i])->getProdMin();
    		//	cout<<"prod Min "<<prodMin<<endl;
    			if (prodMin>0) {

    			for(t=0;t<nbHeures;t++)
    			{
    				for(j=0;j<nbIntervals[i];j++)
    				{
    					IloExpr contrainte2(env);

    					contrainte2 += -prodMin*d[i][j][t][0];
    					for(k=0;k<nbPieces[i];k++)
    					{
    						IloNum prod=systeme->getTurbine(turbines[i])->getProduction(k,j);
    						contrainte2+= prod * x[i][j][t][k];

    					}
    					model.add(contrainte2>=0);
    					contrainte2.end();
    				}
    			}
    			}
    		}
	//	cout<<"contraintes 2 créées"<<endl;
    	    	//Contraintes 3: définitions des tot0:
    	    		for(i=0;i<nbReservoirs;i++)
    	    		{
    	    			if(reservoirs[i]>-1)
    	    			{
    	    			IloExpr contraintes3(env);
    	    			contraintes3 += -tot[reservoirs[i]][0]+ 3600*r[reservoirs[i]][0];
    	    			for(t=0;t<systeme->getReservoir(i)->getNbTurbines();t++)
    	    			{
    	    				int turbine=systeme->getReservoir(i)->getTurbine(t);
    	    				//chercher turbine:
    	    					bool trouve=false;
    	    					l=0;
    	    					while(trouve==false)
    	    					{
    	    						if(turbines[l]==turbine)
    	    						{
    	    							trouve=true;
    	    						}
    	    						else l++;
    	    					}

    	    				for(j=0;j<nbIntervals[l];j++)
    	    				{
    	    					for(k=0;k<nbPieces[l];k++)
    	    					{
    	    						contraintes3+=3600* x[l][j][0][k];

    	    					}
    	    				}

    	    			}
    	    			model.add(contraintes3==0);
    	    			contraintes3.end();
    	    			}
    	    		}
	//		cout<<"contraintes 3 créées"<<endl;
    	    		//Contraintes 4: définitions des tots:
    	    			for(i=0;i<nbReservoirs;i++)
    	    			{
    	    				if(reservoirs[i]>-1)
    	    				{
    	    					for(h=1;h<nbHeures;h++)
    	    					{
    	    					IloExpr contraintes4(env);
    	    					contraintes4+=tot[reservoirs[i]][h-1]-tot[reservoirs[i]][h]+ 3600*r[reservoirs[i]][h];
    	    					for(t=0;t<systeme->getReservoir(i)->getNbTurbines();t++)
    	    					{
    	    						IloInt turbine=systeme->getReservoir(i)->getTurbine(t);
    	    						//chercher turbine:
    	    						bool trouve=false;
    	    						l=0;
    	    						while(trouve==false)
    	    						{
    	    							if(turbines[l]==turbine)
    	    							{
    	    								trouve=true;
    	    							}
    	    							else l++;
    	    						}
    	    						for(j=0;j<nbIntervals[l];j++)
    	    						{
    	    							for(k=0;k<nbPieces[l];k++)
    	    							{
    	    								contraintes4+=3600* x[l][j][h][k];
    	    								//if(turbine==5) cout<<"ici"<<endl;
    	    							}
    	    						}

    	    					}
    	    					model.add(contraintes4==0);
    	    					}
    	    				}
    	    			}
//	cout<<"contraintes 4 créées"<<endl;
    	    			//Contraintes 5; bmax interval R:
    	    				for(i=0;i<nbReservoirs;i++)
    	    				{
    	    					if(reservoirs[i]>-1)
    	    					{

    	    						for(h=1;h<nbH;h++)
    	    						{
    	    							IloNum Vh=V[hDeb+h-1][i]-_genotype.getQuantite(hDeb-1,i)*3600;

    	    						    IloExpr contraintes5(env);
    	    						    for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
    	    						    {
    	    						    	int parent=systeme->getReservoir(i)->getParents()[j];
    	    						    	if(reservoirs[parent]>-1){
    	    						    		contraintes5+=tot[reservoirs[parent]][h-1];
    	    						    		Vh=Vh+_genotype.getQuantite(hDeb-1,parent)*3600;
    	    						    	}
    	    						    	else Vh=Vh+3600*(_genotype.getQuantite(hDeb+h-1,parent));
    	    						    }
    	    						    IloNum borne= systeme->getReservoir(i)->getVmax()-Vh;
    	    						    contraintes5+=-tot[reservoirs[i]][h-1];
    	    						    model.add(contraintes5<=borne);
    	    						}
    	    					}

    	    				}
		//				cout<<"contraintes 5 créées"<<endl;
    	    				//Contraintes 6; bmax interval S:
    	    				for(i=0;i<nbS;i++)
    	    				{
    	    					//calcul Vh ini:

    	    					for(h=1;h<nbH;h++)
    	    					{
    	    						IloNum Vh=V[hDeb+h-1][suc[i]]-_genotype.getQuantite(hDeb-1+h,suc[i])*3600;

    	    						IloExpr contraintes5(env);
    	    						for(j=0;j<systeme->getReservoir(suc[i])->getNbParents();j++)
    	    						{
    	    							int parent=systeme->getReservoir(suc[i])->getParents()[j];
    	    							if(reservoirs[parent]>-1)
    	    							{
    	    								contraintes5+=tot[reservoirs[parent]][h-1];
    	    								Vh=Vh+_genotype.getQuantite(hDeb-1,parent)*3600;
    	    							}
    	    							else Vh=Vh+(_genotype.getQuantite(hDeb+h-1,parent))*3600;
    	    						}
    	    						IloNum borne= systeme->getReservoir(suc[i])->getVmax()-Vh;
    	    					//	if(borne<0)cout<<"borne "<<borne<<" Vmax "<<systeme->getReservoir(suc[i])->getVmax()<<" Vh "
    	    						model.add(contraintes5<=borne);
    	    					}
    	    				}
		//					cout<<"contraintes 6 créées"<<endl;
    	    				//Contraintes 7: bmin intervalle sur R :
    	    				for(i=0;i<nbReservoirs;i++)
    	    				{
    	    					if(reservoirs[i]>-1)
    	    					{


    	    						for(h=1;h<nbHeures;h++)
    	    						{
    	    							IloNum Vh=V[hDeb+h-1][i]-_genotype.getQuantite(hDeb-1,i)*3600;
    	    							for(l=0;l<systeme->getReservoir(i)->getNbParents();l++)
    	    							{
    	    								int parent=systeme->getReservoir(i)->getParents()[l];
    	    								if(reservoirs[parent]==-1)Vh=Vh+ _genotype.getQuantite(h+hDeb-1,parent)*3600;
    	    								else Vh=Vh+ _genotype.getQuantite(hDeb-1,parent)*3600;
    	    							}
    	    							IloNum borne= -systeme->getReservoir(i)->getVmin(h +hDeb-1)+Vh;
    	    							for(t=0;t<systeme->getReservoir(i)->getNbTurbines();t++)
    	    							{
    	    								int turb=systeme->getReservoir(i)->getTurbine(t);
    	    								//chercher turbine:
    	    								bool trouve=false;
    	    								l=0;
    	    								while(trouve==false)
    	    								{
    	    									if(turbines[l]==turb)
    	    									{
    	    										trouve=true;
    	    									}
    	    									else l++;
    	    								}
    	    								IloInt turbine=l;
    	    								for(j=0;j<nbIntervals[turbine];j++)
    	    								{
    	    									//for(k=0;k<nbPieces[turbine];k++)
    	    								//	{
    	    										IloNum coef=systeme->getTurbine(turb)->getBinfIntReservoir(j)-systeme->getReservoir(i)->getVmin(h +hDeb-1);
    	    										IloExpr contraintes6(env);
    	    										contraintes6+=coef*d[turbine][j][h][0];
    	    										contraintes6+=tot[reservoirs[i]][h-1];//bizarre que
    	    										for(l=0;l<systeme->getReservoir(i)->getNbParents();l++)
    	    										{
    	    											int parent=systeme->getReservoir(i)->getParents()[l];
    	    											if(reservoirs[parent]>-1)contraintes6+=-tot[reservoirs[parent]][h-1];

    	    										}
    	    										model.add(contraintes6<=borne);
    	    									//}
    	    								}
    	    							}
    	    							if(systeme->getReservoir(i)->getNbTurbines()==0)
    	    							{
    	    								IloExpr contraintes6(env);
    	    								contraintes6+=tot[reservoirs[i]][h-1];
    	    								for(j=0;j<systeme->getReservoir(i)->getNbParents();j++)
    	    								{
    	    									int parent=systeme->getReservoir(i)->getParents()[j];
    	    									if(reservoirs[parent]>-1)contraintes6+=-tot[reservoirs[parent]][h-1];
    	    									//else borne=borne+( _genotype.getQuantite(h+hDeb-1,parent)-_genotype.getQuantite(hDeb-1,parent))*3600;
    	    								}
    	    								model.add(contraintes6<=borne);
    	    							}
    	    						}
    	    					}
    	    				}
//	cout<<"contraintes 7 créées"<<endl;
    	    				//Contraintes 8: bmin intervalle sur S :
    	    				for(i=0;i<nbS;i++)
    	    				{

    	    					for(h=1;h<nbHeures;h++)
    	    					{
    	    						IloNum Vh=V[hDeb+h-1][suc[i]]-_genotype.getQuantite(hDeb-1+h,suc[i])*3600;
    	    						for(l=0;l<systeme->getReservoir(suc[i])->getNbParents();l++)
    	    						{
    	    							int parent=systeme->getReservoir(suc[i])->getParents()[l];
    	    							if(reservoirs[parent]==-1)Vh=Vh+( _genotype.getQuantite(h+hDeb-1,parent))*3600;
    	    							else Vh=Vh+( _genotype.getQuantite(hDeb-1,parent))*3600;
    	    						}
    	    						IloNum borne= -systeme->getReservoir(suc[i])->getVmin(h +hDeb-1)+Vh;
    	    						for(t=0;t<systeme->getReservoir(suc[i])->getNbTurbines();t++)
    	    						{
    	    							int turb=systeme->getReservoir(suc[i])->getTurbine(t);
    	    							//chercher turbine:
    	    							bool trouve=false;
    	    							l=0;
    	    							while(trouve==false)
    	    							{
    	    								if(turbines[l]==turb)
    	    								{
    	    									trouve=true;
    	    								}
    	    								else l++;
    	    							}
    	    							IloInt turbine=l;
    	    							for(j=0;j<nbIntervals[turbine];j++)
    	    							{
    	    								IloNum coef=systeme->getTurbine(turb)->getBinfIntReservoir(j)-systeme->getReservoir(suc[i])->getVmin(h +hDeb-1);
    	    								IloExpr contraintes6(env);
    	    								contraintes6+=coef*b[turbine-nbT][j][h];
    	    								for(l=0;l<systeme->getReservoir(suc[i])->getNbParents();l++)
    	    								{
    	    									int parent=systeme->getReservoir(suc[i])->getParents()[l];
    	    									if(reservoirs[parent]>-1)contraintes6+= -tot[reservoirs[parent]][h-1];
    	    								}
    	    								//borne=borne-_genotype.getQuantite(h+hDeb-1,suc[i])*3600;

    	    								model.add(contraintes6<=borne);
    	    							}
    	    						}
    	    						if(systeme->getReservoir(suc[i])->getNbTurbines()==0)
    	    						{
    	    							IloExpr contraintes6(env);
    	    							//borne=borne-_genotype.getQuantite(h+hDeb-1,suc[i])*3600;
    	    							for(j=0;j<systeme->getReservoir(suc[i])->getNbParents();j++)
    	    							{
    	    								int parent=systeme->getReservoir(suc[i])->getParents()[j];
    	    								if(reservoirs[parent]>-1)contraintes6+=-tot[reservoirs[parent]][h-1];
    	    								//else borne=borne+ (_genotype.getQuantite(h+hDeb-1,parent)-_genotype.getQuantite(hDeb-1,parent))*3600;
    	    							}
    	    							model.add(contraintes6<=borne);
    	    						}
    	    					}
    	    				}
//	cout<<"contraintes 8 créées"<<endl;
    	    				//Contraintes 9: on est positivement actif dans au plus un interval :
    	    					for(i=0;i<nbT;i++)
    	    					{
    	    						for(t=0;t<nbHeures;t++)
    	    						{
    	    							IloExpr contraintes8(env);
    	    							for(j=0;j<nbIntervals[i];j++)
    	    							{
    	    								contraintes8+=d[i][j][t][0];
    	    							}
    	    							model.add(contraintes8<=1);
    	    						}
    	    					}
//	cout<<"contraintes 9 créées"<<endl;
    	    					//Contraintes 10: on est dans un seul interval:
    	    						for(i=0;i<nbB;i++)
    	    						{
    	    							for(t=0;t<nbHeures;t++)
    	    							{
    	    								IloExpr contraintes8(env);
    	    								for(j=0;j<nbIntervals[i+nbT];j++)
    	    								{
    	    									contraintes8+=b[i][j][t];
    	    								}
    	    								model.add(contraintes8==1);
    	    							}
    	    						}
//	cout<<"contraintes 10 créées"<<endl;
    	    						//Contraintes 11 : Remplissage des pieces ordonné :

    	    							 for(i=0;i<nbT;i++)
    	    							 {
    	    								 for(j=0;j<nbIntervals[i];j++)
    	    								 {
    	    									 for(t=0;t<nbHeures;t++)
    	    									 {
    	    										 for(k=1;k<nbPieces[i];k++)
    	    										 {
    	    											 IloNum coef=1/systeme->getTurbine(turbines[i])->getBmaxMorceau(k-1);
    	    						             			model.add(d[i][j][t][k]<=coef*x[i][j][t][k-1]);
    	    										 }
    	    									 }
    	    								 }
    	    							 }
//	cout<<"contraintes 11 créées"<<endl;
    	    							 //Contraines 12: bornes sup pieces :
    	    							for(i=0;i<nbT;i++)
    	    							 {
    	    								 for(j=0;j<nbIntervals[i];j++)
    	    								 {
    	    									 for(t=0;t<nbHeures;t++)
    	    									 {
    	    										 for(k=0;k<nbPieces[i];k++)
    	    										 {
    	    											 IloNum coef=systeme->getTurbine(turbines[i])->getBmaxMorceau(k);

    	    											 model.add(coef*d[i][j][t][k]-x[i][j][t][k]>=0);
    	    							 			}
    	    							 		}
    	    								 }
    	    							 }
//	cout<<"contraintes 12 créées"<<endl;

    	    							 	//Contraintes 13 : bmin reserves:
    	    							for(i=0;i<nbReservoirs;i++)
    	    							 	{
    	    							 		if(reservoirs[i]>-1)
    	    							 		{
    	    							 			IloNum qmin=systeme->getReservoir(i)->getQmin();
    	    							 			if(qmin>0)
    	    							 			{
    	    							 				for(h=0;h<nbHeures;h++)
    	    							 				{
    	    							 					model.add(r[reservoirs[i]][h]+qmin*R1[reservoirs[i]][h]>=qmin);
    	    							 					//model.add(r[reservoirs[i]][h]+qmin*R1[reservoirs[i]][h]>=qmin);
    	    							 				//if(i==5) cout<<"qmin r"<<qmin<<endl;
    	    							 				}
    	    							 			}
    	    							 		}
    	    							 	}
//	cout<<"contraintes 13 créées"<<endl;
    	    							 	//Contraintes 14: bmax reserves
    	    							 		for(i=0;i<nbReservoirs;i++)
    	    							 		{
    	    							 			if(reservoirs[i]>-1)
    	    							 			{
    	    							 				IloNum qmax=systeme->getReservoir(i)->getQmax();
    	    							 				if(qmax>0)
    	    							 				{
    	    							 					for(h=0;h<nbHeures;h++)
    	    							 					{
    	    							 						model.add(r[reservoirs[i]][h]-10000*R2[reservoirs[i]][h]<=qmax);
    	    							 					}
    	    							 				}
    	    							 			}
    	    							 		}
//	cout<<"contraintes 14 créées"<<endl;
    	    								//Contraintes 15: liquidation

    	    							 			for(i=0;i<nbReservoirs;i++)
    	    							 			{
    	    							 				if(reservoirs[i]>-1)
    	    							 				{
    	    							 					IloNum borne=(_genotype.getQuantite(hDeb+nbH-1,i)-_genotype.getQuantite(hDeb-1,i))*3600;
    	    							 		model.add(tot[reservoirs[i]][nbH-1]==borne);
    	    							 				}
    	    							 			}
//	cout<<"ttes les contraintes créées"<<endl;
    	    							 			//resolution et restitution :
    	    							 			try{
    	    							 				 IloCplex cplex(model);
    	    							 				//cplex.exportModel ("lpex1.lp");
    	    							 				cplex.setOut(env.getNullStream());
    	    							 				 cplex.solve();
    	    							 				//cout << "Solution status = " << cplex.getStatus() << endl;
    	    							 			   // cout<< cplex.getObjValue() << endl;
    	    							 				//Résultats:
    	    							 				for(h=0;h<nbHeures-1;h++)
    	    							 				{
    	    							 					for(i=0;i<nbReservoirs;i++)
    	    							 					{

    	    							 						if(reservoirs[i]>-1){
    	    							 							double qte=cplex.getValue(tot[reservoirs[i]][h])/3600+_genotype.getQuantite(hDeb-1,i);
    	    							 							if(h+hDeb-1>-1&&qte<_genotype.getQuantite(h+hDeb-1,i)+0.00001)qte=_genotype.getQuantite(h+hDeb-1,i);
    	    							 							_genotype.setQuantite(h+hDeb,i,qte);

    	    							 						}

    	    							 					}
    	    							 					_genotype.setModif(h+hDeb,true);

    	    							 				}

    	    							 			  cplex.end();
    	    							 		    							 			  env.end();
    	    							 			

    	    							 			     }
    	    							 			   catch (IloException& e) {
    	    							 			      cerr << "Concert exception caught: " << e << endl;
    	    							 			   }
    	    							 			  

    return isModified;
    // END code for mutation of the _genotype object
  }

private:
// START Private data of an eocascadeMutation object
 	vector< vector<double> > V;
 	Systeme* systeme;
 	int nbH;
 	int nbR;

		   // for example ...
// END   Private data of an eocascadeMutation object
};

#endif

#ifndef _Systeme_h

#define _Systeme_h
#include <iostream>
#include <vector>
#include <fstream>
#include "Reservoir.h"
#include "Turbine.h"
class Systeme
{
	private :
	int nbReservoirs;
	int nbTurbines;
	vector<Reservoir*> reservoirs;
	vector<Turbine> turbines;
	vector<double*> prix;
	int nbPrix;
	int nbH;
	public :
	//Constructeur:
	Systeme()
	{
		nbReservoirs=0;
		nbTurbines=0;
		nbPrix=0;
	}
	Systeme(char* titre1,char* titre2,char*titre3 )
	{	
		//Déclarations
		int i,j,k;
		//Ouverture fichier
		cout<<"là"<<endl;
		ifstream fichier(titre1, ios::in);
		if(fichier){
			cout<<"si fichier"<<endl;
		//création des prix
		
			fichier>>nbPrix;//spot OAcroissants
			fichier>>nbH;
			cout<<"nbH="<<nbH<<endl;
			for(i=0;i<nbPrix;i++)
			{
				double* prixN=new double[nbH];
				for(j=0;j<nbH;j++)
				{
					fichier>>prixN[j];
				}
				prix.push_back(prixN);	
			}
			fichier.close();
		}
		cout<<"prix créés"<<endl;
		ifstream fichier2(titre2,ios::in);	
		//création des turbines
		if(fichier2){
			fichier2 >>nbTurbines;
			
			for(i=0;i<nbTurbines;i++)
			{
				int nbInt;
				fichier2>>nbInt;
				double* Inter=new double[nbInt];
				for(j=0;j<nbInt;j++)
				{
					fichier2>>Inter[j];
				}
				int nbPieces;
				fichier2 >>nbPieces;
				double* pieces=new double[nbPieces];
				for(j=0;j<nbPieces;j++)
				{
					fichier2>>pieces[j];
					
				}
				double* qmax=new double[nbInt];
				for(j=0;j<nbInt;j++)
				{
					fichier2>>qmax[j];
					qmax[j]=qmax[j];
				}
				double** production;
				production=new double*[nbPieces];
				for(j=0;j<nbPieces;j++)
				{
					production[j]=new double[nbInt];
					for(k=0;k<nbInt;k++)
					{
						fichier2>>production[j][k];
						production[j][k]=production[j][k];
					}
				}
				int distance,pr,reservoirParent;
				double prodMin;
				fichier2>>distance;
				fichier2>>prodMin;
				fichier2>>pr;
				fichier2>>reservoirParent;//Ajouter dans fichier données
				Turbine t(nbInt,nbPieces,prodMin, qmax, production,distance,prix[pr],Inter,pieces,reservoirParent);
				turbines.push_back(t);
				//liberation memoire
				for(j=0;j<nbPieces;j++)
				{
				
				
				
					delete [] production[j];
				
				}
				delete[] production;
				delete[] qmax;
				delete[] pieces;
			}
			fichier2.close();
		}
		cout<<"turbines créées"<<endl;
		//création des réservoirs
		ifstream fichier3(titre3,ios::in);
		if(fichier3)
		{
			fichier3>>nbReservoirs;
			fichier3>>nbH;
			cout<<nbH<<endl;
			//Vinit Vmax Vmin
			for(i=0;i<nbReservoirs;i++)
			{
				
				double Vinit;
				fichier3>>Vinit;
				double Vmax;
				fichier3>>Vmax;
				int nbIntVmin;
				fichier3>>nbIntVmin;
				double* Vmin =new double[nbIntVmin];
				int* intVmin=new int[nbIntVmin];
				for(j=0;j<nbIntVmin;j++)
				{
					fichier3>>Vmin[j];
					fichier3>>intVmin[j];//caractérisé par bsup
				}
				double qmin,qmax;
				fichier3>>qmin;
				fichier3>>qmax;
				int distance, deversement,nbParents;
				fichier3 >> distance; fichier3>>deversement; fichier3>>nbParents;
				
				int* parents=new int[nbParents];
				for(j=0;j<nbParents;j++)
				{
					fichier3>>parents[j];
					parents[j]=parents[j]-1;
				}
				
				int nbT;
				fichier3>>nbT;
				int* listeT=new int[nbT];
				for(j=0;j<nbT;j++)
				{
					fichier3>>listeT[j];
					listeT[j]=listeT[j]-1;
				}
				
				double* apports=new double[nbH];
				for(j=0;j<nbH;j++)
				{
					fichier3>> apports[j];
					//apports[j]=apports[j]*3600;
					//cout<<j<<endl;
				}
				//cout<<"là"<<endl;
				//CRÉATION DES Reservoirs :
				
				Reservoir *r = new Reservoir(Vinit,Vmax,nbIntVmin,intVmin, Vmin,apports,deversement-1, distance,nbT, listeT,qmin,qmax,nbParents, parents);
				//cout<<"dgf"<<endl;
				reservoirs.push_back(r);
				cout<<"reservoirs créé"<<endl;
				delete[] parents;
				delete[] listeT;
				delete[] Vmin;
				//cout<<"deletes faits"<<endl;
				//cout<<i<<endl;
			}
			
			fichier3.close();
		}
		//cout<<"réservoirs créés"<<endl;
		cout<<nbH<<endl;
	}
	void afficher()const
	{
		for(unsigned int i=0; i< nbReservoirs; i++)
		{
			cout<< "successeurs du réservoir "<<i<<":";
				cout<<reservoirs[i]->getDeversement()<<endl;
		
		}
		cout << " turbines:"<<endl;
		for(unsigned int i=0; i< nbTurbines; i++)
		{
			cout<<"Intervalles :"<<endl;
			cout<< "premier prix "<<i<<":";
				cout<<turbines[i].getPrix(0)<<endl;
		
		}
		
	}
	//Acsesseurs:
	void adReservoir(Reservoir* R)
	{
		reservoirs.push_back(R);
		nbReservoirs++;
	}
	void adTurbine(Turbine & T)
	{
		turbines.push_back(T);
		nbTurbines++;
	}
	double* getPrix(int cat)
	{
		return prix[cat];
	}
	int getNbPrix()
	{
		return nbPrix;
	}
	void setNbHeures(int n)
	{
		nbH=n;
	}
	int getNbHeures()
	{
		return nbH;
	}
	int getNbReservoirs()const
	{
		return nbReservoirs;
	}
	int getNbTurbines()const
	{
		return nbTurbines;
	}
	 Reservoir* getReservoir(int i)
	{
		return (reservoirs[i]);
	}
	Turbine* getTurbine(int i)
	{
		return &(turbines[i]);
	}
	void adPrix(double* _prix)
	{
		prix.push_back(_prix);
		nbPrix++;
	}
	bool operator ==(Systeme* s2)
	{
		if(s2->getNbHeures()== nbH)
		{
		}
		else
		{
			cout<<" different nombre d'heures "<<nbH <<" v s "<<s2->getNbHeures()<<endl;
			return false;
		}
		if(s2->getNbTurbines()==nbTurbines)
		{
			for(int i=0;i<nbTurbines;i++)
			{
				if (not(turbines[i]==s2->getTurbine(i)))
				{
					cout<<"les turbines numeros "<<i<<" ne sont pas identiques"<<endl;
					return false;
				}
			}
		}
		else
		{
			cout<<" different nombre de turbines "<<nbTurbines<<" v s "<<s2->getNbTurbines()<<endl;
			return false;
		}
		
		if(s2->getNbReservoirs()== nbReservoirs)
		{
			for(int i=0;i<nbReservoirs;i++)
			{
				if (not(*reservoirs[i]==s2->getReservoir(i)))
				{
					cout<<"les reservoirs numéro "<<i<<" ne sont pas identiques"<<endl;
					return false;
				}
			}
		}
		else
		{
			cout<<" different nombre de Reservoirs "<<nbReservoirs<<" v s "<<s2->getNbReservoirs()<<endl;
			return false;
		}
		return true;
	} 
	
};
#endif

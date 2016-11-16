#ifndef _Reservoir_h

#define _Reservoir_h
#include <iostream>
#include <stdlib.h> 
using namespace std;

class Reservoir
{
	private://Atributs
    double Vmax;
	int nbIntVmin;
	int* intVmin;
	double* Vmin;
	double Vinit;
	double* apports;
	double apportAnnuel;
	static int compteur;
	int numero;
	int nbParents;
	int* parents;
	int reservoirDeDeversement;
	int* listeT;
	int nbT;
	int distance;
	double qmin;
	double qmax;
	public:// Méthodes
	//constructeur
	Reservoir(double _Vinit,double _Vmax,int _nbIntVmin,int* _intVmin, double* _Vmin,double* _apports,int _reservoirDeDeversement, int _distance,int _nbT, int* _listeT,double _qmin,double _qmax,int _nbParents, int* _parents)
	{
		Vinit=_Vinit;
		Vmax=_Vmax;
		nbIntVmin=_nbIntVmin;
		Vmin=(double*)malloc(nbIntVmin*sizeof(double));
		for(int i=0;i<nbIntVmin;i++)
		{
			Vmin[i]=_Vmin[i];
		}
		intVmin=(int*)malloc(nbIntVmin*sizeof(int));
	
		for(int i=0;i<nbIntVmin;i++)
		{
			intVmin[i]=_intVmin[i];
		}
		apports=new double[8760];
		 for(int i=0;i<8760;i++)
        	{
            		apports[i]=_apports[i];
        	}
		apportAnnuel=0;
		for(int j=0; j< 8760;j++)
		{
           	 apportAnnuel= apportAnnuel+apports[j];
       		}
		reservoirDeDeversement=_reservoirDeDeversement;
		distance=_distance;
		
		nbT=_nbT;
		listeT=new int[nbT];
	
		for(int i=0;i<nbT;i++)
		{
			listeT[i]=_listeT[i];
		}
		
		 numero= compteur;
		qmin= _qmin;
		qmax=_qmax;
		
		nbParents=_nbParents;
		parents=new int[nbParents];
		
		for(int i=0;i<nbParents;i++)
		{
			parents[i]=_parents[i];
		}
		
       		 ++compteur;	
	
	}
	Reservoir(double _Vinit,double _Vmax,int _nbIntVmin,int* _intVmin, double* _Vmin,int _reservoirDeDeversement, int _distance,int _nbT, int* _listeT,double _qmin,double _qmax,int _nbParents, int* _parents)
	{
		cout<<"const R"<<endl;
		Vinit=_Vinit;
		Vmax=_Vmax;
		nbIntVmin=_nbIntVmin;
		Vmin=(double*)malloc(nbIntVmin*sizeof(double));
		for(int i=0;i<nbIntVmin;i++)
		{
			Vmin[i]=_Vmin[i];
		}
		intVmin=(int*)malloc(nbIntVmin*sizeof(int));
	
		for(int i=0;i<nbIntVmin;i++)
		{
			intVmin[i]=_intVmin[i];
		}
		
		reservoirDeDeversement=_reservoirDeDeversement;
		distance=_distance;
		
		nbT=_nbT;
		listeT=new int[nbT];
	
		for(int i=0;i<nbT;i++)
		{
			listeT[i]=_listeT[i];
		}
		
		 numero= compteur;
		qmin= _qmin;
		qmax=_qmax;
		
		nbParents=_nbParents;
		parents=new int[nbParents];
		
		for(int i=0;i<nbParents;i++)
		{
			parents[i]=_parents[i];
		}
		
       		 ++compteur;
		cout<<"fin const R"<<endl;	
	
	}
	Reservoir()
	{}
	//Accesseurs
	void adApports(double * _apports, int nbH)
	{
		apports=_apports;
		apportAnnuel=0;
		for(int j=0; j< nbH;j++)
		{
           	 apportAnnuel= apportAnnuel+apports[j];
       		}
	}
	double getVinit() 
	{
        	return Vinit;
	}
	int getNbIntVmin()const
	{
		return nbIntVmin;
	}
	double getVminInt(int i)
	{
		return Vmin[i];
	}
	 double getVmin(int i)//faux
    	{
		int j;
		bool b=false;
		j=0;
		while(b==false)
		{
			if(i<intVmin[j]) b=true;
			else j++; 
		}
        	return Vmin[j];
    }
    int getIntVmin(int i)
    {
	return intVmin[i];
    }
	
    double getVmax()
    {
        return Vmax;
    }
	int getDeversement()
    {
        return reservoirDeDeversement;
    }
    int getNumero()
    {
        return numero;
    }
    double getApportAnnuel()
    {
        return apportAnnuel;
    }
	double getApport(int i)const
	{
		return apports[i];
	}
	int getDistance() const
	{
		return distance;
	}
	double getQmin() const
	{
		return qmin;
	}
	double getQmax() const
	{
		return qmax;
	}
	int getNbParents()const
	{
		return nbParents;
	}
	int* getParents()const
	{
		return parents;
	}
	int getNbTurbines()const
	{
		return nbT;
	}
	int getTurbine(int i) const
	{
		return listeT[i];
	}
	bool operator ==(Reservoir* r)
	{
		  if(not (Vmax==r->getVmax()))
		   {
			cout<<"Vmax differents : "<<Vmax<<" vs "<<r->getVmax()<<endl;
			return false;
		   }
		if(not(nbIntVmin==r->getNbIntVmin()))
		{
			cout<<"nbIntVmin différent : "<<nbIntVmin<<" vs "<<r->getNbIntVmin()<<endl;
			return false;
		}
		for(int i=0;i<nbIntVmin;i++)
		{
			if(intVmin[i]!=r->getIntVmin(i))
			{
				cout<<"inVmin "<<i<<" différent : "<<intVmin[i]<<" vs "<<r->getIntVmin(i)<<endl;
				return false;
			}
			if(Vmin[i]!= r->getVminInt(i))
			{
				cout<<"Vmin differents "<<Vmin[i]<<" vs "<<r->getVminInt(i)<<endl;
				return false;
			}

		}
		if(not(Vinit==r->getVinit()))
		{
			cout<<"Vinit différent : "<<Vinit<<" vs "<<r->getVinit()<<endl;
			return false;
		}
		if(not(nbParents==r->getNbParents()))
		{
			cout<<"nbParents différent : "<<nbParents<<" vs "<<r->getNbParents()<<endl;
			return false;
		}
		for(int i=0; i<nbParents; i++)
		{
			if(parents[i]!=r->getParents()[i])
			{
				cout<<" parents differents "<<parents[i]<<" vs "<<r->getParents()[i]<<endl;
				return false;
			 }
		}
		if(not(reservoirDeDeversement==r->getDeversement()))
		{
			cout<<"reservoirDeDeversement différent : "<<reservoirDeDeversement<<" vs "<<r->getDeversement()<<endl;
			return false;
		}
		if(not(nbT==r->getNbTurbines()))
		{
			cout<<"nbT différent : "<<nbT<<" vs "<<r->getNbTurbines()<<endl;
			return false;
		}
		for(int i=0; i<nbT;i++)
		{
			if(listeT[i]!=r->getTurbine(i))
			{
				cout<<"turbine diffreente "<<listeT[i]<<" vs "<<r->getTurbine(i)<<endl;
				return false;
			}
		}
		if(not(distance==r->getDistance()))
		{
			cout<<"distance différent : "<<distance<<" vs "<<r->getDistance()<<endl;
			return false;
		}
		if(not(qmin==r->getQmin()))
		{
			cout<<"qmin différent : "<<qmin<<" vs "<<r->getQmin()<<endl;
			return false;
		}
		if(not(qmax==r->getQmax()))
		{
			cout<<"qmax différent : "<<qmax<<" vs "<<r->getQmax()<<endl;
			return false;
		}
		for(int i=0;i<8760;i++)
		{
			if(apports[i]-r->getApport(i)>0.1 || apports[i]-r->getApport(i)<-0.1)
			{
				cout<<"apport "<<i<<" differents "<<apports[i]<<" vs "<<r->getApport(i)<<endl;
				return false;
			}

		}
		return true;
	
	
	
	
	//double* apports;
	//double apportAnnuel;
	
	
	
	
	
	
	}
};
int Reservoir:: compteur=0;

#endif

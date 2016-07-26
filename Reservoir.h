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
	Reservoir()
	{}
	//Accesseurs
	double getVinit() 
	{
        	return Vinit;
	}
	int getNbIntVmin()const
	{
		return nbIntVmin;
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
};
int Reservoir:: compteur=0;

#endif

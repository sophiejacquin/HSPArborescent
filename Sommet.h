#ifndef _SOMMET
#define _SOMMET
#include <vector>
#include<cstring>
#include <iostream>
#include <fstream>
using namespace std;
class Sommet2{
	public:
	vector<double> contenu;//Faire un essaie là dessu
	double valeur;
	Sommet2* pred;
	//Constructeur par défault :
	Sommet2()
	{
		valeur=-1;
		pred=NULL;
	}
	//Constructeurs specifiques
	Sommet2(vector<double> v)
	{
		contenu=v;
		valeur=-1;
		pred=NULL;
	}
	Sommet2(vector<double> v,double val)
	{
		contenu=v;
		valeur=val;
		pred=NULL;
	}
};
#endif

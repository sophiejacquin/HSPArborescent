#ifndef _eocascade_h
#define _eocascade_h
#include <vector>

template< class FitT>
class eocascade: public EO<FitT> {
public:
	

	
eocascade()
{
    
	nbEtats=0;
	nbReservoirs=0;
	last_fitness=0;
}
double getLast_fitness()
{
	return last_fitness;
}
void setLast_fitness(double _fit)
{
	last_fitness=_fit;
}
void setNbReservoirs(int _nbReservoirs)
{
	nbReservoirs=_nbReservoirs;
}
int getNbReservoirs()
{
	return nbReservoirs;
}
int getNbEtats()
{
	return nbEtats;
}
void adEtat(vector<double> etat)
{
	quantite.push_back(etat);
	eval.push_back(0);
	modif.push_back(true);
	nbEtats++;
}
void setQuantite(int etat,int reservoir,double qte)
{
	quantite[etat][reservoir]=qte;
	
}
double getQuantite(int etat,int reservoir)
{
	return quantite[etat][reservoir];
}
void setEval(int etat,double value)
{
	eval[etat]=value;
	
}
bool getModif(int etat)
{
	return modif[etat];
}
double getEval(int etat)
{
	return eval[etat];
}
void setModif(int etat,bool b)
{
	 modif[etat]= b;
}



virtual ~eocascade()
{

}

virtual string className() const { return "eocascade"; }
void printOn(ostream& os) const
{
	
	os << "time ;";
	for (int i=0;i<nbReservoirs;i++)
	{
		os<<"V_"<<i<<";";
	}
	for (int i=0;i<nbReservoirs;i++)
	{
		os<<"r_"<<i<<";";
	}
	for(int i=0; i<q[0].size(); i++)
	{
		os<<"T_"<<i<<";";
	}
	os<<endl;
	for(int i=0;i<nbEtats+1;i++)
	{
		os<<i<<";";
		for(int j=0;j<nbReservoirs;j++)
		{
			
			os<<V[i][j]<<";";

		}
		for(int j=0;j<nbReservoirs;j++)
		{
			
			if (i<nbEtats) os<<reserve[i][j]<<";";

		}
		for(int j=0;j<q[i].size();j++)
		{
			
			if (i<nbEtats) os<<q[i][j]<<";";

		}
		os<<endl;
	}
	EO<FitT>::printOn(os);
}
void readFrom(istream& is)
{
	EO<FitT>::readFrom(is);
}
void setV(vector < vector<double> > V_)
{
	V=V_;
}
void setQ(vector< vector <double> > q_)
{
	q=q_;
}
void setReserve(vector< vector <double> > reserve_)
{
	reserve=reserve_;
}
private:			  
	vector< vector < double > > quantite;
    	vector<double> eval;
    	vector<bool > modif;
	double last_fitness;
    	int nbEtats;
    	int nbReservoirs;
	vector < vector<double> > V;
	vector< vector <double> > q;
	vector< vector <double> > reserve;
};


#endif

/** -*- mode: c++; c-indent-level: 4; c++-member-init-indent: 8; comment-column: 35; -*-

The above line is usefulin Emacs-like editors
 */

/*
Template for creating a new representation in EO
================================================

Mandatory:
- a default constructor (constructor without any argument)
- the I/O functions (readFrom and printOn)

However, if you are using dynamic memory, there are 2 places
to allocate it: the default constructor (if possible?), or, more in
the EO spirit, the eoInit object, that you will need to write anyway
(template file init.tmpl).

But remember that a COPY CONSTRUCTOR will be used in many places in EO,
so make sure that the default copy constructor works, or, even better,
do write your own if in doubt.
And of course write the corresponding destructor!

*/

#ifndef _eocascade_h
#define _eocascade_h
#include <vector>
/**
 *  Always write a comment in this format before class definition
 *  if you want the class to be documented by Doxygen

 * Note that you MUST derive your structure from EO<fitT>
 * but you MAY use some other already prepared class in the hierarchy
 * like eoVector for instance, if you handle a vector of something....

 * If you create a structure from scratch,
 * the only thing you need to provide are
 *        a default constructor
 *        IO routines printOn and readFrom
 *
 * Note that operator<< and operator>> are defined at EO level
 * using these routines
 */
template< class FitT>
class eocascade: public EO<FitT> {
public:
	
  /** Ctor: you MUST provide a default ctor.
   * though such individuals will generally be processed
   * by some eoInit object
   */
	
  eocascade()
  {
    // START Code of default Ctor of an eocascade object
	  nbEtats=0;
	  nbReservoirs=0;
	last_fitness=0;
	
    // END   Code of default Ctor of an eocascade object
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
  /** Copy Ctor: you MUST provide a copy ctor if the default
   * one is not what you want
   * If this is the case, uncomment and fill the following
   */
  /*
  eocascade(const eocascade &)
  {
    // START Code of copy Ctor of an eocascade object
    // END   Code of copy Ctor of an eocascade object
  }
  */


  virtual ~eocascade()
  {
    // START Code of Destructor of an eoEASEAGenome object
    // END   Code of Destructor of an eoEASEAGenome object
  }

  virtual string className() const { return "eocascade"; }

  /** printing... */
    void printOn(ostream& os) const
    {
      // First write the fitness
      EO<FitT>::printOn(os);
      os << ' ';
    // START Code of default output
	for(int i=0;i<nbEtats;i++)
	{
		
		for(int j=0;j<nbReservoirs;j++)
		{
			if(i>0)os<<quantite[i][j]-quantite[i-1][j]<<"  "; 
			else os<<quantite[i][j]<<"  ";
		}
		os<<endl;
	}
	/** HINTS
	 * in EO we systematically write the sizes of things before the things
	 * so readFrom is easier to code (see below)
	 */

    // END   Code of default output
    }

  /** reading...
   * of course, your readFrom must be able to read what printOn writes!!!
   */
    void readFrom(istream& is)
      {
	// of course you should read the fitness first!
	EO<FitT>::readFrom(is);
    // START Code of input

	/** HINTS
	 * remember the eocascade object will come from the default ctor
	 * this is why having the sizes written out is useful
	 */

    // END   Code of input
      }

private:			   // put all data here
    // START Private data of an eocascade object
    // END   Private data of an eocascade object
    vector< vector < double > > quantite;
    	vector<double> eval;
    	vector<bool > modif;
	double last_fitness;
    	int nbEtats;
    	int nbReservoirs;
};

#endif



#ifndef eocascadeQuadCrossover2_H
#define eocascadeQuadCrossover2_H

#include <eoOp.h>


template<class GenotypeT>
class eocascadeQuadCrossover2: public eoQuadOp<GenotypeT>
{
public:
  /**
   * Ctor - no requirement
   */
// START eventually add or modify the anyVariable argument
  eocascadeQuadCrossover2()
  //  eocascadeQuadCrossover( varType  _anyVariable) : anyVariable(_anyVariable)
// END eventually add or modify the anyVariable argument
  {
   
  }

  /// The class name. Used to display statistics
  string className() const { return "eocascadeQuadCrossover2"; }

  /**
   * eoQuad crossover - _genotype1 and _genotype2 are the (future)
   *       offspring, i.e. _copies_ of the parents, to be modified
   * @param _genotype1 The first parent
   * @param _genotype2 The second parent
   */
  bool operator()(GenotypeT& _genotype1, GenotypeT & _genotype2)
  {
	//cout<<"debut croisement"<<endl;
	int i,j,k;
      bool oneAtLeastIsModified(true);
	int nbHeures=_genotype1.getNbEtats();
	int nbReservoirs=_genotype1.getNbReservoirs();
	int cont=0;
	int der=0;
	bool der2=false;
	double sum1=0;
	double sum2=0;
	bool b;
	for(i=0;i<nbHeures-1;i++)
	{
		sum1=sum1+_genotype1.getEval(i);
		b=false;
		if(cont>3)
		{
			b=true;
			j=0;
			while(b &&j<nbReservoirs)
			{
				if(_genotype1.getQuantite(i,j)>_genotype2.getQuantite(i+1,j)||_genotype2.getQuantite(i,j)>_genotype1.getQuantite(i+1,j))
				{
					b=false;
				}
				j++;
			}
		}
		if(b==true)
		{
			cont=0;
			if(sum2>sum1)
			{
				//if(der2) _genotype1.setModif(der,_genotype2.getModif(der));
				//else
				 
				for(k=der;k<=i;k++)
				{
					_genotype1.setEval(k,_genotype2.getEval(k));
					_genotype1.setModif(k,_genotype2.getModif(k));
					for(j=0;j<nbReservoirs;j++)
					{
						_genotype1.setQuantite(k,j,_genotype2.getQuantite(k,j));
					}
				}
				_genotype1.setModif(der,true);
				_genotype1.setModif(i+1,true);
				der2=true;
			}
			der=i+1;
			sum2=0;
			sum1=0;
			
			
		}
		cont++;
	}
	if(der>0)
	{
		if(sum2>sum1)
			{
				//if(der2) _genotype1.setModif(der,_genotype2.getModif(der));
				//else
				 
				for(k=der;k<nbHeures;k++)
				{
					_genotype1.setEval(k,_genotype2.getEval(k));
					_genotype1.setModif(k,_genotype2.getModif(k));
					for(j=0;j<nbReservoirs;j++)
					{
						_genotype1.setQuantite(k,j,_genotype2.getQuantite(k,j));
					}
				}
				_genotype1.setModif(der,true);
				//_genotype1.setModif(i+1,true);
				
			}
	}
	if(der==0)oneAtLeastIsModified=false;
  //cout<<"fin croisement"<<endl;
    return oneAtLeastIsModified;


    // END code for crossover of _genotype1 and _genotype2 objects
  }

//private:
// START Private data of an eocascadeQuadCrossover object
  //  varType anyVariable;		   // for example ...

// END   Private data of an eocascadeQuadCrossover object
};

#endif

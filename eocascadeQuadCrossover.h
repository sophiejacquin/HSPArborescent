#ifndef eocascadeQuadCrossover_H
#define eocascadeQuadCrossover_H

#include <eoOp.h>


template<class GenotypeT>
class eocascadeQuadCrossover: public eoQuadOp<GenotypeT>
{
public:

  eocascadeQuadCrossover(Systeme* _systeme,vector< vector<double> > _V,int nbHeures)

  {
	int i,h,j;
	systeme=_systeme;
	V= _V;
	int nbReservoirs=systeme->getNbReservoirs();
	//Calcul qTot;
	for(i=0;i<nbReservoirs;i++)
	{
		Reservoir* R=systeme->getReservoir(i);
	  	qTot.push_back(0);
	  	for(h=0;h<nbHeures;h++)
	  	{
	  		qTot[i]=qTot[i]+R->getApport(h);
	  	}
		for(j=0;j<R->getNbParents();j++)
		{
			int parent=R->getParents()[j];
			qTot[i]=qTot[i]+qTot[parent];
		}
	  }
  }

  string className() const { return "eocascadeQuadCrossover"; }

  bool operator()(GenotypeT& _genotype1, GenotypeT & _genotype2)
  {
//cout<<"debcross"<<endl;
	int i,j,h;
	Reservoir* R=systeme->getReservoir(1);
      bool oneAtLeastIsModified(true);
	int nbHeures=_genotype1.getNbEtats();
	int nbReservoirs=_genotype1.getNbReservoirs();
    	eoUniformGenerator<int> choixPt1(0,nbHeures*2/3);
	int pt3,pt31,pt32;
	int pt1=choixPt1();
	int pt2max=pt1+1;;
	int* pt2maxi=new int[nbReservoirs];
	//Recherche pt11 et p21;
	int pt11=pt1+1;
	int pt12=pt1+1;
	bool trouve1=false;
	bool trouve2=false;
	while (trouve1==false||trouve2==false)
	{
		trouve1=true;
		trouve2=true;
		for(i=0;i<nbReservoirs;i++)
		{
			if(_genotype1.getQuantite(pt1,i)>_genotype2.getQuantite(pt11,i))
			{
				trouve1=false;
			}
			if(_genotype2.getQuantite(pt1,i)>_genotype1.getQuantite(pt12,i))
			{
				trouve2=false;
			}
			
		}
		if(trouve1==false)pt11++;
		if(trouve2==false)pt12++;
	}
	if(pt11<pt12)
	{
		//pt1->pt11
		for(h=pt1+1;h<pt11;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qtemin=_genotype1.getQuantite(h-1,i);
				double qmaxC=R->getQmax();
				if(qmaxC<0) qmaxC=R->getQmin();
				if(R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}
				double qtemax;
				if(qmaxC>0) qtemax=qtemin+qmaxC;
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype1.getQuantite(h,p)*3600;
				}
				
				 qtemax=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qtemax)qtemax=(Vh-Vminh)/3600;
					}

				}
				if(qtemin<(Vi-Vmax)/3600)qtemin=(Vi-Vmax)/3600;
				if(qtemax>_genotype2.getQuantite(pt11,i))qtemax=_genotype2.getQuantite(pt11,i);
				if( qtemin<_genotype1.getQuantite(h-1,i))cout<<"pb croisement pt1>-pt11 :qte "<<qtemin<<" "<<h<<" "<<_genotype2.getQuantite(pt11,i)<<" "<<(Vi-Vmin)/3600<<endl;
				eoUniformGenerator <double> random(qtemin,qtemax);
				double qte=random();	
				_genotype1.setQuantite(h,i,qte);
				_genotype1.setModif(h,true);
			}
		}
		_genotype1.setModif(pt11,true);
		//pt11->pt12:
		for(h=pt11;h<pt12;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				_genotype1.setQuantite(h,i,_genotype2.getQuantite(h,i));
				
			}
			_genotype1.setEval(h,_genotype2.getEval(h));
		}
		//pt1->pt12:
		for(h=pt1+1;h<pt12;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qte=_genotype2.getQuantite(h-1,i);
				double qmin=qte;
				double qmax=_genotype1.getQuantite(pt12,i);
				double qmaxC=R->getQmax();
				//if(qmaxC<0) qmaxC=R->getQmin();
				if(qmaxC>0 && R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}

				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype2.getQuantite(h,p)*3600;
				}
				
				if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qmax)qmax=(Vh-Vminh)/3600;
					}

				}
				if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
                                if(qmaxC>0&& qte+qmaxC<qmax&& qte+qmaxC>qmin) qmax=qte+qmaxC;
                                if(qmin>qmax) cout<< qmin<< " "<< qmax<<" "<< (Vi-Vmin)/3600<<" "<<qte<<" "<<(Vi-Vmax)/3600<<endl;
				eoUniformGenerator <double> random(qmin,qmax);
				qte=random();
				_genotype2.setQuantite(h,i,qte);
				_genotype2.setModif(h,true);
			}
		}
		_genotype2.setModif(pt12,true);
		pt2max=pt12;

}

	else
	{
		//pt1->pt12
		for(h=pt1+1;h<pt12;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qte=_genotype2.getQuantite(h-1,i);
				double qmaxC=R->getQmax();
				if(qmaxC<0) qmaxC=R->getQmin();
				if(R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}
				if(qmaxC>0) qte=qte+qmaxC;
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype2.getQuantite(h,p)*3600;
				}
				if(qte>(Vi-Vmin)/3600) qte=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qte)qte=(Vh-Vminh)/3600;
					}

				}
				if(qte<(Vi-Vmax)/3600)qte=(Vi-Vmax)/3600;
				if(qte>_genotype1.getQuantite(pt12,i))qte=_genotype1.getQuantite(pt12,i);
				_genotype2.setQuantite(h,i,qte);
				_genotype2.setModif(h,true);
			}
		}
		_genotype2.setModif(pt12,true);
		//pt12->pt11:
		for(h=pt12;h<pt11;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				_genotype2.setQuantite(h,i,_genotype1.getQuantite(h,i));
				
			}
			_genotype2.setEval(h,_genotype1.getEval(h));
		}
		//pt1->pt11
		for(h=pt1+1;h<pt11;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qte=_genotype1.getQuantite(h-1,i);
				double qmaxC=R->getQmax();
				if(qmaxC<0) qmaxC=R->getQmin();
				if(R->getNbTurbines()>0)
				{
					int turbine=R->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}
				if(qmaxC>0) qte=qte+qmaxC;
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype1.getQuantite(h,p)*3600;
				}
				if(qte>(Vi-Vmin)/3600) qte=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qte)qte=(Vh-Vminh)/3600;
					}

				}
				if(qte<(Vi-Vmax)/3600)qte=(Vi-Vmax)/3600;
				if(qte>_genotype2.getQuantite(pt11,i))qte=_genotype2.getQuantite(pt11,i);
				_genotype1.setQuantite(h,i,qte);
				_genotype1.setModif(h,true);
			}
		}
		_genotype1.setModif(pt11,true);
			pt2max=pt11;


	}
	//choix aléatoire de pt3 :
	if(pt2max<nbHeures-2)
	{
		eoUniformGenerator<int> choixPt3(pt2max,nbHeures-2);
		pt3=choixPt3();
		//Remplissage section du milieu:
		for(j=pt2max;j<=pt3;j++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				double tmp=_genotype1.getQuantite(j,i);
				_genotype1.setQuantite(j,i,_genotype2.getQuantite(j,i));
				_genotype2.setQuantite(j,i,tmp);
				if(j==pt2max)
				{
					_genotype1.setModif(j,true);
					_genotype2.setModif(j,true);
				}
				else
				{
					double ev=_genotype1.getEval(j);
					_genotype1.setEval(j,_genotype2.getEval(j));
					_genotype2.setEval(j,ev);
				}
			}
		}
		//raccordement final:
		//recherche pt31 et pt32
		int pt31=pt3+1;
		int pt32=pt3+1;
		//cout<<"pt31 "<<pt31<<" pt32 "<<pt32<<endl;
		trouve1=false;trouve2=false;
		while (trouve1==false||trouve2==false)
		{
		trouve1=true;
		trouve2=true;
		for(i=0;i<nbReservoirs;i++)
		{
			if(_genotype1.getQuantite(pt3,i)>_genotype1.getQuantite(pt31,i))
			{
				trouve1=false;
			}
			if(_genotype2.getQuantite(pt3,i)>_genotype2.getQuantite(pt32,i))
			{
				trouve2=false;
			}
			
		}
		if(trouve1==false)pt31++;
		if(trouve2==false)pt32++;
		}
		//construction des chemins intermédiares: modifications ds cette boucle :
		for(h=pt3;h<pt31;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qte=_genotype1.getQuantite(h-1,i);
				double qmin=qte;
				double qmax =_genotype1.getQuantite(pt31,i);
				double qmaxC=R->getQmax();
				if(qmaxC<0) qmaxC=R->getQmin();
				if(R->getNbTurbines()>0)
				{
					int turbine=systeme->getReservoir(i)->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}
				
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype1.getQuantite(h,p)*3600;
				}
				if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qmax)qmax=(Vh-Vminh)/3600;
					}

				}
				if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
				if(qmaxC>0 &&qmax>qte+qmaxC&&qmin<qte+qmaxC) qmax=qte+qmaxC;
				eoUniformGenerator <double> random(qmin,qmax);
				 qte=random();
				//if(qte>_genotype1.getQuantite(pt31,i))qte=_genotype1.getQuantite(pt31,i);
				_genotype1.setQuantite(h,i,qte);
				_genotype1.setModif(h,true);
			}
		}
		_genotype1.setModif(pt31,true);
		for(h=pt3;h<pt32;h++)
		{
			for(i=0;i<nbReservoirs;i++)
			{
				
				R=systeme->getReservoir(i);
				//calcul de la quantité max débitable respectant vmin et vmax
				double qte=_genotype2.getQuantite(h-1,i);
				double qmin=qte;
				double qmax =_genotype2.getQuantite(pt32,i);
				double qmaxC=R->getQmax();
				if(qmaxC<0) qmaxC=R->getQmin();
				if(R->getNbTurbines()>0)
				{
					int turbine=systeme->getReservoir(i)->getTurbine(0);
					qmaxC=qmaxC+systeme->getTurbine(turbine)->getQMax(0);
				}
				
				double Vmin=R->getVmin(h);
				double Vmax=R->getVmax();
				//calcul Vi:
				double Vi=V[h][i];
				for(j=0;j<R->getNbParents();j++)
				{
					int p=R->getParents()[j];
					Vi=Vi+_genotype2.getQuantite(h,p)*3600;
				}
				if(qmax>(Vi-Vmin)/3600) qmax=(Vi-Vmin)/3600;
				if(i==1)
				{
					double Vh=Vi;
					for(j=h+1;j<nbHeures;j++)
					{
						Vh=Vh+R->getApport(j)*3600;
						double Vminh=R->getVmin(j);
						if((Vh-Vminh)/3600<qmax)qmax=(Vh-Vminh)/3600;
					}

				}
				if(qmin<(Vi-Vmax)/3600)qmin=(Vi-Vmax)/3600;
				if(qmaxC>0 &&qmax>qte+qmaxC&&qmin<qte+qmaxC) qmax=qte+qmaxC;
				eoUniformGenerator <double> random(qmin,qmax);
				 qte=random();
				_genotype2.setQuantite(h,i,qte);
				_genotype2.setModif(h,true);
			}
		}
		_genotype2.setModif(pt32,true);

	}
	else{
		for(j=pt2max;j<nbHeures;j++)
				{
					for(i=0;i<nbReservoirs;i++)
					{
						double tmp=_genotype1.getQuantite(j,i);
						_genotype1.setQuantite(j,i,_genotype2.getQuantite(j,i));
						_genotype2.setQuantite(j,i,tmp);
						if(j==pt2max)
						{
							_genotype1.setModif(j,true);
							_genotype2.setModif(j,true);
						}
						else
						{
							double ev=_genotype1.getEval(j);
							_genotype1.setEval(j,_genotype2.getEval(j));
							_genotype2.setEval(j,ev);
						}
					}
				}
	}
//cout<<"fin cross"<<endl;
    return oneAtLeastIsModified;

  }

private:

  vector< vector<double> > V;
  Systeme* systeme;
  vector<double> qTot;
};

#endif


#include"Systeme.h"
main()
{
	Systeme systeme("donnePrix","donneTurbines","donneReservoirs3");
	int nbH=8760;
	int nbReservoirs=systeme.getNbReservoirs();
	double quantite[nbH][nbReservoirs];
	int h,i,j;
	vector< vector<double> > V;
       for(i=0;i<NBHEURES;i++)
       {
   	vector<double> vec;
       	for(j=0;j<systeme.getNbReservoirs();j++)
       	{
       		
       		if(i==0)vec.push_back(systeme.getReservoir(j)->getVinit());
       		else vec.push_back(V[i-1][j]);
       		vec[j]=vec[j]+systeme.getReservoir(j)->getApport(i)*3600;
       	}
       	V.push_back(vec);
       }
	for(h=0;h<nbH;h++)
	{
		for(i=0;i<nbReservoirs;i++)
		{
			cin>>quantite[h][i];
			//qminC qmaxC;
			double qminC=systeme.getReservoir(i)->getQmin();
			double qmaxC=systeme.getReservoir(i)->getQmax();
			double qte;
			if(qte<qminC) cout<< "solution non réalisable car qmin non respecté en "<<i<<" à l'h "<<h<<" diff: "<<qte-qminC<<endl;
			else{
				if(h==0) qte =quantite[h][i];
				else qte= quantite[h][i]-quantite[h-1][i];
				//Calcul de Vi
				double Vi;
				if(h==0)Vi=systeme.getReservoir(i)->getVinit();
				else{
				 	Vi=V[h-1][i]-quantite[h-1][i]*3600;
					for(j=0;j<systeme.getReservoir(i)->getNbParents();j++)
					{
						int parent=systeme.getReservoir(i)->getParents()[j];
						Vi=Vi+quantite[h-1][parent]*3600;
					}
			  	}
				double qteT=0;
				if(systeme.getReservoir(i)->getNbTurbines()>0)
				{
					int turbine=systeme.getReservoir(i)->getTurbine(0);
				//calcul qminT et qmaxT:
					int Int=systeme->getTurbine(turbine)->getIntervalle(Vi);
					double qminT=systeme->getTurbine(turbine)->getQmin(Vi);
					if(qminT<0) qminT=0;
					double qmaxT=systeme->getTurbine(turbine)->getQMax(Int);
				//calcul de qteT:
				 	qteT=qte-qminC;
					if(qteT<qminT)qteT=0;
					if(qteT>qmaxT &&qmaxT>0)qteT=qmaxT;
				}
				if(qte-qteT>qmaxC) cout<< "solution non réalisable car qmax non respecté en "<<i<<" à l'h "<<h<<" diff: "<<qte-qteT -qmaxC<<endl;
				
			}
			
		}
	}

}

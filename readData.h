#ifndef ReadDataH
#define ReadDataH
#include <string>
#include <fstream>
#include <iostream>
#include <climits>
#include <cfloat>
#include <cstdlib>
#include "Systeme.h"
#include <sstream>

using namespace std;
class readData{
	private:
	Systeme sys;
	string directory;
	int apportTest;
	vector<string> data;
	int coefProd;


	public:
	readData(string _directory, int _apportTest,int nb,int _coefProd=1000)
	{
	        directory = _directory;
	        apportTest = _apportTest;
	        coefProd=_coefProd;
		sys.setNbHeures(nb);
	}
	Systeme getSysteme()
	{
		return sys;
	}
	void readAllFiles() {
        	readPrix();
        	readReservoirs();
        	readTurbines();
       
    	}
	void readPrix()
    	{
        	string filename = directory+"Prix/dataPrix.csv";
        	readDataFile(filename);
        	initPrix();
	}
	void readReservoirs() {
        	string filename = directory+"Reservoir/reservoirs.csv";
        	readDataFile(filename);
        	initReservoirs();  
    	}
	void readTurbines()
	{
        	string filename = directory+"Turbine/turbines.csv";
        	readDataFile(filename);
		initTurbines();
	}	
	void initPrix()
    	{
        	
        	vector<string> namePrix;
        	vector<int> nbParam;
        	for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); it = it + 2)
        	{
           	 namePrix.push_back(*it);
            	nbParam.push_back(atoi((*(it+1)).c_str()));
            	//sys.adNbParamPrix(atoi((*(it+1)).c_str()));
        	}

        	for (int i=0; i<namePrix.size(); i++)
        	{
            		string filename = directory + "Prix/" + namePrix[i] + ".csv";
            		readDataFile(filename);
            		double * prix = (double*) malloc(8760*sizeof(double));
            		int j=0;
            
            		for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); it = it + nbParam[i] + 1)
            		{
                		if (nbParam[i] == 1)
                   			prix[j] = atof((*(it+1)).c_str());
                		else
                    			prix[j] = atof((*(it+apportTest)).c_str()); /*TODO Revoir apportTest*/
                		j++;
            		}
            		sys.adPrix(prix);
        	}
		
    	}
	void initReservoirs()
    	{
         
        	for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); it = it + 14)
        	{
			//cout<<"deb boucle"<<endl;
			double Vinit=  atof((*(it+1)).c_str());
			//cout<<"'"<<*it<<"'"<<endl;
            		int nbIntVmin = atoi((*(it+3)).c_str());
            		double Vmax;
            		string Vmax_s = *(it+2);
            		if (Vmax_s.empty())
                		Vmax = FLT_MAX;
            		else
                		Vmax = atof(Vmax_s.c_str());
            		int* intVmin = (int*) malloc(nbIntVmin * sizeof(int));
            		double* Vmin = (double*) malloc(nbIntVmin * sizeof(double));
            		string intVmin_s = *(it+4);
            		string Vmin_s = *(it+5);
			
            		if (nbIntVmin>1)
           	 	{
                		std::size_t pos_intVmin_prec = 0;
                		std::size_t pos_Vmin_prec = 0;
                		for (int i=0; i<nbIntVmin; i++)
                		{
                    			size_t pos_intVmin = intVmin_s.find('+',pos_intVmin_prec);
                    			intVmin[i] = atoi((intVmin_s.substr(pos_intVmin_prec,pos_intVmin-pos_intVmin_prec)).c_str());
                    			size_t pos_Vmin = Vmin_s.find('+',pos_Vmin_prec);
                    			Vmin[i] = atof((Vmin_s.substr(pos_Vmin_prec,pos_Vmin-pos_Vmin_prec)).c_str());
                    			pos_intVmin_prec = pos_intVmin+1;
                    			pos_Vmin_prec = pos_Vmin+1;
                		}
            		}
            		else
            		{
                		intVmin[0] = atoi(intVmin_s.c_str());
                		Vmin[0] = atof(Vmin_s.c_str());
            		}
		;
            		double qmin=atof((*(it+6)).c_str());
			double qmax=atof((*(it+7)).c_str());
			int distance= atof((*(it+8)).c_str());
			int successeur=atof((*(it+9)).c_str())-1;
			int nbPred=atof((*(it+10)).c_str());
			int nbTurb=atof((*(it+12)).c_str());
			string pred_s=*(it+11);
			string turb_s=*(it+13);
			
			int* pred=(int*) malloc(nbPred * sizeof(int));
			int* turb=(int*) malloc(nbTurb * sizeof(int));
            		if (nbPred>1)
           	 	{
                		std::size_t pos_pred_prec = 0;
                		for (int i=0; i<nbPred; i++)
                		{
                    			size_t pos_pred = pred_s.find('+',pos_pred_prec);
                    			pred[i] = atoi((pred_s.substr(pos_pred_prec,pos_pred-pos_pred_prec)).c_str())-1;
                    			pos_pred_prec = pos_pred+1;
                		}
            		}
            		else
            		{
                		pred[0] = atoi(pred_s.c_str())-1;
            		}
            		if (nbTurb>1)
           	 	{
                		std::size_t pos_turb_prec = 0;
                		for (int i=0; i<nbTurb; i++)
                		{
                    			size_t pos_turb = turb_s.find('+',pos_turb_prec);
                    			turb[i] = atoi((turb_s.substr(pos_turb_prec,pos_turb-pos_turb_prec)).c_str())-1;
                    			pos_turb_prec = pos_turb+1;
                		}
            		}
            		else if (nbTurb==1)
            		{
				
                		turb[0] = atoi(turb_s.c_str())-1;
            		}

            		Reservoir* R= new Reservoir (Vinit,Vmax,nbIntVmin,intVmin, Vmin,successeur, distance,nbTurb, turb,qmin,qmax,nbPred, pred);
            		sys.adReservoir(R);   

            
        	}

        	
        	ostringstream oss;
       	 	oss << apportTest;
        	string filename = directory+"Apport/apports_serie"+oss.str()+".csv";
        	readDataFile(filename);
        	int nbReservoirs = sys.getNbReservoirs();
        
        	double **apports = (double**) malloc(nbReservoirs * sizeof(double *));
        	for (int i=0; i<nbReservoirs; i++)
            		apports[i] = (double*) malloc(8760 * sizeof(double));

        	int j=0;
        	for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); it = it + nbReservoirs)
        	{
            		for (int i=0; i<nbReservoirs; i++)   
                	apports[i][j] = atof((*(it+i)).c_str());
            		j++;
        	}

        	for (int i=0; i<nbReservoirs; i++)
        	{
            		
            		sys.getReservoir(i)->adApports(apports[i],sys.getNbHeures());
        	}
		
    	}


	void initTurbines()
    	{
		
        	for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); it = it + 7)
        	{
			
            		double prodMin = atof((*(it+2)).c_str());
            		int distance = atoi((*(it+1)).c_str());
            		int resParent = atoi((*(it+6)).c_str())-1;
            		int catPrix = atoi((*(it+5)).c_str());
			int nbPieces = atoi((*(it+3)).c_str());
            		int nbInt = atoi((*(it+4)).c_str());

            		Turbine T(nbInt,nbPieces,  prodMin, distance, sys.getPrix(catPrix), resParent);

            
            		sys.adTurbine(T);
		}
        

        // Production
        int nbTurbines = sys.getNbTurbines();
	
        for (int t=0; t<nbTurbines; t++)
        {
		
            ostringstream oss;
            oss << t+1;
            string filename = directory+"Turbine/Turb"+oss.str()+".csv";
            readDataFile(filename);


            Turbine * T = sys.getTurbine(t);
            int nbInt = T->getNbInt();
            int nbPieces = T->getNbPieces();
            double* intReservoirs = (double *) malloc(nbInt*sizeof(double));
            double* pieces = (double*) malloc((nbPieces-1)*sizeof(double));
            double* qmax = (double*) malloc(nbInt*sizeof(double));
            double** production = (double**) malloc(nbPieces*sizeof(double*));
		
            std::vector<std::string>::iterator it = data.begin();
            for (int i=0; i<nbPieces; i++)
                production[i] = (double*) malloc(nbInt*sizeof(double));
	
            for (int j=0; j<nbInt; j++)
            {
                if (nbInt>1)
                {
                   // cout<<"atof((*(it+j+1)).c_str());="<<atof((*(it+j+1)).c_str())<< " " << j << endl;;
                    intReservoirs[j] = atof((*(it+j+1)).c_str());
                }
                else
                    intReservoirs[j] = 0.0;
                
                qmax[j] = 0.0;
            }

		
            if (nbInt>1)
                it += nbInt + 1;

            int i=0;
		
            for ( ; it != data.end(); it = it + nbInt + 1)
            {
		
                double tmp_d =atof((*it).c_str());
                if (i>0)pieces[i-1] = tmp_d;//*3600;
                
                for (int j=0; j<nbInt; j++)
                {
		    
                   if(t==0 && j==0)cout<<*(it+j+1)<<endl;
                    string tmp_s = *(it+j+1);
			
                    production[i][j] = atof(tmp_s.c_str())*coefProd;
			
                    
                    if (tmp_s != "")
                        qmax[j] = tmp_d;
                }
                i++;
            }
		
            double** coeff = (double**) malloc((nbPieces-1)*sizeof(double*));
            for (int i=0; i<nbPieces; i++)
                coeff[i] = (double*) malloc(nbInt*sizeof(double));
            for (int i=1; i<nbPieces; i++)
            {
                for (int j=0; j<nbInt; j++)
                {
                   
			if(pieces[i-1]>qmax[j]) coeff[i-1][j] =0;
                        else if(i>1) coeff[i-1][j] = (production[i][j]-production[i-1][j]) / (pieces[i-1]-pieces[i-2]);
			 
			else coeff[i-1][j] = (production[i][j]-production[i-1][j]) / (pieces[i-1]);
                    if(t==0 && i==16 && j==15)cout<<production[i][j]<<"   "<<production[i-1][j]<<" coeff="<<coeff[i-1][j]<<" pieces[i]="<<pieces[i-1]<<" pieces[i-1="<<pieces[i-2]<<" qmax="<<qmax[j]<<endl;
                }
            }
	  
            T->setProd(nbInt, nbPieces-1, intReservoirs, pieces, qmax, coeff);
           
        }
        
    }

	bool readDataFile(string fileName) 
	{
		bool ret = true;
        	std::ifstream file(fileName.c_str(), std::ios::in);
        	std::string buff;
        	data.clear();
        	if (file) 
        	{
         		getline(file, buff);
           	 	while (getline(file, buff, ';'))
            		{
               		 	if (buff[0] != '\n')
                  			data.push_back(buff);
                		else if(buff != "\n")
                    			data.push_back(buff.substr(1,buff.size()));
           	 	}
        	}
        	else
            		ret = false;
        	file.close();
        	return ret;
    	}
	
};


#endif

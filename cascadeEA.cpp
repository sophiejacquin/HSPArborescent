#include <iostream>
using namespace std;
#include "eo"
#include "utils/eoRealVectorBounds.h"
#include "eocascade.h"
#include "eocascadeInit.h"
#include "eocascadeInit2.h"
#include "eocascadeInit3.h"
#include "eocascadeInit4.h"
#include "eocascadeInit5.h"
#include "eocascadeInit6.h"
#include "eocascadeEvalFunc.h"
#include "eocascadeQuadCrossover.h"
#include "eocascadeQuadCrossover2.h"
#include "eocascadeMutation.h"
#include "eocascadeMutation2.h"
#include "eocascadeStat.h"

typedef double MyFitT ;	
typedef eocascade<MyFitT> Indi;     

#include <make_pop.h>
eoPop<Indi >&  make_pop(eoParser& _parser, eoState& _state, eoInit<Indi> & _init)
{
  return do_make_pop(_parser, _state, _init);
}


#include <make_continue.h>
eoContinue<Indi>& make_continue(eoParser& _parser, eoState& _state, eoEvalFuncCounter<Indi> & _eval)
{
  return do_make_continue(_parser, _state, _eval);
}


#include <make_checkpoint.h>
eoCheckPoint<Indi>& make_checkpoint(eoParser& _parser, eoState& _state, eoEvalFuncCounter<Indi>& _eval, eoContinue<Indi>& _continue)
{
  return do_make_checkpoint(_parser, _state, _eval, _continue);
}


#include <make_algo_scalar.h>
eoAlgo<Indi>&  make_algo_scalar(eoParser& _parser, eoState& _state, eoEvalFunc<Indi>& _eval, eoContinue<Indi>& _continue, eoGenOp<Indi>& _op, eoDistance<Indi> *_dist = NULL)
{
  return do_make_algo_scalar(_parser, _state, _eval, _continue, _op, _dist);
}


#include <make_run.h>
#include <eoScalarFitness.h>
#include"Systeme.h"
#include <stdio.h>
#include <string.h>
void run_ea(eoAlgo<Indi>& _ga, eoPop<Indi>& _pop)
{
  do_run(_ga, _pop);
}


void make_help(eoParser & _parser);
int main(int argc, char* argv[])
{
	int i,j;
try
{
	eoParser parser(argc, argv);  
	eoState state;    
	int  NBHEURES =parser.createParam(8760,"nbHeures", "nombre d'heures dans la plannification",'N', "My application").value();
	std::string donnePrix = parser.createParam(std::string("./donnePrix"), "donnePrix", "donne pour prix", 'A', "My application").value();
	char Input[600];
	strcpy(Input,donnePrix.c_str());
	std::string donneT = parser.createParam(std::string("./donneTurbines"), "donneTurbines", "donne pour turbines", 'B', "My application").value();
	char InputTurb[600];
	strcpy(InputTurb,donneT.c_str());
	std::string donneR = parser.createParam(std::string("./donneReservoirs3"), "donneReservoirs", "donne pour reservoirs", 'C', "My application").value();
	char InputR[600];
	strcpy(InputR,donneR.c_str());
        Systeme systeme(Input,InputTurb,InputR);
        std::string sortie = parser.createParam(std::string("exit.csv"), "sortie", "sortie résultats", 'x', "My application").value();
        char sortieChar[600];
        strcpy(sortieChar, sortie.c_str());
	cout <<"systeme créé"<<endl;
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
	eocascadeEvalFunc<Indi> plainEval(&systeme,V);
	eoEvalFuncCounter<Indi> eval(plainEval);
	eocascadeInit3<Indi> init(V,&systeme,NBHEURES);
    	eocascadeQuadCrossover<Indi> cross1(&systeme,V,NBHEURES);
 	double cross1Rate = parser.createParam(0.164, "cross1Rate", "Relative rate for crossover 1", 'y', "Variation Operators").value();
	eoPropCombinedQuadOp<Indi> cross(cross1, cross1Rate);
	eocascadeQuadCrossover2<Indi> cross2;//(eoParser parser);
 	double cross2Rate = parser.createParam(1.0, "cross2Rate", "Relative rate for crossover 2", 'x', "Variation Operators").value();
	cross.add(cross2, cross2Rate);
	eocascadeMutation<Indi> mut1(V,&systeme);
	double mut1Rate = parser.createParam(0.402, "mut1Rate", "Relative rate for mutation 1", '1', "Variation Operators").value();
	eoPropCombinedMonOp<Indi> mut(mut1, mut1Rate);
	double lng = parser.createParam(180, "lng", "longueur mut 2", 'L', "Variation Operators" ).value();
 	eocascadeMutation2<Indi> mut2(V,&systeme,lng,100);//(eoParser parser);
	double mut2Rate = parser.createParam(0.012, "mut2Rate", "Relative rate for mutation 2", '2', "Variation Operators").value();
	 mut.add(mut2, mut2Rate);
	vector<double> delta;
	vector<double> pas;
/*delta.push_back(0.5);
pas.push_back(1);
delta.push_back(0.5);
pas.push_back(1);
delta.push_back(0);
pas.push_back(2);
delta.push_back(0);
pas.push_back(1);
delta.push_back(0.5);
pas.push_back(1);
delta.push_back(1.5);
pas.push_back(3);
delta.push_back(1.5);
pas.push_back(3);
 eocascadeMutation3<Indi> mut4(V,&systeme,3,100);//100,delta,pas);//(eoParser parser);
double mut4Rate = parser.createParam(1.0, "mut4Rate", "Relative rate for mutation 2", '4', "Variation Operators").value();
 mut.add(mut4, mut4Rate);*/

	double pCross = parser.createParam(0.718, "pCross", "Probability of Crossover", 'C', "Variation Operators" ).value();

	if ( (pCross < 0) || (pCross > 1) )
		throw runtime_error("Invalid pCross");
	double pMut = parser.createParam(0.089, "pMut", "Probability of Mutation", 'M', "Variation Operators" ).value();
    // minimum check
	if ( (pMut < 0) || (pMut > 1) )
		throw runtime_error("Invalid pMut");
	eoSGAGenOp<Indi> op(cross, pCross, mut, pMut);
	eoPop<Indi>& pop   = make_pop(parser, state, init);
  // stopping criteria
	eoContinue<Indi> & term = make_continue(parser, state, eval);
	eoTimeContinue<Indi> tempsStop(40000);
// output
	eoCheckPoint<Indi> & checkpoint = make_checkpoint(parser, state, eval, term);
	checkpoint.add(tempsStop);
//cout<<"stop"<<endl;
	eocascadeStat<Indi>   myStat;      
	checkpoint.add(myStat);
	eoIncrementorParam<unsigned> generationCounter("Gen.");
	checkpoint.add(generationCounter);
  // need to get the name of the redDir param (if any)
	std::string dirName =  parser.getORcreateParam(std::string("Res"), "resDir", "Directory to store DISK outputs", '\0', "Output - Disk").value() + "/";
// those need to be pointers because of the if's
	eoStdoutMonitor *myStdOutMonitor;
	eoFileMonitor   *myFileMonitor;
#ifdef HAVE_GNUPLOT
	eoGnuplot1DMonitor *myGnuMonitor;
#endif

  // now check how you want to output the stat:
	bool printcascadeStat = parser.createParam(false, "coutcascadeStat", "Prints my stat to screen, one line per generation", '\0', "My application").value();
	bool filecascadeStat = parser.createParam(false, "filecascadeStat", "Saves my stat to file (in resDir", '\0', "My application").value();
	bool plotcascadeStat = parser.createParam(false, "plotcascadeStat", "On-line plots my stat using gnuplot", '\0', "My application").value();
	
  // should we write it on StdOut ?
	if (printcascadeStat)
	{
		myStdOutMonitor = new eoStdoutMonitor(false);
		state.storeFunctor(myStdOutMonitor);
		checkpoint.add(*myStdOutMonitor);
		myStdOutMonitor->add(generationCounter);
		myStdOutMonitor->add(eval);
		myStdOutMonitor->add(myStat);
	}

  // first check the directory (and creates it if not exists already):
	if (filecascadeStat || plotcascadeStat)
		if (! testDirRes(dirName, true) )
			throw runtime_error("Problem with resDir");
	if (filecascadeStat)
	{
		myFileMonitor = new eoFileMonitor(dirName + "myStat.xg");
		state.storeFunctor(myFileMonitor);
		checkpoint.add(*myFileMonitor);
		myFileMonitor->add(generationCounter);
		myFileMonitor->add(eval);
		myFileMonitor->add(myStat);
	}

#ifdef HAVE_GNUPLOT
  // should we PLOT it on StdOut ? (one dot per generation, incremental plot)
	if (plotcascadeStat)
	{
		myGnuMonitor = new eoGnuplot1DMonitor(dirName+"plot_myStat.xg",minimizing_fitness<Indi>());
		state.storeFunctor(myGnuMonitor);
      // and of course to add the monitor to the checkpoint
		checkpoint.add(*myGnuMonitor);
      // and the different fields to the monitor (X = eval, Y = myStat)
		myGnuMonitor->add(eval);
		myGnuMonitor->add(myStat);
	}
#endif

  // algorithm (need the operator!)
	eoAlgo<Indi>& ga = make_algo_scalar(parser, state, eval, checkpoint, op);
	make_help(parser);
	cout<<"go"<<endl;
  //// GO
	cout<<"eval pop ini"<<endl;
	apply<Indi>(eval, pop);

	cout<<"pop ini evaluée"<<endl;

	cout<<pop.best_element().fitness()<<endl;

  	run_ea(ga, pop); // run the ga


  	cout << endl;
	cout << "Meilleur:" <<-1*(double)pop.best_element().fitness();
	cout << endl;
        pop.sort();
        plainEval.details(pop[0]);
        cout << "Final Population\n";
        std::ofstream fichierSortie (sortieChar, ios::out | ios::trunc);
        if(fichierSortie)
        {
            pop.best_element().printOn(fichierSortie);
        }
        else
            pop.best_element().printOn(cout);
  }
  catch(exception& e)
  {
	cout << e.what() << endl;
  }
  return 0;
}

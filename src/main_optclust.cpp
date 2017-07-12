 /**
  * The MIT License (MIT)
  *
  * Copyright (c) 2014 Behrouz Babaki, Tias Guns, Siegfried Nijssen
  *
  * Permission is hereby granted, free of charge, to any person obtaining a copy
  * of this software and associated documentation files (the "Software"), to deal
  * in the Software without restriction, including without limitation the rights
  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  * copies of the Software, and to permit persons to whom the Software is
  * furnished to do so, subject to the following conditions:
  *
  * The above copyright notice and this permission notice shall be included in
  * all copies or substantial portions of the Software.
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  * THE SOFTWARE.
  */

#include "reader.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <getopt.h>
#include <boost/iterator/iterator_concepts.hpp>

#include "init.h"
#include "scip_exception.hpp"
#include "pricer_optclust.h"
#include "branch_ryanfoster.h"
#include "cons_mlcl.h"
#include "SiegBnBCons.h"
#include "cccg_githash.h"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>


using namespace std;
using namespace scip;

bool isMember(SCIP*, SCIP_VAR*, SCIP_CONS*);
void print_usage(char*);
void print_config(Stabilisation stab, bool findAll, unsigned int nclusters, unsigned int ndatapoints, unsigned int nconstraints, unsigned int firstCoefStrategy, char* epsilon_str);

int main(int argc , char* argv[]){

  char* inputFile = NULL;
  char* clusterCount = NULL;
  char* initFile = NULL;
  char* initDir = NULL;
  char* coefFile = NULL;
  char* columnFile = NULL;
  char* consFile = NULL;
  char* epsilonStr = NULL;
  char* muStr = NULL;
  char* firstCoefStr = NULL;
  char* firstCoefStrategyStr = NULL;
  char* timeoutStr = NULL;
  SCIP_Real mu = 0.99, 
    epsilon = 1e-6, 
    firstCoef = 300,
    timeout = -1;
  unsigned int firstCoefStrategy = 0;
  char* stabFlag = NULL;
  Stabilisation stab = ST_None;
  bool findAll = false;
  bool quiet = false;

  int c;
  /* getopt note: an option character followed by a colon takes a required argument */
  while ( (c = getopt(argc, argv, "d:k:i:c:s:m:e:g:v:aqn:r:t:x:")) != -1)
    switch (c){
    case 'd':
      inputFile = strdup(optarg);
      break;

    case 'k':
      clusterCount = strdup(optarg);
      break;

    case 'i':
      initFile = strdup(optarg);
      break;

    case 'r':
      initDir = strdup(optarg);
      break;

    case 'c':
      coefFile = strdup(optarg);
      break;

    case 's':
      columnFile = strdup(optarg);
      break;

    case 'n':
       consFile = strdup(optarg);
       break;
      
    case 'g':
      stabFlag = strdup(optarg);
      if (!strcmp(stabFlag,"bound"))
          stab = ST_Bound;
      else if (!strcmp(stabFlag,"both"))
          stab = ST_Both;
      else if (!strcmp(stabFlag,"smooth"))
          stab = ST_Smooth;
      else if (!strcmp(stabFlag,"smooth2"))
          stab = ST_Smooth2;
      else if (!strcmp(stabFlag,"smoothAuto"))
          stab = ST_SmoothAuto;
      else
	{
	  cout << "invalid argument for parameter 'g'" << endl;
	  exit (EXIT_FAILURE);
	}
      break;
      
    case 'm':
      muStr = strdup(optarg);
      break;

    case 'e':
      epsilonStr = strdup(optarg);
      break;
      
    case 'v':
      firstCoefStr = strdup(optarg);
      break;
      
    case 'a':
      findAll = true;
      break;
      
    case 'q':
      quiet = true;
      break;

    case 't':
      firstCoefStrategyStr = strdup(optarg);
      break;

    case 'x':
      timeoutStr = strdup(optarg);
      break;

    default:
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  
  
  if (inputFile == NULL || clusterCount == NULL)
    {
      print_usage(argv[0]);
      exit (EXIT_FAILURE);
    }
  else
    {
  InstanceReader reader(inputFile);
  vector<vector<double> > distances = reader.get_distances();
  int number_of_elements = reader.get_number_of_instances();
  int number_of_clusters = atoi(clusterCount);
   static const char* OPTCLUST_PRICER_NAME = "MSSC_Pricer";
   SCIP* scip = NULL;

   if (muStr != NULL){
     stringstream ss (muStr);
     ss >> mu;
     if (ss.bad()){
       cout << "The value entered for parameter 'mu' is not valid." << endl;
       exit (EXIT_FAILURE);
     }
   }

   if (epsilonStr != NULL){
     stringstream ss (epsilonStr);
     ss >> epsilon;
     if (ss.bad()){
       cout << "The value entered for parameter 'epsilon' is not valid." << endl;
       exit (EXIT_FAILURE);
     }
   }

   if (firstCoefStr != NULL){
     stringstream ss (firstCoefStr);
     ss >> firstCoef;
     if (ss.bad()){
       cout << "The value entered for parameter 'first coefficients' is not valid." << endl;
       exit (EXIT_FAILURE);
     }
   }

   if (firstCoefStrategyStr != NULL){
     stringstream ss (firstCoefStrategyStr);
     ss >> firstCoefStrategy;
     if (ss.bad()){
       cout << "The value entered for parameter 'first coefficients strategy' is not valid." << endl;
       exit (EXIT_FAILURE);
     }
   }

   if (timeoutStr != NULL){
     stringstream ss (timeoutStr);
     ss >> timeout;
     if (ss.bad()) {
       cout << "The value entered for parameter 'timeout' is not valid." << endl;
       exit (EXIT_FAILURE);
     }
   }


  try
    {

   /**************
    * Setup SCIP *
    **************/

   /* initialize SCIP environment */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* set verbosity parameter */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   // SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) );
   // SCIP_CALL (SCIPsetCharParam (scip, "lp/initalgorithm" , 'c') );
   // SCIP_CALL (SCIPsetCharParam (scip, "lp/resolvealgorithm" , 'c') );

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, "OptClust", 0, 0, 0, 0, 0, 0, 0) );

   /* set the timeout */
   if (timeout > 0.0) 
     SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeout));

   /* remove a variable which has been set to zero for 5 LP iterations   */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/colagelimit", 5) );
 
   /* add setcovering constraints */
   char con_name[255];
   vector<SCIP_CONS*> setcover_con( number_of_elements );
   for (int i = 0; i < number_of_elements; ++i)
   {
         SCIP_CONS* con;
         SCIPsnprintf(con_name, 255, "S%d", i);


         SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 0, NULL, NULL,
                     1.0,                    /* lhs */
                     1.0,                    /* rhs */
                     true,                   /* initial */
                     true,                   /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     true,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
         SCIP_CALL( SCIPaddCons(scip, con) );
         setcover_con[i] = con;
   }

   vector<SCIP_VAR *> pertVars;
   if (stab == ST_Bound || stab == ST_Both){
     /* create perturbation variables */
     char pert_var_name [255];
     for (int counter = 0 ; counter < number_of_elements ; counter++){
       SCIP_VAR* var;
       SCIPsnprintf(pert_var_name, 255, "PV%d", counter);
       SCIP_CALL( SCIPcreateVar(scip, &var, pert_var_name,
				-mu,                     // lower bound
				mu,                      // upper bound
				firstCoef,                       // objective (to be changed later)
				SCIP_VARTYPE_CONTINUOUS,     // variable type
				true,                   /* initial */
				true,                   /* removable */ 
				0, 0, 0, 0, 0) );

       /* add new variable to LP and to the set covering constraints */
       SCIP_CALL( SCIPaddVar(scip, var) );
       SCIP_CALL( SCIPaddCoefLinear (scip, setcover_con[counter] , var, 1.0) );
       pertVars.push_back(var);
       
     }

   }


   /* add side constraint */
      SCIP_CONS* side_con;
      SCIPsnprintf(con_name, 255, "Sd");
      SCIP_CALL( SCIPcreateConsLinear(scip, &side_con, con_name, 0, NULL, NULL,
                  -SCIPinfinity(scip),    /* lhs */
                  number_of_clusters,     /* rhs */
                  true,                   /* initial */
                  true,                   /* separate */
                  true,                   /* enforce */
                  true,                   /* check */
                  true,                   /* propagate */
                  false,                  /* local */
                  true,                   /* modifiable */
                  false,                  /* dynamic */
                  false,                  /* removable */
                  false) );               /* stickingatnode */
      SCIP_CALL( SCIPaddCons(scip, side_con) );

      
      // include the ML/CL handler (used in brancher)
      SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrMLCL(scip, setcover_con), TRUE) );
      
      vector<pair<unsigned int, unsigned int> > mlConstraints;
      vector<pair<unsigned int, unsigned int> > clConstraints;
      if (consFile != NULL)
         SCIP_CALL( init_MlCl_constraints(consFile, scip, mlConstraints, clConstraints) );
      
      vector<double> upper;
      vector<double> lower;
      
      if (initDir != NULL)
         SCIP_CALL( init_sol_dir(initDir, scip,
				 distances, 
				 setcover_con, side_con,
				 mlConstraints, clConstraints,
				 upper, lower
                                       ) );
	
      if (initFile != NULL)
	SCIP_CALL( initialize_solution(initFile, scip,
				       distances, 
				       setcover_con, side_con,
				       mlConstraints, clConstraints,
				       upper, lower
                                       ) );
      
      int verbose = 0;
      if (verbose >= 3){
      for (int counter = 0; counter < number_of_elements; counter++)
	cout << setprecision(5) << "Upper: " << upper[counter] << "\tLower: " 
	     << lower[counter]  << "\tRelative Gap: " 
	     << (upper[counter] - lower[counter])/lower[counter] << endl;
      }
      
      if ((stab == ST_Both || stab == ST_Bound) && firstCoefStr == NULL && firstCoefStrategy)
	for (int counter = 0; counter < number_of_elements; counter++){
	  SCIP_Real new_pertVar_obj;
	  switch (firstCoefStrategy){
	  case 1:
	    new_pertVar_obj = lower[counter];
	    break;
	  case 2:
	    new_pertVar_obj = upper[counter];
	    break;
	  case 3:
	    new_pertVar_obj = (lower[counter] + upper[counter]) / 2;
	    break;
	  default:
	    cout << "The specified strategy for assigning first set of coefficients is not valid" 
		 << endl;
	    exit(EXIT_FAILURE);
	    break;
	  }
	  SCIP_CALL (SCIPchgVarObj(scip, pertVars[counter], new_pertVar_obj));
	}


      unsigned int number_of_constraints = mlConstraints.size() + clConstraints.size();


   /* include clustering pricer */
      ObjPricerOptClust* optclust_pricer_ptr = new
	ObjPricerOptClust(scip,
			  OPTCLUST_PRICER_NAME,
			  number_of_elements,
			  number_of_clusters,
			  distances,
			  setcover_con,
			  side_con,
                          stab,
			  findAll,
			  coefFile,
			  columnFile);

      if (stab == ST_Bound || stab == ST_Both){
	optclust_pricer_ptr -> setPertVars(pertVars, mu, epsilon);
      }
      if ((stab == ST_Smooth || stab == ST_Smooth2 || stab == ST_SmoothAuto) && epsilonStr != NULL){
        // set the alpha parameter fixed
        optclust_pricer_ptr->setEpsilon(epsilon);
      }

   SCIP_CALL( SCIPincludeObjPricer(scip, optclust_pricer_ptr, true) );

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, OPTCLUST_PRICER_NAME)) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeObjBranchrule(scip, new BranchruleRyanFoster(scip, setcover_con), TRUE) );

   /* print the setup */
   print_config(stab, findAll, number_of_clusters, number_of_elements, number_of_constraints, firstCoefStrategy, epsilonStr);

   /* set the message handler to be quiet*/
   if (quiet)
     SCIPsetMessagehdlrQuiet (scip, true);

   cout << "\n\n\nstarting to solve ...\n" << endl;

   /*************
    *  Solve    *
    *************/

   SCIP_CALL( SCIPsolve(scip) );


   /**************
    * Statistics *
    *************/
   //SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   
   

   /* printing the resulting clustering */

   int number_of_variables = SCIPgetNVars(scip);
   SCIP_VAR ** variables = SCIPgetVars(scip);

   SCIP_SOL* sol = SCIPgetBestSol(scip);
  
   // following to check that all elements are covered, and none twice
   vector<bool> cover(number_of_elements, false);
   
   int clusnum = 0;
   for (int counter = 0; counter < number_of_variables; counter++){
     if (SCIPvarGetType(variables[counter]) == SCIP_VARTYPE_BINARY && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, variables[counter])))
     {
       vector<bool> cluster(number_of_elements, false);
       for (int i = 0 ; i < number_of_elements; i++) {
           if (isMember(scip, variables[counter], setcover_con[i]))
               cluster[i] = true;
       }

       // print cluster
       clusnum++;
       char msg[100];
       SCIPsnprintf(msg, 100, "Cluster %i: ", clusnum);
       print_cluster(cluster, msg);
       
       // check constraints NOT branching, ONLY given ones
       // (we could be in a different node of the one where the optimal solution was found)
       if (!cluster_satisfies_mlcl(cluster, mlConstraints, clConstraints)) {
           cout << "ERROR\n";
           cout << "The above cluster does not satisfy the following set of constraints:\n";
           ConshdlrMLCL::print_constraints(scip);
           exit(1);
       }
       
       // check unique cover
       for (unsigned i=0; i!=cluster.size(); i++) {
           if (cluster[i]) {
               if (cover[i]) {
           cout << "ERROR\n";
           cout << "The above cluster has overlap with another cluster for element " << i << " !??\n";
           exit(1);
               }
               cover[i] = true;
           }
       }
                   
     }
   }
   
   // check all examples covered
   if ( find(cover.begin(), cover.end(), false) != cover.end() ) {
           cout << "Not all elements are coverd! \n";
           exit(1);
   }
   
   cout << "Solution verified.\n";

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   } catch(SCIPException& exc)
   {
      cerr << exc.what() << endl;
      exit(exc.getRetcode());
   }

  
   return EXIT_SUCCESS;
    }
}

void print_usage(char* fileName){
  cout << "\nusage: " << fileName << " -d <datafile> -k <number of clusters> [OPTIONS]" << endl;
  cout << "\nOptions:" << endl;
  cout << "\t-i FILE,\t\tUse FILE as initial clustering file." << endl;
  cout << "\t-x [T],\t\t\tTerminate the program after a timeout of T seconds." << endl;
  cout << "\t-r DIRECTORY,\t\tUse DIRECTORY as initial clustering directory." << endl;
  cout << "\t-n FILE,\t\tUse FILE as constraint file." << endl;
  cout << "\t-g bound|both|smooth|smooth2|smoothAuto." << endl;
  cout << "\t\t\t\tUse this type of stabilization." << endl;
  cout << "\t-e [E],\t\t\tUse E for epsilon value." << endl;
  cout << "\t-m [M],\t\t\tUse M for mu value." << endl;
  cout << "\t-v [C]\t\t\tUse C for the value of first coefficients." << endl;
  cout << "\t-a \t\t\tAdd all columns with negative reduced cost in each iteration." << endl;
  cout << "\t-t 0|1|2\t\tUse this strategy for computing the values of first coefficients." << endl;
  cout << "\t-q\t\t\tPrint less verbose output." << endl;

  cout << "\nOptions for storing intermediate values:" << endl;
  cout << "\t-c FILE,\t\tStore coefficients in FILE." << endl;
  cout << "\t-s FILE,\t\tStore columns in FILE." << endl;

}


bool isMember(SCIP* scip, SCIP_VAR* givenVar, SCIP_CONS* givenCons){
  bool ret = false;
  
  SCIP_CONS ** transcons = (SCIP_CONS**) malloc (sizeof(SCIP_CONS*));
  SCIP_CALL (SCIPgetTransformedCons (scip,
				     givenCons,
				     transcons 
				      ) );
  SCIP_VAR** vars;
  int varssize;
  SCIP_Bool success;

  SCIP_CALL( SCIPgetConsNVars ( scip, *transcons, &varssize, &success));
  if (!success) {
      cout << "A problem happened when trying to get the number of variables in a constraint. " << endl;
  } else {
    vars = (SCIP_Var**) malloc (varssize * sizeof(SCIP_VAR*));
  
    SCIP_CALL( SCIPgetConsVars(scip, *transcons, vars, varssize, &success));
    if (success){
      for (int counter = 0 ; counter < varssize; counter++)
        if (SCIPvarGetIndex(vars[counter]) == SCIPvarGetIndex(givenVar)) { //TODO fix this
          ret = true;
          break;
        }
    } else {
      cout << "A problem happened when trying to check the variable membership in a constraint. " << endl;
    }
    free(vars);
  }
  free(transcons);
  
  return ret;
}

void print_config(Stabilisation stab, bool findAll, unsigned int nclusters, unsigned int ndatapoints, unsigned int nconstraints, unsigned int firstCoefStrategy, char* epsilon_str)
{
  cout << endl;
  cout << "Constrained Clustering using Column Generation (CCCG) [GitHash: " << CCCG_GITHASH << "]" << endl;
  cout << endl << "parameters:" << endl;
  switch(stab) {
      case ST_Both:
          cout << "ST_Both\n";
	  cout << "STRATEGY FOR GENERATING FIRST COEFFICIENTS:" << firstCoefStrategy << endl;
          break;
      case ST_Bound:
          cout << "ST_Bound\n";
	  cout << "STRATEGY FOR GENERATING FIRST COEFFICIENTS:" << firstCoefStrategy << endl;
          break;
      case ST_None:
          cout << "ST_None\n";
          break;
      case ST_Smooth:
          cout << "ST_Smooth\n";
          break;
      case ST_Smooth2:
          cout << "ST_Smooth2\n";
          break;
      case ST_SmoothAuto:
          cout << "ST_SmoothAuto\n";
          break;
  }

  if (findAll)
    cout << "FIND_ALL" << endl;
  else 
    cout << "FIND_BEST" << endl;

  if (epsilon_str != NULL)
    cout << "EPSILON:" << epsilon_str << endl;

  cout << "number of clusters:" << nclusters << endl;
  cout << "number of datapoints:" << ndatapoints << endl;
  cout << "number of constraints:" << nconstraints << endl;
  cout << "branch and bound version:" << getBnBVersion() << endl;

}

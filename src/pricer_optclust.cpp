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

#include "pricer_optclust.h"
#include "SiegBnBCons.h"
#include "SiegBnB.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "scip/cons_linear.h"
#include "scip/pub_var.h"

using namespace std;
using namespace scip;
bool map_compare (SCIP* scip, const std::map<const char*, SCIP_Real>& lhs, const std::map<const char*, SCIP_Real>& rhs);

// print a list in ranges, is helper function
void print_subrange(int min, int max, ostream& out) {
    if (min == max)
        out << min << " ";
    else if (min+1 == max)
        out << min << " " << max << " ";
    else
        out << min << ".." << max << " ";
}
void print_range(const vector<int>& vals, ostream& out) {
    if (vals.size() == 0)
        return;
    int min = vals[0];
    int max = min;
    for (unsigned i=1; i!=vals.size(); i++) {
        if (max+1 == vals[i])
            max++;
        else {
            print_subrange(min, max, out);
            min = vals[i];
            max = vals[i];
        }
    }
    print_subrange(min, max, out);
}
void print_cluster(std::vector< bool > cluster, const char* msg) {
     cout << msg;
     vector<int> vals;
     for (unsigned i=0; i!=cluster.size(); i++) {
       if (cluster[i])
         vals.push_back(i);
     }
     print_range(vals, cout);
     cout << endl;
}

/** Constructs the pricer object with the data needed */
ObjPricerOptClust::ObjPricerOptClust(
   SCIP*                                scip,          /**< SCIP pointer */
   const char*                          p_name,        /**< name of pricer */
   const int                            p_num_instances,   /**< number of instances */
   const int                            p_num_clusters, /**< number of clusters in clustering problem*/
   const vector< vector<double> >&       p_distance,    /**< matrix of distances */
   const vector<SCIP_CONS*>&            p_setcover_con,     /**< matrix of arc constraints */
   SCIP_CONS*                           p_side_con,     /**< array of partitioning constraints */
   Stabilisation                        stab,           /**type of stabilisation to apply */
   bool                                 findAll,
   const char*                          coefFileName,
   const char*                          columnFileName
   ):
   ObjPricer(scip, p_name, "Finds the cluster with the smallest reduced cost", 0, TRUE),
   _num_instances(p_num_instances),
   _num_clusters(p_num_clusters),
   _distance(p_distance),
   _setcover_con(p_setcover_con),
   _side_con(p_side_con),
   _coefFileName(coefFileName),
   _columnFileName(columnFileName),
   _stab(stab),
   _findAll (findAll),
   _mu(0),
   _epsilon(0),
   _updateCount (0),
   _prev_colVals(),
   _prev_out_score(-SCIPinfinity(scip)),
   _piHat_score(-SCIPinfinity(scip)),
   _piHat_redcost(-SCIPinfinity(scip)),
   _piHat_lambdas(),
   _piHat_sigma(0),
   _count_mispricings(1) // should be always 1, is used as divisor
{
  this -> _lower_bound = -SCIPinfinity (scip);
  this -> _upper_bound = SCIPinfinity (scip);
  this -> _duality_gap = SCIPinfinity (scip);
  this -> _bestLambSum = SCIPinfinity (scip);
  
  if (_coefFileName != NULL)
    _outStream.open(_coefFileName , ios::out); 

  if (_columnFileName != NULL)
    _colStream.open(_columnFileName, ios::out);
  
  if ((this->_stab == ST_Smooth || this->_stab == ST_Smooth2)) {
      if (SCIPisZero(scip, _epsilon))
          _epsilon = 0.4; // default alpha value
  } else if (this->_stab == ST_SmoothAuto) {
      _epsilon = 0.8; // auto changed
  }
      
}


/** Destructs the pricer object. */
ObjPricerOptClust::~ObjPricerOptClust()
{
   _outStream.close();
   _colStream.close();
}


/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerOptClust::scip_init)
{
   for (int i = 0; i < num_instances(); ++i)
         SCIP_CALL( SCIPgetTransformedCons(scip, _setcover_con[i], &_setcover_con[i]) );

   SCIP_CALL( SCIPgetTransformedCons(scip, _side_con, &_side_con) );

   return SCIP_OKAY;
}


SCIP_RETCODE ObjPricerOptClust::farkasPricing(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
  /* Are we going to check the degeneracy checking also in Farkas pricing? */

  
  /* allocate array for reduced costs */
  vector<SCIP_Real>  lambdas(num_instances() , 0.0);
  SCIP_Real sigma;

  cout << "FARKAS" << endl;
  for (int i = 0; i < num_instances(); ++i)
    {
      assert( setcover_con(i) != 0 );
      lambdas[i] = SCIPgetDualfarkasLinear(scip, setcover_con(i));
    }
  sigma = SCIPgetDualfarkasLinear(scip, side_con());

  /* find the cluster with the minimum reduced cost*/
  (void) find_add_farkas_cluster(scip, lambdas, sigma);
 
  return SCIP_OKAY;
}


/** perform pricing*/

SCIP_RETCODE ObjPricerOptClust::pricing(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   int verbose = 0;

   if (verbose >= 1) {
      int nVars = SCIPgetNVars(scip);
      cout << "Best Sol so far ("<<nVars<<" vars):\n";

      if (verbose >= 6) {
         SCIP_VAR ** vars = SCIPgetVars(scip);
         for (int i=0; i!=nVars; i++) 
	   SCIPprintVar(scip, vars[i], NULL);

         // print best real sol
         SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
      }
   }
   
   // save values of columns
   map<const char*, SCIP_Real> colVals;
   bool was_degenerate = true;
   {  
         int nCols;
         SCIP_COL** cols;
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &nCols) );
         for (int i=0; i!=nCols; i++) {
             SCIP_VAR* colVar = SCIPcolGetVar(cols[i]);
             SCIP_Real val = SCIPcolGetPrimsol(cols[i]);
             if (SCIPvarIsBinary(colVar) and not SCIPisZero(scip, val)) {
                // this ignores the continuous (perturbation) columns
                colVals[SCIPvarGetName(colVar)] = val;
             }
         }

         if (!map_compare(scip, _prev_colVals, colVals)) {
             _prev_colVals = colVals;
             was_degenerate = false;
         }
   }
   
   if (verbose >= 4) {
        if (was_degenerate)
            cout << "YES degen; ";
        else
            cout << "not degen; ";
        cout << "non-zero columns ("<<colVals.size()<<") -- ";
        if (verbose >= 6) {
            for (map<const char*, SCIP_Real>::iterator ii=colVals.begin(); ii!=colVals.end(); ++ii) {
                cout << (*ii).first << "(" << (*ii).second << ")\t";
            }
        }
        cout << endl;  
   }

   /* allocate array for reduced costs */
   vector<SCIP_Real>  lambdas(num_instances() , 0.0);
   SCIP_Real sigma;
   

   /* store the reduced costs for setcover constraints */
   for (int i = 0; i < num_instances(); ++i)
   {
       assert( setcover_con(i) != 0 );
       lambdas[i] = SCIPgetDualsolLinear(scip, setcover_con(i));
       if (SCIPisZero(scip, lambdas[i]))
           lambdas[i] = 0;
   }
   sigma = SCIPgetDualsolLinear(scip, side_con());
   if (SCIPisZero(scip, sigma))
       sigma = 0;

   /* When smoothing, try to cut the speration point */       
   // out point (in point is _piHat_*
   SCIP_Real out_score = getCurSol(scip);
   vector<SCIP_Real> out_lambdas = lambdas;
   SCIP_Real out_sigma = sigma;
   // sep point
   SCIP_Real sep_score = -SCIPinfinity(scip);
   vector<SCIP_Real> sep_lambdas(lambdas.size());
   SCIP_Real sep_sigma;
   if ((this->_stab == ST_Smooth || this->_stab == ST_Smooth2 || this->_stab == ST_SmoothAuto)) {
       // safeguard (unfortunately it happens)
       if (SCIPisGT(scip, _piHat_score, out_score)) {
           _piHat_score = -SCIPinfinity(scip);
           _piHat_redcost  = -SCIPinfinity(scip);
           if (this->_stab == ST_SmoothAuto)
               _epsilon = 0.8; // try to reset eps as well
       }
       
       _prev_out_score = out_score;
       
       // smooth the duals
       SCIP_Real alpha = _epsilon/_count_mispricings; // mispricings counted lower
       if (!SCIPisInfinity(scip, -_piHat_score)) {
       
           if (!SCIPisEQ(scip, out_score, _piHat_score)) {
               assert(_piHat_lambdas.size() == out_lambdas.size());
           
               sep_score = alpha*_piHat_score + (1-alpha)*out_score;
               sep_sigma = alpha*_piHat_sigma + (1-alpha)*out_sigma;
               for (unsigned i=0; i!=out_lambdas.size(); i++)
                   sep_lambdas[i] = alpha*_piHat_lambdas[i] + (1-alpha)*out_lambdas[i];
               // cut sep point
               sigma = sep_sigma;
               lambdas = sep_lambdas;
           
           }
       }
           
   }

   /* find the cluster with the minimum reduced cost*/
   vector<bool> cluster;
   SCIP_Real reduced_cost;
   reduced_cost = find_add_clusters(scip, lambdas, sigma);

     /* update duality gap (only used when stabilizing, but interesting to print anyways)*/
     update_duality_gap(scip, out_lambdas, out_sigma, reduced_cost);
     if (verbose >= 2)
       cout << "Duality gap: "<<this->_duality_gap<< "\t(" << this->_lower_bound << " -- " << this->_upper_bound << ")\n";
     
     /* update the bounds */
     if (this->_stab == ST_Bound || this->_stab == ST_Both){
         if (verbose >= 5)
            cout << "Stabilize (not farkas)\n";

         if (SCIPisLT(scip, (this -> _duality_gap), (this -> _epsilon))) {
	     SCIP_CALL(updatePertVarsBounds(scip));
	     if (this->_stab == ST_Both)
	       SCIP_CALL(updatePertVarsCoefs(scip, lambdas));
         }
     }
     
     /* when smoothing, update in/out point */
     if (SCIPisFeasNegative(scip, reduced_cost) &&
         (this->_stab == ST_Smooth || this->_stab == ST_Smooth2 || this->_stab == ST_SmoothAuto)) {
         SCIP_Real sep_observed_score = this->bigL(out_lambdas,out_sigma,reduced_cost);
         
         if (!SCIPisLT(scip, sep_observed_score, sep_score) || _stab == ST_Smooth2) {         

             if (this->_stab == ST_SmoothAuto && !SCIPisInfinity(scip, -_piHat_score)) {
                 _count_mispricings++;
             }

             
         }
         
             if (SCIPisGT(scip, reduced_cost, _piHat_redcost)) {
                 _piHat_score = sep_observed_score;
                 _piHat_redcost = reduced_cost;
                 _piHat_lambdas = out_lambdas;
                 _piHat_sigma = out_sigma;
             }
     }

   return SCIP_OKAY;
}



/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - find the best cluster with respect to these dual values
 *  - if this cluster has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
 */
SCIP_DECL_PRICERREDCOST(ObjPricerOptClust::scip_redcost)
{
   SCIPdebugMessage("call scip_redcost ...\n");

   /* set result pointer, see above */
   *result = SCIP_SUCCESS;

   /* call pricing routine */
   SCIP_CALL( pricing(scip) );
   
   return SCIP_OKAY;
}


/** Pricing of additional variables if LP is infeasible.
 *
 *  - get the values of the dual Farks multipliers you need
 *  - find the best cluster with respect to these values
 *  - if this tour has negative reduced cost, add it to the LP
 */
SCIP_DECL_PRICERFARKAS(ObjPricerOptClust::scip_farkas)
{
   SCIPdebugMessage("call scip_farkas ...\n");

   /* call pricing routine */
   SCIP_CALL( farkasPricing(scip) );
   

   return SCIP_OKAY;
}


/** add cluster variable to problem */
SCIP_RETCODE ObjPricerOptClust::add_cluster_variable(
   SCIP*                 scip,               /**< SCIP data structure */
   const vector<bool>&      cluster,        /**< the best cluster found by subproblem */
   const vector<SCIP_Real>& lambdas, 
   const SCIP_Real sigma
   )
{
   assert(cluster.size() == (unsigned)num_instances());
   bool verbose = false;

   char var_name[255];
   SCIPsnprintf(var_name, 255, "X_%d", SCIPgetNVars(scip));
   SCIPdebugMessage("new variable <%s>\n", var_name);
   if (verbose)
    cout << "creating new variable " << var_name << endl;

   if (verbose)
     print_cluster(cluster, "\tNew column: ");

   SCIP_Real cluster_cost = ObjPricerOptClust::getClusterCost(this->_distance, cluster);
   
   /* create the new variable */
   SCIP_VAR* var;
   SCIP_CALL( SCIPcreateVar(scip, &var, var_name,
                            0.0,                     // lower bound
                            1.0,      // upper bound
                            cluster_cost,                       // objective
                            SCIP_VARTYPE_BINARY, // variable type
                            false, false, 0, 0, 0, 0, 0) );

   /* add new variable to the list of variables to price into LP */
   SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );

   /* add the coefficients of this variable to setcovering and side constraints */
   for (int i= 0; i < num_instances(); ++i)
      if (cluster[i])
         SCIP_CALL( SCIPaddCoefLinear(scip, setcover_con(i), var, 1.0) );         

   SCIP_CALL( SCIPaddCoefLinear(scip, side_con(), var, 1.0 ) );

   /* cleanup */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** find and add a cluster when in farkas
 *
 *  uses straightforward calculation
 */
SCIP_Real ObjPricerOptClust::find_add_farkas_cluster(
   SCIP* scip, 
   const  vector<SCIP_Real> & lambdas,   /**< matrix of lengths */
   const  SCIP_Real& sigma
   )
{
   int verbose = 0;
   bool found = false;
   SCIP_Real sol = 0;
   vector<bool> cluster;
   SCIP_Real reduced_cost = SCIPinfinity (scip);
  
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "MLCL");
   int nconss = SCIPconshdlrGetNActiveConss(conshdlr);

   /* If there are no constraints, add all data points with a positive lambda */
   if (nconss == 0)
     {
       for (vector<SCIP_Real>::const_iterator itr = lambdas.begin(); itr != lambdas.end() ; itr++)
	 {
	   if (SCIPisPositive(scip,*itr)){
	     found = true;
	     cluster.push_back(true);
	     sol -= *itr;
	   }
	   else
	     cluster.push_back(false);
	 }
     }

   /* If there are constraints, use constrained B&B subsolver with zero distances */
   else{

     vector<unsigned int> mlBlocks;
     vector<pair <unsigned int, unsigned int> > clBlocks;
     getLinkConstriants(scip, mlBlocks, clBlocks);

     if (verbose >= 3){
       //printing constraints
       SCIP_CONSDATA* consdata;
       SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr);
      
       for (int i=0; i!=nconss; i++) {
	 consdata = SCIPconsGetData(conss[i]);
	 assert(consdata != NULL);
	 unsigned int first = SCIPgetMlClFirst(consdata);
	 unsigned int second = SCIPgetMlClSecond(consdata);
	 bool isML = SCIPgetIsML(consdata);
	 cout << (isML ? "ML" : "CL") << " : " << first << "-" << second << endl;
       }

     }

     vector<vector<double> > zeroDistance(_num_instances, vector<double> (_num_instances, 0));
     SiegBnBCons solver(mlBlocks, clBlocks, zeroDistance , lambdas , sigma, lambdas.size()); 
     solver.solve();
     if (!solver.hasSolution())
       solver.solvePositive();

     if (solver.hasSolution()){
       sol = solver.getSolution(cluster);
       found = true;
     }
     
     if (verbose >= 2){
       if (found) {
	 //printing the cluster found by farkas
	 for (auto itr = cluster.begin(), endItr = cluster.end(); itr != endItr; itr++)
	   cout  << (*itr ? "*" : "-");
	 cout << endl;
       }
       else
	 cout << "Farkas could not find an initial solution" << endl;
     }
   }
   
   if (found){
     /* verifying obtained cluster */
     if (!ConshdlrMLCL::clusterSatisfied(scip, cluster)) {
       cout << "ERROR\n";
       cout << "The above cluster does not satisfy the following set of constraints:\n";
       ConshdlrMLCL::print_constraints(scip);
       exit(1);
     }

     reduced_cost = sol - sigma;
   }   

   /* add cluster variable */
   if ( found && SCIPisFeasNegative(scip, reduced_cost) )
   {
     cout << "reduced cost: " << reduced_cost << " sigma: " << sigma  << " sol: " << sol << endl;
     (void) add_cluster_variable(scip, cluster, lambdas, sigma);
   }
   
   else 
     cout << "reduced cost not negative!" << endl;
   
   return reduced_cost;
}

/** find and add clusters with negative reduced cost
 *  uses a subsolver
 */
//TODO check the effect of returning one reduced cost for a set of added columns
SCIP_Real ObjPricerOptClust::find_add_clusters(
   SCIP* scip, 
   const  vector<SCIP_Real> & lambdas,   /**< matrix of lengths */
   const  SCIP_Real& sigma)
{
   // local verbosity
   int verbose = 2;
   bool only_best = !_findAll;

   if (_coefFileName != NULL)
       saveCoef(lambdas, sigma);

   //saveLp(scip, true);

   if (verbose >= 1) {
      int nonzero = 0;
      for (vector<SCIP_Real>::const_iterator itr = lambdas.begin();
            itr != lambdas.end(); itr++){
	if (SCIPisPositive(scip, *itr))
            nonzero++;
      }

      cout << "Sigma: " << sigma;
      cout << " Lambda's ("<<lambdas.size()<<", "<<nonzero<<" non-zero): ";
      vector<int> vals;
      for (unsigned i=0; i != lambdas.size(); i++) {
	if (!SCIPisZero(scip, lambdas[i])) {
            vals.push_back(i);
	    if (verbose >= 3)
	      cout << i << "(" << lambdas[i] << ") ";
	 }
      }
      if (verbose >= 3)
         cout << endl;
      print_range(vals, cout);
      cout << endl;
   }

   boost::posix_time::ptime time1,time2; 
   time1 = boost::posix_time::microsec_clock::local_time();

   if (verbose >= 1)
     cout << endl << "Solving subproblem ..." << flush;

   // initialise solver (with ml/cl blocks from the constraint handler)
   vector<unsigned int> mlblocks;
   vector<pair<unsigned int, unsigned int> > clblocks;
   getLinkConstriants(scip, mlblocks, clblocks);
   SiegBnBCons solver(mlblocks, clblocks, _distance , lambdas , sigma, lambdas.size()); 

   solver.solve();
 
   if (verbose >= 1) {
       time2 = boost::posix_time::microsec_clock::local_time();
       boost::posix_time::time_duration diff = time2 - time1; 
       cout << " runtime: " << boost::posix_time::to_simple_string(diff) << endl;
   }
   SCIP_Real reduced_cost = SCIPinfinity(scip);

   if (only_best) {
      /* get the cluster with the minimum reduced cost*/
      vector<bool> cluster;
      SCIP_Real sol;
      if (!solver.hasSolution()) {
	  if (verbose >= 2)
	    cout << "No solution with negative reduced cost. Now searching for minimum solution with positive reduced cost. " << endl;

	solver.solvePositive();
      }

      sol = solver.getSolution(cluster);
      // sigma was constant in the optimisation and was not yet substracted
      reduced_cost = sol - sigma;
      
     /* verifying obtained cluster */
     if (!ConshdlrMLCL::clusterSatisfied(scip, cluster)) {
       cout << "ERROR\n";
       cout << "The above cluster does not satisfy the following set of constraints:\n";
       ConshdlrMLCL::print_constraints(scip);
       exit(1);
     }

      /* add cluster variable */
      if ( SCIPisFeasNegative(scip, reduced_cost) || !_pertVars.empty() ) {
	cout << "reduced cost: " << reduced_cost << " sigma: " << sigma  << " sol: " << sol << endl;
	(void) add_cluster_variable(scip, cluster, lambdas, sigma);
      } else 
         cout << "reduced cost not negative!" << endl;

   } else {
     /* get and add ALL clusters (if neg. reduced cost) */
     vector< vector<bool> > clusters;
     
     vector< SCIP_Real > quals = solver.getSolutions(clusters);
     assert(quals.size() == clusters.size());
     
     int nrAdded = 0;
     for (unsigned i=0; i!=quals.size(); i++) {
       SCIP_Real redcost = quals[i] - sigma;
       if (SCIPisLT(scip, redcost, reduced_cost))
         reduced_cost = redcost;
       if ( SCIPisFeasNegative(scip, redcost) || !_pertVars.empty()) {
         (void) add_cluster_variable(scip, clusters[i], lambdas, sigma);
	 nrAdded++;
	 if (verbose >= 3)
	    cout << "Sol["<<i<<"] reduced cost: " << redcost << " sigma: " << sigma  << " qual: " << quals[i] << endl;
       } else {
	 if (verbose >= 3)
            cout << "Sol["<<i<<"] reduced cost not negative! (" << redcost << ", sigma: " << sigma  << " qual: " << quals[i] << ")" << endl;
	 if (i+1 == quals.size())
	   cout << "No negative reduced cost clusters found.\n";
       }
     }
     if (verbose >= 1)
	cout << "Added " << nrAdded << " columns with smalles neg. red. cost of " << reduced_cost << " (+sigma="<<reduced_cost+sigma<<")" << endl;
   }
      
   return reduced_cost;
}

// recalculate the duality gap
void ObjPricerOptClust::update_duality_gap(
      SCIP* scip, 
      const vector<SCIP_Real> & lambdas,     /**< matrix of lengths */
      const SCIP_Real& sigma,
      const SCIP_Real& reduced_cost
      )
{
	 int verbose = 0;

	 SCIP_Real lambSum = 0;
	 for (vector<SCIP_Real>::const_iterator itr = lambdas.begin(), endItr = lambdas.end(); itr != endItr; itr++)
	   lambSum += *itr;

	 SCIP_Real ub = (this -> _num_clusters) * sigma + lambSum;
         if (verbose >= 2)
           cout << "ub("<<ub<<") = num_clus("<<this -> _num_clusters<<") * sigma("<<sigma<<") + lambSum("<<lambSum<<") -- was: "<<this -> _upper_bound << endl;
	 if (SCIPisLT(scip, ub, (this -> _upper_bound))) {
	   this -> _upper_bound = ub;
	 }
     
	 SCIP_Real lb = ub + (this -> _num_clusters) * reduced_cost;
         if (verbose >= 2)
           cout << "lb("<<lb<<") = ub("<<ub<<") + num_clus("<<this -> _num_clusters<<") * redCost("<<reduced_cost<<") -- was: "<<this -> _lower_bound << endl;
	 if (SCIPisGT(scip, lb, (this -> _lower_bound))){
	   this -> _lower_bound = lb;
	 }
     
	 if (!SCIPisInfinity(scip, this -> _upper_bound) && 
	     !SCIPisInfinity(scip, - this -> _lower_bound)){
	   SCIP_Real denom = ((this -> _upper_bound) > 1 ? (this -> _upper_bound) : 1);
	   this -> _duality_gap = ((this -> _upper_bound) - (this -> _lower_bound)) / denom;
           if (verbose >= 2)
               cout << "gap("<<this -> _duality_gap<<") = ( ub("<<this -> _upper_bound<<") - lb("<<this -> _lower_bound<<") ) / denom("<<denom<<")"<<endl;
	 }
         // for plotting UB/LB evolution:
         //cout << "oub;olb;"<<ub<<";"<<lb<<"\n";
}


void ObjPricerOptClust::saveCoef(
   const  vector<SCIP_Real> & lambdas,   /**< matrix of lengths */
   const  SCIP_Real& sigma
   )
{
  _outStream << sigma;

   for (vector<SCIP_Real>::const_iterator itr = lambdas.begin();
        itr != lambdas.end(); itr++){
     _outStream << "\t" << *itr;
   }
   _outStream << endl;

}


void ObjPricerOptClust::saveLp(
   SCIP* scip, bool onlyBasic){
  //   _outStream << endl;

   int nrows = SCIPgetNLPRows(scip);
   int ncols = SCIPgetNLPCols(scip);

   SCIP_COL ** cols = SCIPgetLPCols(scip);
   SCIP_ROW ** rows = SCIPgetLPRows(scip);
   
   for (int counter = 0; counter < ncols; counter ++){

      //SCIP_BASESTAT stat = (SCIP_BASESTAT) SCIPcolGetBasisStatus (cols[counter]);
      bool isBasic = (SCIPcolGetPrimsol(cols[counter]) > 0);
      
      if (!onlyBasic || isBasic){
         
      int nNonZero = SCIPcolGetNLPNonz(cols[counter]);
      SCIP_ROW ** currentColRows = SCIPcolGetRows (cols[counter]);
      
      int index = 0;
      for (int counter2 = 0; counter2 < nNonZero; counter2++){
         SCIP_ROW * currentRow = currentColRows[counter2]; 
        int position = SCIProwGetLPPos(currentRow);
        // cout << position << ";";

         while (index < position ){
            cout << "-\t";
            index++;
         }

         cout << position + 1<< "\t";
         index++;
      }

      while (index < nrows - 1){
         cout << "-\t";
         index ++;
      }
      
      cout << "\tObjVal:" << SCIPcolGetObj(cols[counter]);
      cout << "\tSol:" << SCIPcolGetPrimsol(cols[counter]);
      cout << endl;
      }
   }
   
   for (int counter = 0; counter < nrows; counter++)
      cout << SCIProwGetDualsol (rows[counter]) << "\t";
   cout << endl;


   //   _outStream << endl;
}

void ObjPricerOptClust::setPertVars(std::vector<SCIP_VAR*> pertVars, SCIP_Real mu, SCIP_Real epsilon)
{
    this -> _pertVars = pertVars;
    this -> _mu = mu;
    this -> _epsilon = epsilon;
}
void ObjPricerOptClust::setEpsilon(SCIP_Real epsilon)
{
    this -> _epsilon = epsilon;
}


SCIP_RETCODE ObjPricerOptClust::updatePertVarsBounds (SCIP* scip){
  _mu /= pow (2, ++ _updateCount);

  if (SCIPisFeasPositive(scip, _mu))
    cout << "Updating bounds of perturbed vars (gap="<<this -> _duality_gap<<" < "<<this -> _epsilon<<"=eps) -- new bound: " << _mu << "\n";


  SCIP_Bool infeasible = false, fixed = true;
    for (std::vector<SCIP_VAR *>::iterator itr = (this -> _pertVars).begin(), endItr = (this -> _pertVars).end();itr != endItr; itr++){
      /* If the bound is not too tight, shrink it */
      if (SCIPisFeasPositive(scip, _mu)){
	SCIP_CALL (SCIPchgVarLb(scip, *itr, -_mu));
	SCIP_CALL (SCIPchgVarUb(scip, *itr, _mu));
      }
      /* otherwise, fix the variables to zero  */
      else {
        SCIP_CALL (SCIPfixVar (scip,*itr,0,&infeasible,&fixed ) );
	_pertVars.clear();
      }
    }

    this -> _lower_bound = -SCIPinfinity (scip);
    this -> _upper_bound = SCIPinfinity (scip);
    this -> _duality_gap = SCIPinfinity (scip);

    return SCIP_OKAY;
}

SCIP_RETCODE ObjPricerOptClust::updatePertVarsCoefs (SCIP * scip, std::vector<SCIP_Real> lambdas)
{
  /* Verify that the bound is not set to zero */
  if (!SCIPisFeasPositive(scip, _mu)) {
    cout << "updatePerVarsCoefs: _mu("<<_mu<<") set to zero, nothing left to do...\n";
    return SCIP_OKAY;
  }

  // TODO!! Check that updated coefficient is not the same as previous
  // if so: no need to fix/recreate (this happens in experiments)
  
    int verbose=1;

    SCIP_Bool infeasible = false, fixed = true;
    // remove all pertvars for which the lambda changed (and set them to NULL to re-fill later)
    for (unsigned i=0; i!= _pertVars.size(); i++) {
      SCIP_VAR* var = _pertVars[i];
    
      if (!SCIPisEQ(scip, lambdas[i], SCIPvarGetObj(var))) {
          if (verbose >= 3)
            cout << "var " << i << "new lambda: " << lambdas[i] << " old one: " << SCIPvarGetObj(var) << "\n";
          
          SCIP_CALL (SCIPfixVar (scip,var,0,&infeasible,&fixed ) );
          _pertVars[i] = NULL;
          if (!fixed)
            return SCIP_ERROR;
      }
    }

    if (verbose >= 3)
      std::cout << "Continueing (num inst) " << _num_instances << "\n";
    for (int counter = 0; counter < _num_instances; counter++){
      if (_pertVars[counter] != NULL)
        continue;

      if (verbose >= 3)
	std::cout << "Creating a perturb var ("<<counter<<")\n";

      SCIP_VAR* var;
      char pert_var_name [255];
      SCIPsnprintf(pert_var_name, 255, "PV%d", counter);
      SCIP_CALL( SCIPcreateVar(scip, &var, pert_var_name,
			       -_mu,                     // lower bound
			       _mu,      // upper bound
			       lambdas[counter],                       // objective
			       SCIP_VARTYPE_CONTINUOUS, // variable type
			       true, //initial 
			       true, //removable
			       0, 0, 0, 0, 0) );

      /* add new variable to the list of variables to price into LP */
      SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );

      /* add the coefficients of this variable to setcovering and side constraints */
      SCIP_CALL( SCIPaddCoefLinear(scip, setcover_con(counter), var, 1.0) );         

      /* add to the vector of purturbation variables */
      _pertVars[counter] = var;

      /* cleanup */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
    }
    
    this -> _lower_bound = -SCIPinfinity (scip);
    this -> _upper_bound = SCIPinfinity (scip);
    this -> _duality_gap = SCIPinfinity (scip);

    return SCIP_OKAY;
}

SCIP_Real ObjPricerOptClust::getCurSol(SCIP* scip)
{
    SCIP_Real ret = 0;

    int nVars = SCIPgetNVars(scip);

    SCIP_VAR ** vars = SCIPgetVars(scip);
    for (int i=0; i!=nVars; i++) {
        ret += (SCIPgetVarSol(scip,vars[i]) * SCIPvarGetObj(vars[i]));
    }

    return ret;
}


// The bound at each iteration would still be calculated by lambSum + sigmaK + (reduced cost computed for that iteration)
SCIP_Real ObjPricerOptClust::bigL(
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma,
      const SCIP_Real& reduced_cost
                             )
{
    SCIP_Real lambSum = 0;
    for (vector<SCIP_Real>::const_iterator itr = lambdas.begin(), endItr = lambdas.end(); itr != endItr; itr++)
        lambSum += *itr;
    
    SCIP_Real sigmaK = sigma * this->_num_clusters;
    
    SCIP_Real ret = lambSum + sigmaK + reduced_cost;
    
    return ret;
}
SCIP_Real ObjPricerOptClust::bigL(
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma,
      const vector<bool>& cluster
                             )
{
    assert(lambdas.size() == cluster.size());
    
    SCIP_Real lambSum = 0;
    for (vector<SCIP_Real>::const_iterator itr = lambdas.begin(), endItr = lambdas.end(); itr != endItr; itr++)
        lambSum += *itr;
    
    SCIP_Real sigmaK = sigma * this->_num_clusters;
    
    SCIP_Real cx = ObjPricerOptClust::getClusterCost(this->_distance, cluster);
    
    SCIP_Real piAx = sigma;
    for (unsigned i=0; i!=cluster.size(); i++) {
        if (cluster[i])
            piAx += lambdas[i];
    } 
    
    SCIP_Real ret = lambSum + sigmaK + cx - piAx;
    return ret;
}

SCIP_Real ObjPricerOptClust::bigL_t(SCIP* scip, 
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma
                             )
{
         SCIP_Real ret_min = SCIPinfinity(scip);
         
         // calculate for each column
         int nCols;
         SCIP_COL** cols;
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &nCols) );
         for (int i=0; i!=nCols; i++) {
             SCIP_VAR* colVar = SCIPcolGetVar(cols[i]);
             if (SCIPvarIsBinary(colVar)) {
                 vector<bool> cluster = getColCluster(scip, colVar);
                 SCIP_Real valL = bigL(lambdas, sigma, cluster);
                 
                 if (SCIPisLT(scip, valL, ret_min)) {
                     ret_min = valL;
                 }
             }
         }

         return ret_min;
}

std::vector< bool > ObjPricerOptClust::getColCluster(SCIP* scip, SCIP_VAR* colVar)
{
    SCIP_VAR** vars = (SCIP_Var**) malloc (SCIPgetNVars(scip) * sizeof(SCIP_VAR*));
    
    vector<bool> cluster(num_instances());
    int varssize = SCIPgetNVars(scip);
    SCIP_Bool success;
    
    for (int i = 0; i < num_instances(); ++i) {
         SCIPgetConsNVars ( scip, _setcover_con[i], &varssize, &success);
         
         bool inCluster = false;
         if (success) {
             SCIPgetConsVars(scip, _setcover_con[i], vars, varssize, &success);
             if (success) {
                 for (int j=0; j!=varssize; j++) {
                     if (SCIPvarGetIndex(vars[j]) == SCIPvarGetIndex(colVar)) {
                         inCluster = true;
                         break;
                     }
                 }
             }
         }
         cluster[i] = inCluster;
    }
    
    free(vars);
    return cluster;
}


SCIP_Real ObjPricerOptClust::getClusterCost(const vector< std::vector< double > >& distance, const std::vector< bool >& cluster)
{
   /* compute the cost of a cluster */
   SCIP_Real cluster_cost = 0;
      
   int cluster_size = 0;
   /* calculate the size of cluster */
   for (unsigned i = 0 ; i < cluster.size(); i++)
      if (cluster[i])
        cluster_size++;
   
   for (unsigned i = 0 ; i < cluster.size(); i++)
      for (unsigned j = i + 1; j < cluster.size(); j++)
         if (cluster[i] && cluster[j])
            cluster_cost += distance[i][j];

   if (cluster_size > 0)
     cluster_cost /= cluster_size;
   
   return cluster_cost;
}


struct equalityChecker{
  SCIP* scip;
  bool operator() (pair<const char*, SCIP_Real> const &lhs, pair<const char*, SCIP_Real> const &rhs) {
    return (lhs.first == rhs.first) && 
      SCIPisEQ(scip, lhs.second , rhs.second);
      }
};


bool map_compare (SCIP* scip, const std::map<const char*, SCIP_Real>& lhs, const std::map<const char*, SCIP_Real>& rhs) {
  equalityChecker eqCheck;
  eqCheck.scip = scip;
  return lhs.size() == rhs.size()
    && std::equal(lhs.begin(), lhs.end(),
		  rhs.begin(), eqCheck);
}


void ObjPricerOptClust::getLinkConstriants(SCIP* scip, vector<unsigned int>& mlBlocks, vector<pair<unsigned int, unsigned int> >& clBlocks){
  clBlocks.clear();

  vector<int> tempMBs (_num_instances, -1);
  int blockNum = 0;

  SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "MLCL");
    
  int nconss = SCIPconshdlrGetNActiveConss(conshdlr); // only these first ones are active
  SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr);
  assert(nconss == 0 || conss != NULL);
  
  SCIP_CONSDATA* consdata;
  for (int i=0; i!=nconss; i++) {
    consdata = SCIPconsGetData(conss[i]);
    assert(consdata != NULL);
      
    unsigned int first = SCIPgetMlClFirst(consdata);
    unsigned int second = SCIPgetMlClSecond(consdata);
      
    if (!SCIPgetIsML(consdata))
      clBlocks.push_back(make_pair(first, second));
    else{
      int firstInd = tempMBs[first] , secondInd = tempMBs[second];
      if (firstInd == -1 && secondInd == -1){
	tempMBs[first] = blockNum;
	tempMBs[second] = blockNum;
	blockNum++;
      }
      else if (firstInd > -1 && secondInd > -1){
	for (auto itr = tempMBs.begin(), endItr = tempMBs.end(); itr != endItr; itr++)
	  if (*itr == firstInd)
	    *itr = secondInd;
      }
      else{
	if (firstInd == -1)
	  tempMBs[first] = secondInd;
	else
	  tempMBs[second] = firstInd;
      }
    }
  }

  mlBlocks.clear();
  for (auto itr = tempMBs.begin(), endItr = tempMBs.end(); itr != endItr; itr++){
    if (*itr == -1){
      *itr = blockNum++;
    }
    mlBlocks.push_back(*itr);
  }
}

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

#include <scip/scip.h> 
#include <scip/scipdefplugins.h>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include <boost/concept_check.hpp>
#include <sys/types.h>
#include <dirent.h>

#include "init.h"
#include "cons_mlcl.h"
#include "pricer_optclust.h" // for print_cluster()
using namespace std;

/** public-like **/

SCIP_RETCODE initialize_solution(const char* fileName, SCIP* scip,
                                 const vector<vector<double> > & distances,
                                 vector<SCIP_CONS*> setcover_con, SCIP_CONS* side_con,
                                 vector<pair<unsigned int, unsigned int> > mlVec, 
				 vector<pair<unsigned int, unsigned int> > clVec,
				 vector<double>& upper,
				 vector<double>& lower
                                ) {
  vector<vector<bool> > clusters;
  vector<double> costs;

  extract_clusters(fileName, distances, clusters, costs, upper, lower);

  assert (clusters.size() == costs.size());
  int size = costs.size();

  for (int counter = 0; counter < size; counter++) {
    if ( cluster_satisfies_mlcl(clusters[counter], mlVec, clVec) ) {
      SCIP_CALL( add_cluster_variable(scip, distances, setcover_con, side_con, clusters[counter],  costs[counter]) );
    }
    else 
      cout << "cluster " << counter << " does not satisfy the constraints." << endl;
  }
  
  return SCIP_OKAY;
}

SCIP_RETCODE init_sol_dir(const char* dirName, SCIP* scip,
			  const vector<vector<double> > & distances,
			  vector<SCIP_CONS*> setcover_con, SCIP_CONS* side_con,
			  vector<pair<unsigned int, unsigned int> > mlVec, 
			  vector<pair<unsigned int, unsigned int> > clVec,
			  vector<double>& upper,
			  vector<double>& lower
			  ){
  
  priority_queue<struct initClustering, vector<struct initClustering>, compQ> clusteringQ;
  DIR * dir = opendir (dirName);
  struct dirent * entries;

  upper.clear();
  lower.clear();

  while ((entries = readdir(dir)) != NULL)
    if (entries->d_type == 8){
      bool success = false;
      struct initClustering currClustering;
      vector<vector<bool> > clusters;
      vector<double> costs;

      vector<double> upper_bounds;
      vector<double> lower_bounds;

      char fileName [255];
      strcpy(fileName, dirName);
      strcat(fileName, "/");
      strcat(fileName, entries->d_name);
      success = extract_clusters(fileName, distances, clusters, costs, upper_bounds, lower_bounds);

      if (success){
	currClustering.totalCost = 0;
	currClustering.satCount = 0;
	for (int counter = 0, clusteringSize = clusters.size(); counter < clusteringSize; counter++)
	  if (cluster_satisfies_mlcl(clusters[counter], mlVec, clVec)){
	    currClustering.clusters.push_back(clusters[counter]);
	    currClustering.costs.push_back(costs[counter]);
	    currClustering.totalCost += costs[counter];
	    (currClustering.satCount)++;
	  }
	currClustering.upper_bounds = upper_bounds;
	currClustering.lower_bounds = lower_bounds;
	clusteringQ.push(currClustering);
      }

    }
  closedir(dir);

  if (!clusteringQ.empty()){
    vector<vector<bool> > bestClusters = (clusteringQ.top()).clusters;
    vector<double> bClusterCosts = (clusteringQ.top()).costs;
    upper = (clusteringQ.top()).upper_bounds;
    lower = (clusteringQ.top()).lower_bounds;
    for (int counter = 0, clusteringSize = bestClusters.size(); counter < clusteringSize; counter++)
      SCIP_CALL( add_cluster_variable(scip, distances, setcover_con, side_con, bestClusters[counter],  bClusterCosts[counter]) );
  }
 
  return SCIP_OKAY;

}


SCIP_RETCODE init_MlCl_constraints(const char* consFile, SCIP* scip, vector<pair<unsigned int, unsigned int> >& mlVec, vector<pair<unsigned int, unsigned int> >& clVec){

   extract_constraints(consFile, mlVec, clVec);
          
   char con_name[255];

   for (auto itr = mlVec.begin(); itr!= mlVec.end(); itr++){
    SCIP_CONS* consML;
    SCIPsnprintf(con_name, 255, "ML-%i-%i", itr->first, itr->second);
    SCIP_CALL( SCIPcreateConsML(scip, &consML, con_name,
                  itr->first,                       /* the first data point to be linked */
                  itr->second                       /* the second data point to be linked */
                  ) );
    SCIP_CALL( SCIPaddCons(scip, consML) );
   }

   for (auto itr = clVec.begin(); itr!= clVec.end(); itr++){
    SCIP_CONS* consCL;
    SCIPsnprintf(con_name, 255, "CL-%i-%i", itr->first, itr->second);
    SCIP_CALL( SCIPcreateConsCL(scip, &consCL, con_name,
                  itr->first,                       /* the first data point to be linked */
                  itr->second                       /* the second data point to be linked */
                  ) );

    SCIP_CALL( SCIPaddCons(scip, consCL) );
   }
   return SCIP_OKAY;
}

/** private-like **/

bool extract_clusters(const char* fileName, const vector<vector<double> >& distances, vector<vector<bool> >& clusters, vector<double>& costs, vector<double>& upper, vector<double>& lower){ 

  clusters.clear();
  ifstream inputFile (fileName, ios::in);

  if (inputFile.fail()) {
    cerr << "Error! Can't open cluster file '" << fileName << "', quitting...\n";
    return false;
  }
  
  int number_of_elements = distances.size();
  
  int line = -1;
  string s;
  while(getline(inputFile,s)) {
    if (!s.empty()) {
      line += 1;
      int clusIndex = atoi(s.c_str());
      
      if (clusIndex < 0) {
          cerr << "Error reading '" << fileName << "', found negative cluster index '" << clusIndex << "'\n";
	  return false;
      }
      
      if (clusIndex > (int)clusters.size()-1) {
          // a new cluster index
          // create empty clusters upto the index+1 (to compensate offset 0)
          clusters.resize(clusIndex+1, vector<bool>(number_of_elements, false));
      }
    
      clusters[clusIndex][line] = true;
    }
  }
  int number_of_clusters = clusters.size();

  #ifndef NDEBUG
  for (int counter = 0; counter < number_of_clusters; counter++)
    assert (clusters[counter].size() == (unsigned)number_of_elements);
  #endif

  /* computing costs of clusters*/
  costs.clear();
  for (int counter = 0; counter < number_of_clusters; counter++)
    costs.push_back(ObjPricerOptClust::getClusterCost(distances, clusters[counter]));

  /* print the sum of all cluster costs (cost of the whole clustering) */
  double clusteringCost = 0;
  for (auto itr = costs.begin(), itrEnd = costs.end(); itr != itrEnd; itr++)
    clusteringCost += *itr;
  cout << "\ncost of the initial clustering: " << setprecision(15) << clusteringCost << endl;

  /* compute an estimate for first coefficients */
  getFirstCoefs (distances, clusters, costs, upper, lower);

return true;

}

SCIP_RETCODE add_cluster_variable(
				  
   SCIP*                 scip,
   const vector<vector<double> >& distances,
   vector<SCIP_CONS*> setcover_con,
   SCIP_CONS* side_con,
   const vector<bool>&      cluster,
   double cluster_cost
   )
{
   char var_name[255];
   SCIPsnprintf(var_name, 255, "X_%d", SCIPgetNVars(scip));
   SCIPdebugMessage("new variable <%s>\n", var_name);

   /* create the new variable */

   SCIP_VAR* var;
   SCIP_CALL( SCIPcreateVar(scip, &var, var_name,
                            0.0,                     // lower bound
                            1.0,      // upper bound
                            cluster_cost,                       // objective
                            SCIP_VARTYPE_BINARY, // variable type
                            false, false, 0, 0, 0, 0, 0) );

   /* add new variable to the list of variables to price into LP */
   SCIP_CALL( SCIPaddVar(scip, var) );

   /* add the coefficients of this variable to setcovering and side constraints */
   for (unsigned i= 0; i < cluster.size(); ++i)
      if (cluster[i])
         SCIP_CALL( SCIPaddCoefLinear(scip, setcover_con[i], var, 1.0) );         

   SCIP_CALL( SCIPaddCoefLinear(scip, side_con, var, 1.0 ) );

   // /* cleanup */
   // SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

void extract_constraints(const char* fileName, vector<pair<unsigned int, unsigned int> >& mlVec, vector<pair< unsigned int, unsigned int> >& clVec){
   mlVec.clear();
   clVec.clear();
   ifstream inputFile (fileName, ios::in);

  if (inputFile.fail()) {
    cerr << "Error! Can't open constraint file '" << fileName << "', quitting...\n";
    exit(-1);
  }

  string line;
  while(getline(inputFile, line)){
     istringstream ss(line);
     unsigned int first, second;
     int consType;
     ss >> first >> second >> consType;
     if (!ss.fail()){
        pair<unsigned int, unsigned int> p = 
           make_pair(first, second);
        if (consType == 1)
           mlVec.push_back(p);
        else if (consType == -1)
           clVec.push_back(p);
     }
  }
}

bool cluster_satisfies_mlcl(std::vector< bool > cluster, vector< std::pair< unsigned int, unsigned int > > mlVec, vector< std::pair< unsigned int, unsigned int > > clVec)
{
    // check ML violation
    for (unsigned i=0; i!=mlVec.size(); i++) {
        unsigned first = mlVec[i].first;
        unsigned second = mlVec[i].second;
        if ( ( cluster[first] == true && cluster[second] == false )
           ||( cluster[first] == false && cluster[second] == true ) )
            return false;
    }
    // check CL violation
    for (unsigned i=0; i!=clVec.size(); i++) {
        if ( cluster[clVec[i].first] == true && cluster[clVec[i].second] == true )
            return false;
    }
    return true;
}

void getFirstCoefs (const vector<vector<double> >& distances, const vector<vector<bool> >& clusters, const vector<double>& costs, vector<double>& upper, vector<double>& lower){

  unsigned int num_instances = distances.size();
  unsigned int num_clusetrs = clusters.size();
  vector<int> cluster_sizes (num_clusetrs , -1);
  vector<int> assignments (num_instances, -1);

  upper.clear();
  upper.resize(num_instances);
  lower.clear();
  lower.resize(num_instances);

  for (size_t i = 0; i < num_clusetrs; i++){
    unsigned int clusterSize = 0;
    for (size_t j = 0; j < num_instances; j++)
      if (clusters[i][j]){
	assignments[j] = i;
	clusterSize++;
      }

    cluster_sizes[i] = clusterSize;
  }

  #ifndef NDEBUG
  for (auto itr = assignments.begin(), iend = assignments.end(); itr!=iend; itr++)
    if (*itr >= num_clusetrs || *itr < 0){
      cout << "something went wrong when assigning instances to clusters" << endl;
      exit (EXIT_FAILURE);
    }

  for (auto itr = cluster_sizes.begin(), iend = cluster_sizes.end(); itr!=iend; itr++)
    if (*itr < 0){
      cout << "cluster with negative size?" << endl;
      exit (EXIT_FAILURE);
    }
  #endif

  for (size_t curr_inst = 0; curr_inst < num_instances; curr_inst++){
    double curr_dist_sum = 0;

    for (size_t other_inst = 0; other_inst < num_instances; other_inst++)
      if (assignments[curr_inst] == assignments[other_inst])
	curr_dist_sum += distances[curr_inst][other_inst];
    
    unsigned int curr_assigned_cluster = assignments[curr_inst];
    double curr_new_cost = (costs[curr_assigned_cluster] * cluster_sizes[curr_assigned_cluster] -
			    curr_dist_sum ) / (cluster_sizes[curr_assigned_cluster] - 1);
    lower[curr_inst] = costs[curr_assigned_cluster] - curr_new_cost;
  }

  for (size_t curr_inst = 0; curr_inst < num_instances; curr_inst++){
    bool firstRound = true;
    double min_diff = 0;

    for (size_t curr_clust = 0; curr_clust < num_clusetrs; curr_clust++)
      if (curr_clust != static_cast <unsigned int> (assignments[curr_inst])){
	double dist_sum = 0;
	for (size_t other_inst = 0; other_inst < num_instances; other_inst++)
	  if (clusters[curr_clust][other_inst])
	    dist_sum += distances[curr_inst][other_inst];

	double curr_clust_diff = ((costs[curr_clust] * cluster_sizes[curr_clust] 
				   + dist_sum) / (cluster_sizes[curr_clust] + 1)) - costs[curr_clust];
	if (firstRound){
	  min_diff = curr_clust_diff;
	  firstRound = false;
	}
	else if (curr_clust_diff < min_diff)
	  min_diff = curr_clust_diff;
      }

    assert (firstRound == false);
    upper[curr_inst] = min_diff;
  }

}

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

#ifndef __MSSCINIT_H__
#define __MSSCINIT_H__

#include <vector>
#include <string>

#include <scip/scip.h>

using namespace std;

struct initClustering {
  vector<vector<bool> > clusters;
  vector<double> costs;
  double totalCost;
  unsigned int satCount;
  vector<double> upper_bounds;
  vector<double> lower_bounds;
};

// public-like
SCIP_RETCODE init_MlCl_constraints(const char* fileName, SCIP*, vector<pair<unsigned int, unsigned int> >&, vector<pair<unsigned int, unsigned int> >&);
SCIP_RETCODE initialize_solution(const char* fileName, SCIP* scip,
                                 const vector<vector<double> > & distances,
                                 vector<SCIP_CONS*> setcover_con, SCIP_CONS* side_con,
                                 vector<pair<unsigned int, unsigned int> > mlVec, 
				 vector<pair<unsigned int, unsigned int> > clVec,
				 vector<double> &upper,
				 vector<double> &lower
                                );

SCIP_RETCODE init_sol_dir(const char* dirName, SCIP* scip,
			  const vector<vector<double> > & distances,
			  vector<SCIP_CONS*> setcover_con, SCIP_CONS* side_con,
			  vector<pair<unsigned int, unsigned int> > mlVec, 
			  vector<pair<unsigned int, unsigned int> > clVec,
			  vector<double>& upper,
			  vector<double>& lower
                                );

// private-like
bool extract_clusters(const char* fileName, const vector< std::vector< double > >& distances, vector< std::vector< bool > >& clusters, std::vector< double >& costs, vector<double>& upper, vector<double>& lower);
SCIP_RETCODE add_cluster_variable(SCIP*, const vector<vector<double> >&, vector<SCIP_CONS*>, SCIP_CONS*, const vector<bool>&, double);
void extract_constraints(const char* fileName, vector< std::pair< unsigned int, unsigned int > >& mlVec, vector< std::pair< unsigned int, unsigned int > >& clVec);
bool cluster_satisfies_mlcl(vector<bool> cluster, vector<pair<unsigned int, unsigned int> > mlVec, vector<pair<unsigned int, unsigned int> > clVec);
void getFirstCoefs (const vector<vector<double> >& distances, const vector<vector<bool> >& clusters, const vector<double>& costs, vector<double>& upper, vector<double>& lower);
class compQ{
 public:
  bool operator() ( const struct initClustering &a, const struct initClustering &b ) const
  { 
    if (a.satCount > b.satCount)
      return false;
    else if (a.satCount == b.satCount)
      return (a.totalCost > b.totalCost);
    
    return true;
  } 
};

#endif

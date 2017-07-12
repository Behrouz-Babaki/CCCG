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

#ifndef SIEGBNBCONS_H
#define SIEGBNBCONS_H

#include "subsolver.h"
#include "SiegBnB.h"
#include <utility>
#include <vector>

using std::vector;
using std::pair;

struct CandidateBlock {
  int block;
  double quality;
  double distance_sum;
};

struct cSolution {
  double quality;
  vector<int> cluster;
  
  cSolution(double _quality, vector<int> _cluster) : quality(_quality), cluster(_cluster) {};
  
  void print();
};

class SiegBnBCons : public SubSolver {

public:
  
  SiegBnBCons (const vector<unsigned int>&, const vector<pair<unsigned int, unsigned int> >&, const vector<vector<double> >&, const vector<double> & , const double);
  SiegBnBCons (const vector<unsigned int>&, const vector<pair<unsigned int, unsigned int> >&, const vector<vector<double> >&, const vector<double> & , const double, int);
  virtual void solve(void);
  void solvePositive(void);
  bool hasSolution(void);
  virtual double getSolution(vector<bool>&);
  virtual vector<double> getSolutions(vector< vector<bool> >&);
  
protected:

  void check_dist();
  void pruning_power_lambdas();
  void set_compatibility();
  void print_best_cluster();
  void check_best(vector<int>&, const CandidateBlock&);
  void calc_dist_bound (const vector<CandidateBlock>&);
  bool prune0 (vector<int>& , const vector<CandidateBlock>&, double, double);
  bool prune1 (vector<int>& , const vector<CandidateBlock>&, double, double);
  bool prune2 (vector<int>& , const vector<CandidateBlock>&, double, double);
  bool prune3 (vector<int>& , const vector<CandidateBlock>&, double, double);
  void search ( vector<int> &, const vector<CandidateBlock> &, int, double, double);
  void search ( vector<int> &, const vector<CandidateBlock> &, double, double);
  int getSize();
  double getDist(int, int);
  double getPointDist(int, int);
  double getBlockDist(int, int);
  void computeBlockDist ();
  void initialize_lambSum();
  void initialize_within_distances ();
  int getIndex(int, int, int);
  void initialize_mlBlocks(void);
  bool checkCompatibility ( int i, int j );

  

  vector<vector<int> > _compatibility;
  vector<double> _dist_bound;
  vector<cSolution> _best_clusters;
  vector<unsigned int> _inputMlBlocks;
  vector<unsigned int> _mlBlockIndices;
  vector<pair<unsigned int, unsigned int> > _inputClBlocks;
  vector<double> _within_dist;
  vector<double> _block_lambSum;
  vector<vector<unsigned int> > _mlBlocks;
  vector<double> _block_distances;
  vector<vector<double> > _point_block_distances; 


  double _best;
  bool _solved;
  bool _mlBlocksInitialized;
  bool _onlyNegative;
  int _size;
  unsigned int _numberOfMlBlocks;

};

bool comp1 ( const pair<int,double> &a, const pair<int,double> &b );
bool comp2 ( const CandidateBlock &a, const CandidateBlock &b );
bool comp3 ( const CandidateBlock &a, const CandidateBlock &b ); 

inline unsigned int getBnBVersion(){
  return 1;
}

#endif

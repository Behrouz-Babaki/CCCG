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

#ifndef SIEGBNB_H
#define SIEGBNB_H

#include "subsolver.h"
#include <utility>
#include <vector>

struct Candidate {
  int point;
  double quality;
  double distance_sum;
};

struct Solution {
  double quality;
  std::vector<int> cluster;
  
  Solution(double _quality, std::vector<int> _cluster) : quality(_quality), cluster(_cluster) {};
  
  void print();
};

class SiegBnB : public SubSolver {

public:
  
  SiegBnB (const std::vector<std::vector<double> >&, const std::vector<double> & , const double);
  SiegBnB (const std::vector<std::vector<double> >&, const std::vector<double> & , const double, int);
  virtual void solve(void);
  void solvePositive(void);
  bool hasSolution(void);
  virtual double getSolution(std::vector<bool>&);
  virtual std::vector<double> getSolutions(std::vector< std::vector<bool> >&);
  
protected:

  void check_dist();
  void pruning_power_lambdas();
  void set_compatibility();
  void print_best_cluster();
  void check_best(std::vector<int>&, const Candidate&);
  void calc_dist_bound (const std::vector<Candidate>&);
  bool prune0 (std::vector<int>& , const std::vector<Candidate>&, double, double);
  bool prune1 (std::vector<int>& , const std::vector<Candidate>&, double, double);
  bool prune2 (std::vector<int>& , const std::vector<Candidate>&, double, double);
  bool prune3 (std::vector<int>& , const std::vector<Candidate>&, double, double);
  void search ( std::vector<int> &, const std::vector<Candidate> &, int, double, double);
  void search ( std::vector<int> &, const std::vector<Candidate> &, double, double);
  int getSize();
  double getDist(int, int);
  

  std::vector<std::vector<int> > _compatibility;
  std::vector<double> _dist_bound;
  // best clusters, last one is best so far
  std::vector<Solution> _best_clusters;

  double _best;
  bool _solved;
  bool _onlyNegative;
  int _size;

};

#endif

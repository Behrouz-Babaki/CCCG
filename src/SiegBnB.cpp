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

#include "SiegBnB.h"
#include <iostream>
#include <math.h>
#include <algorithm>
using namespace std;

bool comp1 ( const std::pair<int,double> &a, const std::pair<int,double> &b ) 
{ return a.second < b.second; } 

bool comp2 ( const Candidate &a, const Candidate &b )
{ return a.distance_sum > b.distance_sum; } 

bool comp3 ( const Candidate &a, const Candidate &b ) 
{ return a.quality < b.quality ; } 


void SiegBnB::check_dist () {
  for ( int i = 0, size = getSize(); i < size; ++i )
    for ( int j = 0; j < i; ++j ) {
      if ( sqrt(getDist ( i, j )) > sqrt(_lambdas[i]) + sqrt(_lambdas[j]) )
	cout << "PAIR " << i << " AND " << j << ": " << getDist ( i, j ) << ">" << sqrt(_lambdas[i]) + sqrt(_lambdas[j]) << endl;
    }
}

void SiegBnB::pruning_power_lambdas() {
  int prune = 0;
  int total = 0;
  for ( int i = 0, size = getSize (); i < size; ++i )
    for ( int j = 0; j < i; ++j ) 
      if ( _lambdas[i] > 0 && _lambdas[j] > 0 ) {
	total++;
        if ( sqrt(getDist ( i, j )) > sqrt(_lambdas[i]) + sqrt(_lambdas[j]) )
	  prune++;
      }
  cout << total << " edges " << prune << " pruned" << endl;
}

void SiegBnB::set_compatibility () {
  _compatibility.resize ( getSize () );
  for ( int i = 0, size = getSize (); i < size; ++i ) {
    _compatibility[i].resize ( getSize () );
    for ( int j = 0; j < i; ++j ) {
      if ( sqrt(getDist ( i, j )) > sqrt(_lambdas[i]) + sqrt(_lambdas[j]) )
        _compatibility[i][j] = _compatibility[j][i] = 0;
      else
        _compatibility[i][j] = _compatibility[j][i] = 1;
    }
  }
}

void SiegBnB::solvePositive() {
  cout << "SiegBnB: solving also with positives\n";
  this -> _onlyNegative = false;
  this -> solve();
  this -> _onlyNegative = true;
}

bool SiegBnB::hasSolution(){
  return (!_best_clusters.empty());
}


void SiegBnB::print_best_cluster () {
  _best_clusters.back().print();
}

void Solution::print()
{
  cout << this->quality << " ";
  int size2 = this->cluster.size ();
  cout << "(" << size2 << ") ";
  for ( int i = 0; i < size2; ++i )
    cout << this->cluster[i] << " ";
  cout << endl;
}


void SiegBnB::check_best ( vector<int> &cluster, const Candidate &candidate ) {
  if ( candidate.quality < _best ) {
    _best = candidate.quality;
    cluster.push_back ( candidate.point );
    _best_clusters.push_back(Solution(_best, cluster));
    cluster.pop_back ();
    //print_best_cluster ();
  }
}

void SiegBnB::calc_dist_bound ( const vector<Candidate> &candidates ) {
  _dist_bound.resize ( 0 );
  vector<double> dist_list;
  for ( vector<Candidate>::const_iterator x = candidates.begin (); x != candidates.end (); ++x ) 
    for ( vector<Candidate>::const_iterator y = candidates.begin (); y != x; ++y )
      dist_list.push_back ( getDist ( x->point, y->point ) );
  sort ( dist_list.begin (), dist_list.end () );
  int pos = 0;
  for ( int i = 0, size = candidates.size (); i < size; ++i ) {
    double sum = 0;
    for ( int j = 0; j < i; ++j ) {
      sum += dist_list[pos];
      ++pos;
    }
    _dist_bound.push_back ( sum );
  }
}

// no additional bound
bool SiegBnB::prune0 ( vector<int> &in_cluster, const vector<Candidate> &candidates, double distance_sum, double lambda_sum ) {
  return false;
}

// this bound seems easy to calculate and reasonably effective in practice;
// its complexity is O(n log n), where n is the number of remaining candidate
// points.
bool SiegBnB::prune1 ( vector<int> &in_cluster, const vector<Candidate> &candidates, double distance_sum, double lambda_sum ) {
  static vector<double> lambda_values;
  static vector<double> distance_values;
  lambda_values.resize ( 0 );
  distance_values.resize ( 0 );
  double new_lambda_sum = lambda_sum;
  double new_distance_sum = distance_sum;
  int new_size = in_cluster.size ();
  for ( vector<Candidate>::const_iterator x = candidates.begin (); x != candidates.end (); ++x ) {
    lambda_values.push_back ( _lambdas[x->point] );
    distance_values.push_back ( x->distance_sum );
  }
  sort ( lambda_values.rbegin (), lambda_values.rend () );
  sort ( distance_values.begin (), distance_values.end () );
  for ( int i = 0, size = candidates.size (); i < size; ++i ) {
    new_lambda_sum += lambda_values[i];
    new_distance_sum += distance_values[i];
    new_size++;
    if ( ( ( new_distance_sum + _dist_bound[i] ) / new_size ) - new_lambda_sum < _best )
      return false;
  }
  return true;
}

// modification of prune1, which takes more time as it recalculates within-remaining
// points-distances. Although the bound is better, it seems not worth the effort 
// in practice. O(n^2 log n) to calculate.
bool SiegBnB::prune2 ( vector<int> &in_cluster, const vector<Candidate> &candidates, double distance_sum, double lambda_sum ) {
  calc_dist_bound ( candidates );
  return prune1 ( in_cluster, candidates, distance_sum, lambda_sum );
}

// another tighter bound, which is harder to calculate and does not 
// prune effectively in practice, it seems.
// the problem is that this bound is O(n^2log n) to calculate.
// I guess it should be O(n log n) to be acceptable.
bool SiegBnB::prune3 ( vector<int> &in_cluster, const vector<Candidate> &candidates, double distance_sum, double lambda_sum ) {
  static vector<pair<int,double> > cs;
  cs.resize ( 0 );
  int size = candidates.size ();
  int is = in_cluster.size () + 1;
  for ( int i = 0; i < size; ++i ) 
    cs.push_back ( make_pair ( candidates[i].point, candidates[i].distance_sum - is * _lambdas[candidates[i].point] ) );
  for ( int i = 1; i <= size; ++i ) {
    sort (cs.begin (), cs.end (), comp1);
    double val = distance_sum - is * lambda_sum;
    for ( int j = 0; j < i; ++j ) 
      val += cs[j].second;
    if ( val / is < _best )
      return false;
    for ( int j = 0; j < size; ++j )
      cs[j].second -= _lambdas[cs[j].first];
    ++is;
  }
  return true;
}

void SiegBnB::search ( vector<int> &in_cluster, const vector<Candidate> &candidates, int i, double distance_sum, double lambda_sum ) {
  const Candidate &x = candidates[i];
  in_cluster.push_back ( x.point );
  vector<Candidate> new_candidates;
  Candidate c;
  int new_size = in_cluster.size () + 1;

  lambda_sum += _lambdas[x.point];
  distance_sum += x.distance_sum;
  
  // scores for all candidates, exclude candidates that do not score well
  for ( int j = i + 1, size = candidates.size (); j < size; ++j ) {
    c.point = candidates[j].point;
    if ( _compatibility[x.point][c.point] == 0 && _onlyNegative)
      continue;
    c.distance_sum = candidates[j].distance_sum + getDist ( x.point, c.point );
    c.quality = ( distance_sum + c.distance_sum ) / new_size - lambda_sum - _lambdas[c.point];
    if ( c.quality < 0 || !_onlyNegative ) {
      new_candidates.push_back ( c );
      check_best ( in_cluster, c );      
    }
  }
  
  sort ( new_candidates.begin (), new_candidates.end (), comp2);
  search ( in_cluster, new_candidates, distance_sum, lambda_sum );
  
  in_cluster.pop_back ();
}

void SiegBnB::search ( vector<int> &in_cluster, const vector<Candidate> &candidates, double distance_sum, double lambda_sum ) {
  if ( prune1 ( in_cluster, candidates, distance_sum, lambda_sum ) )
    return;
  int size = candidates.size ();
  double bound = 0.0;
  // simple idea: can we still use the remaining  lambdas to improve the current
  // best value, assuming the error will not get worse?
  for ( int i = 0; i < size; ++i )
    bound += _lambdas[candidates[i].point];
  for ( int i = 0; i < size; ++i ) {
    bound -= _lambdas[candidates[i].point];
    // additional condition does not seem to help much in practice.
    if ( candidates[i].quality - bound < _best ) 
      search ( in_cluster, candidates, i, distance_sum, lambda_sum );
  }
}

void SiegBnB::solve () {

  vector<int> in_cluster;
  vector<Candidate> candidates;
  set_compatibility ();
  Candidate c;
  c.distance_sum = 0;

  for ( int i = 0; i < getSize (); ++i )
    if ( _lambdas[i] > 0 ) {
      c.point = i;
      c.quality = -_lambdas[i];
      candidates.push_back ( c );
      check_best ( in_cluster, c );
    }
      
  sort ( candidates.begin (), candidates.end (), comp3);
  calc_dist_bound ( candidates );
  
  search ( in_cluster, candidates, 0, 0 ) ;
  
  this -> _solved = true;
}

SiegBnB::SiegBnB (const std::vector<std::vector<double> > &distances, const std::vector<double> & lambdas , const double sigma) : SubSolver(distances, lambdas, sigma), _best(sigma), _solved(false), _onlyNegative (true) { 
   
  this -> _size = lambdas.size();
}

SiegBnB::SiegBnB (const std::vector<std::vector<double> > &distances, const std::vector<double> & lambdas , const double sigma, int size) : SubSolver(distances, lambdas, sigma), _best(sigma), _solved(false), _onlyNegative (true), _size(size) { }

double SiegBnB::getDist(int first, int second){
  return this -> _distances[first][second];
}


double SiegBnB::getSolution(vector<bool>& outSolution){

  if (! _solved)
    this -> solve();

  outSolution.clear();
  outSolution.resize(getSize(), false);

  Solution &best = _best_clusters.back();

  for (vector<int>::iterator itr = best.cluster.begin(); itr != best.cluster.end(); itr++)
    outSolution[*itr] = true;

  return best.quality;
}

std::vector<double> SiegBnB::getSolutions(vector< std::vector<bool> >& outSolutions)
{

  if (! _solved)
    this -> solve();

  unsigned nSols = _best_clusters.size();
  
  outSolutions.clear();
  outSolutions.reserve(nSols);
  vector<double> outQuals;
  outQuals.reserve(nSols);
  
  for (unsigned i=0; i!=nSols; i++) {
    Solution &cur = _best_clusters[i];
    vector<bool> clus(getSize(), false);
    for (vector<int>::iterator itr = cur.cluster.begin(); itr != cur.cluster.end(); itr++)
      clus[*itr] = true;    
    outSolutions.push_back(clus);
    outQuals.push_back(cur.quality);
  }
  
  return outQuals;
}


int SiegBnB::getSize(void){
  return this -> _size;
}



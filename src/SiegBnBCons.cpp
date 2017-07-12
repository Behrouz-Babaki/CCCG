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

#include "SiegBnBCons.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <map>

#define SAFE_GT(a,b) ((a) > (b))
#define SAFE_LT(a,b) ((a) < (b))

using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::make_pair;
using std::swap;
using std::map;
using std::isfinite;

int verbosity = 0;

bool comp2 ( const CandidateBlock &a, const CandidateBlock &b )
{ return a.distance_sum > b.distance_sum; } 

bool comp3 ( const CandidateBlock &a, const CandidateBlock &b ) 
{ return a.quality < b.quality ; } 


void SiegBnBCons::initialize_mlBlocks() {
  _mlBlocks.clear();
  _mlBlockIndices.resize(_size, -1);
  map<unsigned int, unsigned int> indexMap;
  for (int counter = 0, size = _inputMlBlocks.size(); counter < size; counter++){
    map<unsigned int, unsigned int>::iterator itr = indexMap.find(_inputMlBlocks[counter]);
    if (itr == indexMap.end()){
      vector<unsigned int> vec;
      vec.push_back(counter);
      _mlBlocks.push_back(vec);
      indexMap[_inputMlBlocks[counter]] = _mlBlocks.size() - 1;
      _mlBlockIndices[counter] = _mlBlocks.size() - 1;
    }
    else
      {
	unsigned int index = itr->second;
	_mlBlocks[index].push_back(counter);
	_mlBlockIndices[counter] = index;
      }
  }
  this -> _numberOfMlBlocks = _mlBlocks.size();
  this -> _mlBlocksInitialized = true;

  //printing ml_Blocks
  if (verbosity >= 2){
    cout << "printing _mlBlocks: " << endl;
    for (size_t blockCoutner = 0; blockCoutner < _numberOfMlBlocks; blockCoutner++){
      cout << "data points in block " << blockCoutner << ": " << endl;
      for (size_t counter = 0; counter < _mlBlocks[blockCoutner].size(); counter++)
	cout << _mlBlocks[blockCoutner][counter] << " ";
      cout << endl;
    }

    cout << "checking indices:" << endl;
    for (int instanceCounter = 0; instanceCounter < _size; instanceCounter++)
      cout << "instance " << instanceCounter << ": block " << _mlBlockIndices[instanceCounter] << endl;
  }
}

bool SiegBnBCons::checkCompatibility ( int i, int j ) {
  if ( ( _within_dist[i] + _within_dist[j] + getBlockDist(i,j) ) / (double) ( _mlBlocks[i].size () + _mlBlocks[j].size () ) > _block_lambSum[i] + _block_lambSum[j] )
    return false;
  else
    return true;
}

void SiegBnBCons::set_compatibility () {
  if (!_mlBlocksInitialized)
    this->initialize_mlBlocks();

  _compatibility.clear();
  _compatibility.resize (_numberOfMlBlocks);
  for (size_t counter = 0; counter < _numberOfMlBlocks; counter++)
    _compatibility[counter].resize(_numberOfMlBlocks, 1);

  for (auto itr = _inputClBlocks.begin(), endItr = _inputClBlocks.end(); itr != endItr; itr++){
    unsigned int firstBlock, secondBlock;
    firstBlock = _mlBlockIndices[itr->first];
    secondBlock = _mlBlockIndices[itr->second];
    _compatibility[firstBlock][secondBlock] = _compatibility[secondBlock][firstBlock] = 0;
  }

  for (size_t counter1 = 0; counter1 < _numberOfMlBlocks; counter1++) {
    for (size_t counter2 = 0; counter2 < _numberOfMlBlocks; counter2++) {
      _compatibility[counter1][counter2] = _compatibility[counter2][counter1] =
        _compatibility[counter2][counter1] && checkCompatibility ( counter1, counter2 ) 
	&& !( ( _mlBlocks[counter1].size () == 1 && _mlBlocks[counter2].size () == 1  &&
	// special case, original Hansen can be applied
	sqrt(getBlockDist(counter1,counter2)) > sqrt(_block_lambSum[counter1]) + sqrt(_block_lambSum[counter2]) ) );
    }
  }
  
  if (verbosity >= 2){
    cout << "checking correct enforcement of can-not-links at the block level" << endl;
    int cl = 0;
    for (size_t counter1 = 0; counter1 < _numberOfMlBlocks; counter1++)
      for (size_t counter2 = 0; counter2 < _numberOfMlBlocks; counter2++, cl++)
	if (!_compatibility[counter1][counter2])
	  cout << " " << counter1 << "-" << counter2 
	       << (cl % 10 == 0 ?"\n" : " ");
  }

}

void SiegBnBCons::solvePositive() {
  
  cout << "SiegBnBCons: solving also with positives\n";
  
  double maxLamb = _block_lambSum[0];
  size_t maxIndex = 0;
  for (size_t counter = 1, size = _block_lambSum.size(); counter < size; counter++)
    if (maxLamb > _block_lambSum[counter]){
      maxIndex = counter;
      maxLamb = _block_lambSum[counter];
    }

  assert (maxLamb <= 0);

  vector<int> bestPosBlock;
  bestPosBlock.resize(_mlBlocks[maxIndex].size());
  for (size_t counter = 0, size = _mlBlocks[maxIndex].size();
       counter < size; counter++)
    bestPosBlock[counter] = (int) _mlBlocks[maxIndex][counter];

  cSolution positiveSol(maxLamb, bestPosBlock);
  _best_clusters.clear();
  _best_clusters.push_back(positiveSol);

}

bool SiegBnBCons::hasSolution(){
  return (!_best_clusters.empty());
}

//#define INTERRUPT_HACK
#ifdef INTERRUPT_HACK
long recursion_counter = 0;
bool count_recursion = false;
long recursion_timeout = 10000;
#endif


void SiegBnBCons::print_best_cluster () {
}

void SiegBnBCons::check_best ( vector<int> &cluster, const CandidateBlock &candidate ) {

  if (verbosity >= 3){
    cout << "quality: " << candidate.quality << endl;
    if (isfinite(_best))
      cout << "_best is:" << _best << endl;
    else 
      cout << "_best is infinity" << endl;
  }

  if ( _best - candidate.quality > 0.000001 ) { // NEED epsilon: SCIP ignores columns with too small redcost
#ifdef INTERRUPT_HACK
    // in case of INTERRUPT_HACK, epsilon is uber important as we may not start the count for to-be-ignored columns
    count_recursion = true;
    //recursion_counter = 0; // start over again, only stop if no progress for a long time
#endif
    _best = candidate.quality;
    for (auto itr = _mlBlocks[candidate.block].begin(), endItr = _mlBlocks[candidate.block].end(); itr != endItr; itr++)
      cluster.push_back ( *itr );
    _best_clusters.push_back(cSolution(_best, cluster));

    if (verbosity >= 3){
      for (auto itr = cluster.begin(), endItr = cluster.end(); itr != endItr; itr++)
	cout << (*itr) << " ";
      cout << endl;
    }
    for (size_t counter = 0; counter < _mlBlocks[candidate.block].size(); counter++)
      cluster.pop_back ();
  }
}

void SiegBnBCons::calc_dist_bound ( const vector<CandidateBlock> &candidates ) {
  vector<double> dist_list;
  for ( auto  x = candidates.begin (); x != candidates.end (); ++x ) 
    for ( auto y = candidates.begin (); y < x; ++y ) {
        int size = ( _mlBlocks[x->block].size () *  _mlBlocks[y->block].size () );
        double score = _point_block_distances[x->block][y->block] / size;
	for ( int i = 0; i < size; ++i )
          dist_list.push_back ( score );
    }
  sort ( dist_list.begin (), dist_list.end () );

  int pos = 0;
  //TODO keep track of the total size of candidate clusters somewhere else
  int size = 0;
  for (auto itr = candidates.begin(), endItr = candidates.end(); itr != endItr; itr++)
    size += _mlBlocks[itr->block].size();

  for ( int i = 0; i < size; ++i ) {
    double sum = 0;
    for ( int j = 0; j < i; ++j ) {
      sum += dist_list[pos];
      ++pos;
    }
    _dist_bound.push_back ( sum );
  }
}

bool SiegBnBCons::prune1 ( vector<int> &in_cluster, const vector<CandidateBlock> &candidates, double distance_sum, double lambda_sum ) {
  static vector<double> lambda_values;
  static vector<double> distance_values;

  lambda_values.resize ( 0 );
  distance_values.resize ( 0 );

  double new_lambda_sum = lambda_sum;
  double new_distance_sum = distance_sum;
  int new_size = in_cluster.size ();
  
  for ( auto  x = candidates.begin (); x != candidates.end (); ++x ) 
    {
      int blockSize = _mlBlocks[x->block].size();
      double avgLamb = _block_lambSum[x->block] / blockSize;
      double avgDist = x->distance_sum / blockSize;

      for (int blockpointCount = 0; blockpointCount < blockSize; blockpointCount++){
	distance_values.push_back(avgDist);
	lambda_values.push_back(avgLamb);
      }
    }

  sort(lambda_values.rbegin(), lambda_values.rend());
  sort(distance_values.begin(), distance_values.end());

  int pointCount = 0;
  for (auto itr = candidates.begin(), endItr = candidates.end(); itr != endItr; itr++){
    int blockNumber = itr->block;

    for (int blockpointCount = 0, nblockpoints = _mlBlocks[blockNumber].size(); blockpointCount < nblockpoints; blockpointCount++, pointCount++){
      new_lambda_sum += lambda_values[pointCount];
      new_distance_sum += distance_values[pointCount];
      new_size++;
      if ( SAFE_LT ( ( ( new_distance_sum /*+ _dist_bound[pointCount]  + std::max(0.0,std::max((_dist_bound[pointCount+in_cluster.size ()]-distance_sum ),_dist_bound[pointCount] ) )*/  ) / new_size ) - new_lambda_sum, _best ) )
	return false;
    }
  }
  if (verbosity > 5){
	cout << "branch pruned!\tsum of lambdas:" << new_lambda_sum << "\tsum of distances:" 
	     << new_distance_sum << "\tsize of cluster:" << new_size << "\tscore: " 
	     << (( new_distance_sum /* + _dist_bound[pointCount]*/ ) / new_size ) - new_lambda_sum 
	     << " which is larger than best:"  << _best << endl;
  }
  return true;
}


void SiegBnBCons::search ( vector<int> &in_cluster, const vector<CandidateBlock> &candidates, int i, double distance_sum, double lambda_sum ) {
  const CandidateBlock &x = candidates[i];
  /* Isn't it better if we represent in_cluster as a set of blocks? */
  for (auto itr = _mlBlocks[x.block].begin(), endItr = _mlBlocks[x.block].end(); itr != endItr; itr++)
  in_cluster.push_back (*itr);

  vector<CandidateBlock> new_candidates;
  new_candidates.reserve ( candidates.size () ); // a lot of time turns out to be spent allocating in this loop... we can avoid this a little.
  CandidateBlock c;
  int new_size;
  lambda_sum += _block_lambSum[x.block];
  distance_sum += x.distance_sum;
  
  for ( int j = i + 1, size = candidates.size (); j < size; ++j ) {
    c.block = candidates[j].block;
    /* We do not have a notion of block compatibility yet.
     * Instead, we use _compatibility for storing can_not_link constraints */
    if ( _compatibility[x.block][c.block] == 0 /*&& _onlyNegative*/)
      continue;
    c.distance_sum = candidates[j].distance_sum + getBlockDist ( x.block, c.block );
    new_size = in_cluster.size () + _mlBlocks[c.block].size();
    c.quality = ( distance_sum + c.distance_sum ) / new_size - lambda_sum - _block_lambSum[c.block];

        if (verbosity >= 3)
      {
	cout << "computing quality: distance_sum: " << distance_sum << " c.distance_sum is: " << c.distance_sum << " _within_dist[c.block]: " << _within_dist[c.block] << " new_size: " << new_size << " lambda_sum: " << lambda_sum << " _block_lambSum[c.block]: " << _block_lambSum[c.block] << endl;
	cout << "c.quality is: " << c.quality << endl;
      }

    /* In the unconstrained B&B, here we exclude the candidates
     * with non-negative qualities. */
    new_candidates.push_back ( c );
    check_best(in_cluster, c);
  }
  
  /* In the unconstrained case, using this order turned out to be VERY important.
   * Here, we again use 'quality' of blocks for sorting them.  */
  sort ( new_candidates.begin (), new_candidates.end (), comp2);
  
  search ( in_cluster, new_candidates, distance_sum, lambda_sum );

  for (size_t counter = 0, blockSize = _mlBlocks[x.block].size(); counter < blockSize; counter++)
    in_cluster.pop_back ();
}


void SiegBnBCons::search ( vector<int> &in_cluster, const vector<CandidateBlock> &candidates, double distance_sum, double lambda_sum ) {
#ifdef INTERRUPT_HACK
  if ( count_recursion ) {
    recursion_counter++;
    if ( recursion_counter == recursion_timeout )
      return;
  }
#endif
  int size = candidates.size ();
  double bound = 0.0;
  // simple idea: can we still use the remaining  lambdas to improve the current
  // best value, assuming the error will not get worse?
  for ( int i = 0; i < size; ++i )
    bound += _block_lambSum[candidates[i].block];
  if ( SAFE_GT(distance_sum / in_cluster.size () - lambda_sum - bound, _best ) ) {
    return;
  }


  if ( prune1 ( in_cluster, candidates, distance_sum, lambda_sum ) )
    return;
  for ( int i = 0; i < size; ++i ) {
    bound -= _block_lambSum[candidates[i].block];
    // additional condition does not seem to help much in practice.
    if ( SAFE_LT(candidates[i].quality /* + _dist_bound[size - i - 1]/(in_cluster.size () + size - i - 1) */ - bound, _best )) 
      search ( in_cluster, candidates, i, distance_sum, lambda_sum );
   
    if (verbosity >= 2){
      cout << "partial cluster: " << endl;
      for (auto itr = in_cluster.begin(), endItr = in_cluster.end(); itr != endItr; itr++)
	cout << *itr << " ";
      cout << endl;
      
      cout << "candidate block to be added: " << i << "(";
      for (auto itr = _mlBlocks[candidates[i].block].begin(), endItr =  _mlBlocks[candidates[i].block].end(); itr != endItr; itr++)
	cout << *itr << " ";
      cout << ")" << endl;
      cout << "distance_sum: " << distance_sum << endl;
      cout << "lambda_sum: " << lambda_sum << endl;
    }

#ifdef INTERRUPT_HACK
    if ( recursion_counter == recursion_timeout )
      return;
#endif
  }
}

void SiegBnBCons::solve () {
#ifdef INTERRUPT_HACK
  count_recursion = false;
  recursion_counter = 0;
  recursion_timeout = _numberOfMlBlocks * 10000; // parameter can have big influence: too small is slow convergence (too high has no effect)
#endif

  vector<int> in_cluster;
  vector<CandidateBlock> candidates;
  /* set_compatibility could not be applied to the constrained setting
   * Could we come up with an alternative version? (Currently, we're
   * using it for enforcing can_not_link constraints.
   */
  set_compatibility ();
  CandidateBlock c;
  c.distance_sum = 0;
  
  #ifndef NDEBUG
  size_t wd_size = _within_dist.size();
  size_t mb_size = _mlBlocks.size();
  size_t ls_size = _block_lambSum.size();
  assert (wd_size == mb_size);
  assert (mb_size == ls_size);
  #endif

  for ( size_t counter = 0; counter < _numberOfMlBlocks; ++counter ){
    double blockQuality = _within_dist[counter] / _mlBlocks[counter].size() - _block_lambSum[counter];

    if (verbosity >= 3)
      cout << "Block quality for block " << counter << ": " << blockQuality << endl;
      if ( blockQuality < 0 || !_onlyNegative) {
	c.block = counter;
	c.quality = blockQuality;
	c.distance_sum = _within_dist[counter];
	candidates.push_back ( c );
	check_best ( in_cluster, c );
      }
    }
  
  sort ( candidates.begin (), candidates.end (), comp3);
  
  search ( in_cluster, candidates, 0, 0 ) ;
  
  this -> _solved = true;
}

SiegBnBCons::SiegBnBCons (const std::vector<unsigned int>& inputMlBlocks, const std::vector<std::pair<unsigned int, unsigned int> >& inputClBlocks, const std::vector<std::vector<double> > &distances, const std::vector<double> & lambdas , const double sigma) : SubSolver(distances, lambdas, sigma), _inputMlBlocks(inputMlBlocks), _inputClBlocks(inputClBlocks), _best(sigma), _solved(false), _mlBlocksInitialized(false), _onlyNegative (true), _numberOfMlBlocks(-1) { 
  this -> _size = lambdas.size();
  this -> initialize_mlBlocks();
  this -> initialize_within_distances();
  this -> initialize_lambSum();
  this -> computeBlockDist();
}

SiegBnBCons::SiegBnBCons (const std::vector<unsigned int>& inputMlBlocks, const std::vector<std::pair<unsigned int, unsigned int> >& inputClBlocks, const std::vector<std::vector<double> > &distances, const std::vector<double> & lambdas , const double sigma, int size) : SubSolver(distances, lambdas, sigma), _inputMlBlocks(inputMlBlocks), _inputClBlocks(inputClBlocks), _best(sigma), _solved(false), _mlBlocksInitialized(false), _onlyNegative (true), _size(size), _numberOfMlBlocks(-1) { 
  this -> initialize_mlBlocks();
  this -> initialize_within_distances();
  this -> initialize_lambSum();
  this -> computeBlockDist();
}

double SiegBnBCons::getPointDist(int first, int second){
  return this -> _distances[first][second];
}


double SiegBnBCons::getSolution(vector<bool>& outSolution){

  if (! _solved)
    this -> solve();

  outSolution.clear();
  outSolution.resize(getSize(), false);

  cSolution &best = _best_clusters.back();

  for (vector<int>::iterator itr = best.cluster.begin(); itr != best.cluster.end(); itr++)
    outSolution[*itr] = true;

  return best.quality;
}

std::vector<double> SiegBnBCons::getSolutions(vector< std::vector<bool> >& outSolutions)
{

  if (! _solved)
    this -> solve();

  unsigned nSols = _best_clusters.size();
  
  outSolutions.clear();
  outSolutions.reserve(nSols);
  vector<double> outQuals;
  outQuals.reserve(nSols);
  
  for (unsigned i=0; i!=nSols; i++) {
   cSolution &cur = _best_clusters[i];
    vector<bool> clus(getSize(), false);
    for (vector<int>::iterator itr = cur.cluster.begin(); itr != cur.cluster.end(); itr++)
      clus[*itr] = true;    
    outSolutions.push_back(clus);
    outQuals.push_back(cur.quality);
  }
  
  return outQuals;
}


int SiegBnBCons::getSize(void){
  return this -> _size;
}


double SiegBnBCons::getBlockDist (int firstBlock , int secondBlock){

  if (firstBlock == secondBlock)
    return 0;

  if (firstBlock > secondBlock)
    swap (firstBlock, secondBlock);

  return _block_distances[getIndex(firstBlock, secondBlock, _numberOfMlBlocks)];
  
}

void SiegBnBCons::computeBlockDist (){
  if (!_mlBlocksInitialized)
    initialize_mlBlocks();

   _block_distances.clear();
   _block_distances.resize(_numberOfMlBlocks * (_numberOfMlBlocks-1) / 2);
   _point_block_distances.resize(_size, vector<double> (_numberOfMlBlocks, 0));
  for (size_t i = 0; i < _numberOfMlBlocks; i++)
      for (size_t j = i+1; j < _numberOfMlBlocks; j++){
	double total_dist =  0;
	for (auto itr1 = _mlBlocks[i].begin(), endItr1 = _mlBlocks[i].end(); itr1 != endItr1; itr1++)
	  for (auto itr2 = _mlBlocks[j].begin(), endItr2 = _mlBlocks[j].end(); itr2 != endItr2; itr2++){
	    _point_block_distances[*itr1][j] += getPointDist(*itr1, *itr2);
	    total_dist += this-> getPointDist(*itr1, *itr2);
	}
        _block_distances [getIndex(i,j,_numberOfMlBlocks)] = total_dist;
	if (verbosity >= 3)
	  cout << "The distance between blocks " << i << " and " << j << ": " << total_dist << endl;
    }
}

void SiegBnBCons::initialize_lambSum(){
  if (!_mlBlocksInitialized)
    initialize_mlBlocks();
  
  _block_lambSum.clear();
  _block_lambSum.reserve(_numberOfMlBlocks);
  for (size_t i = 0; i < _numberOfMlBlocks; i++)
    {
      double lambSum = 0;
      for (auto itr = _mlBlocks[i].begin(), endItr = _mlBlocks[i].end(); itr != endItr; itr++)
	lambSum += _lambdas[*itr];
      _block_lambSum.push_back (lambSum);
      if (verbosity >= 3)
	cout << "sum of lambdas for block " << i << ": " << lambSum << endl;
    }
}


void SiegBnBCons::initialize_within_distances () {
  if (!_mlBlocksInitialized)
    initialize_mlBlocks();

  _within_dist.clear();
  _within_dist.reserve(_numberOfMlBlocks);

  int counter = 0;
  for (auto itr1 = _mlBlocks.begin(), endItr1 = _mlBlocks.end(); itr1 != endItr1; itr1++, counter++)
    {
      double within = 0;
      for (auto itr2 = itr1->begin(), endItr2 = itr1->end(); itr2 != endItr2; itr2++)
	for (auto itr3 = itr1->begin(); itr3 != itr2; itr3++)
	  within += this->getPointDist(*itr2, *itr3);
      
      _within_dist.push_back(within);
      if (verbosity >= 3)
	cout << "sum of within-block distances for block " << counter << ": " << within << endl;
    }
}

int SiegBnBCons::getIndex(int row, int column, int n){
   return - ((row)*(row + 3) /2) + row * n + column - 1;
}

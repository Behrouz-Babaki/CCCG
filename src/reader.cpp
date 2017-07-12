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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
using namespace std;

InstanceReader::InstanceReader(string fileName){

  ifstream f;
  istringstream ss;
  string str;

  _coordinates.clear();

  f.open(fileName.c_str());
  if (f.is_open()){
    getline(f, str);
    ss.clear();
    ss.str(str);
    
    vector<double> values;
    double value;
    while (ss >> value)
      values.push_back(value);
    
    _coordinates.push_back(values);
    _dimension = values.size();
    
    ss.clear();
    getline(f,str);
    
    while (f.good())
      {
	values.clear();
	ss.clear();
	ss.str(str);
	
	while (ss >> value)
	  values.push_back(value);
	
	if (!values.empty())
	_coordinates.push_back(values);
	
	getline(f,str);
	
      }
  } else {
    cerr << "Error! Can't open instance file '" << fileName << "', quitting...\n";
    exit(-1);
  }
  
  _number_of_elements = _coordinates.size();

}


vector<vector<double> > InstanceReader::get_distances(void){
  if (_distances.empty())
    this -> compute_distances();

  return this -> _distances;
}

void InstanceReader::compute_distances(void){
  _distances.clear();
  
  for (vector<vector<double> >::size_type i = 0; i < _number_of_elements; i++)
    {
      vector<double> current_distances;

      for (vector<vector<double> >::size_type j = 0; j < _number_of_elements; j++)
	current_distances.push_back(compute_distance(_coordinates[i], _coordinates[j]));
      
      _distances.push_back(current_distances);
	
    }

}

int InstanceReader::get_number_of_instances(void){
    return this -> _number_of_elements;
 }

int InstanceReader::get_dimension(void){
  return this -> _dimension;
}

vector<vector<double> > InstanceReader::get_coordinates(void){
  return this -> _coordinates;
}
  

double InstanceReader::compute_distance(const vector<double>& first , const vector<double>& second){

  double distance = 0;

  if (first.size() != second.size() )
    {
      cerr << "different dimensions" << endl;
      exit (EXIT_FAILURE);
    }

  for (vector<double>::size_type i = 0; i < first.size(); i++)
    distance += (first[i] - second[i]) * (first[i] - second[i]);

  return distance;
}

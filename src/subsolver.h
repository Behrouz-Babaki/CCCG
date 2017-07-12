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

#ifndef SUBSOLVER_H
#define SUBSOLVER_H

#include "reader.hpp"

#include <vector>
#include <iostream>

   /**@class MIPSubSolver
    * @brief solver class for the optimal clustering problem
    *
    *  this class implements a solver for the optimal clustering problem as an mip model, which will be solved using SCIP
    */
class SubSolver
   {
   protected:

     /** @brief distances between instances */
     std::vector<std::vector<double> > _distances;

      /** lambda values */
      std::vector<double> _lambdas;

      double _sigma;

      

   public:
     /** @brief constructs the BP model for the optimal clustering subproblem */
      SubSolver(const std::vector<std::vector<double> > &, const std::vector<double> & , const double);
      virtual ~SubSolver();
      virtual void solve(void) = 0; ///< solves the subproblem
      virtual double getSolution(std::vector<bool> &) = 0;
      
      
   };

#endif


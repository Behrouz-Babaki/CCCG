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

#ifndef __MSSCBRANCHRYANFOSTER_H__
#define __MSSCBRANCHRYANFOSTER_H__

#include <vector>
#include <unordered_map>

#include "scip/scip.h"
#include "objscip/objbranchrule.h"

#include "cons_mlcl.h"

using namespace std;
using namespace scip;

/**
 *  @brief Ryan Foster branching for MSSC
 */
class BranchruleRyanFoster : public ObjBranchrule
{
private:
    unordered_map<string,int> _conname_map;
    
public:
   /** constructor */
   BranchruleRyanFoster(SCIP* scip, const vector<SCIP_CONS*>& setcover_con)
   : ObjBranchrule(scip,                /**< SCIP data structure */
                   "branch_RF",         /**< name of branching rule */
                   "Ryan Foster branching for MSSC",    /**< description of branching rule */
                   50000,               /**< priority of the branching rule */
                   -1,                  /**< maximal depth level, up to which this branching rule should be used (or -1) */
                   1.0)                 /**< maximal relative distance from current node's dual bound to primal bound
                                          *   compared to best node's dual bound for applying branching rule
                                          *   (0.0: only on current best node, 1.0: on all nodes) */
   {

       for (unsigned i=0; i != setcover_con.size(); i++) {
           string name = SCIPconsGetName(setcover_con[i]);
           _conname_map[name] = i;
       }
   }

   /** destructor */
   virtual ~BranchruleRyanFoster()
   {
   }

   /** branching execution method for fractional LP solutions
    */
   virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);
   
   SCIP_COL* getColOfVar(SCIP* scip, SCIP_VAR* var);
   
   int getBestCand(SCIP* scip, int nlpcands, SCIP_VAR** lpcands, SCIP_Real* lpcandsfrac);
   
   std::vector<bool> getClusterOfCol(SCIP* scip, SCIP_COL* col);
   
   // finds two entries in clus1/2: first one where both are 1, second one where one is 1 and other is 0 (or inverse)
   bool findClusterConflict(vector<bool>& cluster1, vector<bool>& cluster2, int *first, int *second);
   
   int getMLBlock (SCIP*, int);
};

#endif

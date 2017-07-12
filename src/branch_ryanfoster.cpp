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

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <algorithm>

#include "branch_ryanfoster.h"
#include "wrappers.h"

using namespace scip;
using namespace std;

SCIP_DECL_BRANCHEXECLP(BranchruleRyanFoster::scip_execlp)
{
    assert(scip != NULL);
    assert(result != NULL);
    
    int verbose = 2;
    if (verbose >= 1)
        printf("\nBranching rule 'Ryan-Foster' started\n");

    /* get branching candidates */
    SCIP_VAR** lpcands;
    SCIP_Real* lpcandsfrac;
    int nlpcands;
    SCIP_CALL( WrapperSCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands) );
    assert(nlpcands > 0);
    
    typedef pair<SCIP_Real, SCIP_VAR*> BranchCands; // fracvalue, lpcand
    vector<BranchCands> cands(nlpcands);
    for (int i = 0; i != nlpcands; i++) {
        cands[i].first = lpcandsfrac[i];
        cands[i].second = lpcands[i];
    }
    // reverse sorts according to fracvalue (the SCIP_Real)
    sort(cands.rbegin(), cands.rend());
    
    if (verbose >= 4) {
        for (unsigned i = 0; i != cands.size(); i++)
            cout << cands[i].first << " -- " << SCIPvarGetName(cands[i].second) << "\n";
    }
    
    
    int first = -1;
    int second = -1;
    bool match = false;
    for (unsigned c1=0; c1!=cands.size(); c1++) {
        if (verbose >= 3)
            printf(" -> c1 candidate: variable <%s> (frac=%g)\n", SCIPvarGetName(cands[c1].second), cands[c1].first);
        SCIP_COL* col1 = BranchruleRyanFoster::getColOfVar(scip, cands[c1].second);
        vector<bool> cluster1 = BranchruleRyanFoster::getClusterOfCol(scip, col1);
        if (verbose >= 3)
            cout << "Above is comparing cluster\n";
        
        // TODO: check that c2 not already constrained with c1
        for (unsigned c2=c1+1; c2!=cands.size(); c2++) {
            if (verbose >= 3)
                printf(" -> c2 candidate: variable <%s> (frac=%g)\n", SCIPvarGetName(cands[c2].second), cands[c2].first);
            SCIP_COL* col2 = BranchruleRyanFoster::getColOfVar(scip, cands[c2].second);
            vector<bool> cluster2 = BranchruleRyanFoster::getClusterOfCol(scip, col2);
        
            match = BranchruleRyanFoster::findClusterConflict(cluster1, cluster2, &first, &second);
            if (match) {
                if (verbose >= 2)
                    printf("Found cluster conflict match: %i -- %i\n", first, second);
                break;
            }   
        }
        if (match)
            break;
    }
    
    if (!match) {
        printf("Err... no match found : (\n");
        *result = SCIP_DIDNOTRUN;
        return SCIP_OKAY;
    }

    /* create the MIP-b&b-tree child nodes of the current node */
    // TODO INV what are other alternatives to SCIPgetLocalTransEstimate(scip) ?
    SCIP_NODE* childML;
    SCIP_CALL( SCIPcreateChild(scip, &childML, 0.0, SCIPgetLocalTransEstimate(scip)) );
    SCIP_NODE* childCL;
    SCIP_CALL( SCIPcreateChild(scip, &childCL, 0.0, SCIPgetLocalTransEstimate(scip)) );
    
    /* create corresponding constraints */
    char con_name[255];
    SCIP_CONS* consML;
    SCIPsnprintf(con_name, 255, "ML-%i-%i", first, second);
    SCIP_CALL( SCIPcreateConsML(scip, &consML, con_name,
                  first,                       /* the first data point to be linked */
                  second                       /* the second data point to be linked */
                  ) );
    SCIP_CONS* consCL;
    SCIPsnprintf(con_name, 255, "CL-%i-%i", first, second);
    SCIP_CALL( SCIPcreateConsCL(scip, &consCL, con_name,
                  first,                       /* the first data point to be linked */
                  second                       /* the second data point to be linked */
                  ) );
    
    /* add constraints to nodes */
    SCIP_CALL( SCIPaddConsNode(scip, childML, consML, NULL) );
    SCIP_CALL( SCIPaddConsNode(scip, childCL, consCL, NULL) );

    /* release constraints */
    SCIP_CALL( SCIPreleaseCons(scip, &consML) );
    SCIP_CALL( SCIPreleaseCons(scip, &consCL) );

    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}

/**
 * iterate over cluster1 and cluster2
 * find two entries:
 * - one with true and true in both clusters
 * - one with XOR true false // false true in the clusters
 */
bool BranchruleRyanFoster::findClusterConflict(std::vector<bool>& cluster1, std::vector<bool>& cluster2, int *first, int *second) {
    assert(cluster1.size() == cluster2.size());
    
    int match_and = -1;
    int match_xor = -1;
    
    unsigned length = cluster1.size();
    for (unsigned i = 0; i != length; i++) {
        if (match_and == -1
           && (cluster1[i] == true && cluster2[i] == true)) {
            match_and = i;
        
            if (match_xor != -1)
                break; // both found
        }
        if (match_xor == -1) {
            if ( (cluster1[i] == false && cluster2[i] == true)
               || (cluster1[i] == true && cluster2[i] == false) ) {
                match_xor = i;
                
                if (match_and != -1)
                    break; // both found
            }
        }
    }
    
    if (match_and != -1 and match_xor != -1) {
        *first = match_and;
        *second = match_xor;
        return true;
    }
    
    return false;
}

SCIP_COL* BranchruleRyanFoster::getColOfVar(SCIP* scip, SCIP_VAR* var){
    SCIP_COL * col;
    int nCols;
    SCIP_COL ** cols;
    int counter;
    bool found;
    
    col = NULL;
    nCols = SCIPgetNLPCols (scip);
    cols = SCIPgetLPCols (scip);
    found = false;
    for (counter = 0; !found && counter < nCols; counter++) {
        if (SCIPvarGetIndex(SCIPcolGetVar (cols[counter])) == SCIPvarGetIndex(var)){
            col = cols[counter];
            found = true;
        }
    }
        
    assert (found == true);
    return col;
}

int BranchruleRyanFoster::getBestCand(SCIP* scip, int nlpcands, SCIP_VAR** lpcands, SCIP_Real* lpcandsfrac) {
    int bestcand;
    SCIP_Real fractionality;
    SCIP_Real bestfractionality;
    int i;
    
    bestcand = -1;
    bestfractionality = 1;
    
    /* search the least fractional candidate */
    for(i = 0; i < nlpcands; ++i ) {
        assert(lpcands[i] != NULL);
        fractionality = lpcandsfrac[i];
        fractionality = MIN( fractionality, 1.0-fractionality );
        if ( fractionality < bestfractionality )
        {
            bestfractionality = fractionality;
            bestcand = i;
        }
    }
    
    assert(bestcand >= 0);
    assert(SCIPisFeasPositive(scip, bestfractionality));
    
    return bestcand;
}

vector<bool> BranchruleRyanFoster::getClusterOfCol(SCIP* scip, SCIP_COL* col) {
    vector<bool> rowIndices(_conname_map.size(), false);
    
    int verbose = 0;
    
    int nRows = SCIPcolGetNNonz(col);
    SCIP_ROW** rows = SCIPcolGetRows(col);
    
    for (int i=0; i!= nRows; i++) {
        try {
            string name = SCIProwGetName(rows[i]);
            int index = _conname_map.at(name);
            
            rowIndices[index] = true;
        } catch (const std::out_of_range& oor) {
            // OK, non-setcover constraint (e.g. not a datapoint)
        }
    }
    
    if (verbose >= 1) {
        cout << "cluster of col '" << SCIPvarGetName(SCIPcolGetVar(col)) << ": { ";
        for (unsigned i = 0; i != rowIndices.size(); i++)
            if (rowIndices[i])
                cout << i << " ";
        cout << "}\n";
    }
    
    return rowIndices;
}

/* This function is supposed to return the representative datapoint for a block of must-linked points
 * Currently, we assume that there are no must-link constraints, and we simply return the datapoint itself. 
 */
//TODO update this so that it returns the must-link block number
int BranchruleRyanFoster::getMLBlock (SCIP* scip, int rowId){
    return rowId;
}

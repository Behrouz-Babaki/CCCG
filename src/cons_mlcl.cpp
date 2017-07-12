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
#include <string>
#include <iostream>
#include <boost/iterator/iterator_concepts.hpp>
#include "cons_mlcl.h"


#include "objscip/objscip.h"

using namespace scip;
using namespace std;

struct SCIP_ConsData
{
  unsigned first;
  unsigned second;
  bool ML;
  bool propagated;
};


/** frees specific constraint data */
SCIP_DECL_CONSDELETE(ConshdlrMLCL::scip_delete)
{
   assert(consdata != NULL);
   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrMLCL::scip_trans)
{
   SCIP_CONSDATA* sourcedata = NULL;
   SCIP_CONSDATA* targetdata = NULL;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocMemory(scip, &targetdata) );
   targetdata->first = sourcedata->first;
   targetdata->second = sourcedata->second;
   targetdata->ML = sourcedata->ML;
   targetdata->propagated = sourcedata->propagated;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

SCIP_DECL_CONSENFOLP(ConshdlrMLCL::scip_enfolp)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

SCIP_DECL_CONSENFOPS(ConshdlrMLCL::scip_enfops)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

SCIP_DECL_CONSCHECK(ConshdlrMLCL::scip_check)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

SCIP_DECL_CONSLOCK(ConshdlrMLCL::scip_lock)
{
   return SCIP_OKAY;
}

SCIP_DECL_CONSPROP(ConshdlrMLCL::scip_prop)
{
    int verbose = 0;
    
    assert(nconss == 0 || conss != NULL);
    
    if (nconss == 0) {
        *result = SCIP_DIDNOTRUN;   
        return SCIP_OKAY;
    }
    *result = SCIP_DIDNOTFIND;
    
    if (verbose >= 1)
        cout << "Running MLCL propagation, "<<nconss<<" constraints\n";
    
    // for *result:
    int nfixedvars = 0;
    bool cutoff = false;
    
    SCIP_CONSDATA* consdata;
    for (int i=0; i!=nconss; i++) {
        consdata = SCIPconsGetData(conss[i]);
        assert(consdata != NULL);
        
        if (consdata->propagated) {
            if (verbose >= 3)
                cout << "Already propagated " << SCIPconsGetName(conss[i]) << "\n";
            continue;
        }
        
        if (verbose >= 2)
            cout << "Doing propagation for " << SCIPconsGetName(conss[i]) << "\n";
        
        string name_first = setcover_names[consdata->first];
        string name_second = setcover_names[consdata->second];
        
        // loop over all columns
        bool fixedvars = false;
        bool checked = false;
        int nVars = SCIPgetNVars(scip);
        SCIP_VAR** vars = SCIPgetVars(scip);
        for (int v=0; v!=nVars; v++) {
            if (SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY || SCIPvarGetStatus(vars[v]) != SCIP_VARSTATUS_COLUMN) {
                // scip non-binary vars (e.g. perturbation vars)
                if (verbose >= 5)
                    cout << "Skipping non-binary (or non-column) var " << SCIPvarGetName(vars[v]) << "\n";
                continue;
            }
            checked = true;
            SCIP_COL* col = SCIPvarGetCol(vars[v]);
            
            bool match = ConshdlrMLCL::colHasConflict(scip, col, name_first, name_second, consdata->ML);
            if (match) {
                SCIP_Bool fixed;
                SCIP_Bool infeasible;
                
                // fix the var
                SCIP_CALL( SCIPfixVar(scip, vars[v], 0.0, &infeasible, &fixed) );
                
                if (verbose >= 2)
                    cout << "Fixing var " << SCIPvarGetName(vars[v]) << " ("<<SCIPvarGetLbLocal(vars[v])<<") infeas: " << infeasible << " fixed: " << fixed << "\n";
                
                if (infeasible) {
                    assert( SCIPvarGetLbLocal(vars[v]) > 0.5 );
                    SCIPdebugMessage("-> cutoff\n");
                    cutoff = true;
                }
                if (fixed) {
                    fixedvars = true;
                    nfixedvars++;
                }
            }
        }
        
        if (verbose >= 3) {
            if (!fixedvars)
                cout << " -> did not fix any vars...\n";
            else 
                cout << "Propping this constraint done\n";
        }
        
        if (checked)
           consdata->propagated = true;
    }
    
    if (cutoff) {
        cout << "Hmmm.... did a cutoff in MLCL propagation !!! ----------------------\n";
        *result = SCIP_CUTOFF;
    }
    else if (nfixedvars > 0) {
        if (verbose >= 1)
            cout << "Propping done, fixed " << nfixedvars << " vars\n";
        *result = SCIP_REDUCEDDOM;
    }
    
    return SCIP_OKAY;
}

void ConshdlrMLCL::print_constraints(SCIP* scip, const char* msg)
{
    SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "MLCL");
    
    int nconss = SCIPconshdlrGetNActiveConss(conshdlr); // only these first ones are active
    SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr);
    assert(nconss == 0 || conss != NULL);
 
    if (nconss != 0)
        cout << msg << "\n";
    
    SCIP_CONSDATA* consdata;
    for (int i=0; i!=nconss; i++) {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      
      if (consdata->ML) {
          cout << "ML ";
      } else {
          cout << "CL ";
      }
      cout << consdata->first << " -- " << consdata->second << "\n";
    }
}

bool ConshdlrMLCL::colHasConflict(SCIP* scip, SCIP_COL* col, string name_first, string name_second, bool ML)
{
    int verbose = 0;
    
    int nRows = SCIPcolGetNNonz(col);
    SCIP_ROW** rows = SCIPcolGetRows(col);
    
    bool has_first = false;
    bool has_second = false;
    for (int i=0; i!= nRows; i++) {
        string name = SCIProwGetName(rows[i]);
        if (name_first == name)
            has_first = true;
        if (name_second == name)
            has_second = true;
    }
    
    if ( ML &&
         ( (has_first == true && has_second == false) ||
           (has_first == false && has_second == true) ) ) {
        if (verbose >= 1) {
            cout << "ML and " << has_first << " sec: " << has_second << " -- conflict: " << SCIPvarGetName(SCIPcolGetVar(col)) << "\n";
            
            cout << "cluster of col '" << SCIPvarGetName(SCIPcolGetVar(col)) << "':";
            for (int i=0; i!= nRows; i++)
                cout << SCIProwGetName(rows[i]) << " ";
            cout << "\n";
        }
            
        return true;
    }
    if ( !ML &&
         (has_first == true && has_second == true) ) {
        if (verbose >= 1) {
            cout << "CL and " << has_first << " sec: " << has_second << " -- conflict: " << SCIPvarGetName(SCIPcolGetVar(col)) << "\n";
            
            cout << "cluster of col '" << SCIPvarGetName(SCIPcolGetVar(col)) << "':";
            for (int i=0; i!= nRows; i++)
                cout << SCIProwGetName(rows[i]) << " ";
            cout << "\n";
        }
            
        return true;
    }
    
    return false;
}

bool ConshdlrMLCL::clusterSatisfied(SCIP* scip, std::vector< bool > cluster)
{
    SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "MLCL");
    
    int nconss = SCIPconshdlrGetNActiveConss(conshdlr); // only these first ones are active
    SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr);
    assert(nconss == 0 || conss != NULL);
    
    SCIP_CONSDATA* consdata;
    for (int i=0; i!=nconss; i++) {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      
      bool has_first = cluster[consdata->first];
      bool has_second = cluster[consdata->second];
      
      if (consdata->ML) {
          if ( (has_first == true && has_second == false) ||
               (has_first == false && has_second == true) )
              return false;
      } else { // CL
          if (has_first == true && has_second == true)
              return false;
      }
    }
    return true;
}


/** creates and captures an MSSC Must-Link constraint */
SCIP_RETCODE SCIPcreateConsML(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                  first,               /**< the first data point to be linked */
   int                  second               /**< the second data point to be linked */
   )
{
   SCIP_CONSDATA* consdata = NULL;

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   consdata->first = first;
   consdata->second = second;
   consdata->ML = true;
   consdata->propagated = false;
   
   return SCIPcreateConsMLCL(scip, cons, name, consdata);
}
/** creates and captures an MSSC Cannot-Link constraint */
SCIP_RETCODE SCIPcreateConsCL(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                  first,               /**< the first data point to be linked */
   int                  second               /**< the second data point to be linked */
   )
{
   SCIP_CONSDATA* consdata = NULL;

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   consdata->first = first;
   consdata->second = second;
   consdata->ML = false;
   consdata->propagated = false;
   
   return SCIPcreateConsMLCL(scip, cons, name, consdata);
}

/** creates and captures an MSSC Must-Link or Cannot-Link constraint (shared code) */
SCIP_RETCODE SCIPcreateConsMLCL(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;

   /* find the ML/CL constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "MLCL");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Must-Link constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, 
			     FALSE,            /**< should the LP relaxation of constraint be in the initial LP? */
			     FALSE,            /**< should the constraint be separated during LP processing? */
			     FALSE,            /**< should the constraint be enforced during node processing? */
			     FALSE,            /**< should the constraint be checked for feasibility?
						*   TODO: Check the constraint after the solution is generated */
			     TRUE,             /**< should the constraint be propagated during node processing?
						*   This is set to true to fix the variables based on ML constraints */
					       // TODO: check which routine is called by this...
			     FALSE,            /**< is constraint only valid locally? */
			     FALSE,            /**< is constraint modifiable (subject to column generation)? */
			     FALSE,            /**< is constraint dynamic? */
			     FALSE,            /**< should the constraint be removed from the LP due to aging or cleanup? */
			     FALSE             /**< should the constraint always be kept at the node where it was added? */
			     ) );

   return SCIP_OKAY;
}

SCIP_Bool SCIPgetIsML(SCIP_CONSDATA* consdata){
  return consdata -> ML;
}

unsigned int SCIPgetMlClFirst(SCIP_CONSDATA* consdata){
  return consdata -> first;
}

unsigned int SCIPgetMlClSecond(SCIP_CONSDATA* consdata){
  return consdata -> second;
}

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

#ifndef __MSSCCONSHDLRMLCL_H__
#define __MSSCCONSHDLRMLCL_H__

#include "objscip/objscip.h"
#include "wrappers.h"
#include <vector>
#include <string>
using namespace std;


/** C++ constraint handler for MSSC Must-Link and Cannot-Link constraints */
class ConshdlrMLCL : public scip::ObjConshdlr
{
private:
  vector<string> setcover_names;

public:
   /** default constructor */
   ConshdlrMLCL(
      SCIP* scip,
      const vector<SCIP_CONS*>& setcover_con
      )
     : WrapperObjConshdlr(scip,	                   /**< SCIP data structure */
		   "MLCL",	                   /**< name of connstraint handler */
		   "MSSC Must-Link/Cannot-Link constraints",   /**< description of constraint handler */
		   1000000,	                   /**< priority of the constraint handler for separation  
						        NOT APPLICABLE */
		   0,   	                   /**< priority of the constraint handler for constraint enforcing  
						        NOT APPLICABLE */
		   -2000000,	                   /**< priority of the constraint handler for checking infeasibility (and propagation)  
						        NOT APPLICABLE */
		   1,		                   /**< frequency for separating cuts; zero means to separate only in the root node  
						        NOT APPLICABLE */
		   1,		                   /**< frequency for propagating domains; zero means only preprocessing */
		   1,		                   /**< frequency for using all instead of only the useful constraints 
				                        in separation, propagation and enforcement, -1 for no eager evaluations, 
							0 for first only  
						        NOT APPLICABLE */
		   -1,		                   /**< maximal number of presolving rounds the constraint handler participates in 
						        (-1: no limit)  */
		   FALSE,	                   /**< should separation method be delayed, if other separators found cuts?  */
		   FALSE,	                   /**< should propagation method be delayed, if other propagators found reductions?  */
		   FALSE,	                   /**< should presolving method be delayed, if other presolvers found reductions?  */
		   TRUE,	                   /**< should the constraint handler be skipped, if no constraints are available?  */
		   SCIP_PROPTIMING_BEFORELP,       /**< positions in the node solving loop where propagation method of constraint 
						        handler should be executed  
						        We are assuming that ML constraints are automatically enforced by the 
						        pricing algorithm, so we do not execute the handler after each LP solving
						        during the pricing loop */
		   SCIP_PRESOLTIMING_FAST          /**<  timing mask of the constraint handler's presolving method */
    
		   ), setcover_names()
    {
        for (unsigned i=0; i != setcover_con.size(); i++) {
            setcover_names.push_back(SCIPconsGetName(setcover_con[i]));
        }
    }

   /** destructor */
   virtual ~ConshdlrMLCL()
   {
   }

   /** frees specific constraint data */
   virtual SCIP_DECL_CONSDELETE(scip_delete);

   /** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_DECL_CONSTRANS(scip_trans);


   /** constraint enforcing method of constraint handler for LP solutions */
   virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

   /** constraint enforcing method of constraint handler for pseudo solutions */
   virtual SCIP_DECL_CONSENFOPS(scip_enfops);

   /** feasibility check method of constraint handler for primal solutions */
   virtual SCIP_DECL_CONSCHECK(scip_check);


   /** variable rounding lock method of constraint handler */
   virtual SCIP_DECL_CONSLOCK(scip_lock);


   /** domain propagation method of constraint handler */
   virtual SCIP_DECL_CONSPROP(scip_prop);
   
   static void print_constraints(SCIP* scip, const char* msg="");
   
   static bool colHasConflict(SCIP* scip, SCIP_COL* col, std::string name_first, std::string name_second, bool ML);
   
   static bool clusterSatisfied(SCIP* scip, vector<bool> cluster);
};

/** creates and captures an MSSC Must-Link constraint */
SCIP_RETCODE SCIPcreateConsML(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                  first,               /**< the first data point to be linked */
   int                  second               /**< the second data point to be linked */
   );
/** creates and captures an MSSC Cannot-Link constraint */
SCIP_RETCODE SCIPcreateConsCL(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                  first,               /**< the first data point to be linked */
   int                  second               /**< the second data point to be linked */
   );

/** creates and captures an MSSC Must-Link or Cannot-Link constraint (shared code) */
SCIP_RETCODE SCIPcreateConsMLCL(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   );

/** returns whether a constraint is must-link or not */
SCIP_Bool SCIPgetIsML(SCIP_CONSDATA* consdata);

/** returns the id of first element of a link constraint */
unsigned int SCIPgetMlClFirst(SCIP_CONSDATA* consdata);

/** returns the id of second element of a link constraint */
unsigned int SCIPgetMlClSecond(SCIP_CONSDATA* consdata);

#endif

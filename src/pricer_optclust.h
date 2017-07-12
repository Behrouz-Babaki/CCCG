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

#ifndef __SCIP_PRICER_OPTCLUST_H__
#define __SCIP_PRICER_OPTCLUST_H__

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/pub_var.h>
#include <objscip/objscip.h>

#include "cons_mlcl.h"

#include <vector>
#include <map>
#include <fstream>

using namespace std;
using namespace scip;

enum Stabilisation { ST_None, ST_Bound, ST_Both, ST_Smooth, ST_Smooth2, ST_SmoothAuto};

void print_cluster(std::vector< bool > cluster, const char* msg="");

/** pricer class */
class ObjPricerOptClust : public ObjPricer
{
public:

      /** Constructs the pricer object with the data needed */
   ObjPricerOptClust(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,      /**< name of pricer */
      const int                           p_num_instances, /**< number of instances */
      const int                           p_num_clusters, /**< number of clusters */
      const vector< vector<double> >&      p_distance,  /**< matrix of distances */
      const vector<SCIP_CONS*>&           p_setcover_con,   /**< matrix of arc constraints */
      SCIP_CONS*                          p_side_con,   /**< array of partitioning constraints */
      Stabilisation                       stab,
      bool                                findAll,
      const char*                         coefFileName = NULL ,
      const char*                         columnFileName = NULL
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerOptClust();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** farkas pricing method of variable pricer for infeasible LPs */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas);

   /** perform (non-farkas) pricing */
   SCIP_RETCODE pricing(
      SCIP*              scip               /**< SCIP data structure */
      );

   /** perform pricing */
   SCIP_RETCODE farkasPricing(
      SCIP*              scip               /**< SCIP data structure */
      );


   /** add variable (column) to problem */
   SCIP_RETCODE add_cluster_variable(
      SCIP*              scip,               /**< SCIP data structure */
      const vector<bool>& cluster,           /**< cluster indicator */
      const vector<SCIP_Real>& lambdas,
      const SCIP_Real sigma
      );

   /** find+add cluster when in farkas */
   SCIP_Real find_add_farkas_cluster(
      SCIP* scip, 
      const vector<SCIP_Real> & lambdas,     /**< matrix of lengths */
      const SCIP_Real& sigma
      );
   
   /** find+add  clusters with negative reduced cost */
   SCIP_Real find_add_clusters(
      SCIP* scip, 
      const vector<SCIP_Real> & lambdas,     /**< matrix of lengths */
      const SCIP_Real& sigma
      );

   /** sets the vector of perturbation variables */
   void setPertVars (std::vector<SCIP_VAR *> vars, SCIP_Real mu, SCIP_Real epsilon);
   void setEpsilon(SCIP_Real epsilon);

   static SCIP_Real getClusterCost(const vector<vector<double> > &distance, const std::vector<bool> &cluster);


protected:

   /** return number of nodes */
   inline int num_instances() const
   {
      return _num_instances;
   }

   /** return number of nodes */
   inline double distance(int i , int j) const
   {
      return _distance[i][j];
   }


   /** return constraint corresponding to setcover constraint for node i */
   inline SCIP_CONS* setcover_con(
      const int          i                  /**< instance number */
      ) const
   {
      return ( _setcover_con[i] );
   }

   /** return partitioning constraint for node i */
   inline SCIP_CONS* side_con(void) const
   {
      return _side_con;
   }


   /** saves the coefficient values in a file */
   void saveCoef(const vector<SCIP_Real>&, const SCIP_Real&);
   
   /** stores the coefficient matrix of the master problem in a file */
   void saveLp(SCIP* scip, bool);

   
   /** modifying the bounds of perturbation variables */
   SCIP_RETCODE updatePertVarsBounds (SCIP*);

   /** modifying the objective values of perturbation variables */
   SCIP_RETCODE updatePertVarsCoefs (SCIP* , std::vector<SCIP_Real>);

   /** gets the constraints from constraint handler and stores them in vectors */
   void getLinkConstriants(SCIP*, std::vector<unsigned int>&, std::vector<std::pair<unsigned int, unsigned int> >&);

   // updates the duality gap
   void update_duality_gap(
      SCIP* scip, 
      const vector<SCIP_Real> & lambdas,     /**< matrix of lengths */
      const SCIP_Real& sigma,
      const SCIP_Real& reduced_cost
      );
      
   // get the current (fractional) solution value
   SCIP_Real getCurSol(SCIP* scip);
   
   // get the cluster corresponding to the column (variable)
   vector<bool> getColCluster(SCIP* scip, SCIP_VAR* colVar);
   
   // see paper "In-Out separation and col gen stab by dual price smoothing"
   SCIP_Real bigL(
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma,
      const SCIP_Real& reduced_cost
                             );
   SCIP_Real bigL(
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma,
      const vector<bool>& cluster
                             );
   SCIP_Real bigL_t(
      SCIP* scip,
      const vector<SCIP_Real> & lambdas,
      const SCIP_Real& sigma
   );
   
protected:

   const int _num_instances;
   const int _num_clusters;
   const vector<vector<double> > _distance;
   std::vector<SCIP_CONS*>  _setcover_con;
   SCIP_CONS*  _side_con;
   ofstream _outStream, _colStream;
   const char* _coefFileName;
   const char* _columnFileName;
   std::vector<SCIP_VAR *> _pertVars;

   Stabilisation _stab;
   bool _findAll;
   SCIP_Real _mu, _epsilon;
   SCIP_Real _lower_bound;
   SCIP_Real _upper_bound;
   SCIP_Real _duality_gap;
   SCIP_Real _bestLambSum;
   int _updateCount;

   // values of non-zero columns of previous pricer iteration
   // (used to check degeneracy)
   std::map<const char*, SCIP_Real> _prev_colVals;

   SCIP_Real _prev_out_score;
   SCIP_Real _piHat_score;
   SCIP_Real _piHat_redcost;
   std::vector<SCIP_Real> _piHat_lambdas;
   SCIP_Real _piHat_sigma;
   int _count_mispricings;
};

#endif


















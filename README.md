
 CCCG -- Constrained Clustering using Column Generation 


 What is it?
 -----------

 The problem  of clustering in presence of  additional constraints has
 been extensively  studied already, and a number  of efficient methods
 to solve these problems (approximately) have been developed. However,
 these  methods  are  often  limited  to only  a  number  of  specific
 constraints. CCCG is an attempt to build a more general framework for
 constrained clustering. It is  based on an integer linear programming
 formulation  of  clustering  problem,  and  hence  obtains  an  exact
 solution for this  problem. The clustering criterion used  in CCCG is
 similar to that  of k-means algorithm, namely to  minimize the sum of
 squared distances from cluster centers.

 Installation
 ------------

 CCCG depends on Boost C++ libraries, which could be obtained from:

 <http://www.boost.org/>

 CCCG also depends on SCIP which is a Mixed Integer programming
 solver. CCCG has been tested with versions 3.0.1 and 3.1.0. You can
 obtain SCIP for free from here:
 
 <http://scip.zib.de/>

 We recommend that you install not just the SCIP MIP solver, but
 rather the whole SCIP optimization suite, which also contains LP
 solver SOPLEX. SCIP will need such an LP solver for solving the LP
 relaxations. Alternatively, you could only install SCIP, and later
 configure it to use your favorite LP solver.

 Once you have installed SCIP, define the environment variable
 SCIP_HOME such that it points to the directory where you have
 installed SCIP. For example:

 $ export SCIP_HOME=/usr/local/scipoptsuite-3.0.1/scip-3.0.1

 Then go to the CCCG home directory and run:
 
 $ make depend
 $ make

 There are a number of example datasets and example sets of
 constraints in /examples directory. To verify that CCCG is working
 correctly, you could run it on these examples. Perhaps the easiest
 way to do so is to go to the CCCG home directory and run

 $ make test

 Usage
 -----

 Like many other clustering algorithms, CCCG assumes that the number
 of clusters are specified in advance by user. So the minimum input
 for CCCG are information about objects, and the number of
 clusters. After going to the CCCG home directory, one could run the
 following command:

 $ ./bin/cccg -d <DATA> -k <#CLUSTERS>

 In the above command, DATA is a file in which the dimensions of data
 objects are stroed, and #CLUSTERS is the number of clusters.

 If solving a problem takes too long, you can opt for the best
 solution found within a time limit. This can be done by using flag
 -x:

 $ ./bin/cccg -d <DATA> -k <#CLUSTERS> -x <TIME>

 In the above command, TIME is the timeout in seconds. To see if any
 solution has been found within the time limit and whether it is the
 optimal solution or not, you should check the value of "Primal
 Bound" and "SCIP Status" (reported near the end of output).

 To have CCCG consider constraints when generating the optimal
 clustering, they should be specified in an input file which is
 communicated to CCCG code using flag -n:
 
 $ ./bin/cccg -d <DATA> -k <#CLUSTERS> -n <FILE>

 In the input constraint file, the Must-Link (ML) and Can-Not-Link
 (CL) constraints are specified in the usual way: A pair of
 must-linked objects are followed by +1, and a pair of can-not-linked
 objects are followed by -1. For an example of constraint file, look
 into the /examples directory in CCCG home directory.

 There are several other options that could be used with CCCG. For a
 list of example commands using these options, look into the /examples
 directory. Some of these options are described below:

 	     * The efficiency of CCCG could be improved by providing an initial
	     clustering. This could be done by using the flag -i:

	     $ ./bin/cccg -d <DATA> -k <#CLUSTERS> -i <INIT>

	     * Instead of specifying a single clustering as an initial
               solution for CCCG, one could alternatively specify a
               directory containing several clusterings. In this case,
               CCCG code will use collects all clusters from all these
               clusterings, and will add each and every one of them as
               a column to the problem. This could be done by using
               the flag -r:

	     $ ./bin/cccg -d <DATA> -k <#CLUSTERS> -r <DIR>

	     * In each call to subproblem-solver, several columns with
               negative reduced cost are found. One could decide all
               these columns or only the most negative one. The
               default behavior is to do the latter. If one prefers to
               add all columns, this could be done via flag -a:

	     $ ./bin/cccg -d <DATA> -k <#CLUSTERS> -a

	     For a list of all available options, run ./bin/cccg with
	     no arguments.

 Availability
 ------------

 CCCG software source can be found on the CCCG page under
 <http://dtai.cs.kuleuven.be/CP4IM/cccg/>

 Authors
 -------

 Behrouz Babaki 
 Tias Guns 
 Siegfried Nijssen

 Licensing
 ---------

 Please see the file named LICENSE.

 Publications
 ------------

 Babaki,  B., Guns, T.,  & Nijjsen,  S. (2014)  Constrained Clustering
 using  Column   Generation.  Eleventh  International   Conference  on
 Integration of  Artificial Intelligence (AI)  and Operations Research
 (OR) techniques in Constraint Programming (CPAIOR 2014)

 Acknowledgements
 ----------------

 We are thankful to Khanh-Chuong Duong and Mohadeseh Ganji for their
 feedback.

 Contact
 -------
 
 For questions, bug reports, and any other matters, you could contact
 Behrouz Babaki via this email address:

 <Behrouz.Babaki@cs.kuleuven.be>

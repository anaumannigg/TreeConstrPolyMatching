#ifndef _linear_program_included_
#define _linear_program_included_

#include "gurobi_c++.h"
#include "cgal_includes.h"
#include "graph_computations.h"
#include "solution.h"
#include "tree_computations.h"

using namespace std;

class LinearProgram {
	

public:
    //sets up and solves ILP based on a bipartite graph G.
    //maximizes sum of chosen edges in constrained bipartite matching problem
    //the result is written into the provided solution object
	static void solveILP(GRBEnv& env,const Graph& g, int num_polys1, int num_polys2, Solution* solution);

    //solves ILP considering tree constraints
    static void solveILP_trees(GRBEnv& env,TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2, Solution* solution);

    //solves the tree constraint matching problem approximately using Canzar's Algorithm
    static void solveViaCanzar(GRBEnv& env,TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2, Solution& solution);
};

#endif
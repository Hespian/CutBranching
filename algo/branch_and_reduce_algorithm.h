/******************************************************************************
 * branch_and_reduce_algorithm.h
 *
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef BRANCH_AND_REDUCE_SOLVER_H
#define BRANCH_AND_REDUCE_SOLVER_H

// local includes
#include "fast_set.h"
#include "modified.h"

//#include "definitions.h"
#include "../VieCut/lib/data_structure/graph_access.h"
#include "../VieCut/lib/algorithms/global_mincut/viecut.h"
#include "../VieCut/lib/algorithms/flow/push_relabel.h"
#include "../VieCut/lib/algorithms/flow/preflow_push.h"

#include "../VieCut/lib/tools/timer.h"


// system includes
#include <vector>
#include <string>
#include <memory>

class branch_and_reduce_algorithm
{
	friend class modified;
	friend class fold;
	friend class alternative;

public:
	static int REDUCTION;
	static int LOWER_BOUND;
	static int BRANCHING;
	static bool outputLP;

	// 1 for extra decompose step
	static int EXTRA_DECOMP;

	static long defaultBranchings;
	static bool defaultBranch;

	std::vector<std::vector<int>> adj;
	static long nBranchings;
	static int debug;
	double SHRINK;
	int depth;
	int maxDepth;
	int rootDepth;
	int n;
	int N;

	/**
	 * current best solution
	 */
	int opt;
	std::vector<int> y;

	/**
	 * current solution (-1: not determined, 0: not in the vc, 1: in the vc, 2: removed by foldings)
	 */
	int crt;
	std::vector<int> x;

	/**
	 * #remaining vertices
	 */
	int rn;

	/**
	 * max flow
	 */
	std::vector<int> in;
	std::vector<int> out;

	/**
	 * lower bound
	 */
	int lb;

	std::vector<int> que, level, iter;

	std::vector<int> modTmp;

	std::vector<std::shared_ptr<modified>> modifieds;
	int modifiedN;

	/**
	 * Packing constraints
	 */
	////	std::list<std::vector<int>> packing;
	std::vector<std::vector<int>> packing;

	std::vector<int> vRestore;

	int reductionSnapshotSize;
	std::vector<int> snapshotX;

	branch_and_reduce_algorithm(std::vector<std::vector<int>>& _adj, int const _N);

	int deg(int v);
	void set(int v, int a);

	fast_set used;

	// helpers for modifying the graph
	void compute_fold(std::vector<int> const &S, std::vector<int> const &NS);
	void compute_alternative(std::vector<int> const &A, std::vector<int> const &B);
	void restore(int n);
	void reverse();

	// helpers for lpReduction
	bool dinicDFS(int v);
	void updateLP();

	// reduction methods
	bool lpReduction();
	bool deg1Reduction();
	bool dominateReduction();
	bool fold2Reduction();
	bool twinReduction();
	bool funnelReduction();
	bool deskReduction();
	bool unconfinedReduction();
	int packingReduction();

	// lower bounds for pruning
	int lpLowerBound();
	int cycleLowerBound();
	int cliqueLowerBound();
	int lowerBound();

	// recursive methods
	void branching(timer &t, double time_limit);
	bool decompose(timer &t, double time_limit);
	bool reduce();
	void rec(timer &t, double time_limit);

	// Track how much the starting solution helped
	bool startingSolutionIsBest = false;
	int numBranchesPrunedByStartingSolution = 0;

	// vestiges of original Java code
#if 0
	
	void debug(String str, Object...os) {
		StringBuilder sb = new StringBuilder();
		Calendar c = Calendar.getInstance();
		sb.append(String.format("%02d:%02d:%02d  ", c.get(Calendar.HOUR_OF_DAY), c.get(Calendar.MINUTE), c.get(Calendar.SECOND)));
		for (int i = 0; i < depth && i <= maxDepth; i++) sb.append(' ');
		System.err.print(sb);
		System.err.printf(str, os);
	}
#endif // 0

	int solve(timer &t, double time_limit);

	void initial_reduce_graph();
	void reduce_graph();

	void restore_to_snapshot();

	std::string debugString() const;
	void PrintState() const;

	size_t get_current_is_size() const;
	size_t get_current_is_size_with_folds() const;
	bool folded_vertices_exist() const;
	std::vector<int> compute_maximal_is();
	size_t compute_alternative_maximal_is_size();
	size_t number_of_nodes_remaining() const;
	size_t number_of_edges_remaining() const;
	//void force_into_independent_set(std::vector<NodeID> const &nodes);
	void extend_finer_is(std::vector<bool> &independent_set);
	void get_solved_is(std::vector<bool> &independent_set);

	void convert_adj_lists(graph_access &G, std::vector<NodeID> &reverse_mapping) const;

	void convert_to_ga(std::shared_ptr<graph_access> G, std::vector<NodeID> &reverse_mapping, std::vector<NodeID> &mapping);
	void convert_to_metis(int32_t* nNodes, std::vector<int32_t> &xadj, std::vector<int32_t> &adjncy, std::vector<NodeID> &reverse_mapping);

	void addStartingSolution(std::vector<int> solution, int solutionSize);


	// NEW BRANCHING RULES      

	// articulation points
    std::vector<int> articulation_points;
	void get_articulation_points();

	std::vector<int> visited;
    std::vector<int> minNr;      
    int current_dfs_num = 0;

	int get_articulation_point();
	void dfs_root(int s);
	void dfs(int v, int in);
	// ----------

	// global mincuts
	viecut cut_algo;
	int get_mincut_vertex();
	// ----------


	// st cuts
	NodeID s, t;
	int branch_t = 10;
	push_relabel flow;
	std::vector<NodeID> cut;

	void get_stcut_vertices();
	void find_st_vtcs(std::shared_ptr<graph_access> graph, NodeID ss = -1, NodeID tt = -1);
	//void generate_flow_graph(std::shared_ptr<graph_access> graph, std::shared_ptr<flow_graph> flow);
	// ------------

	// utility
	inline int get_max_deg_vtx();
	inline int get_min_deg_vtx(std::shared_ptr<graph_access> g);
	inline bool is_neighbour_of(std::shared_ptr<graph_access> graph, NodeID v, NodeID u);



#if 0
}
#endif // 0
};

#endif //BRANCH_AND_REDUCE_SOLVER_H
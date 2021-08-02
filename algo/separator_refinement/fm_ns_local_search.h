//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#ifndef FM_NS_LOCAL_SEARCH_P621XWW8
#define FM_NS_LOCAL_SEARCH_P621XWW8

#include "partial_boundary.h"
#include "../data_structures/maxNodeHeap.h"

#include "../definitions.h"

// information to perform undo 
struct change_set {
        NodeID node;
        PartitionID block;
};

class fm_ns_local_search {
public:
        fm_ns_local_search(std::vector<int> &part_ind, std::vector<int> &act, AdjList &adj);
        virtual ~fm_ns_local_search();

        EdgeWeight perform_refinement(AdjList & G, std::vector<int> &sep, int ub, bool balance = false, PartitionID to = 4);
        EdgeWeight perform_refinement(AdjList & G, PartialBoundary & separator, int ub, bool balance = false, PartitionID to = 4);
        
        EdgeWeight perform_local_search(bool balance = false, PartitionID to = 4);
        
        void set_target_separator(std::vector<int> &sep) 
        {
                std::vector<NodeID> start_nodes(sep);

                //random_functions::permutate_vector_good(start_nodes, false);
        
                for(NodeID node : start_nodes) {
                        Gain toLHS = 0;
                        Gain toRHS = 0;
                        compute_gain( G, node, toLHS, toRHS);

                        queues[0].insert(node, toLHS);
                        queues[1].insert(node, toRHS);
                }


                queues.resize(0);
                queues.resize(2);

                // initialize block weights/sizes
                block_weights.resize(0);
                block_weights.resize(3, 0);
                for (NodeID node = 0; node < G.size(); node++) {
                        if (x[node] < 0)
                        {
                                if (partition_index[node] == 0) {
                                        block_weights[0] += 1;
                                } else if (partition_index[node] == 1) {
                                        block_weights[1] += 1;
                                } else {
                                        block_weights[2] += 1;
                                }
                        }
                }
        }

        std::vector<int> &partition_index;
        std::vector<int> &x;
        AdjList &G;

        std::vector<NodeWeight> block_weights;
        std::vector<maxNodeHeap> queues;

private: 
        void compute_gain( AdjList & G, NodeID node, Gain & toLHS, Gain & toRHS);
        void move_node( AdjList & G,  NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                        std::vector< NodeWeight > & block_weights,
                        std::vector< bool > & moved_out_of_S,
                        std::vector< maxNodeHeap > & heaps,
                        std::vector< change_set > & rollback_info);

        void move_node( AdjList & G,  NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                        std::vector< NodeWeight > & block_weights,
                        std::vector< bool > & moved_out_of_S,
                        std::vector< maxNodeHeap > & heaps,
                        std::vector< change_set > & rollback_info,
                        PartialBoundary & separator);

        std::vector<int> moved_nodes;
};


inline
void fm_ns_local_search::compute_gain( AdjList & G, NodeID node, Gain & toLHS, Gain & toRHS) {
        toLHS = 1;
        toRHS = 1;

        for (int neighbor : G[node])
        {
                if (x[neighbor] < 0) {
                        int part_ind = partition_index[neighbor];
                        if (part_ind == 0) {
                                toRHS -= 1;
                        } else if (part_ind == 1) {
                                toLHS -= 1;
                        }
                }
        }
}

inline
void fm_ns_local_search::move_node( AdjList & G, NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                                    std::vector< NodeWeight > & block_weights,
                                    std::vector< bool > & moved_out_of_S, 
                                    std::vector< maxNodeHeap > & queues,
                                    std::vector< change_set > & rollback_info) {

        change_set cur_move;
        cur_move.node = node;
        cur_move.block = partition_index[node];
        rollback_info.push_back(cur_move);

        partition_index[node] = to_block;
        block_weights[to_block] += 1;
        block_weights[2] -= 1;
        moved_out_of_S[node] = true;

        std::vector<NodeID> to_be_added;
        std::vector<NodeID> to_be_updated; // replace by hashmap?
        Gain gain_achieved = 1;

        for (int neighbor : G[node])
        {
                if (x[neighbor] < 0) {
                        if (partition_index[neighbor] == other_block) {
                                change_set cur_move;
                                cur_move.node = neighbor;
                                cur_move.block = partition_index[neighbor];
                                rollback_info.push_back(cur_move);

                                partition_index[neighbor] = 2;
                                block_weights[other_block] -= 1;
                                block_weights[2]           += 1;
                                gain_achieved              -= 1;

                                if(!moved_out_of_S[neighbor]) {
                                        to_be_added.push_back(neighbor);
                                }

                                for (int v : G[neighbor]) {
                                        if(x[v] < 0) {
                                                if (queues[0].contains(v)) {
                                                        to_be_updated.push_back(v);
                                                }
                                        }
                                }

                        } else if(partition_index[neighbor] == 2) {
                                to_be_updated.push_back(neighbor);
                        }
                }
        }

        Gain toLHS = 0;
        Gain toRHS = 0;

        for(NodeID node : to_be_added) {
                compute_gain(G, node, toLHS, toRHS);
                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }

        for(NodeID node : to_be_updated) {
                compute_gain(G, node, toLHS, toRHS);
                queues[0].changeKey(node, toLHS);
                queues[1].changeKey(node, toRHS);
        }
}

inline
void fm_ns_local_search::move_node( AdjList & G, NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                                    std::vector< NodeWeight > & block_weights,
                                    std::vector< bool > & moved_out_of_S, 
                                    std::vector< maxNodeHeap > & queues,
                                    std::vector< change_set > & rollback_info,
                                    PartialBoundary & separator) {

        change_set cur_move;
        cur_move.node = node;
        cur_move.block = partition_index[node];
        rollback_info.push_back(cur_move);
        separator.deleteNode(node);

        partition_index[node] = to_block;
        block_weights[to_block] += 1;
        block_weights[2] -= 1;
        moved_out_of_S[node] = true;
        moved_nodes.push_back(node);

        std::vector<NodeID> to_be_added;
        std::vector<NodeID> to_be_updated; // replace by hashmap?
        Gain gain_achieved = 1;

        for (int neighbor : G[node]) 
        {
                if (x[neighbor] < 0) {
                        if(partition_index[node] == other_block) {
                                change_set cur_move;
                                cur_move.node = neighbor;
                                cur_move.block = partition_index[neighbor];
                                rollback_info.push_back(cur_move);

                                partition_index[neighbor] = 2;
                                separator.insert(neighbor);

                                block_weights[other_block] -= 1;
                                block_weights[2]           += 1;
                                gain_achieved              -= 1;

                                if(!moved_out_of_S[neighbor]) {
                                        to_be_added.push_back(neighbor);
                                }

                                for (int v : G[neighbor]) {
                                        if(x[v] < 0) {
                                                if(queues[0].contains(v)) {
                                                        to_be_updated.push_back(v);
                                                }   
                                        }
                                }

                        } else if(partition_index[neighbor] == 2) {
                                to_be_updated.push_back(neighbor);
                        }
                }
        }        

        Gain toLHS = 0;
        Gain toRHS = 0;

        for(NodeID node : to_be_added) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }

        for(NodeID node : to_be_updated) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].changeKey(node, toLHS);
                queues[1].changeKey(node, toRHS);
        }
}

#endif /* end of include guard: FM_NS_LOCAL_SEARCH_P621XWW8 */

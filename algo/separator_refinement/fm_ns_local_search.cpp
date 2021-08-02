//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#include <algorithm>
#include "fm_ns_local_search.h"
#include "../data_structures/maxNodeHeap.h"
#include "../tools/random_functions.h"

fm_ns_local_search::fm_ns_local_search(std::vector<int> &part_ind, std::vector<int> &act, AdjList &adj) :
        partition_index(part_ind), x(act), G(adj)
{
}

fm_ns_local_search::~fm_ns_local_search() {
                
}

#define sep_fm_unsucc_steps 25

inline
EdgeWeight fm_ns_local_search::perform_local_search(bool balance, int to) {

        std::vector<change_set> rollback_info;
        std::vector<bool> moved_out_of_separator(G.size(), false);
        int upper_bound_partition = 1000000000;
        
        NodeWeight best_separator  = block_weights[2];
        NodeWeight input_separator = block_weights[2];
        int best_diff              = abs((int)block_weights[1]-(int)block_weights[0]);
        int undo_idx = 0;        
        
        int steps_till_last_improvement = 0;
        //roll forwards
        while( steps_till_last_improvement < sep_fm_unsucc_steps) {
                Gain gainToA = queues[0].maxValue();
                Gain gainToB = queues[1].maxValue();

                Gain top_gain        = 0;
                PartitionID to_block = 0;
               
                if(balance) {
                        top_gain = queues[to].maxValue();
                        to_block = to;
                } else {
                        if( gainToA == gainToB ) {
                                top_gain = gainToA;
                                to_block = random_functions::nextInt(0,1); 
                        } else {
                                top_gain = gainToA > gainToB ? gainToA : gainToB;
                                to_block = top_gain == gainToA ? 0 : 1;
                        }
                }

                Gain other_gain = gainToA > gainToB ? gainToB : gainToA;
                PartitionID other_block = to_block == 0 ? 1 : 0;

                NodeID nodeToBlock = queues[to_block].maxElement();
                if( block_weights[to_block] + 1 < upper_bound_partition) {
                        queues[to_block].deleteMax();
                        queues[other_block].deleteNode(nodeToBlock);
                        move_node(G, nodeToBlock, to_block, other_block, block_weights, moved_out_of_separator, queues, rollback_info);
                } else {
                        NodeID nodeOtherBlock = queues[other_block].maxElement();
                        if( other_gain >= 0 && block_weights[other_block] + 1 < upper_bound_partition) {
                                queues[other_block].deleteMax();
                                queues[to_block].deleteNode(nodeOtherBlock);
                                move_node(G, nodeOtherBlock, other_block, to_block, block_weights, moved_out_of_separator, queues, rollback_info);
                        } else {
                                // need to make progress (remove a random node from the queues)
                                if( nodeOtherBlock == nodeToBlock ) {
                                        queues[0].deleteMax();
                                        queues[1].deleteMax();
                                } else {
                                        int block = random_functions::nextInt(0,1);
                                        queues[block].deleteMax();
                                }
                        }
                }

                int cur_diff = abs((int)block_weights[1]-(int)block_weights[0]);
                if( block_weights[2] < best_separator || (block_weights[2] == best_separator && cur_diff < best_diff)  ) {
                        best_separator = block_weights[2];
                        undo_idx = rollback_info.size();
                        steps_till_last_improvement = 0;
                }  else {
                        steps_till_last_improvement++;
                }

                if( queues[0].empty() || queues[1].empty() ) {
                        break;
                }
        }

        // roll back 
        for( int i = rollback_info.size()-1; i >= undo_idx; i--) {
                partition_index[rollback_info[i].node] = rollback_info[i].block;
        }
                  
        return input_separator - best_separator;        
}

EdgeWeight fm_ns_local_search::perform_refinement(AdjList & G, std::vector<int> &sep, int ub, bool balance, PartitionID to) {

        std::vector<maxNodeHeap> queues; queues.resize(2);
        std::vector<bool> moved_out_of_separator(G.size(), false);
        std::vector<change_set> rollback_info;
        int upper_bound_partition = ub;

        std::vector<NodeID> start_nodes;
        
        // Initialize Todo-List with separator
        for (NodeID node : sep) {
                if (x[node] < 0) {
                        start_nodes.push_back(node);
                }
        }

        if(start_nodes.empty()) return 0;

        random_functions::permutate_vector_good(start_nodes, false);
        
        // insert nodes into left and right pq's
        for(NodeID node : start_nodes) {
                Gain toLHS = 0;
                Gain toRHS = 0;
                compute_gain( G, node, toLHS, toRHS);

                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }
        
        // initialize block weights/sizes
        std::vector< NodeWeight > block_weights(3,0);
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

        NodeWeight best_separator  = block_weights[2];
        NodeWeight input_separator = block_weights[2];
        int best_diff              = abs((int)block_weights[1]-(int)block_weights[0]);
        int undo_idx = 0;

        int steps_till_last_improvement = 0;
        //roll forwards
        while( steps_till_last_improvement < sep_fm_unsucc_steps) {
                Gain gainToA = queues[0].maxValue();
                Gain gainToB = queues[1].maxValue();

                Gain top_gain        = 0;
                PartitionID to_block = 0;
               
                if(balance) {
                        top_gain = queues[to].maxValue();
                        to_block = to;
                } else {
                        if( gainToA == gainToB ) {
                                top_gain = gainToA;
                                to_block = random_functions::nextInt(0,1); 
                        } else {
                                top_gain = gainToA > gainToB ? gainToA : gainToB;
                                to_block = top_gain == gainToA ? 0 : 1;
                        }
                }

                Gain other_gain = gainToA > gainToB ? gainToB : gainToA;
                PartitionID other_block = to_block == 0 ? 1 : 0;

                NodeID nodeToBlock = queues[to_block].maxElement();
                if( block_weights[to_block] + 1 < upper_bound_partition) {
                        queues[to_block].deleteMax();
                        queues[other_block].deleteNode(nodeToBlock);
                        move_node(G, nodeToBlock, to_block, other_block, block_weights, moved_out_of_separator, queues, rollback_info);
                } else {
                        NodeID nodeOtherBlock = queues[other_block].maxElement();
                        if( other_gain >= 0 && block_weights[other_block] + 1 < upper_bound_partition) {
                                queues[other_block].deleteMax();
                                queues[to_block].deleteNode(nodeOtherBlock);
                                move_node(G, nodeOtherBlock, other_block, to_block, block_weights, moved_out_of_separator, queues, rollback_info);
                        } else {
                                // need to make progress (remove a random node from the queues)
                                if( nodeOtherBlock == nodeToBlock ) {
                                        queues[0].deleteMax();
                                        queues[1].deleteMax();
                                } else {
                                        int block = random_functions::nextInt(0,1);
                                        queues[block].deleteMax();
                                }
                        }
                }

                int cur_diff = abs((int)block_weights[1]-(int)block_weights[0]);
                if( block_weights[2] < best_separator || (block_weights[2] == best_separator && cur_diff < best_diff)  ) {
                        best_separator = block_weights[2];
                        undo_idx = rollback_info.size();
                        steps_till_last_improvement = 0;
                }  else {
                        steps_till_last_improvement++;
                }

                if( queues[0].empty() || queues[1].empty() ) {
                        break;
                }
        }

        // roll back 
        for( int i = rollback_info.size()-1; i >= undo_idx; i--) {
                partition_index[rollback_info[i].node] = rollback_info[i].block;
        }

        return input_separator - best_separator;

}

EdgeWeight fm_ns_local_search::perform_refinement(AdjList & G, PartialBoundary & separator, int ub, bool balance, PartitionID to) {

        std::vector<maxNodeHeap> queues; queues.resize(2);
        std::vector<change_set> rollback_info;
        int upper_bound_partition = ub;
        std::vector<bool> moved_out_of_separator(G.size(), false);

        std::vector<NodeID> start_nodes;
        forall_boundary_nodes( separator, node ) {
                if(x[node] < 0)
                        start_nodes.push_back(node);
        } endfor 

        random_functions::permutate_vector_good(start_nodes, false);
        
        for( NodeID node : start_nodes ) {
                Gain toLHS = 0;
                Gain toRHS = 0;
                compute_gain( G, node, toLHS, toRHS);

                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }

        // initialize block weights/sizes
        std::vector< NodeWeight > block_weights(3,0);
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

        
        std::vector< NodeWeight > best_block_weights(3,0);
        best_block_weights = block_weights;
        NodeWeight best_separator  = block_weights[2];
        NodeWeight input_separator = block_weights[2];
        int best_diff              = abs((int)block_weights[1]-(int)block_weights[0]);
        int undo_idx = 0;

        int steps_till_last_improvement = 0;
        //roll forwards
        while( steps_till_last_improvement < sep_fm_unsucc_steps) {
                Gain gainToA = queues[0].maxValue();
                Gain gainToB = queues[1].maxValue();

                Gain top_gain        = 0;
                PartitionID to_block = 0;
               
                if(balance) {
                        top_gain = queues[to].maxValue();
                        to_block = to;
                } else {
                        if( gainToA == gainToB ) {
                                top_gain = gainToA;
                                to_block = random_functions::nextInt(0,1); 
                        } else {
                                top_gain = gainToA > gainToB ? gainToA : gainToB;
                                to_block = top_gain == gainToA ? 0 : 1;
                        }
                }
                
                Gain other_gain = gainToA > gainToB ? gainToB : gainToA;
                PartitionID other_block = to_block == 0 ? 1 : 0;

                NodeID nodeToBlock = queues[to_block].maxElement();
                if( block_weights[to_block] + 1 < upper_bound_partition ) {
                        queues[to_block].deleteMax();
                        queues[other_block].deleteNode(nodeToBlock);
                        move_node(G, nodeToBlock, to_block, other_block, block_weights, moved_out_of_separator, queues, rollback_info, separator);
                } else {
                        NodeID nodeOtherBlock = queues[other_block].maxElement();
                        if( other_gain >= 0 && block_weights[other_block] + 1 < upper_bound_partition) {
                                queues[other_block].deleteMax();
                                queues[to_block].deleteNode(nodeOtherBlock);
                                move_node(G, nodeOtherBlock, other_block, to_block, block_weights, moved_out_of_separator, queues, rollback_info, separator);
                        } else {
                                // need to make progress (remove a random node from the queues)
                                if( nodeOtherBlock == nodeToBlock ) {
                                        queues[0].deleteMax();
                                        queues[1].deleteMax();
                                } else {
                                        int block = random_functions::nextInt(0,1);
                                        queues[block].deleteMax();
                                }
                        }
                }

                int cur_diff = abs((int)block_weights[1]-(int)block_weights[0]);
                if( block_weights[2] < best_separator || (block_weights[2] == best_separator && cur_diff < best_diff)  ) {
                        best_separator              = block_weights[2];
                        undo_idx                    = rollback_info.size();
                        steps_till_last_improvement = 0;
                        best_block_weights          = block_weights;
                }  else {
                        steps_till_last_improvement++;
                }

                if( queues[0].empty() || queues[1].empty() ) {
                        break;
                }
        }

        // roll back 
        for( int i = rollback_info.size()-1; i >= undo_idx; i--) {
                if( partition_index[rollback_info[i].node] == 2 ) separator.deleteNode( rollback_info[i].node );
                partition_index[rollback_info[i].node] = rollback_info[i].block;
                if( partition_index[rollback_info[i].node] == 2 ) separator.insert( rollback_info[i].node );
        }
        block_weights = best_block_weights;      
        return input_separator - best_separator;
}


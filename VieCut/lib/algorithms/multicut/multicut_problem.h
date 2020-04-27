/******************************************************************************
 * multiterminal_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../../data_structure/mutable_graph.h"
#include "../../io/graph_io.h"

struct terminal {
    terminal() { }

    terminal(NodeID position, NodeID original_id)
        : position(position), original_id(original_id), invalid_flow(true) { }

    terminal(NodeID position, NodeID original_id, bool invalid_flow)
        : position(position),
          original_id(original_id),
          invalid_flow(invalid_flow) { }

    NodeID position;
    NodeID original_id;
    bool   invalid_flow;
};

struct multicut_problem {
    multicut_problem() { }

    explicit multicut_problem(std::shared_ptr<mutable_graph> G)
        : multicut_problem(G, std::vector<terminal>()) { }

    multicut_problem(std::shared_ptr<mutable_graph> G,
                     std::vector<terminal> term)
        : multicut_problem(G,
                           term,
                           std::vector<
                               std::shared_ptr<std::vector<NodeID> > >(),
                           -1,
                           std::numeric_limits<FlowType>::max(),
                           0,
                           { UNDEFINED_NODE, UNDEFINED_EDGE },
                           std::unordered_set<NodeID>()) { }

    multicut_problem(std::shared_ptr<mutable_graph> G,
                     std::vector<terminal> term,
                     std::vector<std::shared_ptr<std::vector<NodeID> > >
                     mappings,
                     FlowType lower,
                     FlowType upper,
                     EdgeWeight deleted,
                     std::pair<NodeID, EdgeID> prio,
                     std::unordered_set<NodeID> finished_bp)
        : graph(G),
          terminals(term),
          mappings(mappings),
          lower_bound(lower),
          upper_bound(upper),
          deleted_weight(deleted),
          priority_edge(prio),
          finished_blockpairs(finished_bp) { }

    NodeID mapped(NodeID n) const {
        NodeID n_coarse = n;
        for (const auto& map : mappings) {
            n_coarse = (*map)[n_coarse];
        }
        return n_coarse;
    }

    void addFinishedPair(NodeID a, NodeID b, NodeID numOriginalTerminals) {
        if (a == b) {
            //LOG1 << "Error. Pair between " << a << " and itself!";
            exit(1);
        }
        NodeID min = std::min(a, b);
        NodeID max = std::max(a, b);

        NodeID combined = min * numOriginalTerminals + max;
        finished_blockpairs.insert(combined);
    }

    void removeFinishedPair(NodeID a, NodeID b, NodeID numOriginalTerminals) {
        if (a == b) {
            //LOG1 << "Error. Pair between " << a << " and itself!";
            exit(1);
        }
        NodeID min = std::min(a, b);
        NodeID max = std::max(a, b);

        NodeID combined = min * numOriginalTerminals + max;
        finished_blockpairs.extract(combined);
    }

    bool isPairFinished(NodeID a, NodeID b, NodeID numOriginalTerminals) {
        if (a == b) {
            //LOG1 << "Error. Searching pair between " << a << " and itself!";
            exit(1);
        }
        NodeID min = std::min(a, b);
        NodeID max = std::max(a, b);

        NodeID combined = min * numOriginalTerminals + max;
        return (finished_blockpairs.count(combined) > 0);
    }

    static void writeGraph(std::shared_ptr<multicut_problem> problem,
                           std::string path) {
        std::shared_ptr<mutable_graph> g = std::make_shared<mutable_graph>();

        // LOG1 << problem->graph->n() << " nodes and "
        //      << problem->graph->m() << " edges";
        // LOG1 << "currently deleted " << problem->deleted_weight;

        // bfs around all terminals, print all edges between first 3000 nodes
        std::queue<NodeID> Q;
        std::unordered_set<NodeID> S;
        std::unordered_set<NodeID> terms;
        std::unordered_map<NodeID, NodeID> gMapping;

        for (auto p : problem->terminals) {
            Q.push(p.position);
            S.insert(p.position);
            terms.insert(p.position);
        }

        while (S.size() < UNDEFINED_NODE && !Q.empty()) {
            NodeID f = Q.front();
            Q.pop();
            for (EdgeID e : problem->graph->edges_of(f)) {
                NodeID tgt = problem->graph->getEdgeTarget(f, e);
                if (S.count(tgt) == 0) {
                    S.insert(tgt);
                    Q.push(tgt);
                }
            }
        }

        g->start_construction(S.size());

        size_t gMapIndex = 0;
        for (NodeID n : problem->graph->nodes()) {
            if (S.count(n) > 0) {
                if (gMapping.count(n) == 0) {
                    gMapping.emplace(n, gMapIndex);
                    gMapIndex++;
                }

                if (terms.count(n) > 0) {
                    //LOG1 << "terminal in " << gMapping[n];
                }

                for (EdgeID e : problem->graph->edges_of(n)) {
                    auto [t, w] = problem->graph->getEdge(n, e);
                    if (t > n && S.count(t) > 0) {
                        if (gMapping.count(t) == 0) {
                            gMapping.emplace(t, gMapIndex);
                            gMapIndex++;
                        }

                        g->new_edge_order(gMapping[n], gMapping[t], w);
                    }
                }
            }
        }

        if (gMapIndex < problem->graph->n()) {
            //LOG1 << "Warning: Graph has more than 5000 vertices, only printing "
             //    << gMapIndex << " out of " << problem->graph->n();
        }

        //LOG1 << "Writing...";
        graph_io::writeGraphWeighted(g->to_graph_access(), path);
        //LOG1 << "...done!";
    }

    std::shared_ptr<mutable_graph>                      graph;
    std::vector<terminal>                               terminals;
    std::vector<std::shared_ptr<std::vector<NodeID> > > mappings;
    FlowType                                            lower_bound;
    FlowType                                            upper_bound;
    EdgeWeight                                          deleted_weight;
    std::pair<NodeID, EdgeID>                           priority_edge;
    std::unordered_set<NodeID>                          finished_blockpairs;
};

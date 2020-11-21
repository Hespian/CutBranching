#include <vector>
#include <queue>

using namespace std;

typedef int NodeID;
typedef int FlowType;


class push_relabel {

const int    WORK_OP_RELABEL    = 9;
const double GLOBAL_UPDATE_FRQ  = 0.51;
const int    WORK_NODE_TO_EDGES = 4;

public:
        push_relabel(std::vector<std::vector<int>> & adjj, std::vector<int> & x) 
            : adj( adjj), active(x)
        {
            flow.resize(adj.size());

            for (int i = 0; i < adj.size(); i++) {
                flow[i].resize(adj[i].size(), 0);
            }
        };

        void init(int rn, NodeID source, NodeID sink ) {
                this->rn = rn;
                n = adj.size();
                m = 0;

                for (NodeID node = 0; node < adj.size(); ++node)
                    if (isActive(node))
                    {
                        for (int const neighbor : adj[node])
                            if (isActive(neighbor))
                                m++;
                    }

                m_excess.resize(n,0);
                m_distance.resize(n,0);
                m_active.resize(n, false);
                m_count.resize(2*rn,0);
                m_bfstouched.resize(n);

                m_count[0] = rn-1;
                m_count[rn] = 1;
                
                m_distance[source] = rn;
                m_active[source]   = true;
                m_active[sink]     = true;

                for (int i = 0; i < adj[source].size(); i++)
                {
                    int u = adj[source][i];
                    if (isActive(u)) {
                        m_excess[source] += getEdgeCapacity(source, u);
                        push(source, u, i);
                    }
                }    
        }

        // perform a backward bfs in the residual starting at the sink
        // to update distance labels
        void global_relabeling( NodeID source, NodeID sink ) {
                std::queue< NodeID > Q;                 
                for (NodeID node = 0; node < adj.size(); node++)
                {
                    m_distance[node]   = std::max(m_distance[node], rn);
                    m_bfstouched[node] = false;
                }
                Q.push(sink);
                m_bfstouched[sink]   = true;
                m_bfstouched[source] = true;
                m_distance[sink]     = 0;

                NodeID node = 0;
                while( !Q.empty() ) {
                        node = Q.front();
                        Q.pop();  

                        for (int i = 0; i < adj[node].size(); i++)
                        {
                            NodeID target = adj[node][i];
                            if (isActive(target)) 
                            {
                                if(m_bfstouched[target]) continue;
                                int reverseIndex = getReverseIndex(target, node);
                                if (getEdgeCapacity(target, node) - getEdgeFlow(target, reverseIndex) > 0) 
                                {
                                    m_count[ m_distance[target] ] --;
                                    m_distance[target] = m_distance[node]+1;
                                    m_count[ m_distance[target] ] ++;
                                    Q.push(target);
                                    m_bfstouched[target] = true;
                                }
                            }
                        }
                }
        }

        // push flow from source to target if possible
        void push( NodeID source, NodeID target, int targetIndex) {
                m_pushes++;
                FlowType capacity = getEdgeCapacity(source, target);
                FlowType flow     = getEdgeFlow(source, targetIndex);
                FlowType amount   = std::min((int)(capacity-flow), m_excess[source]);

                if( m_distance[source] <= m_distance[target] || amount == 0) return;

                setEdgeFlow(source, targetIndex, flow+amount);

                int rev_index = getReverseIndex(target, source);
                FlowType rev_flow = getEdgeFlow(target, rev_index);
                setEdgeFlow(target, rev_index, rev_flow-amount);

                m_excess[source] -= amount;
                m_excess[target] += amount;

                enqueue(target);
        }

        // put a vertex in the FIFO queue of the algorithm 
        void enqueue( NodeID target ) {
                if( m_active[target] ) return;
                if( m_excess[target] > 0) {
                        m_active[target] = true;
                        //m_Q.push(target, m_distance[target]);
                        m_Q.push(target);
                }
        }

        // try to push as much excess as possible out of the node node
        void discharge( NodeID node ) {

                for (int i = 0; i < adj[node].size() && m_excess[node] > 0; i++)
                {
                    NodeID target = adj[node][i];
                    if (isActive(target))
                        push(node, target, i);
                }

                if( m_excess[node] > 0 ) {
                        if( m_count[ m_distance[node] ] == 1 && m_distance[node] < rn) {
                                // hence this layer will be empty after the relabel step
                                gap_heuristic(m_distance[node]);
                        } else {
                                relabel(node);
                        }
                }
        }


        // gap heuristic
        void gap_heuristic( NodeID level ) {
                m_gaps++;
                for (NodeID node = 0; node < adj.size(); node++)
                {
                    if (isActive(node))
                    {                        
                        if(m_distance[node] < level ) continue;
                        m_count[m_distance[node]]--;
                        m_distance[node] = std::max(m_distance[node], rn);
                        m_count[m_distance[node]]++;
                        enqueue(node);
                    }
                }
        }

        // relabel a node with respect to its 
        // neighboring nodes
        void relabel( NodeID node ) {
                m_work += WORK_OP_RELABEL;
                m_num_relabels++;

                m_count[m_distance[node]]--;
                m_distance[node] = 2*rn;
              

                for (int i = 0; i < adj[node].size(); i++)
                {
                    NodeID target = adj[node][i];
                    if (isActive(target))
                    {
                        if (getEdgeCapacity(node, target) - getEdgeFlow(node, i) > 0)
                            m_distance[node] = std::min(m_distance[node], m_distance[target]+1);
                        m_work++;
                    }
                }

                m_count[m_distance[node]]++;
                enqueue(node);
        }

        FlowType solve_max_flow_min_cut(int rn, 
                                         NodeID source, 
                                         NodeID sink, 
                                         bool compute_source_set, 
                                         std::vector<NodeID> & source_set,
                                         bool compute_cut, 
                                         std::vector<std::pair<int, int>> & cut) {
                m_work               = 0;
                m_num_relabels       = 0;
                m_gaps               = 0;
                m_pushes             = 0;
                m_global_updates     = 1;


                flow.resize(0);
                flow.resize(adj.size());


                for (int i = 0; i < adj.size(); i++) {
                    flow[i].resize(adj[i].size(), 0);
                }

                init(rn, source, sink);
                global_relabeling( source, sink );
         
                int work_todo = WORK_NODE_TO_EDGES*rn + m;
                // main loop
                while(!m_Q.empty()) {
                        NodeID v    = m_Q.front(); m_Q.pop();
                        m_active[v] = false;
                        discharge(v);
                        if( m_work > GLOBAL_UPDATE_FRQ*work_todo) {
                                global_relabeling( source,  sink );
                                m_work = 0;
                                m_global_updates++;
                        }
                }

                std::vector<bool> n_r(n, false);

                if(compute_source_set) {
                        // perform bfs starting from source set 
                        source_set.clear();

                        for (NodeID node = 0; node < n; node++)
                                m_bfstouched[node] = false;

                        std::queue< NodeID > Q;
                        Q.push(source);
                        m_bfstouched[source] = true;

                        while( !Q.empty() ) {
                                NodeID node = Q.front();
                                Q.pop();
                                source_set.push_back(node);

                                for (int i = 0; i < adj[node].size(); i++)
                                {
                                    NodeID target = adj[node][i];
                                    if (isActive(target))
                                    {
                                        FlowType resCap = getEdgeCapacity(node, target) - getEdgeFlow(node, i);
                                        if (resCap > 0 && !m_bfstouched[target])
                                        {
                                            Q.push(target);
                                            m_bfstouched[target] = true;
                                            n_r[target] = false;
                                        }
                                        else if (resCap == 0 && !m_bfstouched[target])
                                        {
                                            n_r[target] = true;
                                        }
                                    }
                                }
                        }
                }
                if (compute_cut)
                {
                    for (NodeID node = 0; node < n; node++)
                    {
                        if (n_r[node])
                        {
                            for (NodeID neigh : adj[node])
                            {
                                if (m_bfstouched[neigh])
                                {
                                    cut.emplace_back(neigh, node);
                                }
                            }
                        }
                    }
                }

                //return value of flow
                return m_excess[sink];
        }
private:
        std::vector<int> m_excess;
        std::vector<NodeID>    m_distance;
	    std::vector<bool>      m_active; // store which nodes are in the queue already
	    std::vector<int>       m_count;
        std::queue<NodeID>     m_Q;
	    std::vector<bool>      m_bfstouched; 
        //highest_label_queue    m_Q;
        int m_num_relabels;
        int m_gaps;
        int m_global_updates;
        int m_pushes;
        int m_work;

        // graph
        int rn;
        int n;
        int m;
        std::vector<std::vector<FlowType>> flow;
        std::vector<std::vector<int>> & adj;

        std::vector<int> & active;

        inline int getReverseIndex(NodeID source, NodeID target)
        {
            for (int i = 0; i < adj[source].size(); i++)
                if (adj[source][i] == target)
                    return i;
        }

        inline int getEdgeCapacity(NodeID v, NodeID u) 
        {
            return 1;
        }

        inline int getEdgeFlow(NodeID source, int targetIndex)
        {
            return flow[source][targetIndex];
        }

        inline void setEdgeFlow(NodeID source, int targetIndex, int amount)
        {
            flow[source][targetIndex] = amount;
        }

        inline bool isActive(NodeID node) 
        {
            return active[node] < 0;
        }
};
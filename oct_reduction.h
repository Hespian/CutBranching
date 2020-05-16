#include <vector>

std::vector<std::vector<int>> oct_to_vc(std::vector<std::vector<int>>& adjj)
{
    int n = adjj.size();
    std::vector<std::vector<int>>adj (n * 2);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < adjj[i].size(); j++)
        {
            adj[i].push_back(adjj[i][j]);
            adj[i+n].push_back(adjj[i][j] + n);
        }
        adj[i].push_back(i+n);
        adj[i+n].push_back(i);
    }

    for (int i = 0; i < adj.size(); ++i)
    {
        std::sort(adj[i].begin(), adj[i].end());
    }
    
    return adj;
}
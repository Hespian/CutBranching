#include <vector>

bool validate_solution(std::vector<bool> &sol, std::vector<std::vector<int>> &adj)
{
    for (int i = 0; i < sol.size(); i++)
    {
        if (!sol[i]) continue;
        for (int j = 0; j < adj[i].size(); j++)
        {
            if (sol[adj[i][j]]) return false;
        }
    }
    return true;
}
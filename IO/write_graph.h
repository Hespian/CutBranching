#include <string>
#include <vector>
#include <algorithm>
#include <iostream>


void writeGraph(std::vector<std::vector<int>>& adj, std::ostream &f) 
{
    for (int i = 0; i < adj.size(); i++)
    {
        for(int j : adj[i])
            f << i + 1 << " " << j + 1 << std::endl;
    }

}


void writeGraphToFile(std::vector<std::vector<int>>& adj, std::string filename)
{
    std::ofstream f(filename.c_str(), std::ofstream::app);
    writeGraph(adj, f) ;
    f.close();
}
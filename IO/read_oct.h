#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <algorithm>


std::vector<std::vector<int>> readOctGraph(std::istream &infile)
{
    std::string line;
    while (getline(infile, line))
    {
        if (line.compare(0, 14, "# Vertex names") == 0)
            break;
    }

    int next_index = 0;
    std::map<std::string, int> node_map;
    std::vector<std::vector<int>> graph(0);

    // skip vertices
    while (getline(infile, line))
    {
        if (line.empty())
            continue;

        if (line.compare(0, 7, "# Edges") == 0)
            break;

    }

    while (getline(infile, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (line.compare(0, 5, "# EOF") == 0)
            break;

        std::istringstream iss(line);

        std::string u, v;
        int uu, vv ;
        iss >> u >> v;

        std::map<std::string,int>::iterator itu;
        std::map<std::string,int>::iterator itv;

        itu = node_map.find(u);
        itv = node_map.find(v);

        if (itu == node_map.end())
        {
            uu = next_index++;
            node_map[u] = uu;
            graph.push_back(std::vector<int>(0));
        }
        else uu = itu->second;

        if (itv == node_map.end())
        {
            vv = next_index++;
            node_map[v] = vv;
            graph.push_back(std::vector<int>(0));
        }
        else vv = itv->second;


        if (u == v)
        {
            continue;
        }

        if (std::find(graph[uu].begin(), graph[uu].end(), vv) == graph[uu].end())
        {
            graph[uu].push_back(vv);
        }

        if (std::find(graph[vv].begin(), graph[vv].end(), uu) == graph[vv].end())
        {
            graph[vv].push_back(uu);
        }
    }

    for (int i = 0; i < graph.size(); ++i)
    {
        std::sort(graph[i].begin(), graph[i].end());
    }

    return graph;
}

std::vector<std::vector<int>> readOctGraphFromFile(std::string fileName)
{
    std::ifstream fin(fileName);
    return readOctGraph(fin);
}

std::vector<std::vector<int>> readOctGraphFromCin()
{
    return readOctGraph(std::cin);
}
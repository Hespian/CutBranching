#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <algorithm>


std::vector<std::vector<int>> readWebGraph(std::istream &infile)
{
    int numVertices;
    int numEdges;
    std::string line;
    while (getline(infile, line))
    {
        std::istringstream iss(line);
        char firstSymbol;
        if (!(iss >> firstSymbol))
        {
            break;
        } // error

        if (firstSymbol == 'p')
        {
            std::string td;
            iss >> td >> numVertices >> numEdges;
            break;
        }
    }

    int next_index = 0;
    std::map<int, int> node_map;
    std::vector<std::vector<int>> graph(0);

    while (getline(infile, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (line[0] == 'c')
        {
            continue;
        }

        std::istringstream iss(line);

        int u, v, uu, vv ;
        iss >> u >> v;

        std::map<int,int>::iterator itu;
        std::map<int,int>::iterator itv;

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

std::vector<std::vector<int>> readWebGraphFromFile(std::string fileName)
{
    std::ifstream fin(fileName);
    return readWebGraph(fin);
}

std::vector<std::vector<int>> readWebGraphFromCin()
{
    return readWebGraph(std::cin);
}
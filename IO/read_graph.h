#include <iostream>
#include <vector>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <algorithm>


std::vector<std::vector<int>> readGraph(std::istream &infile)
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
            iss >> numVertices >> numEdges;
            break;
        }
    }

    std::vector<std::vector<int>> graph(numVertices);

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

        int u, v;
        iss >> u >> v;

        if (u < 0 || u >= numVertices)
        {
            std::cout << "Invalid node ID: " << u << std::endl;
            exit(1);
        }
        if (v < 0 || v >= numVertices)
        {
            std::cout << "Invalid node ID: " << v << std::endl;
            exit(1);
        }

        if (u == v)
        {
            continue;
        }

        graph[u].push_back(v);
    }

    for (int i = 0; i < numVertices; ++i)
    {
        std::sort(graph[i].begin(), graph[i].end());
    }

    return graph;
}

std::vector<std::vector<int>> readGraphFromFile(std::string fileName)
{
    std::ifstream fin(fileName);
    return readGraph(fin);
}

std::vector<std::vector<int>> readGraphFromCin()
{
    return readGraph(std::cin);
}
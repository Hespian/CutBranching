#include <iostream>
#include <vector>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <algorithm>


std::vector<NodeID> readBO(std::istream &infile)
{
    int numVertices;
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
            iss >> td >> numVertices;
            break;
        }
    }

    std::vector<NodeID> bo(numVertices);
    int crnt = 0;

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

        int u;
        iss >> u;
        
        bo[crnt] = u;
        crnt++;
    }

    return bo;
}

std::vector<NodeID> readBOFromFile(std::string fileName)
{
    std::ifstream fin(fileName);
    return readBO(fin);
}

std::vector<NodeID> readBOFromCin()
{
    return readBO(std::cin);
}
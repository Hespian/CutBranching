#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

void writeBranchOrderToFile(std::vector<int>& order, std::string filename)
{
    std::ofstream f(filename.c_str(), std::ofstream::app);

    for (int v : order)
        f << v << std::endl;

    f.close();
}

void writeBranchOrderLevelToFile(std::vector<std::vector<int>>& order, std::string filename)
{
    std::ofstream f(filename.c_str(), std::ofstream::app);

    for (std::vector<int> v : order)
        f << v[0] << " " << v[1] << std::endl;

    f.close();
}
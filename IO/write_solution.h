#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

struct Result {

public:

    std::string instance;
    long nBranchings;
    double time;
    long misSize;

    long nDefaultBranchings = -1;
    long nDefaultPicks = -1;
    long nStratPicks = -1;
    long nDecomps = -1;  
    long maxDepth = -1;

    Result(std::string i, long nBrnch, double t, long mis) : instance {i}, nBranchings {nBrnch}, time {t}, misSize {mis}
    {
    }

};


void writeResult(Result& sol, std::ostream &f) {
    f << "Instance: " << sol.instance
      << " branches: " << sol.nBranchings 
      << " time: " << sol.time 
      << " MIS size: " << sol.misSize
      << " nDefPicks: " << sol.nDefaultPicks
      << " nStratPicks: " << sol.nStratPicks
      << " nDefBranchings: " << sol.nDefaultBranchings
      << std::endl;
}


void writeResultToFile(Result& sol, std::string filename) {
    std::ofstream f(filename.c_str(), std::ofstream::app);
    writeResult(sol, f) ;
    f.close();
}

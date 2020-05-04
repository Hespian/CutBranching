#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

class Result {

public:

    std::string instance;
    long nBranchings;
    double time;

    long mis_size;
    long nDefaultBranchings;

    Result(std::string i, long nBrnch, double t, long mis, long nDB) : instance {i}, nBranchings {nBrnch}, time {t}, mis_size {mis}, nDefaultBranchings {nDB}
    {
    }

    ~Result() {}
};


void writeResult(Result& sol, std::ostream &f) {
    f << "Inst: " << sol.instance
      << " nBranch: " << sol.nBranchings 
      << " time: " << sol.time 
      << " MIS size: " << sol.mis_size
      << " nDefaultBranchings: " << sol.nDefaultBranchings
      << std::endl;
}


void writeResultToFile(Result& sol, std::string filename) {
    std::ofstream f(filename.c_str(), std::ofstream::app);
    writeResult(sol, f) ;
    f.close();
}
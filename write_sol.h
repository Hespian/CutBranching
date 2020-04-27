#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

class Result {

public:

    int instanceNr;
    int nBranchings;
    double time;

    Result(int iNr, int nBrnch, double t) : instanceNr {iNr}, nBranchings {nBrnch}, time{t}
    {
    }

    ~Result() {}
};


void writeResult(Result& sol, std::ostream &f) {
    f << "Inst: " << sol.instanceNr
      << " nBranch: " << sol.nBranchings 
      << " time: " << sol.time 
      << std::endl;
}


void writeResultToFile(Result& sol, std::string filename) {
    std::ofstream f(filename.c_str(), std::ofstream::app);
    writeResult(sol, f) ;
    f.close();
}
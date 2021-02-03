#include <vector>
#include <string>
#include <filesystem>

#include "algo/branch_and_reduce_algorithm.h"
#include "algo/timer.h"

#include "IO/write_solution.h"
#include "IO/read_graph.h"

int branching_strat = 2;

long tuningParam1;
double tuningParam2;
long tuningParam3;

void setParams(int argc)
{
    if (argc < 5)
    {
        if (branching_strat == 4)
            tuningParam1 = 25;
        else if (branching_strat == 5)
            tuningParam1 = 50;
        else if (branching_strat >= 6)
            tuningParam1 = 2;
    }

    if (argc < 6)
    {
        if (branching_strat == 4)
            tuningParam2 = 0.1;
        else if (branching_strat == 5)
            tuningParam2 = 0.4;
    }

    if (argc < 7)
    {
        if (branching_strat == 4)
            tuningParam3 = 10;
        else if (branching_strat == 5)
            tuningParam3 = 3;
    }

    branch_and_reduce_algorithm::BRANCHING = branching_strat;
    branch_and_reduce_algorithm::TUNING_PARAM1 = tuningParam1;
    branch_and_reduce_algorithm::TUNING_PARAM2 = tuningParam2;
    branch_and_reduce_algorithm::TUNING_PARAM3 = tuningParam3;   
}

int main(int argc, char **argv)
{
    if (argc < 3)
        exit(-1);

    string instances_path = argv[1];
    string out_path = argv[2];

    if (argc > 3)
        branching_strat = atoi(argv[3]);
    if (argc > 4)
        tuningParam1 = atoi(argv[4]);
    if (argc > 5)
        tuningParam2 = atof(argv[5]);
    if (argc > 6)
        tuningParam3 = atoi(argv[6]); 

    setParams(argc);


    out_path += "/resuts_" + to_string(branching_strat) + "_" 
                          + to_string(tuningParam1) + "_" 
                          + to_string(tuningParam2) + "_" 
                          + to_string(tuningParam3) + ".txt";

    #ifdef USE_IFC
    std::cout << "HIER" << endl;
    #endif
    
    for (const auto &entry : std::filesystem::directory_iterator(instances_path))
    {
        if (!entry.is_directory())
        {
            std::vector<std::vector<int>> adj = readGraphFromFile(entry.path());
             
            branch_and_reduce_algorithm::resetStatistics();
            branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());
            timer t;
            t.restart();

            int vcSize = algo.solve(t, 6000);
            double secs = t.elapsed();

            Result res(entry.path().filename(), algo.nBranchings, secs, vcSize);
            // Timeout
            if (vcSize == -1)
            {
                res.nBranchings = -1;
                res.time = -1;
            }
            else
                res.misSize = adj.size() - vcSize;

            writeResultToFile(res, out_path);
        }
    }

    return 0;
}
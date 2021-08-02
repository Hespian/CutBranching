#include <vector>
#include <string>
#include <filesystem>

#include "algo/branch_and_reduce_algorithm.h"
#include "algo/timer.h"

#include "algo/tools/debug_info_logger.h"

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

    
    string out_path_log = out_path + "/log_" + to_string(branching_strat) + "_" 
                          + to_string(tuningParam1) + "_" 
                          + to_string(tuningParam2) + "_" 
                          + to_string(tuningParam3) + ".txt";

    out_path += "/results_" + to_string(branching_strat) + "_" 
                          + to_string(tuningParam1) + "_" 
                          + to_string(tuningParam2) + "_" 
                          + to_string(tuningParam3) + ".txt";
                          
    branch_and_reduce_algorithm::USE_DEPENDENCY_CHECKING = true;
    
    if (branch_and_reduce_algorithm::BRANCHING == 2 && branch_and_reduce_algorithm::TUNING_PARAM3 == 1) 
    {
        branch_and_reduce_algorithm::USE_DEPENDENCY_CHECKING = false;
    }
    
    
    for (const auto &entry : std::filesystem::directory_iterator(instances_path))
    {
        if (!entry.is_directory())
        {
            std::cout << entry.path() << std::endl;
            std::vector<std::vector<int>> adj = readGraphFromFile(entry.path());
            int N = adj.size();
            int M = 0;
            for (auto entr: adj)
                M += entr.size();

            debug_info_logger logger(entry.path().filename(), N, M/2);
             
            branch_and_reduce_algorithm::resetStatistics();
            branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());
            algo.logger = &logger;
            timer t;
            std::cout << "start" << std::endl;
            t.restart();

            int vcSize = algo.solve(t, 28800);
            double secs = t.elapsed();

            Result res(entry.path().filename(), algo.nBranchings, secs, N - vcSize);
            // Timeout
            if (vcSize == -1)
            {
                res.nBranchings = -1;
                res.time = -1;
            }
            else
                res.misSize = N - vcSize;

            writeResultToFile(res, out_path);
            logger.write_log(out_path_log);
        }
    }

    return 0;
}

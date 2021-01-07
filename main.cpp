#include <vector>
#include <string>
#include <algorithm>

#include "algo/branch_and_reduce_algorithm.h"
#include "algo/timer.h"

#include "IO/read_pace.h"
#include "IO/read_web_graph.h"
#include "IO/read_dimacs.h"
#include "IO/write_sol.h"
#include "IO/read_oct.h"

#include "validate_sol.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 5)
        exit(0);

    branch_and_reduce_algorithm::BRANCHING = atoi(argv[1]);
    branch_and_reduce_algorithm::EXTRA_DECOMP = atoi(argv[2]);
    branch_and_reduce_algorithm::TUNING_PARAM1 = atoi(argv[3]);
    branch_and_reduce_algorithm::TUNING_PARAM2 = atof(argv[4]);
    branch_and_reduce_algorithm::TUNING_PARAM3 = atoi(argv[5]);    

    std::string of = "../Benchmarks/" + to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_" + "tuning" + "_" + to_string(branch_and_reduce_algorithm::TUNING_PARAM1) + "::" + to_string(branch_and_reduce_algorithm::TUNING_PARAM2) + "::" + to_string(branch_and_reduce_algorithm::TUNING_PARAM3) + ".txt";

    int instances[] = {38,
                       41,42,43,44,45,
                       51,52,53,54,55,
                       61,62,63,64,65,
                       71,72,73,74,75,
                       36,37,38,39,40,
                       46,47,48,49,50,
                       56,57,58,59,60,
                       66,67,68,69,70, };

    double ts = 0;
    for (int i : instances)
    {
        std::string fn = "../demo_instances/vc-exact_";
        if (i < 10)
            fn += "00";
        else if (i < 100)
            fn += "0";

        std::vector<std::vector<int>> adj = readPaceGraphFromFile(fn + to_string(i) + ".gr");
        cout << adj.size() << endl;

        branch_and_reduce_algorithm::resetStatistics();
        branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());
        timer t;
        t.restart();
        int sol = algo.solve(t, 1800);
        double secs = t.elapsed();
        ts += secs;

        Result res(to_string(i), algo.nBranchings, secs, sol, algo.defaultBranchings);
        res.nDefaultPicks = branch_and_reduce_algorithm::defaultPicks;
        res.nStratPicks = branch_and_reduce_algorithm::stratPicks;
        res.nDecomps = branch_and_reduce_algorithm::nDecomps;
        res.max_depth = algo.max_depth;

        if (sol == -1)
        {
            res.time = -1;
        }

        std::vector<bool> solution(algo.adj.size(), false);
        algo.get_solved_is(solution);
        if (!validate_solution(solution, algo.adj))
        {
            res.instance += "--failed";
        }

        writeResultToFile(res, of);
    }

    return 0;
}
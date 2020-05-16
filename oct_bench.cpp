#include <vector>
#include <string>
#include <algorithm>

#include "algo/branch_and_reduce_algorithm.h"
#include "VieCut/lib/tools/timer.h"
#include "VieCut/lib/data_structure/graph_access.h"

#include "read_pace.h"
#include "read_dimacs.h"
#include "read_web_graph.h"
#include "read_oct.h"
#include "write_sol.h"
#include "validate_sol.h"
#include "oct_reduction.h"

int main(int argc, char **argv)
{
    if (argc < 5)
        exit(-1);

    branch_and_reduce_algorithm::BRANCHING = atoi(argv[1]);
    branch_and_reduce_algorithm::EXTRA_DECOMP = atoi(argv[2]);
    branch_and_reduce_algorithm::ND_LEVEL = atoi(argv[4]);
    int start = atoi(argv[3]);
    int end = atoi(argv[4]) + 1;

    std::string basePath = "../demo_instances/OCT/afro-americans/";
    std::string of = "../Benchmarks/";
    std::string conf = to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_";

    timer tt;
    tt.restart();
    for (int i = 10; i < 55; i++)
    {
        std::string path  = basePath + to_string(i) + ".graph";
        std::vector<std::vector<int>> adjj = readOctGraphFromFile(path);
        std::vector<std::vector<int>> adj = oct_to_vc(adjj);
              
        branch_and_reduce_algorithm::resetStatistics();
        branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());

        timer t;
        t.restart();
        int sol = algo.solve(t, 14400);
        double secs = t.elapsed();

        Result res("afro-" + to_string(i), algo.nBranchings, secs, sol, algo.defaultBranchings);
        res.nDefaultPicks = branch_and_reduce_algorithm::defaultPicks;
        res.nStratPicks = branch_and_reduce_algorithm::stratPicks;
        res.nDecomps = branch_and_reduce_algorithm::nDecomps;

        if (sol == -1)
        {
            res.nBranchings = -1;
            res.time = -1;
        }

        std::vector<bool> solution(algo.adj.size(), false);
        algo.get_solved_is(solution);
        if (!validate_solution(solution, algo.adj))
        {
            res.instance += "--failed";
        }

        writeResultToFile(res, of + conf + "_afro.txt");
    }

    double secs = tt.elapsed();
    cout << "ready" << endl;
    cout << secs << endl;

    return 0;
}
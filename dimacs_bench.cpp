#include <vector>
#include <string>
#include <algorithm>

#include "algo/branch_and_reduce_algorithm.h"
#include "VieCut/lib/tools/timer.h"
#include "VieCut/lib/data_structure/graph_access.h"

#include "read_pace.h"
#include "read_dimacs.h"
#include "read_web_graph.h"
#include "write_sol.h"

#include "validate_sol.h"

int main(int argc, char **argv)
{
    if (argc < 5)
        exit(-1);

    branch_and_reduce_algorithm::BRANCHING = atoi(argv[1]);
    branch_and_reduce_algorithm::EXTRA_DECOMP = atoi(argv[2]);
    branch_and_reduce_algorithm::ND_LEVEL = atoi(argv[4]);
    int start = atoi(argv[3]);
    int end = atoi(argv[4]) + 1;

    std::string basePath = "../demo_instances/DIMACS/";
    std::string of = "../Benchmarks/";
    std::string conf = to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_";

    std::string instances[] = {"brock200_1.clq-compliment.txt", "brock200_2.clq-compliment.txt",
                               "brock200_3.clq-compliment.txt", "brock200_4.clq-compliment.txt",
                               "gen200_p0.9_44.clq-compliment.txt", "gen200_p0.9_55.clq-compliment.txt",
                               "hamming6-4.clq-compliment.txt", "hamming8-4.clq-compliment.txt",
                               "johnson8-2-4.clq-compliment.txt", "johnson8-4-4.clq-compliment.txt", "johnson16-2-4.clq-compliment.txt",
                               "p_hat300-1.clq-compliment.txt", "p_hat300-2.clq-compliment.txt", "p_hat500-1.clq-compliment.txt", "p_hat500-2.clq-compliment.txt"};

    for (auto str : instances)
    {
        std::vector<std::vector<int>> adj = readDimacsGraphFromFile(basePath + str);
       
        branch_and_reduce_algorithm::resetStatistics();
        branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());

        timer t;
        t.restart();
        int sol = algo.solve(t, 14400);
        double secs = t.elapsed();

        Result res(str, algo.nBranchings, secs, sol, algo.defaultBranchings);
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

        writeResultToFile(res, of + conf + str);
    }

    cout << "ready" << endl;

    return 0;
}
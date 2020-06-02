#include <vector>
#include <string>
#include <algorithm>

#include "read_pace.h"
#include "read_web_graph.h"
#include "read_dimacs.h"
#include "write_sol.h"
#include "algo/branch_and_reduce_algorithm.h"
#include "timer.h"
#include "data_structure/graph_access.h"
//#include "Metis/include/metis.h"

#include "validate_sol.h"
#include "read_branch_order.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 6)
        exit(-1);

    branch_and_reduce_algorithm::BRANCHING = atoi(argv[1]);
    branch_and_reduce_algorithm::EXTRA_DECOMP = atoi(argv[2]);
    int start = atoi(argv[3]);
    int end = atoi(argv[4]) + 1;
    int mode = atoi(argv[5]);

    branch_and_reduce_algorithm::ND_LEVEL = 64;
    branch_and_reduce_algorithm::ro = start > 0 ? 1 : 0;

    if (mode == 1) // pace
    {
        std::string of = "../Benchmarks/" + to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_" + to_string(start) + "_" + to_string(end) + "-centr" + ".txt";
        for (int i = start; i < end; i++)
        {
            std::string fn = "../demo_instances/vc-exact_";
            if (i < 10)
                fn += "00";
            else if (i < 100)
                fn += "0";

            std::vector<std::vector<int>> adj = readPaceGraphFromFile(fn + to_string(i) + ".gr");
            std::vector<NodeID> bo = readBOFromFile("../branch_order/" + to_string(i) + "-order.txt");
            cout << i << endl;
            branch_and_reduce_algorithm::resetStatistics();
            branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());

            timer t;
            t.restart();
            // algo.reduce();

            // std::vector<NodeID> m, rm;
            // algo.convert_to_adj(adj, rm, m);
            // for (int i = 0; i < bo.size(); i++)
            //     bo[i] = rm[bo[i]];

            // algo.nd_computed = true;
            // algo.nd_order.swap(bo);

            int sol = algo.solve(t, 900);
            double secs = t.elapsed();

            Result res(to_string(i), algo.nBranchings, secs, sol, algo.defaultBranchings);
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

            writeResultToFile(res, of);

            cout << "nbranch: " << algo.nBranchings << endl;
            cout << sol << " found in: " << secs << " secs" << endl;
        }
    }
    else if (mode == 2) // sparse networks
    {
        std::string of = "../Benchmarks/";
        std::string conf = to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_" + to_string(branch_and_reduce_algorithm::ND_LEVEL);
        std::string basePath = "../demo_instances/web/";
        std::string instances[] = {"web-Stanford", "web-NotreDame", "web-BerkStan", "as-skitter", "out.libimseti"};
        for (auto str : instances)
        {
            std::vector<std::vector<int>> adj = readWebGraphFromFile(basePath + str + ".txt");
            branch_and_reduce_algorithm::resetStatistics();
            branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());
            std::vector<NodeID> bo = readBOFromFile("../branch_order/" + str + "-order.txt");
            cout << "Hier" << endl;
            timer t;
            t.restart();
            // algo.reduce();

            // std::vector<NodeID> m, rm;
            // algo.convert_to_adj(adj, rm, m);
            // for (int i = 0; i < bo.size(); i++)
            //     bo[i] = rm[bo[i]];

            // algo.nd_computed = true;
            // algo.nd_order.swap(bo);

            int sol = algo.solve(t, 86400);

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

            writeResultToFile(res, of + conf + str + "-centr.txt");
        }
    }
    else if (mode == 3) // dimacs
    {
        std::string of = "../Benchmarks/";
        std::string conf = to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP) + "_" + to_string(branch_and_reduce_algorithm::ND_LEVEL);
        std::string basePath = "../demo_instances/DIMACS/";

        std::string instances[] = {/*"brock200_1.clq-compliment",*/ "brock200_2.clq-compliment",
                                   "brock200_3.clq-compliment", "brock200_4.clq-compliment",
                                   "gen200_p0.9_44.clq-compliment", "gen200_p0.9_55.clq-compliment",
                                   "hamming6-4.clq-compliment", "hamming8-4.clq-compliment",
                                   "johnson8-2-4.clq-compliment", "johnson8-4-4.clq-compliment", "johnson16-2-4.clq-compliment",
                                   "p_hat300-1.clq-compliment", "p_hat300-2.clq-compliment", "p_hat500-1.clq-compliment",
                                   "p_hat500-2.clq-compliment"};

        for (auto str : instances)
        {
            std::vector<std::vector<int>> adj = readDimacsGraphFromFile(basePath + str + ".txt");
            branch_and_reduce_algorithm::resetStatistics();
            branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());
            std::vector<NodeID> bo = readBOFromFile("../branch_order/dimacs/" + str + "-order.txt");
            cout << "Hier" << endl;
            timer t;
            t.restart();
            // algo.reduce();

            // std::vector<NodeID> m, rm;
            // algo.convert_to_adj(adj, rm, m);
            // for (int i = 0; i < bo.size(); i++)
            //     bo[i] = rm[bo[i]];

            // algo.nd_computed = true;
            // algo.nd_order.swap(bo);

            int sol = algo.solve(t, 86400);

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

            writeResultToFile(res, of + conf + str + "-centr.txt");
        }
    }
    return 0;
}
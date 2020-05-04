#include <vector>
#include <string>
#include <algorithm>

#include "read_pace.h"
#include "write_sol.h"
#include "algo/branch_and_reduce_algorithm.h"
#include "VieCut/lib/tools/timer.h"
#include "VieCut/lib/data_structure/graph_access.h"
//#include "Metis/include/metis.h"

#include "validate_sol.h"

using namespace std;

std::vector<std::vector<int>> readPaceGraphFromFile(std::string fileName);

int main(int argc,  char** argv)
{
    if (argc < 5) exit(-1);

    branch_and_reduce_algorithm::BRANCHING = atoi(argv[1]);
    branch_and_reduce_algorithm::EXTRA_DECOMP = atoi(argv[2]);
    int start = atoi(argv[3]);
    int end = atoi(argv[4]) + 1;

    std::string of = "../Benchmarks/" + to_string(branch_and_reduce_algorithm::BRANCHING) + "_" + to_string(branch_and_reduce_algorithm::EXTRA_DECOMP)
                    + "_" + to_string(start) + "_" + to_string(end) + ".txt";

    for (int i = start; i < end; i++)
    {
        std::string fn = "../demo_instances/vc-exact_";
        if (i < 10) fn += "00";
        else if (i < 100) fn += "0";

        std::vector<std::vector<int>> adj = readPaceGraphFromFile(fn + to_string(i) + ".gr");
        int size = adj.size();
        cout << "Instance: " << i << " of size: " << adj.size() << endl;
        
        branch_and_reduce_algorithm::nBranchings = 0;
        branch_and_reduce_algorithm algo = branch_and_reduce_algorithm(adj, adj.size());

        timer t;
        t.restart();
        int sol = algo.solve(t, 900);
        double secs = t.elapsed();

        Result res(to_string(i), algo.nBranchings, secs, sol, algo.defaultBranchings);
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

    return 0;
}
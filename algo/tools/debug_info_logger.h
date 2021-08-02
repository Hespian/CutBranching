#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#pragma once

enum info_type { Decompose, Branching, Reductions };

struct debug_info {
    info_type type;
    int depth;
    int rn;

    std::vector<std::pair<string, int>> add_stats;
};


class debug_info_logger {

private:

    int num_nodes, num_edges;
    string graph_name;
    std::vector<debug_info> infos;

public:

    debug_info_logger(string name, int n, int m)  {
        graph_name = name;
        num_nodes = n;
        num_edges = m;
    }

    void add_info(debug_info info) {
        infos.push_back(info);
    }

    void write_log(std::string file_path) {
        std::ofstream f(file_path.c_str(), std::ofstream::app);
        f << "-------------" << std::endl;
        f << "###: " << graph_name << " N: " << num_nodes << " M: " << num_edges << std::endl;
        
        for (auto info : infos) {
            switch (info.type)
            {
                case Decompose:
                    f << "[Decompose]:";
                    break;
                case Reductions:
                    f << "[Reductions]:";
                    break;
                case Branching:
                    f << "[Branching]:";
                    break;

                default:
                    f << "[Info]:";
                    break;
            }

            f << " depth: " << info.depth;
            f << " num. rn: " << info.rn;

            for (auto add_stat : info.add_stats)
                f << " " << add_stat.first << to_string(add_stat.second);

            f << std::endl;
        }

        f.close();
    }
};
#include <iostream>
#include <fstream>
#include "Graph.h"
#include <ctime>
#include <string>
#include <array>

///
/// \brief init
/// \param in
/// \return
///
Graph* init(std::istream &in)
{
    char pr_type[4]; ////problem type;

    std::vector<int> edge_tails, edge_heads;
    std::vector<double> edge_weights;

    int n = 0, m = 0; //no of nodes, no of arcs;
    std::string line_inf;
    std::getline(in, line_inf);

    sscanf(line_inf.c_str(), "%*c %3s %d %d", pr_type, &n, &m);

    double en_weight = 0;
    double ex_weight = 0;

    auto *resG = new Graph(n, m, 0, n - 1, en_weight, ex_weight);
    int edges = 0;
    int edge_id = 0;
    std::cout << "Parsing edges: " << std::endl;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        switch (line[0]) {
            case 'c':                  /* skip lines with comments */
            case '\n':                 /* skip empty lines   */
            case 'n':
            case '\0':                 /* skip empty lines at the end of file */
                break;
            case 'p':
            case 'a':
        {
            int tail = 0;
            int head = 0;
            double weight = 0;
                sscanf(line.c_str(), "%*c %d %d %lf", &tail, &head, &weight);
                edges++;
                if (tail>n || head>n)
                    std::cout<<"ERROR: more nodes than expected!"<<std::endl;
                resG->add_edge(tail - 1, head - 1, edge_id, weight);
                edge_id++;
                if (edges % 10000 == 0)
                    std::cout << edges << std::endl;
                break;
            }
            default:
                break;
        }
    }
    std::cout << "Parsing done!" << std::endl;
    return resG;
}

///
/// \brief print_solution
///        output the min-cost flow results
/// \param resG
/// \param path_set
/// \param outputstream
///
void print_solution(Graph* resG, const std::vector<std::vector<int>>& path_set, std::ostream &outputstream)
{
    std::vector<bool> edge_visited_flag(resG->num_edges_);
    for (size_t i = 0; i < path_set.size(); i++) {
        for (size_t j = 0; j < path_set[i].size() - 1; j++) {
            int tail = path_set[i][j];
            int head = path_set[i][j + 1];
            int edge_idx = resG->node_id2edge_id[Graph::node_key(tail, head)];
            edge_visited_flag[edge_idx] = !edge_visited_flag[edge_idx];
        }
    }

    outputstream << "path,tail,head\n";
    int total_paths = (int) path_set.size();
    for (int i = 0; i < total_paths; i++) {
        for (size_t j = 0; j < path_set[i].size() - 1; j++) {
            int tail = path_set[i][j];
            int head = path_set[i][j+1];
            outputstream << i << "," << tail << "," << head << "\n";
        }
    }
}

///
/// \brief main
/// \param argc
/// \param argv
/// \return
///
int main(int argc, char *argv[])
{
    std::unique_ptr<Graph> org_graph;
    if (argc > 1) {
        std::ifstream fin(argv[1]);
        if (!fin) {
            std::cerr << "Error opening input file: " << argv[1] << "\n";
            return 1;
        }
        org_graph.reset(init(fin));
    }
    else {
        org_graph.reset(init(std::cin));
    }

    //// reading data
    clock_t t_start = clock();
    clock_t t_end = clock();
    long double parsing_time = t_end - t_start;

    std::array<long double, 10> duration;
    duration.fill(0);
    //// 1: remove dummy edges
    t_start = clock();
    org_graph->invalid_edge_rm();
    t_end = clock();
    duration[0] += t_end - t_start;

    ////save path and path cost
    std::vector<double> path_cost;
    std::vector<std::vector<int>> path_set;
    int path_num = 0;
    t_start = clock();
    //// 2: initialize shortest path tree from the DAG
    org_graph->shortest_path_dag();
    t_end = clock();
    duration[1] += t_end - t_start;

    path_cost.push_back(org_graph->distance2src[org_graph->sink_id_]);
    org_graph->cur_path_max_cost = -org_graph->distance2src[org_graph->sink_id_]; // the largest cost we can accept

    //// 3: convert edge cost (make all weights positive)
    t_start = clock();
    org_graph->update_allgraph_weights();
    t_end = clock();
    duration[2] += t_end - t_start;

    //// 8: extract shortest path
    t_start = clock();
    org_graph->extract_shortest_path();
    t_end = clock();
    duration[7] += t_end - t_start;

    path_set.push_back(org_graph->shortest_path);
    path_num++;

    std::vector<unsigned long> update_node_num;

    //// 4: find nodes for updating based on branch node
    std::vector<int> node_id4updating;
    t_start = clock();
    org_graph->find_node_set4update(node_id4updating);
    t_end = clock();
    duration[3] += t_end - t_start;

    //// 10: rebuild residual graph by flipping paths
    t_start = clock();
    org_graph->flip_path();//also erase the top sinker
    t_end = clock();
    duration[9] += t_end - t_start;
    while (true) {
        //// 6: update shortest path tree based on the selected sub-graph
        t_start = clock();
        org_graph->update_shortest_path_tree_recursive(node_id4updating);
        printf("Iteration #%d, updated node number  %ld \n", path_num, org_graph->upt_node_num);
        t_end = clock();
        duration[5] += t_end - t_start;

        //// 7: update sink node (heap)
        t_start = clock();
        org_graph->update_sink_info(node_id4updating);
        t_end = clock();
        duration[6] += t_end - t_start;

        update_node_num.push_back(node_id4updating.size());

        //// 8: extract shortest path
        t_start = clock();
        org_graph->extract_shortest_path();
        t_end = clock();
        duration[7] += t_end - t_start;

        // test if stop
        double cur_path_cost = path_cost[path_num - 1] + org_graph->distance2src[org_graph->sink_id_];

        if (cur_path_cost > -0.0000001) {
            break;
        }

        path_cost.push_back(cur_path_cost);
        org_graph->cur_path_max_cost = -cur_path_cost;
        path_set.push_back(org_graph->shortest_path);
        path_num++;

        //// 9: update weights
        t_start = clock();
        org_graph->update_subgraph_weights(node_id4updating);
        t_end = clock();
        duration[8] += t_end - t_start;

        //// 4: find nodes for updating
        t_start = clock();
        org_graph->find_node_set4update(node_id4updating);
        t_end = clock();
        duration[3] += t_end - t_start;
        //// 10: rebuild the graph
        t_start = clock();
        org_graph->flip_path();
        t_end = clock();
        duration[9] += t_end - t_start;
    }

    //// out put results and time consuming
    std::cout << "Parsing time is: " << parsing_time / CLOCKS_PER_SEC << " s" << std::endl;

    long double all_cpu_time = 0;
    for (int i = 0; i < 10; i++) {
        auto time_elapsed_ms = 1000.0 * duration[i] / CLOCKS_PER_SEC;
        all_cpu_time += time_elapsed_ms;
        std::cout << "the " << i + 1 << " step used: " << time_elapsed_ms / 1000.0 << " s\n";
    }
    std::cout << "The overall time is " << all_cpu_time / 1000.0 << " s\n\n";

    double cost_sum = 0;
    for (auto &&i : path_cost) {
        cost_sum += i;
    }
    printf("Number of paths: %ld, total cost is %.7f, final path cost is: %.7f.\n",
           path_cost.size(), cost_sum, path_cost[path_cost.size() - 1]);

    for (int i = 0; i < path_set.size(); i++) {
        std::cout << "  Path #" << i << " length: " << path_set[i].size() << std::endl;
    }

    /********* write paths to csv ********/
    if (argc == 3) {
        std::string outputFileName = argv[2];
        std::ofstream fout(outputFileName);
        if (!fout) {
            std::cerr << "Error opening output file: " << outputFileName << "\n";
            return 1;
        }
        print_solution(org_graph.get(), path_set, fout);
        fout.close();
    } else {
        std::cout << "\nSolution:\n" << std::endl;

        print_solution(org_graph.get(), path_set, std::cout);
    }

    return 0;
}

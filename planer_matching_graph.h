#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array> 
#include <random>

#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>


struct vertex_property {
    int                 name;
    bool                is_boundary;
    std::vector<int>    nodes;      // 統合したnodeのidをスタックしておく配列
};

struct edge_property {
    int                 fault_id;
    float               error_probability;
    double              weight;
    std::pair<int,int>  original_nodes;
    bool                erasure;
};


typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    vertex_property,        // 頂点のBundleプロパティ
    edge_property         // 辺のBundleプロパティ
> Matching_Graph;
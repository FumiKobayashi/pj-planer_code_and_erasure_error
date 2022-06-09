#include "planer_matching_graph.h"
#include <boost/python.hpp>

// #define _GLIBCXX_DEBUG

//add_edge用のboolを宣言
bool inserted = false;


// ---------- matching graphの作成関数 -----------------
Matching_Graph planer_code_matching_graph(const int L){

//グラフの宣言
Matching_Graph g;

const int N_V = L*(L-1)+1;          //vertexの数
const int N_E = L*L+(L-1)*(L-1);    //edgeの数

//vertexの配列の宣言
std::vector<Matching_Graph::vertex_descriptor> v(N_V);
//edgeの配列を宣言
std::vector<Matching_Graph::edge_descriptor> e(N_E);

//nodeをグラフに追加
for (int i=0; i< N_V; i++){
    v.at(i) = add_vertex(g);
    g[v.at(i)].name     = i;
    // g[v.at(i)].nodes.reserve(L*(L-1)/2);     高速化するならメモリを確保しといたほうがいいかも
    g[v.at(i)].nodes.push_back(i);
if (i!=N_V-1){
    g[ v.at(i)].is_boundary = false;
    }
else {
    g[v.at(N_V-1)].is_boundary = true;
    }
}

//edgeをグラフに追加
//下のboundary nodeに接続するedge
for (int i=0; i< L; i++){
    boost::tie(e.at(i), inserted) = add_edge(v.at(i), v.at(N_V-1), g);
    g[e.at(i)].fault_id = i;
    g[e.at(i)].error_probability =0.0;
    g[e.at(i)].weight = 0.0;
    g[e.at(i)].erasure = false;
}

//間の縦のedge
for (int i=L; i< L+ L*(L-2); i++){
    boost::tie(e.at(i), inserted) = add_edge(v.at(i-L), v.at(i), g);
    g[e.at(i)].fault_id = i;
    g[e.at(i)].error_probability =0.0;
    g[e.at(i)].weight = 0.0;
    g[e.at(i)].erasure = false;
}

//上のboundary nodeに接続するedge
for (int i=0; i< L; i++){
    boost::tie(e.at(i), inserted) = add_edge(v.at(N_V-1-L+i), v.at(N_V-1), g);
    g[e.at(i)].fault_id = L*L-L + i;
    g[e.at(i)].error_probability =0.0;
    g[e.at(i)].weight = 0.0;
    g[e.at(i)].erasure = false;
}

int id;
//間の横のedge
for (int j=0 ; j<L-1; j++){
    for (int i=0; i<L-1; i++){
        id = L*L+j*(L-1)+i;
        boost::tie(e.at(id), inserted) = add_edge(v.at(i+j*L), v.at(i+j*L+1), g);
        g[e.at(id)].fault_id = id;
        g[e.at(id)].error_probability =0.0;
        g[e.at(id)].weight = 0.0;
        g[e.at(id)].erasure = false;
}
}

return g;
}
// --------------------------------------------------------------------

// ------------- graph中の特定のe属性を持ったedegeの検索 -------------------
int find_specific_id_edge(Matching_Graph g, int num_id, Matching_Graph::edge_descriptor &e_discript){

    auto edge_range = edges(g);
    for (auto first = edge_range.first, last = edge_range.second; first != last; ++first) {
            Matching_Graph::edge_descriptor e = *first;
            if (g[e].fault_id == num_id){
                e_discript = e;
                break;
                }
            else{
             continue;
            }
            std::cout << 
            "ERROR:Number of the id do not exist in the graph" 
            << std::endl;
            return -1;
            } 
    return 0;
    }
// --------------------------------------------------------------------
// -----------Erasured_matching_graph_creater内で使う関数 ---------------
//頂点u,vに対して、u_id > v_id の時、uにvの辺を接合してvを消す関数
void Merge_vertex(Matching_Graph& g_erasure,Matching_Graph::vertex_descriptor u,Matching_Graph::vertex_descriptor v, int u_id, int v_id) 
{
    // 消去する頂点が今までどの頂点と合体してきたかの履歴の更新
    for (const auto& node : g_erasure[v].nodes){
                      g_erasure[u].nodes.push_back(node);
                    }
    // 新しいedgeを追加して、頂点を合体
    std::pair<Matching_Graph::adjacency_iterator, Matching_Graph::adjacency_iterator>
    neighbors = boost::adjacent_vertices(v, g_erasure);     //近隣のnodeを検索
    // vの近隣のu以外の頂点ごとに操作するためのfor文
    for (auto first = neighbors.first, last = neighbors.second; first != last; ++first){
        Matching_Graph::vertex_descriptor neighbor = *first;    //近隣の頂点のdescripter
        if (g_erasure[neighbor].name != u_id){              //u以外であることの検査
            // Matching_Graph::edge_descriptor neighbor_edge;
            //     boost::tie(neighbor_edge, inserted) = edge(v,  neighbor, g_erasure);
            std::pair<Matching_Graph::out_edge_iterator, Matching_Graph::out_edge_iterator>
                neighbor_edges = boost::out_edges(neighbor, g_erasure);     //vの近隣頂点から出るedgeのitrator
            
            std::vector<edge_property> edge_property_list;  //一時的にedgeを格納しておく配列
            for(auto first = neighbor_edges.first, last = neighbor_edges.second; first != last; ++first){
                if (boost::target(*first, g_erasure) == v){     //neighbor-vのedgeを抽出するためのif
                    // edgeのpropertyをメモリ
                    edge_property edge_datas;
                    edge_datas.fault_id             = g_erasure[*first].fault_id;
                    edge_datas.error_probability    = g_erasure[*first].error_probability;
                    edge_datas.weight               = g_erasure[*first].weight;
                    edge_property_list.push_back(edge_datas);
                }
                else continue;
            }
            //新しいedgeの追加
            for (auto first = edge_property_list.begin(), last  = edge_property_list.end()
                ;first!=last; ++first){
                Matching_Graph::edge_descriptor new_e;
                boost::tie(new_e, inserted) = add_edge(u, neighbor, g_erasure);
                edge_property neighbor_v_edge_property = *first;
                // neighbor-vのedgeのedge_propertyを新しいedgeに継承
                g_erasure[new_e].fault_id           = neighbor_v_edge_property.fault_id;
                g_erasure[new_e].error_probability  = neighbor_v_edge_property.error_probability;
                g_erasure[new_e].weight             = neighbor_v_edge_property.weight;
                g_erasure[new_e].erasure            = false; 
                }
        }
        else
            continue;
    }
    boost::clear_vertex(v, g_erasure);
}
// --------------------------------------------------------------------

// ---------------------- Erasure graphの作成関数 -----------------------
Matching_Graph Erasured_matching_graph_creater(Matching_Graph *g_matching, const int L, std::vector<int> lossqubits){
    // g_matchingのコピー
    Matching_Graph g_erasure;
    boost::copy_graph(*g_matching, g_erasure);

    // edgeにErasureエラー属性を付加
    auto first = std::begin(lossqubits) ;
    auto last = std::end(lossqubits) ;
    // auto edge_range = edges(g_erasure);

    for (const auto& loss_id : lossqubits){
        std::cout << "loss id: "<< loss_id << std::endl;
        auto edge_range = edges(g_erasure);
        for (auto first = edge_range.first, last = edge_range.second; first != last; ++first) {
            Matching_Graph::edge_descriptor e = *first;
            if (g_erasure[e].fault_id == loss_id){
                g_erasure[e].erasure = true;
                // loss edgeに接続している頂点の検索
                Matching_Graph::vertex_descriptor 
                    u = boost::source(e, g_erasure);
                    int u_id = g_erasure[u].name;
                Matching_Graph::vertex_descriptor 
                    v = boost::target(e, g_erasure);
                    int v_id = g_erasure[v].name;
                
                // 頂点のidが大きい方にidが小さい方のidをnodesに追加
                if (u_id > v_id){
                    Merge_vertex(g_erasure, u, v, u_id, v_id);              
                    }
                else if (u_id < v_id){
                    Merge_vertex(g_erasure, v, u, v_id, u_id);
                    }
                else{
                std::cout << 
                "ERROR:I cannot merge vertecies because node id is same."
                << std::endl;
                }
            break;
            }}}
        // Check_the_graph(g_erasure);
    return g_erasure;
    }    
// --------------------------------------------------------------------

// ------------------ edgeのダブりをカウントする関数 -----------------------
int number_of_edges(Matching_Graph::vertex_descriptor n, Matching_Graph::vertex_descriptor m, Matching_Graph g_erasure){
    int l;
    std::pair<Matching_Graph::out_edge_iterator, Matching_Graph::out_edge_iterator>
        edge_range = boost::out_edges(n, g_erasure);

    for(auto edge = edge_range.first, last = edge_range.second; edge != last; ++edge){
        if (target(*edge, g_erasure) == m)
            l++;
        else continue;
    }
    return l;
}


//上のカウントを出力する関数
std::map<int, int> get_number_of_edge_degenarate(Matching_Graph g_erasure){
    
    std::map<int, int>      degraded_edges;    // {fault_id : number of degraded}のmap(辞書型)

    std::pair<Matching_Graph::edge_iterator, Matching_Graph::edge_iterator> 
        edge_range = boost::edges(g_erasure);

    for (auto e = edge_range.first, last = edge_range.second; e != last; ++e) {
        Matching_Graph::vertex_descriptor n = source(*e, g_erasure);
        Matching_Graph::vertex_descriptor m = target(*e, g_erasure);
        int l        = number_of_edges(n, m, g_erasure);
        int fault_id = g_erasure[*e].fault_id;

        degraded_edges[fault_id] = l;
    }                       
    return degraded_edges;
}
// -------- loss_edgeの乱択関数 -----------------
std::vector<int> rand_lossqubits(const int L, const float p_loss){
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::bernoulli_distribution dist(p_loss);
    std::vector<int> lossqubits;
    const int N_E = L*L+(L-1)*(L-1); 
    lossqubits.reserve(N_E);

    for (size_t i = 0; i < N_E; ++i) {
        if (dist(engine)){
            lossqubits.push_back(i);
        }
    }
    return lossqubits;
}


std::vector<double> spacelike_weight_generator(double p_comp, const int L, const float p_loss){
    //グラフの宣言
    Matching_Graph g_matching;
    Matching_Graph g_erasure;
    std::vector<int> lossqubits = rand_lossqubits(L, p_loss);

    g_matching = planer_code_matching_graph(L);
    g_erasure = Erasured_matching_graph_creater(&g_matching, L, lossqubits);

    std::map<int, int>
    degraded_edges = get_number_of_edge_degenarate(g_erasure);  //((node_id, node_id), number_of_degenarate)のpair

    const int N_E = L*L+(L-1)*(L-1);
    std::vector<double> spacelike_weights(N_E);

    for (size_t i = 0; i < N_E; ++i) {
        auto itr = degraded_edges.find(i);
        if (itr != degraded_edges.end()){
            float n = float(degraded_edges[i]);
            float val = std::pow(1-2*p_comp, n);
            float pn = (1-val)*0.5;
            spacelike_weights.at(i) = log((1-pn)/pn);
        }
        else{
            double pn = 1.0E-308;
            spacelike_weights.at(i) = log((1-pn)/pn);
        }
    }
    // //出力検査用のprint
    // std::cout << "edge_num \t"<< "weight"<<std::endl;
    // for (int i=0; i<N_E ; ++i){
    //     std::cout << i << "\t\t" << space_like_weights.at(i) <<std::endl;
    // } 

return spacelike_weights;
}

BOOST_PYTHON_MODULE(weight_culcurator)
{
    using namespace boost::python;
    def("spacelike_weight_generator", spacelike_weight_generator);
}
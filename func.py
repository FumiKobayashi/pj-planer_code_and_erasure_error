import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag
import random
from pymatching import Matching
import networkx as nx
import copy

# 高速化まわり
import joblib



def planer_code_x_stabilisers(L):                               # planer codeのstabilizerを表す行列を生成
    row_ind = []
    col_ind = []

    for i in range(L**2-L):
        if i%L == 0:
            row_ind.append(i)
            col_ind.append(i)
            row_ind.append(i)
            col_ind.append(L+i)
            row_ind.append(i)
            col_ind.append(L**2+(L-1)*i/L)
        elif i%L == L-1:
            row_ind.append(i)
            col_ind.append(i)
            row_ind.append(i)
            col_ind.append(L+i)
            row_ind.append(i)
            col_ind.append(L**2+(L-1)*(i+1)/L-1)
        else:
            row_ind.append(i)
            col_ind.append(i)
            row_ind.append(i)
            col_ind.append(L+i)
            row_ind.append(i)
            col_ind.append(L**2+(L-1)*(i//L) + i%L-1)
            row_ind.append(i)
            col_ind.append(L**2+(L-1)*(i//L) + i%L)
        
    row_ind = [int(i) for i in row_ind]
    col_ind = [int(i) for i in col_ind]
    row_ind = tuple(row_ind)
    col_ind = tuple(col_ind)
    
    data = np.ones((L-1)*(4*L-2), dtype=np.uint8)                                         # 長さ 2nの全ての要素が1であるベクトル
    H = csr_matrix((data, (row_ind, col_ind)))  

    H.data = H.data % 2                                 # csr_matrix().dataで行列が持つ非零要素の列を取得。ここではmod 2をしている
    H.eliminate_zeros()                                 # 行列のzero要素を取り除く
    return csr_matrix(H)                                # CSRtypeの疎行列に変換


def planer_code_x_logicals(L):
    """
    Sparse binary matrix with each row corresponding to an X logical operator
    of a toric code with lattice size L. Constructed from the
    homology groups of the repetition codes using the Kunneth
    theorem.
    """
    x_logicals = [1 for i in range(L)]
    x_logicals = x_logicals + [0 for i in range(L**2+(L-1)**2 -L)]
    return csr_matrix(x_logicals)


def num_decoding_failures(H,logicals, L, p_comp, num_trials, ploss=0.0, multi_thread=False):
    node_num = np.array([i for i in range(L**2+(L-1)**2)])           # lossしたqubitを見つけるための数列
    num_errors = 0                                              # エラーカウンターの設定
    if multi_thread:
        print('p=',p_comp) 
        result_list = joblib.Parallel(n_jobs=-1, verbose=3)([joblib.delayed(onestep)(H,logicals, L, p_comp, ploss, node_num) for i in range(num_trials)])
        num_errors = sum(result_list)                                          
    else:
        for i in range(num_trials):
            flag_loss = np.random.binomial(1, ploss, L**2+(L-1)**2)
            lossqubits = node_num  * flag_loss
            lossqubits = lossqubits[lossqubits.nonzero()]
            g_matching = multi_graph_generator(H, L)

            # weightの定義
            spacelike_weights = spacelike_weight_generator(g_matching, p_comp, L, lossqubits)
            # spacelike_weights=[np.log((1-1e-300)/1e-300) if i in lossqubits else np.log((1-p_comp)/p_comp) for i in range(L**2+(L-1)**2)]

            matching = Matching(H, spacelike_weights=spacelike_weights) 
            noise = np.random.binomial(1, p_comp, H.shape[1])            # nioseの生成
            syndrome = H@noise % 2                                  # syndrome値の検出
            correction = matching.decode(syndrome, num_neighbours=16)                  # デコード結果の出力
            error = (noise + correction) % 2                        # 訂正後の残留エラーの計算
            if np.any(error@logicals.T % 2):                        # 訂正の可否を確認
                num_errors += 1
    return num_errors

# 上の関数の並列処理するパートの関数

def onestep(H,logicals, L, p_comp, ploss, node_num):
    num_errors = 0 
    flag_loss = np.random.binomial(1, ploss, L**2+(L-1)**2)
    lossqubits = node_num  * flag_loss
    lossqubits = lossqubits[lossqubits.nonzero()]
    g_matching = multi_graph_generator(H, L)

    # weightの定義
    spacelike_weights = spacelike_weight_generator(g_matching, p_comp, L, lossqubits)
    # spacelike_weights=[np.log((1-1e-300)/1e-300) if i in lossqubits else np.log((1-p_comp)/p_comp) for i in range(L**2+(L-1)**2)]

    matching = Matching(H, spacelike_weights=spacelike_weights) 
    noise = np.random.binomial(1, p_comp, H.shape[1])            # nioseの生成
    syndrome = H@noise % 2                                  # syndrome値の検出
    correction = matching.decode(syndrome)                  # デコード結果の出力
    error = (noise + correction) % 2                        # 訂正後の残留エラーの計算
    if np.any(error@logicals.T % 2):                        # 訂正の可否を確認
        num_errors += 1
    return num_errors



def spacelike_weight_generator(graph_matching, p, L, lossqubits=None):
    dict_loss_degraded_edges = get_edge_numbers(Erasured_matching_graph_creater(graph_matching, L,lossqubits))
    graph_matching = add_degrade_attribute(graph_matching, dict_loss_degraded_edges, lossqubits)
    spacelike_weights_dict = {}
    for e in graph_matching.edges(keys = True):
        if graph_matching.edges[e]['n'] == 0:
            spacelike_weights_dict[graph_matching.edges[e]['fault_id']] = np.log((1-1e-300)/1e-300)
        else: 
            n  = graph_matching.edges[e]['n']
            pn = (1-(1-2*p)**n)*0.5
            spacelike_weights_dict[graph_matching.edges[e]['fault_id']]= np.log((1-pn)/pn)
    spacelike_weights_sort = sorted(spacelike_weights_dict.items(), key=lambda i: i[0])
    spacelike_weights = [item[1] for item in spacelike_weights_sort]
    return spacelike_weights



def multi_graph_generator(H, L):
    matching = Matching(H)
    g = matching.to_networkx()
    y=0
    i=-1
    for n in g.nodes():
        if not g.nodes[n]['is_boundary']:
            i +=1
            if i !=0:
                if i%L == 0:
                    y+=1
                else:
                    pass
            g.nodes[n]['pos'] = (i%L, y)
        else:
            g.nodes[n]['pos'] = ((L-1)/2, -2)
    for e in g.edges():
        g.edges[e]['fault_id'] = list(g.edges[e]['fault_ids'])[0]
    
    g_multi = nx.MultiGraph(g)
    return g_multi
    
def draw_matching_graph(self, edge_labels=True, ax = None):
    plt.rcParams['figure.figsize'] = (10.0, 10.0)
    # グラフ描画 ============================
    ngrid = self.ngrid
    g = self.matching_graph()
    node_pos={n: g.nodes[n]['pos'] for n in g.nodes()}
    # node_size=[]
    # for n in g.nodes():
    #     if  g.nodes[n]['is_boundary']:
    #         node_size.append(100/ngrid)
    #     else:
    #         node_size.append(0)
    nx.draw(g,
            pos=node_pos, 
            node_color='white', 
            # edge_color='black',
            node_size=700/ngrid, 
            font_size=40/ngrid,
            ax=ax
            )
    if edge_labels:
        # node_posから座標の指定
        nx.draw_networkx_edge_labels(g, node_pos, 
                {e: g.edges[e]['fault_ids'] for e in g.edges()}, 
                font_size=8, 
                rotate=False, 
                clip_on=True,
                ax=ax
                )
    plt.axis('equal')
    #=======================================

#属性を検索する関数
def find_specific_attribute_edge(G, attr, value):
    d = nx.get_edge_attributes(G, attr)
    return  [key for key, v in d.items() if v==value]


" matchingグラフのedgeに対する\"ダブり\"の数カウントを返す関数"

def add_degrade_attribute(g_matching, dict_loss_degraded_edges, lossqubits):

    for e in g_matching.edges(keys=True):
        if g_matching.edges[e]['fault_id'] in lossqubits:
            g_matching.edges[e]['n'] = 0
        elif g_matching.edges[e]['fault_id'] in dict_loss_degraded_edges.keys():
            g_matching.edges[e]['n'] = dict_loss_degraded_edges[g_matching.edges[e]['fault_id']]
        else:
            g_matching.edges[e]['n'] = 1

    return g_matching


def Erasured_matching_graph_creater(g, L, erasure_errors): # g:multi graphのmatching_graph
    g_matching = copy.deepcopy(g)

    # edgeにエラー属性を付加
    for e in g_matching.edges(keys=True):
        g_matching.edges[e]['erasure'] = False
    for n in g_matching.nodes():
        g_matching.nodes[n]['node_ids'] = [n]
  
    for id in erasure_errors:
            l = find_specific_attribute_edge(g_matching, 'fault_id', id)
            if len(l)>0:
                e = l[0]
                # 消失したedgeに関するnodeを統合 
                node_id0= e[0]
                node_id1= e[1]
                nodeset = {node_id0,node_id1}
                if L**2-L in nodeset:
                    node_id0 = max(nodeset)
                    node_id1 = min(nodeset)
                new_pos = ((g_matching.nodes[node_id0]['pos'][0]+g_matching.nodes[node_id1]['pos'][0])/2, (g_matching.nodes[node_id0]['pos'][1]+g_matching.nodes[node_id1]['pos'][1])/2)
                g_matching.nodes[node_id0]['pos'] = new_pos
                
                for e1 in g_matching.edges(node_id1, keys=True):
                    if node_id0!=e1[1]:
                        g_matching.add_edge(node_id0, e1[1], **g_matching.edges[e1])
                g_matching.remove_node(node_id1)
            else:
                pass

    return  g_matching

def get_edge_numbers(g_matching):
    return {g_matching.edges[e]['fault_id']: g_matching.number_of_edges(e[0], e[1]) for e in g_matching.edges(keys = True)}
#########################################

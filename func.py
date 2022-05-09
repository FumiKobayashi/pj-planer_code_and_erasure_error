import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag
import random
from pymatching import Matching
import networkx as nx
import copy

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


def num_decoding_failures(H,logicals, L, p, num_trials, ploss=0.0):
    num = np.array([i for i in range(L**2+(L-1)**2)])           # lossしたqubitを見つけるための数列
    num_errors = 0                                              # エラーカウンターの設定
    for i in range(num_trials):
        flag_loss = np.random.binomial(2, ploss, L**2+(L-1)**2)
        lossqubits = num  * flag_loss
        lossqubits = lossqubits[lossqubits.nonzero()]
        g_matching = multi_graph_generator(H, L)


        spacelike_weights=[np.log((1-1e-300)/1e-300) if i in lossqubits else np.log((1-p)/p) for i in range(L**2+(L-1)**2)]

        matching = Matching(H, spacelike_weights=spacelike_weights) 
        noise = np.random.binomial(1, p, H.shape[1])            # nioseの生成
        syndrome = H@noise % 2                                  # syndrome値の検出
        correction = matching.decode(syndrome)                  # デコード結果の出力
        error = (noise + correction) % 2                        # 訂正後の残留エラーの計算
        if np.any(error@logicals.T % 2):                        # 訂正の可否を確認
            num_errors += 1
    return num_errors


# def weight_generator(p, L, lossqubits=None):
#     spacelike_weights = []
#     for i in range(L**2+(L-1)**2):
#         if i in lossqubits:
#             spacelike_weights.append(np.log((1-1e-300)/1e-300))
#         elif i in *** :
#             n  = ***
#             pn = (1-(1-2*p)**n)*0.5
#             spacelike_weights.append(np.log((1-pn)/pn))
#         else:
#             spacelike_weights.append(np.log((1-p)/p))
#     return spacelike_weights

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


def contenious_loss_detecter(lossqubits, L, g_matching):
    superweight_qubits = []
    already = []
    n_dict = {}
    for lossqubit in lossqubits:
        if lossqubit < L**2:
            if (lossqubit+1 in lossqubits) and (not (lossqubit+1 in already)):
                already.append(lossqubit+1)
                
            else:
                pass
        elif lossqubit > (L**2-1):
            """ ここ作る """
            pass
        else:
            print("There is some error in the contenious_loss_detecter.")


#属性を検索する関数
def find_specific_attribute_edge(G, attr, value):
    result = []
    d = nx.get_edge_attributes(G, attr)
    for key, v in d.items():
        if(v == value):
            result.append(key)
    return result


def find_specific_node_id(G, attr, id):
    result=None
    d = nx.get_node_attributes(G, attr)
    for key, v in  d.items():
        if(id in v):
            result=key
            break
        else:
            pass
    if result:
        print("find_specific_node_id could not find {0} on attribute {1}".format(id, attr))
    return result

def erasure_error(g_matching, error_prob={"erasure":0.0}):
    # エラーが起きるedgeの乱択
    # 二項分布からエラーが起こるedgeの個数を選択
    num_eerrer = np.random.binomial(g_matching.number_of_edges(), error_prob["erasure"], 1)

    # エラーの起きるedgeの選択
    erasure_errors = random.sample([i for i in range(g_matching.number_of_edges())], k=num_eerrer[0])

    # edegeにエラー属性を付加
    for e in g_matching.edges():
        g_matching.edges[e]['erasure'] = False
    erasure_edges = []
    for id in erasure_errors:
        temp=find_specific_attribute_edge(g_matching, "fault_ids", {id})
        erasure_edges.extend(temp)
    for e in erasure_edges:
        pass     
        
    return g_matching

def Multigraph_creater(g, error_prob={"erasure":0.5}): # g:multi graphのmatching_graph
    g_matching = copy.deepcopy(g)
    num_eerror = np.random.binomial(g_matching.number_of_edges(), error_prob["erasure"], 1)

    # エラーの起きるedgeの選択
    erasure_errors = random.sample([i for i in range(g_matching.number_of_edges())], k=num_eerror[0])
    print("erasure errors:",erasure_errors)

    # edegeにエラー属性を付加
    for e in g_matching.edges(keys=True):
        g_matching.edges[e]['erasure'] = False
    for n in g_matching.nodes():
        g_matching.nodes[n]['node_ids'] = [n]

    erasure_edges = []
    for id in erasure_errors:
        temp=find_specific_attribute_edge(g_matching, "fault_ids", {id})
        erasure_edges.extend(temp)

    print("erasure_edges",erasure_edges)

    # 消失したedgeに関するnodeを統合
    for e0 in erasure_edges:
        node_id0=find_specific_node_id(g_matching, 'node_ids', e0[0])
        node_id1=find_specific_node_id(g_matching, 'node_ids', e0[1])
        print("node_id:", node_id0, node_id1)
        if node_id0 !=node_id1:
            new_node_ids = g_matching.nodes[node_id0]['node_ids']
            new_node_ids.extend(g_matching.nodes[node_id1]['node_ids'])
            g_matching.nodes[node_id0]['node_ids'] = new_node_ids

            new_pos = ((g_matching.nodes[node_id0]['pos'][0]+g_matching.nodes[node_id1]['pos'][0])/2, (g_matching.nodes[node_id0]['pos'][1]+g_matching.nodes[node_id1]['pos'][1])/2)
            g_matching.nodes[node_id0]['pos'] = new_pos

            for e1 in g_matching.edges(node_id1, keys=True):
                print("e1:",e1)
                if node_id0!=e1[1]:
                    g_matching.add_edge(node_id0, e1[1], **g_matching.edges[e1])
            g_matching.remove_node(node_id1)
    
    return  g_matching


def get_edge_number(g_matching, L):
    multi_loss_edges = {}
    for n in g_matching.nodes():
        for m in range(n+1, L*(L-1)):
            l = g_matching.number_of_edges(n, m)
            if l >1:
                for node0 in g_matching.nodes[n]['node_ids']:
                    for node1 in g_matching.nodes[m]['node_ids']:
                        if node0<node1:
                            multi_loss_edges[(node0, node1)]=l
                        else:
                            pass
    return multi_loss_edges
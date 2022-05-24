import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag
import pickle
from time import time
import os

from func import planer_code_x_stabilisers, planer_code_x_logicals, num_decoding_failures

# 実行開始時間の入力
import datetime
t_delta = datetime.timedelta(hours=9)
JST = datetime.timezone(t_delta, 'JST')
now = datetime.datetime.now(JST)
day = now.strftime('%Y-%m-%d-%H%M')

# configの設定
plosses = np.linspace(0.0, 0.3, 5)                           # 線形で値を決めるとき
plosses = np.hstack((plosses, np.linspace(0.35, 0.5, 7)))      # 0.5付近だけ細かく調べる
num_trials = 5000
Ls = range(24,37,4)                                                 # code distanceを指定するiteraterの作成
# Ls = range(8,21,4)
 
ploss_th = 0.35 # 計算エラーの刻みを変える閾値の指定   
ps1 = np.linspace(0.001, 0.2, 10)                                      # 計算に用いる確率の列1
ps2 = np.logspace(-3, -0.3010299956639812, 10, base=10)                                     # 計算に用いる確率の列2


np.random.seed(2)

# dict_config = {'plosses':plosses, 'Ls':Ls, 'ps':ps, "num_trials":num_trials}
dict_config = {'plosses':plosses, 'Ls':Ls, 'ps':list(ps1).extend(list(ps2)), "num_trials":num_trials}
output_dict_name = './{}/'.format(day)
os.mkdir("./"+output_dict_name)
f = open(output_dict_name+'config.pickle', 'wb')
pickle.dump(dict_config, f)
f.close()

start = time()
dict_log_errors_all_L = {}
for ploss in plosses:
    if ploss >ploss_th:
        ps = ps1
    else:
        ps = ps2
    print("Simulating p_loss={}...".format(ploss))
    log_errors_all_L = []                                               # code distanceごとのエラー確率を格納するlistを格納するlist
    for L in Ls:
        print("Simulating L={}...".format(L))
        Hx = planer_code_x_stabilisers(L)                                # toric codeのstabilizer生成
        logX = planer_code_x_logicals(L)                                 # toric codeの論理qubit空間のstabilzier
        log_errors = [num_decoding_failures(Hx, logX, L, p, num_trials, ploss=ploss, multi_thread=True)/num_trials for p in ps]                                                 # code distanceごとのエラー確率を格納するlist
        # for p in ps:
        #     num_errors = num_decoding_failures(Hx, logX, L, p, num_trials, ploss=ploss, multi_thread=False)
        #     log_errors.append(num_errors/num_trials)
        log_errors_all_L.append(np.array(log_errors))

        dict_log_errors_all_L[ploss] = log_errors_all_L
        output_file_name = output_dict_name+'dict_log_errors_all_L.pickle'
        f = open(output_file_name, 'wb')
        pickle.dump(dict_log_errors_all_L, f)
        f.close()

print('Calculation time: {}(sec)'.format(time() - start))
    
    
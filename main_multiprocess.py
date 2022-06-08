import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag
import pickle
from time import time
import joblib
import os

from func import planer_code_x_stabilisers, planer_code_x_logicals, num_decoding_failures

# 実行開始時間の指定
import datetime
t_delta = datetime.timedelta(hours=9)
JST = datetime.timezone(t_delta, 'JST')
now = datetime.datetime.now(JST)
day = now.strftime('%Y-%m-%d-%H%M')


plosses = np.linspace(0.001, 0.5, 3)
num_trials = 500    
Ls = range(8,21,4)                                                  # code distanceを指定するiteraterの作成
ps = np.linspace(0.001, 0.2, 3)                                      # 計算に用いる確率の列
np.random.seed(2)

dict_config = {'plosses':plosses, 'Ls':Ls, 'ps':ps, "num_trials":num_trials}
output_dict_name = './{}/'.format(day)
os.mkdir("./"+output_dict_name)
f = open(output_dict_name+'config.pickle', 'wb')
pickle.dump(dict_config, f)
f.close()

'''
並列化する部分
'''
def onestep_culc(ploss, Ls, ps):
    print("Simulating p_loss={}...".format(ploss))
    log_errors_all_L = []                                               # code distanceごとのエラー確率を格納するlistを格納するlist
    for L in Ls:
        print("Simulating L={}...".format(L))
        Hx = planer_code_x_stabilisers(L)                                # toric codeのstabilizer生成
        logX = planer_code_x_logicals(L)                                 # toric codeの論理qubit空間のstabilzier
        log_errors = []                                                 # code distanceごとのエラー確率を格納するlist
        for p in ps:
            num_errors = num_decoding_failures(Hx, logX, L, p, num_trials, ploss=ploss)
            log_errors.append(num_errors/num_trials)
        # output_file_name = output_dict_name+'./results/dict_log_errors_L={}.pickle'.format(L)
        # f = open(output_file_name, 'wb')
        # pickle.dump(log_errors, f)
        # f.close()
        log_errors_all_L.append(np.array(log_errors))
    return (ploss, log_errors_all_L)


start = time()
result_list = joblib.Parallel(n_jobs=-2, verbose=2)( [joblib.delayed(onestep_culc)(ploss, Ls, ps) for ploss in plosses] )

dict_log_errors_all_L = {}
for l in result_list:
    dict_log_errors_all_L[l[0]] = l[1]
output_file_name = output_dict_name+'dict_log_errors_all_L.pickle'
f = open(output_file_name, 'wb')
pickle.dump(dict_log_errors_all_L, f)
f.close()

print('Calculation time: {}(sec)'.format(time() - start))
    
    
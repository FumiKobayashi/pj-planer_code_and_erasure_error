import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import hstack, kron, eye, csr_matrix, block_diag
import pickle
from time import time

from func import planer_code_x_stabilisers, planer_code_x_logicals, num_decoding_failures


plosses = np.linspace(0.001, 0.5, 5)
num_trials = 2000
Ls = range(12,24,4)                                                  # code distanceを指定するiteraterの作成
ps = np.linspace(0.001, 0.2, 10)                                      # 計算に用いる確率の列
np.random.seed(2)

dict_config = {'plosses':plosses, 'Ls':Ls, 'ps':ps, "num_trials":num_trials}
output_file_name = 'config.pickle'
f = open(output_file_name, 'wb')
pickle.dump(dict_config, f)
f.close()

start = time()
dict_log_errors_all_L = {}
for ploss in plosses:
    log_errors_all_L = []                                               # code distanceごとのエラー確率を格納するlistを格納するlist
    for L in Ls:
        print("Simulating L={}...".format(L))
        Hx = planer_code_x_stabilisers(L)                                # toric codeのstabilizer生成
        logX = planer_code_x_logicals(L)                                 # toric codeの論理qubit空間のstabilzier
        log_errors = []                                                 # code distanceごとのエラー確率を格納するlist
        for p in ps:
            num_errors = num_decoding_failures(Hx, logX, L, p, num_trials, ploss=ploss)
            log_errors.append(num_errors/num_trials)
        log_errors_all_L.append(np.array(log_errors))

    dict_log_errors_all_L[ploss] = log_errors_all_L
    output_file_name = 'dict_log_errors_all_L.pickle'
    f = open(output_file_name, 'wb')
    pickle.dump(dict_log_errors_all_L, f)
    f.close()

print('Calculation time: {}(sec)'.format(time() - start))
    
    
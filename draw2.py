import networkx as nx
import matplotlib.pyplot as plt
import pickle
import datetime

import os


# GUIでファイルを選択
# import tkinter, tkinter.filedialog, tkinter.messagebox
# root = tkinter.Tk()
# root.withdraw()
# fTyp = [("","*")]
# iDir = os.path.abspath(os.path.dirname(__file__))
# tkinter.messagebox.showinfo('draw2.py','処理ディレクトリを選択してください')
# input_dir_name = tkinter.filedialog.askdirectory(initialdir = iDir)
# # 処理ファイル名の出力
# tkinter.messagebox.showinfo('draw2.py',input_dir_name)

# 計算結果のdirの取得
path = './'
files = os.listdir(path)
files_dir = [f for f in files if os.path.isdir(os.path.join(path, f))]

# 描画するdirの指定
for i, d in enumerate(files_dir):
    print('{0}:{1}'.format(i,d))
input_dir_number = int(input('描画するdirectryの番号を選択してください:'))
input_dir_name = files_dir[input_dir_number]
input_file_name = 'dict_log_errors_all_L.pickle'
f = open('{}/{}'.format(input_dir_name,input_file_name), 'rb')
dict_log_errors_all_L = pickle.load(f)
f.close()
# print(type(dict_log_errors_all_L))
input_file_name = 'config.pickle'
f = open('{}/{}'.format(input_dir_name,input_file_name), 'rb')
dict_config = pickle.load(f)
f.close()

t_delta = datetime.timedelta(hours=9)
JST = datetime.timezone(t_delta, 'JST')
now = datetime.datetime.now(JST)
day = now.strftime('%Y-%m-%d-%H%M')
os.mkdir(input_dir_name+"/fig_{0}".format(day))

num_trials = dict_config['num_trials']
Ls = dict_config['Ls']            # code distanceを指定するiteraterの作成
ps = dict_config['ps'] 

for key in dict_log_errors_all_L.keys():
    log_errors_all_L = dict_log_errors_all_L[key]
    plt.figure()
    for L, logical_errors in zip(Ls, log_errors_all_L):
        std_err = (logical_errors*(1-logical_errors)/num_trials)**0.5
        plt.errorbar(ps, logical_errors, yerr=std_err, label="L={}".format(L))
    plt.xlabel("Physical error rate")
    plt.ylabel("Logical error rate")
    plt.legend(loc=0);
    plt.savefig(input_dir_name+"/fig_{0}/{1}.png".format(day,key))

    # plt.figure()
    # for L, logical_errors in zip(Ls, log_errors_all_L):
    #     std_err = (logical_errors*(1-logical_errors)/num_trials)**0.5
    #     plt.errorbar(ps, logical_errors, yerr=std_err, label="L={}".format(L))
    # plt.xlabel("Physical error rate")
    # plt.ylabel("Logical error rate")
    # plt.xlim(0.075, 0.125)
    # plt.ylim(0.1, 0.3)
    # plt.legend(loc=0);
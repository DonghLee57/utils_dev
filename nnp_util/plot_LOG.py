# LOG from SIMPLE_NN v2.0.0
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    info = readLOG(sys.argv[1])
    fig, ax = plt.subplots(1,3) 
    ax[0].semilogy(info['epoch'],info['e_t'],'r')
    ax[0].semilogy(info['epoch'],info['e_v'],'r:')
    ax[0].set_ylabel('E RMSE',fontsize=15)
    ax[0].set_ylim([1E-3,1E+0])
    ax[1].semilogy(info['epoch'],info['f_t'],'b')
    ax[1].semilogy(info['epoch'],info['f_v'],'b:')
    ax[1].set_ylabel('F RMSE',fontsize=15)
    ax[1].set_ylim([1E-1,1E+1])
    ax[2].semilogy(info['epoch'],info['s_t'],'g')
    ax[2].semilogy(info['epoch'],info['s_v'],'g:')
    ax[2].set_ylabel('S RMSE',fontsize=15)
    ax[2].set_ylim([1E+0,1E+2])
    for i in range(3):
        ax[i].tick_params(axis='both', which='major', labelsize=15)
        ax[i].set_xlabel('Epoch',fontsize=15)

    plt.tight_layout()
    plt.show()
    return 1

#--------------------------------------------------------------------------------------
def readLOG(filename):
    info={'epoch':[],\
          'e_t':[], 'e_v':[],\
          'f_t':[], 'f_v':[],\
          's_t':[], 's_v':[]}
    o = enumerate(open(filename,'r'))
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            if 'SIMPLE_NN' in tmp:
                if len(info['epoch'])==0: last_epo = 0
                else: last_epo = info['epoch'][-1]
           
            if 'Epoch' in tmp:
                info['epoch'].append(int(tmp[1])+last_epo)
                info['e_t'].append(float(tmp[5]))
                info['e_v'].append(float(tmp[6]))
                if len(tmp) > 7:
                    info['f_t'].append(float(tmp[10]))
                    info['f_v'].append(float(tmp[11]))
                if len(tmp) > 12:
                    info['s_t'].append(float(tmp[15]))
                    info['s_v'].append(float(tmp[16]))
    return info

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()

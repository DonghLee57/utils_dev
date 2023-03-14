# python 3.X
import sys
import numpy as np
import matplotlib.pyplot as plt
import ase
import ase.io
from ase import Atoms

tu = 1 # fs


#--------------------------------------------------------------------------------------
def main(time_unit):
    info = readLOG('vacf.dat')
    vacf = np.array(info['vacf'])/info['vacf'][0]

    # Original unit: fs = 10**(-15) sec
    # Unit: ps = 10**(-12) sec
    time = np.array(info['step'])*time_unit/1000
    dt = time[1]-time[0]
    
    fft_vacf = np.fft.fft(vacf)
    # Unit: THz = 10**(12) Hz
    full_axis = np.fft.fftfreq(len(fft_vacf), dt)
    full_pdos = np.absolute(fft_vacf)

    fig, ax = plt.subplots(2,1)
    ax[0].plot(time,vacf,'ro')
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlabel('Time (ps)',fontsize=15)
    ax[0].set_ylabel('VACF',fontsize=15)
    plt.tight_layout()

    ax[1].plot(full_axis,full_pdos,'b-')
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_xlim([0,20])
    ax[1].set_ylim(bottom=0)
    # 1 THz = 4.136 meV = 33.356 cm-1
    ax[1].set_xlabel('frequency (THz)',fontsize=15)
    ax[1].set_ylabel('PDOS',fontsize=15)

    plt.tight_layout()
    plt.show()
    return 1
    
def test():
    # sin function
    Fs = 1000
    T = 1/Fs
    time = np.linspace(0,1,Fs)
    s1 = 2*np.sin(10*2*np.pi*time)
    s2 = 1*np.sin(20*2*np.pi*time)
    s3 = 0.5*np.sin(30*2*np.pi*time)
    s4 = 0.2*np.sin(40*2*np.pi*time)
    signal = s1+s2+s3+s4
    sfft=np.fft.fft(signal)
    amp = np.absolute(sfft)*(2/len(sfft))
    freq = np.fft.fftfreq(len(sfft),T)

    fig, ax = plt.subplots(2,1)
    ax[0].plot(time,signal,'ro')
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlabel('Time',fontsize=15)
    ax[0].set_ylabel('Signal',fontsize=15)
    plt.tight_layout()

    ax[1].stem(freq,amp,'bo')
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_xlim([0,50])
    ax[1].set_xlabel('frequency',fontsize=15)
    ax[1].set_ylabel('amplitude',fontsize=15)

    plt.tight_layout()
    plt.show()
    return 1

def check_vacf():
    from ase.units import Ang,fs,eV

    info = readLOG('vacf.dat')
    trj = ase.io.read('dump.lammps',index=':',format='lammps-dump-text')
    tot_frame = len(trj)
    NIONS = len(trj[0].positions)
    # velocity unit: Ang/fs
    v_0 = trj[0].get_velocities()*fs
    v = 0
    for n in range(NIONS):
        v += np.dot(v_0[n],v_0[n])
    my_v = [v]
    for f in range(1,tot_frame):
        v_t = trj[f].get_velocities()*fs
        v = 0 
        for n in range(NIONS):
            v += np.dot(v_0[n],v_t[n])
        my_v.append(v)

    # Normalization
    log_vacf = np.array(info['vacf'])/info['vacf'][0]
    cal_vacf = np.array(my_v)/my_v[0]

    # Plot
    fig, ax = plt.subplots() 
    ax.plot(log_vacf[:100],'ro',label='Calc. in LAMMPS')
    ax.plot(cal_vacf[:100],'b*',label='Direct calc.')
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel('Time (fs)',fontsize=15)
    ax.set_ylabel('VACF',fontsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()
    return 1

def readLOG(filename):
    info={'step':[],\
          'temperature':[],\
          'vacf':[]}
    o = enumerate(open(filename,'r'))
    next(o)
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            info['step'].append(int(tmp[0]))
            info['temperature'].append(float(tmp[1]))
            info['vacf'].append(float(tmp[2]))
    return info

#--------------------------------------------------------------------------------------
if __name__ == "__main__":
    #check_vacf()
    #test()
    main(tu)

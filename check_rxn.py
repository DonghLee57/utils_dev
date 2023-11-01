import numpy as np
from scipy import constants
q = constants.e
k_eV = constants.Boltzmann/q
import matplotlib.pyplot as plt

def calc(MAT, T_K):
    dH = []
    S = []
    if MAT == 'H2':
        for i in range(len(T_K)):
            t = T_K[i]/1000
            if T_K[i] <1000:
                A = 33.066178
                B = -11.363417
                C = 11.432816
                D = -2.772874
                E = -0.158558
                F = -9.980797
                G = 172.707974
                H = 0.0
            elif T_K[i] <2500:
                A = 18.563083
                B = 12.257357
                C = -2.859786
                D = 0.268238
                E = 1.977990
                F = -1.147438
                G = 156.288133
                H = 0.0
            dH.append(A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 -E/t + F - H)
            S.append(A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 -E/(2*t**2) + G)
    elif MAT == 'N2':
        for i in range(len(T_K)):
            t = T_K[i]/1000
            if T_K[i] <500:
                A = 28.98641
                B = 1.853978
                C = -9.647459
                D = 16.63537
                E = 0.000117
                F = -8.671914
                G = 226.4168
                H = 0.0
            elif T_K[i] < 2100:
                A = 19.50583
                B = 19.88705
                C = -8.598535
                D = 1.369784
                E = 0.527601
                F = -4.935202
                G = 212.3900
                H = 0.0
            dH.append(A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 -E/t + F - H)
            S.append(A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 -E/(2*t**2) + G)
    elif MAT == 'NH3':
        for i in range(len(T_K)):
            t = T_K[i]/1000
            if T_K[i] <1400:
                A = 19.99563
                B = 49.77119
                C = -15.37599
                D = 1.921168
                E = 0.189174
                F = -53.30667
                G = 203.8591
                H = -45.89806
            elif T_K[i] <6000:
                A = 52.02427
                B = 18.48801
                C = -3.765128
                D = 0.248541
                E = -12.45799
                F = -85.53895
                G = 223.8022
                H = -45.89806
            dH.append(A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 -E/t + F - H -45.90)
            S.append(A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 -E/(2*t**2) + G)
    return [np.array(dH),np.array(S)/1000]

T = np.arange(300,2001,100)

fig,ax = plt.subplots(1,1)
H2 = calc('H2',T)
N2 = calc('N2',T)
NH3 = calc('NH3',T)
H = (2*(NH3[0]) - ((N2[0]) + 3*(H2[0])))
S = (2*(NH3[1]) - ((N2[1]) + 3*(H2[1])))
G = (2*(NH3[0]-T*NH3[1]) - ((N2[0]-T*N2[1]) + 3*(H2[0]-T*H2[1])))
print(np.mean(H/S))
ax.plot(T-298, H,'ro-')
ax.plot(T-298, G,'ro-')
plt.show()

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:51:12 2022

@author: pvrma
"""

from TDMA_main import LBL_TDMA_dif as LBL_TDMA_diff
from TDMA_main import LBL_TDMA as LBL_TDMA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
matplotlib.use('agg')

domain_dimension = 1
grid_resolution = 31
def u(t):
    return(5)
def v(t):
    return(5)
T = 1
rho = 6

# def u(t):
#     return(0)
# def v(t):
#     return(0)
# T = 0.4
# rho = 6
gr = grid_resolution
CM = np.zeros([grid_resolution,grid_resolution,12])
del_x = domain_dimension/grid_resolution
del_y = domain_dimension/grid_resolution

# C = [0:phi1 , 1:phi0, 2:D(constant), 3:Fn, 4:Fe, 5:Fs, 6:Fw,7:Fno,8:Feo,9:Fso,10:Fwo, 11:boundary_status]
# boundary status:
    #north face = 1
    #east face = 2
    #south face = 3
    #west face = 4
    #NE = 12
    #NW = 14
    #SE = 32
    #SW = 34
# populate C initially


for i in range(gr):
    for j in range(gr):
        CM[i,j,2] = T
        CM[i,j,7] = rho*v(0)*del_x
        CM[i,j,8] = rho*u(0)*del_y
        CM[i,j,9] = rho*v(0)*del_x
        CM[i,j,10] = rho*u(0)*del_y
        if i == 0:
            CM[i,j,11] = 4
        if i == gr - 1:
            CM[i,j,11] = 2
        if j == 0:
            CM[i,j,11] = 3
        if j == gr-1:
            CM[i,j,11] = 1
        if i == 0 and j == 0:
            CM[i,j,11] = 34
        if i == 0 and j == gr - 1:
            CM[i,j,11] = 14
        if i == gr - 1 and j == 0:
            CM[i,j,11] = 32
        if i == gr - 1 and j == gr-1:
            CM[i,j,11] = 12

#Define_initial_phi:
#lower one fourt quad phi 1
# CM[0:int(gr/4),0:int(gr/4),1] = 1
CM[15,15,1] = 1


plt.rcParams["font.family"] = "Times New Roman"
text_size = 20
r = 0.8
plt.figure()
plt.imshow(np.transpose(CM[:,:,1]),cmap = 'plasma',origin = 'lower',extent = [0,domain_dimension,0,domain_dimension])
bar = plt.colorbar()
bar.ax.tick_params(labelsize=text_size)
plt.xlabel('x',fontsize = text_size)
plt.ylabel('y',fontsize = text_size)
plt.xticks(fontsize=r*text_size)
plt.yticks(fontsize=r*text_size)
plt.savefig('at_%i.png'%0,dpi = 400, pad_inches = 0, bbox_inches='tight')
np.savetxt('at_%i.csv'%0,CM[:,:,1])


phi_out = 0
f = 0.5
ite = 1
dt = 0.001
total_time = 0.1
t = 0
while t < total_time:
    print(t)
    A = np.zeros([gr**2,gr**2])
    B = np.zeros([gr**2])
    for i in range(gr):
        for j in range(gr):
            CM[i,j,2] = T
            CM[i,j,3] = rho*v(t)*del_x
            CM[i,j,4] = rho*u(t)*del_y
            CM[i,j,5] = rho*v(t)*del_x
            CM[i,j,6] = rho*u(t)*del_y
            if CM[i,j,11] == 0:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = CM[i+1,j,1]
                phi_S_o = CM[i,j-1,1]
                phi_W_o = CM[i-1,j,1]
                Dn = CM[i,j+1,2]
                De = CM[i+1,j,2]
                Ds = CM[i,j-1,2]
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + Ae + As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = -f*Aw
                A[A_i,A_i+1] = -f*Ae
                A[A_i,A_i-gr] = -f*As
                A[A_i,A_i+gr] = -f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                 Aeo*phi_E_o - Aso*phi_S_o -
                                 Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt
                # print(i,j,A_i)
                # print(B[A_i])
                # if B[A_i] != 0:
                #     print(i,j,B[A_i])
            elif CM[i,j,11] == 1:
                phi_P_o = CM[i,j,1]
                phi_N_o = phi_out
                phi_E_o = CM[i+1,j,1]
                phi_S_o = CM[i,j-1,1]
                phi_W_o = CM[i-1,j,1]
                Dn = 2*T
                De = CM[i+1,j,2]
                Ds = CM[i,j-1,2]
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = Ae + As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = -f*Aw
                A[A_i,A_i+1] = -f*Ae
                A[A_i,A_i-gr] = -f*As
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*An*phi_out
            elif CM[i,j,11] == 2:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = phi_out
                phi_S_o = CM[i,j-1,1]
                phi_W_o = CM[i-1,j,1]
                Dn = CM[i,j+1,2]
                De = 2*T
                Ds = CM[i,j-1,2]
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = -f*Aw
                A[A_i,A_i-gr] = -f*As
                A[A_i,A_i+gr] = -f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*Ae*phi_out
            elif CM[i,j,11] == 3:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = CM[i+1,j,1]
                phi_S_o = phi_out
                phi_W_o = CM[i-1,j,1]
                Dn = CM[i,j+1,2]
                De = CM[i+1,j,2]
                Ds = 2*T
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + Ae + As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = f*Aw
                A[A_i,A_i+1] = f*Ae
                A[A_i,A_i+gr] = f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o + Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*As*phi_out
            elif CM[i,j,11] == 4:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = CM[i+1,j,1]
                phi_S_o = CM[i,j-1,1]
                phi_W_o = phi_out
                Dn = CM[i,j+1,2]
                De = CM[i+1,j,2]
                Ds = CM[i,j-1,2]
                Dw = 2*T
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + Ae + As + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i+1] = -f*Ae
                A[A_i,A_i-gr] = -f*As
                A[A_i,A_i+gr] = -f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*Aw*phi_out
            elif CM[i,j,11] == 12:
                phi_P_o = CM[i,j,1]
                phi_N_o = phi_out
                phi_E_o = phi_out
                phi_S_o = CM[i,j-1,1]
                phi_W_o = CM[i-1,j,1]
                Dn = 2*T
                De = 2*T
                Ds = CM[i,j-1,2]
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = -f*Aw
                A[A_i,A_i-gr] = -f*As
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*An*phi_out + f*Ae*phi_out
            elif CM[i,j,11] == 32:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = phi_out
                phi_S_o = phi_out
                phi_W_o = CM[i-1,j,1]
                Dn = CM[i,j+1,2]
                De = 2*T
                Ds = 2*T
                Dw = CM[i-1,j,2]
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i-1] = -f*Aw
                A[A_i,A_i+gr] = -f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*As*phi_out + f*Ae*phi_out
            elif CM[i,j,11] == 34:
                phi_P_o = CM[i,j,1]
                phi_N_o = CM[i,j+1,1]
                phi_E_o = CM[i+1,j,1]
                phi_S_o = phi_out
                phi_W_o = phi_out
                Dn = CM[i,j+1,2]
                De = CM[i+1,j,2]
                Ds = 2*T
                Dw = 2*T
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = An + Ae + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i+1] = -f*Ae
                A[A_i,A_i+gr] = -f*An
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*As*phi_out + f*Aw*phi_out
            elif CM[i,j,11] == 14:

                phi_P_o = CM[i,j,1]
                phi_N_o = phi_out
                phi_E_o = CM[i+1,j,1]
                phi_S_o = CM[i,j-1,1]
                phi_W_o = phi_out
                Dn = 2*T
                De = CM[i+1,j,2]
                Ds = CM[i,j-1,2]
                Dw = 2*T
                Fnc = CM[i,j,3]
                Fec = CM[i,j,4]
                Fsc = CM[i,j,5]
                Fwc = CM[i,j,6]
                An = (Dn + max(-Fnc,0))
                Ae = (De + max(-Fec,0))
                As = (Ds + max(Fsc,0))
                Aw = (Dw + max(Fwc,0))
                Ap = As + Aw + (Fnc + Fec - Fsc - Fwc)
                Fno = CM[i,j,7]
                Feo = CM[i,j,8]
                Fso = CM[i,j,9]
                Fwo = CM[i,j,10]
                Ano = (Dn + max(-Fno,0))
                Aeo = (De + max(-Feo,0))
                Aso = (Ds + max(Fso,0))
                Awo = (Dw + max(Fwo,0))
                Apo = Ano + Aeo + Aso + Awo + (Fno + Feo - Fso - Fwo)
                A_i = j*gr + i
                A[A_i,A_i] = rho*del_x*del_y/dt + f*Ap
                A[A_i,A_i+1] = -f*Ae
                A[A_i,A_i-gr] = -f*As
                B[A_i] = -(1-f)*(Apo*phi_P_o - Ano*phi_N_o -
                                  Aeo*phi_E_o - Aso*phi_S_o -
                                  Awo*phi_W_o)  + phi_P_o*rho*del_x*del_y/dt + f*An*phi_out + f*Aw*phi_out
                
    phi_old = CM[:,:,1]
    # CM[:,:,0] = LBL_TDMA(A, B, [gr, gr], 0.1, phi_old,1000, del_x, rho)
    CM[:,:,0] = LBL_TDMA(A,  B,  [gr, gr],  0.9, phi_old, 1000)
    CM[:,:,1] = CM[:,:,0]
    CM[:,:,7:11] = CM[:,:,3:7]
    
    if ite == int(total_time/(dt*4)) or ite == int(total_time/(dt*2)) or ite == int(3*total_time/(dt*4)) or ite == int(total_time/(dt)):
        plt.figure()
        plt.imshow(np.transpose(CM[:,:,0]),cmap = 'plasma',origin = 'lower',extent = [0,domain_dimension,0,domain_dimension])
        bar = plt.colorbar()
        bar.ax.tick_params(labelsize=text_size)
        plt.xlabel('x',fontsize = text_size)
        plt.ylabel('y',fontsize = text_size)
        plt.xticks(fontsize=r*text_size)
        plt.yticks(fontsize=r*text_size)
        plt.savefig('at_%i.png'%ite,dpi = 400, pad_inches = 0, bbox_inches='tight')
        np.savetxt('at_%i.csv'%ite,CM[:,:,0])
    # np.savetxt('at_%i.csv'%ite, CM[:,:,0])
    t += dt
    ite += 1
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:31:11 2022

@author: pvrma
"""

# -*- coding: utf-8 -*-



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

BC = np.zeros([gr,gr])
for i in range(gr):
    for j in range(gr):
        if i == 0:
            BC[i,j] = 4
        if i == gr - 1:
            BC[i,j] = 2
        if j == 0:
            BC[i,j] = 3
        if j == gr-1:
            BC[i,j] = 1
        if i == 0 and j == 0:
            BC[i,j] = 34
        if i == 0 and j == gr - 1:
            BC[i,j] = 14
        if i == gr - 1 and j == 0:
            BC[i,j] = 32
        if i == gr - 1 and j == gr-1:
            BC[i,j] = 12

#Define_initial_phi:
#lower one fourt quad phi 1
# CM[0:int(gr/4),0:int(gr/4),1] = 1


phi_out = 0
f = 0
ite = 1
dt = 0.00001
total_time = 0.1
t = 0
k1 = np.zeros([gr,gr])
k0 = np.zeros([gr,gr])
phi = np.zeros([gr,gr])
phi_1 = np.zeros([gr,gr])
phi_0 = np.zeros([gr,gr])
phi_0[int(gr/2),int(gr/2)] = 1

plt.rcParams["font.family"] = "Times New Roman"
text_size = 20
r = 0.8
plt.figure()
plt.imshow(np.transpose(phi_0),cmap = 'plasma',origin = 'lower',extent = [0,domain_dimension,0,domain_dimension])
bar = plt.colorbar()
bar.ax.tick_params(labelsize=text_size)
plt.xlabel('x',fontsize = text_size)
plt.ylabel('y',fontsize = text_size)
plt.xticks(fontsize=r*text_size)
plt.yticks(fontsize=r*text_size)
plt.savefig('at_%i.png'%0,dpi = 400, pad_inches = 0, bbox_inches='tight')
np.savetxt('at_%i.csv'%0,phi)


def k(phi_p_p,phi_p_n,phi_p_e,phi_p_s,phi_p_w,u,v,rho,gamma,bc,phi_out,dt):
    if bc == 0:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    elif bc == 1:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_out*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_out*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    elif bc == 2:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2 -
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    elif bc == 3:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    if bc == 4:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_out*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_out*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    elif bc == 12:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_out*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_out*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    elif bc == 32:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_p_w*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_w*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_out-phi_p_p)*2-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_p_w))-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_out*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))                                                                      
    elif bc == 34:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_out*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_out*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_out*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_p_n-phi_p_p)+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_out)*2-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_n*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))
    if bc == 14:
        if u >= 0:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_out*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_out*del_x +
                                                 rho*u*phi_p_p*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_out*del_x))
        else:
            if v>= 0:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_p_p*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_s*del_x -
                                                 rho*u*phi_p_p*del_x))
            else:
                k_r = dt*(1/(rho*del_x*del_y))*((gamma*(phi_out-phi_p_p)*2+
                                                gamma*(phi_p_e-phi_p_p)-
                                                gamma*(phi_p_p-phi_p_s)-
                                                gamma*(phi_p_p-phi_out)*2)-(rho*v*phi_out*del_x +
                                                 rho*u*phi_p_e*del_y -
                                                 rho*v*phi_p_p*del_x -
                                                 rho*u*phi_p_p*del_x))                                                                     
    return(k_r)

while t <= total_time:
    print(t)
    for i in range(gr):
        for j in range(gr):
            if BC[i,j] == 0:
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],phi_0[i+1,j],
                            phi_0[i,j-1],phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 1:     
                k0[i,j] = k(phi_0[i,j],0,phi_0[i+1,j],
                            phi_0[i,j-1],phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 2:     
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],0,
                            phi_0[i,j-1],phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 3:     
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],phi_0[i+1,j],
                            0,phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 4:     
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],phi_0[i+1,j],
                            phi_0[i,j-1],0,u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 12:     
                k0[i,j] = k(phi_0[i,j],0,0,
                            phi_0[i,j-1],phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 32:     
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],0,
                            0,phi_0[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 34:     
                k0[i,j] = k(phi_0[i,j],phi_0[i,j+1],phi_0[i+1,j],
                            0,0,u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 14:     
                k0[i,j] = k(phi_0[i,j],0,phi_0[i+1,j],
                            phi_0[i,j-1],0,u(t),v(t),rho,T,0,phi_out,dt)
    phi_1 = phi_0 + k0
    for i in range(gr):
        for j in range(gr):
            if BC[i,j] == 0:
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],phi_1[i+1,j],
                            phi_1[i,j-1],phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 1:     
                k1[i,j] = k(phi_1[i,j],0,phi_1[i+1,j],
                            phi_1[i,j-1],phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 2:     
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],0,
                            phi_1[i,j-1],phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 3:     
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],phi_1[i+1,j],
                            0,phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 4:     
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],phi_1[i+1,j],
                            phi_1[i,j-1],0,u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 12:     
                k1[i,j] = k(phi_1[i,j],0,0,
                            phi_1[i,j-1],phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 32:     
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],0,
                            0,phi_1[i-1,j],u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 34:     
                k1[i,j] = k(phi_1[i,j],phi_1[i,j+1],phi_1[i+1,j],
                            0,0,u(t),v(t),rho,T,0,phi_out,dt)
            elif BC[i,j] == 14:     
                k1[i,j] = k(phi_1[i,j],0,phi_1[i+1,j],
                            phi_1[i,j-1],0,u(t),v(t),rho,T,0,phi_out,dt)
        phi = 0.5 * (k0 + k1) + phi_0
        phi_0 = phi
    
    if ite == int(total_time/(dt*4)) or ite == int(total_time/(dt*2)) or ite == int(3*total_time/(dt*4)) or ite == int(total_time/(dt)):
        plt.figure()
        plt.imshow(np.transpose(phi),cmap = 'plasma',origin = 'lower',extent = [0,domain_dimension,0,domain_dimension])
        bar = plt.colorbar()
        bar.ax.tick_params(labelsize=text_size)
        plt.xlabel('x',fontsize = text_size)
        plt.ylabel('y',fontsize = text_size)
        plt.xticks(fontsize=r*text_size)
        plt.yticks(fontsize=r*text_size)
        plt.savefig('at_%i.png'%ite,dpi = 400, pad_inches = 0, bbox_inches='tight')        
        
        np.savetxt('at_%i.csv'%ite,phi)
    t += dt
    ite += 1
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import copy

# C_test_tdma = np.array([[1,  3,  0,  0],
#                   [2,  3,  4,  0],
#                   [0,  9,  2,  2],
#                   [0,  0,  5,  7]])
# D_test_tdma = np.array([9,  23,  2,  21])


def TDMA(C,  D):
    n = np.shape(C)[0]  # Obtaining number of nodes
    # Forward elimination
    CN = np.zeros(n-1)
    DN = np.zeros(n)
    X = np.zeros(n)
    CN[0] = C[0,  1]/C[0,  0]
    DN[0] = D[0]/C[0,  0]

    i = 1
    while i < (n - 1):
        CN[i] = C[i,  i+1]/(C[i,  i]-C[i,  i-1]*CN[i-1])
        DN[i] = (- C[i,  i-1]*DN[i-1] + D[i])/(C[i,  i]-C[i,  i-1]*CN[i-1])
        i += 1

    DN[n-1] = (- C[n-1,  n-2]*DN[n-2] + D[n-1])/(C[n-1,  n-1]-C[n-1,  n-2]*CN[n-2])

    # Backward substitution
    X[n-1] = DN[n-1]
    i = n-2
    while i >= 0:
        X[i] = - CN[i]*X[i+1] + DN[i]
        i -= 1

    return(X)


# X_test = TDMA(C_test_tdma, D_test_tdma)
# X_direct_inv = np.matmul(np.linalg.inv(C_test_tdma), np.transpose(D_test_tdma))

# C_test = np.array([[4, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
#                     [-1, 4, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0],
#                     [0, -1, 4, 0, 0, -1, 0, 0, 0, 0, 0, 0],
#                     [-1, 0, 0, 4, -1, 0, -1, 0, 0, 0, 0, 0],
#                     [0, -1, 0, -1, 4, -1, 0, -1, 0, 0, 0, 0],
#                     [0, 0, -1, 0, -1, 4, 0, 0, -1, 0, 0, 0],
#                     [0, 0, 0, -1, 0, 0, 4, -1, 0, -1, 0, 0],
#                     [0, 0, 0, 0, -1, 0, -1, 4, -1, 0, -1, 0],
#                     [0, 0, 0, 0, 0, -1, 0, -1, 4, 0, 0, -1],
#                     [0, 0, 0, 0, 0, 0, -1, 0, 0, 4, -1, 0],
#                     [0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 4, -1],
#                     [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 4]])

# D_test = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# grid_dim_test = [3, 4]

# relaxation_test = 1

# phi_init_test = 0

# phi_test = np.zeros(12)


def rearrange_order(C,  D,  phi,  N_native,  N_new):
    n = np.shape(D)[0]
    old_order = np.arange(0,  n)
    new_order = np.zeros(n)
    i = 0
    while i < N_native:
        new_order[N_new * i: N_new * (i + 1)] = old_order[i::N_native]
        i += 1
    new_order = new_order.astype(int)
    DN = np.array([D[new_order[i]] for i in old_order])
    phiN = np.array([phi[new_order[i]] for i in old_order])
    CN = np.zeros(np.shape(C))
    CNN = np.zeros(np.shape(C))
    for i in old_order:
        CN[:,  i] = C[:,  new_order[i]]
    for i in old_order:
        CNN[i, :] = CN[new_order[i], :]
    return(CNN,  DN,  phiN)

# CNT, DNT, PNT = rearrange_order(C_test, D_test, phi_test, 3, 4)


def LBL_TDMA(C,  D,  grid_dimensions,  relaxation,  phi_init, n_ite):
    ite = 0
    DT = D
    CT = C
    gdx = grid_dimensions[0]
    gdy = grid_dimensions[1]
    phi = np.reshape(phi_init,[grid_dimensions[0]*grid_dimensions[1]])
    while ite < n_ite:
        # print(ite)
        # y_sweep
        for i in range(gdy):
            D_local = np.zeros(gdx)
            C_local = np.zeros([gdx,gdx])
            D_local[0:gdx] = DT[gdx*i:gdx*(i+1)]
            C_local[0:gdx,0:gdx] = CT[gdx*i:gdx*(i+1),  gdx*i:gdx*(i+1)]
            if 0 < i < gdy - 1:
                ii = 0
                while ii < gdx:
                    D_local[ii] = D_local[ii] - (phi[gdx*(i+1)+ii]*CT[gdx*i+ii,  gdx*(i+1)+ii]
                                                            + phi[gdx*(i-1)+ii]*CT[gdx*i+ii,  gdx*(i-1)+ii])
                    ii += 1
                
            elif i == 0:
                ii = 0
                while ii < gdx:
                    
                    D_local[ii] = D_local[ii] - (phi[gdx*(i+1)+ii]*CT[gdx*i+ii,  gdx*(i+1)+ii])
                    ii += 1

            elif i == gdy - 1:
                ii = 0
                while ii < gdx:
                    D_local[ii] = D_local[ii] - (phi[gdx*(i-1)+ii]*CT[gdx*i+ii,  gdx*(i-1)+ii])
                    ii += 1
            # print('ft')
            phi[gdx*i:gdx*(i+1)] = phi[gdx*i:gdx*(i+1)] + relaxation*(TDMA(C_local,  D_local)-phi[gdx*i:gdx*(i+1)])
            # print(DT[3]) 
        #reorient for x sweep
        CT,  DT,  phi = rearrange_order(CT,  DT,  phi,  gdx,  gdy)
        #x_sweep
        for i in range(gdx):
            D_local = np.zeros(gdy)
            C_local = np.zeros([gdy,gdy])
            D_local[0:gdy] = DT[gdy*i:gdy*(i+1)]
            C_local[0:gdy,0:gdy] = CT[gdy*i:gdy*(i+1),  gdy*i:gdy*(i+1)]
            if 0 < i < gdx - 1:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i+1)+ii]*CT[gdy*i+ii,  gdy*(i+1)+ii]
                                                            + phi[gdy*(i-1)+ii]*C[gdy*i+ii,  gdy*(i-1)+ii])
                    ii += 1
            elif i == 0:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i+1)+ii]*CT[gdy*i+ii,  gdy*(i+1)+ii])
                    ii += 1

            elif i == gdx - 1:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i-1)+ii]*CT[gdy*i+ii,  gdy*(i-1)+ii])
                    ii += 1
            # print('st')
            phi[gdy*i:gdy*(i+1)] = phi[gdy*i:gdy*(i+1)] + relaxation*(TDMA(C_local,  D_local)-phi[gdy*i:gdy*(i+1)])
            # print(DT[3]) 
            
        # reorient for y sweep
        CT,  DT,  phi = rearrange_order(CT,  DT,  phi,  gdy,  gdx)
        ite += 1
    CT,  DT,  phi = rearrange_order(CT,  DT,  phi,  gdy,  gdx)
    phi = np.reshape(phi,  [gdy,  gdx])
    return(np.array(phi))

def face_vel(x,y,cs):
    # return ue, uw, un, us
    return([(x+cs/2)**2+1,(x-cs/2)**2+1,(y+cs/2)**2+1,(y-cs/2)**2+1])

def LBL_TDMA_dif(C,  DOO,  grid_dimensions,  relaxation,  phi_init, n_ite, csc, rho):
    ite = 0
    gdx = grid_dimensions[0]
    gdy = grid_dimensions[1]
    phi = np.reshape(phi_init,[grid_dimensions[0]*grid_dimensions[1]])
    D = np.zeros(np.shape(DOO)[0])
    while ite < n_ite:
        D[0:np.shape(D)[0]] = DOO[0:np.shape(DOO)[0]]
        for ii in range(2,gdx-2):
            x_pos = ii*csc+csc/2
            for jj in range(2,gdy-2):
                
                y_pos = jj*csc+csc/2
                
                n = ((jj)*gdx)+ii
                # print(n)
                defE = face_vel(x_pos, y_pos, csc)[0]*csc*rho*(((phi[n+1] + phi[n])/2) - ((phi[n+1] + phi[n-1] - 2*phi[n])/8) - phi[n])
                defN = face_vel(x_pos, y_pos, csc)[2]*csc*rho*(((phi[n+gdx] + phi[n])/2) - ((phi[n+gdx] + phi[n-gdx] - 2*phi[n])/8) - phi[n])
                defW = face_vel(x_pos, y_pos, csc)[1]*csc*rho*(((phi[n] + phi[n-1])/2) - ((phi[n] + phi[n-2] - 2*phi[n-1])/8) - phi[n-1])
                defS = face_vel(x_pos, y_pos, csc)[3]*csc*rho*(((phi[n] + phi[n-gdx])/2) - ((phi[n] + phi[n-2*gdx] - 2*phi[n-gdx])/8) - phi[n-gdx])
                defT = - defE - defN + defW + defS
                # print(defT)
                D[n] += defT
        # print(ite)
        # y_sweep
        for i in range(gdy):
            D_local = np.zeros(gdx)
            C_local = np.zeros([gdx,gdx])
            D_local[0:gdx] = D[gdx*i:gdx*(i+1)]
            C_local[0:gdx,0:gdx] = C[gdx*i:gdx*(i+1),  gdx*i:gdx*(i+1)]
            if 0 < i < gdy - 1:
                ii = 0
                while ii < gdx:
                    D_local[ii] = D_local[ii] - (phi[gdx*(i+1)+ii]*C[gdx*i+ii,  gdx*(i+1)+ii]
                                                            + phi[gdx*(i-1)+ii]*C[gdx*i+ii,  gdx*(i-1)+ii])
                    ii += 1
            elif i == 0:
                ii = 0
                while ii < gdx:
                    D_local[ii] = D_local[ii] - (phi[gdx*(i+1)+ii]*C[gdx*i+ii,  gdx*(i+1)+ii])
                    ii += 1

            elif i == gdy - 1:
                ii = 0
                while ii < gdx:
                    D_local[ii] = D_local[ii] - (phi[gdx*(i-1)+ii]*C[gdx*i+ii,  gdx*(i-1)+ii])
                    ii += 1
            phi[gdx*i:gdx*(i+1)] = phi[gdx*i:gdx*(i+1)] + relaxation*(TDMA(C_local,  D_local)-phi[gdx*i:gdx*(i+1)])
        # reorient for x sweep
        C,  D,  phi = rearrange_order(C,  D,  phi,  gdx,  gdy)
        # x_sweep
        for i in range(gdx):
            D_local = np.zeros(gdy)
            C_local = np.zeros([gdy,gdy])
            D_local[0:gdy] = D[gdy*i:gdy*(i+1)]
            C_local[0:gdy,0:gdy] = C[gdy*i:gdy*(i+1),  gdy*i:gdy*(i+1)]
            if 0 < i < gdx - 1:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i+1)+ii]*C[gdy*i+ii,  gdy*(i+1)+ii]
                                                            + phi[gdy*(i-1)+ii]*C[gdy*i+ii,  gdy*(i-1)+ii])
                    ii += 1
            elif i == 0:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i+1)+ii]*C[gdy*i+ii,  gdy*(i+1)+ii])
                    ii += 1

            elif i == gdx - 1:
                ii = 0
                while ii < gdy:
                    D_local[ii] = D_local[ii] - (phi[gdy*(i-1)+ii]*C[gdy*i+ii,  gdy*(i-1)+ii])
                    ii += 1

            phi[gdy*i:gdy*(i+1)] = phi[gdy*i:gdy*(i+1)] + relaxation*(TDMA(C_local,  D_local)-phi[gdy*i:gdy*(i+1)])
        # reorient for y sweep
        C,  D,  phi = rearrange_order(C,  D,  phi,  gdy,  gdx)
        ite += 1
        
        
    phi = np.reshape(phi,  [gdy,  gdx])
    return(np.array(phi))


# phi = LBL_TDMA(C_test,  D_test,  grid_dim_test,  relaxation_test,  phi_init_test, 10)
# phi_inv = np.matmul(np.linalg.inv(C_test), np.transpose(D_test))
# phi_inv = np.reshape(phi_inv,  [4,  3])
# print(phi, phi_inv)
# X = np.arange(0, 3)
# Y = np.arange(0, 4)
# X, Y = np.meshgrid(X, Y)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(X, Y, phi)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(X, Y, phi_inv)

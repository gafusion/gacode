import scipy.io as sio
import scipy.linalg
import os,sys
import numpy as np
import cmath,math
import matplotlib.pyplot as plt

def dodmd(rhot, tlist, ncut, ns, ng, index):
    # # Preparation
    # plt.rcParams['lines.linewidth'] = 3
    # plt.rcParams['font.size'] = 16 
    # plt.rcParams['lines.markersize'] = 10
    
    [nk, nt] = rhot.shape
    nt_list = tlist.size

    
    if ncut > nt:
        ncut = nt
        
    nh = nk*ns # the height of Y
    nl = math.floor(ncut/ns) # the length of Y
    

    # DMD
    # construct the data matrices        
    Y = np.zeros((nh, nl-1), dtype = 'complex_')
    Yt = np.zeros((nh, nl-1), dtype = 'complex_')
    
    
    for l in range(nl-1):
        for s in range(ns):
            Y[s*nk:((s+1))*nk, l] = rhot[:, l*ng+s]
            Yt[s*nk:((s+1))*nk, l] = rhot[:, l*ng+s+1]
    

    ############# DMD ##############

    # Perform Singular Value Decomposition (SVD) on X
    U, S, Vh = scipy.linalg.svd(Y, full_matrices=False)

    V = Vh.conj().T
    U_out = U
    
    # Truncate SVD to remove small singular values
    # r = np.sum(abs(Sigma) > 1.0e-6)
    # r = np.sum(abs(Sigma) > 9.0e3)
    r = np.sum(abs(S) > 4.0e4)
    r = 7
    # r = np.sum(abs(Sigma) > 5.33e25)
    U = U[:, :index]
    S = S[:index]
    Vh = Vh[:index, :]
    
    
    # Compute A_tilde
    A_tilde = np.dot(np.dot(U.conj().T, Yt), np.dot(Vh.conj().T, np.diag(1/S)))

    

    # Compute eigenvalues and eigenvectors of A_tilde
    D, W = scipy.linalg.eig(A_tilde)
    
    # Compute DMD modes
    Psi = np.dot(np.dot(np.dot(Yt, Vh.conj().T), np.diag(1/S)), W)
    
    # Compute DMD mode amplitudes
    b = np.dot(np.linalg.pinv(Psi), Y[:, 0])

    # lmda = np.diag(D) # discrete-time eigenvalues
    dt = tlist[1] - tlist[0] # the time step size
    omega = np.log(D)/dt # continuous-time eigenvalues

    return Psi, omega, b, U, S, V



if __name__ == "__main__":


    # Arguments passed
    print("\nName of Python script:", sys.argv[0])

    print("\naky: ",str(sys.argv[1]))
    print("delt:",str(sys.argv[2]))     # temporal resolution
    print("ntheta:",str(sys.argv[3]))   # spatial resolution
    print("nperiod:",str(sys.argv[4]))
    print("\nDMD-order:",str(sys.argv[5]))
    print("# of SVD-truncation:",str(sys.argv[6]))
    print("snapshot_start_ratio:",str(sys.argv[7])) # means the start_point = len(t)//snapshot_start_ratio
    print("snapshot_end_ratio:",str(sys.argv[8]))   # means the end_point = len(t)*(100-snapshot_end_ratio)//100
    # print("\n")


    DTYPE='float64'

    aky     = sys.argv[1]
    delt    = np.float64(sys.argv[2])    # temporal resolution
    ntheta  = sys.argv[3]                # spatial resolution
    nperiod = sys.argv[4]                # another sparial resolution parameter, but 3 is sufficient, 6/8/12 do not seem to be helpful

    ns      = int(sys.argv[5])  # DMD-order
    index   = int(sys.argv[6])  # number of SVD-truncation
    snapshot_start_ratio = int(sys.argv[7])    # means the start_point = len(t)//snapshot_start_ratio
    snapshot_end_ratio   = int(sys.argv[8])    # means the end_point = len(t)*(100-snapshot_end_ratio)//100



    t = np.loadtxt("data/time.txt--aky="+str(aky)+"-delt="+str(delt)+"-ntheta="+str(ntheta)+"_nperiod="+str(nperiod), dtype=DTYPE)
    theta = np.loadtxt("data/theta.txt--aky="+str(aky)+"-delt="+str(delt)+"-ntheta="+str(ntheta)+"_nperiod="+str(nperiod), dtype=DTYPE)
    phi_re = np.loadtxt("data/phi_real_time.txt--aky="+str(aky)+"-delt="+str(delt)+"-ntheta="+str(ntheta)+"_nperiod="+str(nperiod), dtype=DTYPE)
    phi_im = np.loadtxt("data/phi_imag_time.txt--aky="+str(aky)+"-delt="+str(delt)+"-ntheta="+str(ntheta)+"_nperiod="+str(nperiod), dtype=DTYPE)
    # dens_amp = np.loadtxt("data/dens_time.txt--aky="+str(aky)+"-delt="+str(delt)+"-ntheta="+str(ntheta)+"_nperiod="+str(nperiod), dtype=DTYPE)

    data0 = np.zeros(len(phi_re), dtype=complex)

    data0.real = phi_re
    data0.imag = phi_im


    Nt = len(t)
    Ntheta = len(theta)
    Nspecies = 2 # number of species, currently only main ions and electrons, hence 2 
    j = Ntheta/4

    delta_t = delt*10 # 1.0 # 0.25
    tt = delta_t - (t[-1]-t[-2]) # 0.2


    # tmp = np.zeros([Nt, Nspecies, Ntheta], dtype=DTYPE)
    tmp = np.zeros([Nt, Ntheta], dtype=DTYPE)

    # data = dens_amp.reshape(tmp.shape)

    # data1 = np.transpose(data[:,0,:])
    # data2 = np.transpose(data[:,1,:])

    # data = data1

    data = data0.reshape(tmp.shape)


    # ############# POD ################


    # snapshot_start = len(t)//snapshot_start_ratio # 150 # 10
    # # snapshot_end   = len(t)-snapshot_end_point # 2100 # 40 # 90 # 140
    # snapshot_end   = len(t)*(100-snapshot_end_point)//100

    # POD_U, POD_S, POD_Vh = scipy.linalg.svd(data[snapshot_start:snapshot_end,:], full_matrices=False)
    # # POD_U, POD_S, POD_Vh = scipy.linalg.svd(data[:,:], full_matrices=False)

    # for i in range(10):
    #     print(POD_S[i])
    
    #     plt.title("omega (real part) = "+str(POD_S[i])) 
    #     # plt.xlabel("Theta") 
    #     # plt.ylabel("dens") 
    #     # plt.ylim(0,7E-4)
    #     # plt.plot   (theta_set2, cgyro[:], color ="purple") 
    #     # plt.scatter(theta, X[:], color ="orange") 
    #     plt.scatter(theta, np.absolute(POD_Vh[i,:].T), color ="orange") 
    #     plt.show()

    # # C = np.dot(data[snapshot_start:snapshot_end,:].T,data[snapshot_start:snapshot_end,:])#*1.0/(float(Nt)-1.0)

    # # POD_lambda, POD_phi = np.linalg.eig(C)

    # # idx = POD_lambda.argsort()[::-1]   
    # # POD_lambda = POD_lambda[idx]
    # # POD_phi = POD_phi[:,idx]


    # # # for i in range(len(POD_lambda)-10,len(POD_lambda)):
    # # for i in range(20):
    # #     print(POD_lambda[i])
    

    # #     plt.title("omega (real part) = "+str(POD_lambda[i])) 
    # #     # plt.xlabel("Theta") 
    # #     # plt.ylabel("dens") 
    # #     # plt.ylim(0,7E-4)
    # #     # plt.plot   (theta_set2, cgyro[:], color ="purple") 
    # #     # plt.scatter(theta, X[:], color ="orange") 
    # #     plt.scatter(theta, np.absolute(POD_phi[:,i]), color ="orange") 
    # #     plt.show()

    # exit()








    data = np.transpose(data)

    print("\n")
    print("\n+++++++ DMD input data +++++++")

    print("\nN_theta",len(data),", N_time",len(data[0]),"\n")



    ############ DMD ############

    # ns = 1      # DMD order
    ng = ns # 2 # ns
    # index = 10  # number of SVD truncation

    snapshot_start = len(t)//snapshot_start_ratio # 150 # 10
    # snapshot_end   = len(t)-snapshot_end_point # 2100 # 40 # 90 # 140
    snapshot_end   = len(t)*(100-snapshot_end_ratio)//100
    # snapshot_end   = -(snapshot_end_point//100)

    print("t is ","[",t[0],",",t[-1],"] seconds")
    print("T is ","[",0,",",len(t),"]\n")

    print("snapshot_start")
    print("t =",t[snapshot_start],"second ")
    print("T =",snapshot_start)
    print("snapshot_end")
    print("t =",t[snapshot_end],"second")
    print("T =",snapshot_end,"\n")
    print("delta_t =",delta_t,"second","\n")

    Phi, omega, b, U, Sigma, V = dodmd(data[:,snapshot_start:], t[snapshot_start:], snapshot_end-snapshot_start, ns, ng, index)



    print("\n+++++++ DMD results +++++++\n")
    print("Sigma")

    for i in range(len(Sigma)):
        print(Sigma[i])
    
    print("=======================")


    print("Omega (real/imag)")

    
    # omega=np.log(omega)/delta_t

    for i in range(len(omega)):
        print(omega.real[i],omega.imag[i])
    
    print("=======================")


    extrapolation = len(t)-snapshot_end # 8
    tt = delta_t - (t[-1]-t[-2]) # The last time step is slightly smaller than the normal time step, so it needs to be adjusted using "tt".

    X_mostunstable = np.zeros(len(theta))

    for i in range(len(omega)):

        if omega[i]>0.0:

            X = Phi[:,i]*cmath.exp(omega[i]*(delta_t*(snapshot_end-snapshot_start+extrapolation-1)-tt))*b[i]

            plt.title("omega (real part) = "+str(omega.real[i])) 
            plt.xlabel("Theta") 
            plt.ylabel("abs|Phi|") 
            # plt.ylim(0,7E-4)
            # plt.plot   (theta_set2, cgyro[:], color ="purple") 
            # plt.scatter(theta, X[:], color ="orange") 
            plt.scatter(theta, np.absolute(Phi[0:len(theta),i]), color ="orange") 
            plt.show()

            X_mostunstable = X_mostunstable + X[0:len(theta)]

            # plt.title("omega (real part) = "+str(omega.real[mode])) 
            # plt.xlabel("Theta") 
            # plt.ylabel("dens") 
            # # plt.ylim(0,7E-4)
            # # plt.plot   (theta_set2, cgyro[:], color ="purple") 
            # plt.scatter(theta, X_mostunstable[:], color ="orange") 
            # # plt.scatter(theta, Phi.real[:,mode], color ="orange") 
            # plt.show()



    # X_mostunstable = np.dot(Phi[:,:],np.exp(omega[:]*(delta_t*(snapshot_end-snapshot_start+extrapolation-1)-tt))*b[:])

    cgyro = data[:,snapshot_end+extrapolation-1]

    plt.title("Trajectory at t = "+str(t[-1])+" second (the last one)") 
    plt.xlabel("Theta") 
    plt.ylabel("abs|Phi|") 
    plt.plot   (theta, np.absolute(cgyro[:]), color ="blue") 
    plt.scatter(theta, np.absolute(X_mostunstable[:]), color ="orange") 
    # plt.plot   (theta, cgyro[:], color ="blue") 
    # plt.scatter(theta, X_mostunstable[:], color ="orange") 
    plt.show()

    # # plt.title(filename) 
    # # plt.xlabel("X axis") 
    # # plt.ylabel("n") 
    # # plt.plot   (xn, cgyro_imag[:], color ="blue") 
    # # plt.scatter(xn, X_mostunstable.imag[:], color ="orange") 
    # # plt.show()

    exit()

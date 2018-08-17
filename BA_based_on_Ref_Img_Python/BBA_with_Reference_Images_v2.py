# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 00:42:44 2018

@author: admins
"""
import numpy as np
from numpy import linalg as LA
import AT

data_path = "data_BBA/"

std_IP = 0.00452646 * 0.1          # 0.1px, for ref.
std_IP2 = 0.00122331 * 0.1      # 0.1px, for smart
std_GPS = 0.0001                # 0.0001m, for ref.
std_INS = np.deg2rad(0.0001)    # 0.0001deg, for ref
std_GPS2 = 1                    # 1m, for smart
std_INS2 = np.deg2rad(1.)       # 1deg, for smart

## Import Data Files
print("(1) Import Data Files")
# Import IOP of Reference & Target Images
IO_ref = AT.import_IOP_with_Name(data_path + "IO_ref.txt", "Reference Images")   # Reference - UAV Images
IO_tar = AT.import_IOP_with_Name(data_path + "IO_tar.txt", "Target Images")   # Target - Smartphone Images
IO = np.asarray([IO_ref[0:3], IO_tar])

# Import EOP of Reference Images to be used for EOP Constraints
EO_c = AT.import_EOP_with_Name (data_path + "EO_c_ref.txt", "Reference Images")
no_EO_c = EO_c.shape[0]

# Import Initial EOP of Target Images
EO_i_tar = AT.import_EOP_with_Name (data_path + "EO_i_tar.txt", "Target Images")
EO_i_tar[:,0] = no_EO_c + 1
EO_i_tar[:,1:4] = EO_i_tar[:,1:4] #+ np.random.randn(1,3) * std_GPS2
EO_i_tar[:,4:7] = EO_i_tar[:,4:7] #+ np.random.randn(1,3) * std_INS2
no_EO_tar = EO_i_tar.shape[0]

# Set Initial EOP for All Images
EO_i = np.vstack((EO_c, EO_i_tar))

# Import Image Points
IP = AT.import_IP (data_path + "IP_m.txt")

## Preprocessing - checking, selecting, indexing
print("(2) Preprocessing - checking, selecting, indexing")

# Find Image Indexes and No. Images
id_IM = np.unique(IP[:,0])
no_IM = len(id_IM)
I2C = {}
for i in id_IM:
    if np.count_nonzero(EO_i_tar[:,0]==i)==0:
        I2C[int(i)] = 0
    else:
        I2C[int(i)] = 1

# Find GP Indexes and No. GP
id_GP = np.unique(IP[:,1])
no_GP = len(id_GP)

# Compute No. Images of Each GP Appearing
cnt_GP = np.zeros(no_GP)
for i in range(no_GP):
    cnt_GP[i] = np.count_nonzero(IP[:,1]==id_GP[i])

# Identify GPs inclusive and exclusive
id_GP_inc = id_GP[np.nonzero(cnt_GP>=2)]
id_GP_inc = np.sort(id_GP_inc)
no_GP_inc = len(id_GP_inc)
id_GP_exc = id_GP[np.nonzero(cnt_GP<2)]
print(" No. GP inc. : {0:d}".format(no_GP_inc))

# Identify IPs inclusive
id_IP_inc = np.nonzero(IP[:,1]==id_GP_inc[0])
for i in range(1,no_GP_inc):
    id_IP_inc = np.append(id_IP_inc, np.nonzero(IP[:,1]==id_GP_inc[i]))
IP_inc = IP[id_IP_inc,:]
no_IP_inc = IP_inc.shape[0]
print(" No. IP inc. : {0:d}".format(no_IP_inc))

# Identify Images Inclusive
id_IM_inc = np.unique(IP_inc[:,0])
no_IM_inc = len(id_IM_inc)
print(" No. IM inc. : {0:d}".format(no_IM_inc))

# Identify Corresponding EO Index
I2E = {}
for i in range(no_IM_inc):
    ix = np.flatnonzero(EO_i[:,0]==id_IM_inc[i])
    if len(ix)==1:
        I2E[int(id_IM_inc[i])] = ix[0]
    else:
        print("Image ID ({0:d}) is not found from EO_i".format(int(id_IM_inc[i])))

## Computing Initial Approximation for GP
print("(3) Computing Initial Approximation for GP")

# Compute Rotational Matrix
Re = np.zeros([3,3,no_IM_inc])
for i in range(no_IM_inc):
    Re[:,:,i] = AT.Rot3D(EO_i[I2E[id_IM_inc[i]],4:7])

# Compute GP Estimates
#for i in range(no_GP_inc):
pid = [0,0]
iid = [0,0]
cid = [0,0]
IPc = np.zeros([2,3])
GPc = np.zeros([2,3])
GP_i = np.zeros([no_GP_inc,4])
Est_stat = np.zeros([no_GP_inc,5])
vct = np.zeros((no_GP_inc))
cnd = np.zeros((no_GP_inc))
for i in range(no_GP_inc):
    idx = np.flatnonzero(IP_inc[:,1]==id_GP_inc[i])

# Find two tie points of the maximum baseline
    no_IM_cur = len(idx)
    max_BL = 0.
    for n in range(no_IM_cur-1):
        I1 = IP_inc[idx[n],0]
        for m in range(n+1,no_IM_cur):
            I2 = IP_inc[idx[m],0]
            BL = LA.norm( EO_i[I2E[int(I2)],1:4] - EO_i[I2E[int(I1)],1:4] )
            if BL > max_BL:
                max_BL = BL;
                max_BL_idx = np.array([n,m], dtype=np.int64)
 
# Determine the ground point corresponding to the tie points    
    for n in range(2):
        pid[n] = idx[max_BL_idx[n]]
        iid[n] = int(IP_inc[pid[n],0])
        cid[n] = I2C[iid[n]]
        IPc[n,0:2] = IP_inc[pid[n],2:4]-IO[cid[n],0:2]
        IPc[n,2] = -IO[cid[n],2]

    A = np.zeros([3,2])
    A[:,0] = np.matmul(np.transpose(Re[:,:,I2E[iid[0]]]),IPc[0,:])
    A[:,1] = np.matmul(np.transpose(-Re[:,:,I2E[iid[1]]]),IPc[1,:])
    y = np.zeros([3,1])
    y[:,0] = np.transpose(EO_i[I2E[iid[1]],1:4] - EO_i[I2E[iid[0]],1:4])

#    N = np.matmul(np.transpose(A), A)
#    cnd[i] = LA.cond(N)
    
    [kt,vct[i],et] = AT.LSE(A,y)  
    
    for n in range(2):
        GPc[n,:] = kt[n] * np.matmul(np.transpose(Re[:,:,I2E[iid[n]]]), IPc[n,:]) + np.transpose(EO_i[I2E[iid[n]],1:4])
    GP_i[i,0] = id_GP_inc[i]
    GP_i[i,1:4] = (GPc[0,:]+GPc[1,:])/2.
    Est_stat[i,0] = id_GP_inc[i]
    Est_stat[i,1] = vct[i]
    Est_stat[i,2:5] = np.transpose(et)

AT.export_GP_with_Name (GP_i, data_path + "GP_i.txt", "Initial")
AT.export_IA_Summary(Est_stat, data_path + "IA_summary.txt")

# Selecting Effective GPs

ix = np.flatnonzero(vct<1)
id_GP_inc2 = id_GP_inc[ix]
no_GP_inc2 = len(id_GP_inc2)
GP_i2 = GP_i[ix,:]
print(" No. GP inc filtered : {0:d}".format(no_GP_inc2))

ix = np.empty((0), dtype=np.int)
for i in range(no_GP_inc2):
    ix = np.append(ix, np.flatnonzero(IP_inc[:,1]==id_GP_inc2[i]))
IP_inc2 = IP_inc[ix,:]
no_IP_inc2 = IP_inc2.shape[0]

# Find Image Indexes and No. Images
id_IM_inc2 = np.unique(IP_inc2[:,0])
no_IM_inc2 = len(id_IM_inc2)
I2C2 = {}
for i in id_IM_inc2:
    if np.count_nonzero(EO_i_tar[:,0]==i)==0:
        I2C2[int(i)] = 0
    else:
        I2C2[int(i)] = 1

ix = np.lexsort((IP_inc2[:,1], IP_inc2[:,0]));
temp = IP_inc2[ix,:]
IP_inc2 = temp

AT.export_IP_with_Name(IP_inc2, data_path + 'IP_mf.txt', 'Filtered')

ix = np.empty((0), dtype=np.int)
for i in range(no_IM_inc2):
    ix = np.append(ix, np.flatnonzero(EO_c[:,0]==id_IM_inc2[i]))
EO_c2 = EO_c[ix,:]
no_EO_c2 = EO_c2.shape[0]

AT.export_EO_with_Name(EO_c2, data_path + 'EO_cf.txt', 'Filtered')

## Performing Bundle Block Adjustment
print("(4) Performing Bundle Block Adjustment")

# Kt_e : initial approximations of the EO of the images
Kt_e = np.zeros((no_IM_inc2*6,1));
for n in range(no_IM_inc2):
    ix = np.flatnonzero(EO_i[:,0]==id_IM_inc2[n])
    if len(ix)==1:
        Kt_e[n*6:(n+1)*6,0] = EO_i[ix,1:7]

# Kt_g : initial approximations of the ground points
Kt_g = np.zeros((no_GP_inc2*3,1))
for n in range(no_GP_inc2):
    ix = np.flatnonzero(GP_i2[:,0]==id_GP_inc2[n])
    if len(ix)==1:
        Kt_g[n*3:(n+1)*3,0] = GP_i2[ix,1:4]

kt = np.zeros((no_IM_inc2*6+no_GP_inc2*3,1))

## Perform AT
delta = 1e-6
cnst_ws = 1

for k in range(100):
    print( " Iteration {0:d}".format(k) );

# The rotational matrix and its derivatives
    Re = np.zeros((3,3,no_IM_inc2))
    dRe_om = np.zeros((3,3,no_IM_inc2))
    dRe_ph = np.zeros((3,3,no_IM_inc2))
    dRe_kp = np.zeros((3,3,no_IM_inc2))
    for i in range(no_IM_inc2):
        Re[:,:,i] = AT.Rot3D(Kt_e[i*6+3:(i+1)*6,0])
        [dRe_om[:,:,i], dRe_ph[:,:,i], dRe_kp[:,:,i]] = AT.dRot3D(Kt_e[i*6+3:(i+1)*6,0])

# Initialize the design matrix and observations
    Aei = np.zeros((2,6))
    Agi = np.zeros((2,3))
    Nee = np.zeros((no_IM_inc2*6, no_IM_inc2*6))
    Neg = np.zeros((no_IM_inc2*6, no_GP_inc2*3))
    Ngg = np.zeros((no_GP_inc2*3, no_GP_inc2*3))
    Ce = np.zeros((no_IM_inc2*6,1))
    Cg = np.zeros((no_GP_inc2*3,1))
    y = np.zeros((no_IP_inc2*2,1))
    ye = np.zeros((no_EO_c2*6,1))
    Om_y = 0.
    Om_ye = 0.;
    
# Set the design matrix and observations
    imi = -1
    gpi = -1
    for n in range(no_IP_inc2):
#    for n in range(0,50):
        imi = np.flatnonzero(id_IM_inc2==IP_inc2[n,0])[0] 
        gpi = np.flatnonzero(id_GP_inc2==IP_inc2[n,1])[0] 
        
        GC = Kt_g[gpi*3:(gpi+1)*3,0] - Kt_e[imi*6:imi*6+3,0]
        ND = np.matmul(Re[:,:,imi], np.transpose(GC))
        F0 = np.zeros([2,1])
        if I2C2[IP_inc2[n,0]]==0:    # Reference Images    
            P_nu = - ND[0:2] / ND[2]
            r = np.sqrt(P_nu[0]**2 + P_nu[1]**2)
            radial = 1 + IO_ref[5]*r**2 + IO_ref[6]*r**4 + IO_ref[7]*r**6 + IO_ref[8]*r**8;
            x_d = P_nu[0] * radial + (IO_ref[10] * (r**2 + 2*P_nu[0]**2) + 2*IO_ref[9]*P_nu[0]*P_nu[1]) * (1 + IO_ref[11]*r**2 + IO_ref[12]*r**4)
            y_d = P_nu[1] * radial + (IO_ref[9] * (r**2 + 2*P_nu[1]**2) + 2*IO_ref[10]*P_nu[0]*P_nu[1]) * (1 + IO_ref[11]*r**2 + IO_ref[12]*r**4)
            F0[0] = IO_ref[0] + (IO_ref[2] + IO_ref[3]) * x_d + IO_ref[4] * y_d
            F0[1] = IO_ref[1] + IO_ref[2] * y_d
            FL = IO_ref[2]
            Pi = 1./std_IP**2
        else:   # Target Images (Smartphone)
            F0 = IO_tar[0:2] - IO_tar[2] / ND[2] * ND[0:2]
            F0 = np.expand_dims(F0,axis=1)
            FL = IO_tar[2]
            Pi = 1./std_IP2**2
        
        dND = np.zeros([3,9])
        dND[:,0:3] = -Re[:,:,imi]
        dND[:,3] = np.matmul(dRe_om[:,:,imi],GC)
        dND[:,4] = np.matmul(dRe_ph[:,:,imi],GC)
        dND[:,5] = np.matmul(dRe_kp[:,:,imi],GC)
        dND[:,6:9] = Re[:,:,imi]

        CM = np.asarray([[-ND[2], 0, ND[0]], [0, -ND[2], ND[1]]])
        Aei = FL / ND[2]**2 * np.matmul(CM, dND[:,0:6])
        Agi = FL / ND[2]**2 * np.matmul(CM, dND[:,6:9])
        yi = np.expand_dims(IP_inc2[n,2:4],axis=1) - F0
        y[n*2:(n+1)*2,0] = np.squeeze(yi)

        Nee[imi*6:(imi+1)*6,imi*6:(imi+1)*6] = Nee[imi*6:(imi+1)*6,imi*6:(imi+1)*6] + np.matmul(np.transpose(Aei),Aei)*Pi
        Ngg[gpi*3:(gpi+1)*3,gpi*3:(gpi+1)*3] = Ngg[gpi*3:(gpi+1)*3,gpi*3:(gpi+1)*3] + np.matmul(np.transpose(Agi),Agi)*Pi
        Neg[imi*6:(imi+1)*6,gpi*3:(gpi+1)*3] = Neg[imi*6:(imi+1)*6,gpi*3:(gpi+1)*3] + np.matmul(np.transpose(Aei),Agi)*Pi
        Ce[imi*6:(imi+1)*6,0] = Ce[imi*6:(imi+1)*6,0] + np.squeeze(np.matmul(np.transpose(Aei),yi))*Pi
        Cg[gpi*3:(gpi+1)*3,0] = Cg[gpi*3:(gpi+1)*3,0] + np.squeeze(np.matmul(np.transpose(Agi),yi))*Pi
        Om_y = Om_y + np.squeeze(np.matmul(np.transpose(yi),yi))*Pi

    Pe = np.zeros((6,6))
    Pe[0:3,0:3] = np.eye(3) * cnst_ws / std_GPS**2;
    Pe[3:6,3:6] = np.eye(3) * cnst_ws / std_INS**2;

    Kei = np.zeros((6,6))
    yei = np.zeros((6,1))
    for ni in range(no_EO_c2):
        ix = np.flatnonzero(id_IM_inc2==EO_c2[ni,0])
        if len(ix)==1:
            imi = ix[0]
            yei = EO_c2[ni,1:7] - Kt_e[imi*6:(imi+1)*6,0]
            ye[ni*6:(ni+1)*6,0] = np.squeeze(yei)
            Nee[imi*6:(imi+1)*6,imi*6:(imi+1)*6] = Nee[imi*6:(imi+1)*6,imi*6:(imi+1)*6] + Pe
            Ce[imi*6:(imi+1)*6,0] = Ce[imi*6:(imi+1)*6,0] + np.matmul(Pe,yei)
            Om_ye = Om_ye + np.squeeze(np.matmul(np.transpose(yei),np.matmul(Pe,yei)))
        else:
            print("EO Constraints of ID {0:d} is not found".format(int(EO_c2[ni,0])))
    iNgg = np.zeros(Ngg.shape)
    for n in range(no_GP_inc2):
        iNgg[n*3:(n+1)*3,n*3:(n+1)*3] = LA.inv(Ngg[n*3:(n+1)*3,n*3:(n+1)*3])

    Nr = (Nee - np.matmul(Neg, np.matmul(iNgg, np.transpose(Neg))))
    iNr = LA.inv(Nr)
    kt_e = np.matmul(iNr, (Ce - np.matmul(Neg, np.matmul(iNgg, Cg))))
    kt_g = np.matmul(iNgg, (Cg - np.matmul(np.transpose(Neg), kt_e)))
    Kt_e = Kt_e + kt_e
    Kt_g = Kt_g + kt_g
    kt[0:no_IM_inc2*6,0] = np.squeeze(kt_e)
    kt[no_IM_inc2*6:no_IM_inc2*6+no_GP_inc2*3,0] = np.squeeze(kt_g)

    if LA.norm(kt)<delta:
        break

no_IT = k

EO_e = np.zeros([no_IM_inc2,7])
for n in range(no_IM_inc2):
    EO_e[n,0] = id_IM_inc2[n]
    EO_e[n,1:7] = Kt_e[n*6:(n+1)*6,0]
GP_e = np.zeros([no_GP_inc2,4])
for n in range(no_GP_inc2):
    GP_e[n,0] = id_GP_inc2[n]
    GP_e[n,1:4] = Kt_g[n*3:(n+1)*3,0]

AT.export_EO_with_Name (EO_e, data_path + "EO_e.txt", "Estimated")
AT.export_GP_with_Name (GP_e, data_path + "GP_e.txt", "Estimated")

Om = Om_y + Om_ye
vct = Om / (no_IP_inc2*2 + no_EO_c*6 - no_IM_inc2*6 - no_GP_inc2*3)

iNee = iNr
iNeg = np.matmul(-iNr, np.matmul(Neg, iNgg))
iNge = np.transpose(iNeg)
iNgg = iNgg - np.matmul(iNge, np.matmul(Neg, iNgg))

Dkt_e = vct * iNee
Skt_e = np.sqrt(np.diag(Dkt_e))
Dkt_g = vct * iNgg
Skt_g = np.sqrt(np.diag(Dkt_g))
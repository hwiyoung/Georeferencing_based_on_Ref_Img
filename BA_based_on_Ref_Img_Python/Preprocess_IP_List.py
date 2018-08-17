# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 00:14:37 2018

@author: admins
"""

# Preprocess IP, EO for AT
# 1. Converting coordinate system from ICS to CCS
# 2. Arranging index of image and IPs
# 3. Generate EO corresponding to IP
# 13/12/2017
# Hwiyoung Kim / Impyeong Lee

# Input data
#   Tie-points
#   Image list
# Output data
#   Arranged tie-points

import numpy as np
import AT

data_path = "data_Pre/"

# Reference Images
pixCnt_ref = np.asarray([5304, 7952])   # px
pixSize_ref = 0.00452646                # mm/px

## Read TP & Image list & EO
# iPhone 7
pixCnt_acq = np.asarray([3024, 4032])   # px
pixSize_acq = 0.00122331                # mm/px

TP = AT.import_IP(data_path + 'TP.txt')
#TP = AT.import_IP('ResultText/' + 'TP.txt')
imgList = AT.import_ImageList(data_path + 'ImgList.txt')
for n in range(len(imgList)):
    imgList[n] = imgList[n].replace('\n', '')
idxAcqImg = len(imgList)  # Index of the acquired image

idxRef = AT.import_IndexRef (data_path + 'indexRef.txt')

EO_ref = AT.import_EOP(data_path + 'EO_ref_true.txt')
EO_ref[:,4:7] = np.deg2rad(EO_ref[:,4:7])

EO_acq = AT.import_EOP(data_path + 'EO_smart_true.txt')
EO_acq[:,4:7] = np.deg2rad(EO_acq[:,4:7])

## 1. Converting coordinate system from ICS to CCS
# Find the row which is not the acquired image
A = TP[TP[:,0].argsort(),:]
flag = np.min(np.nonzero(A[:,0] != 0));

# Transform from ICS to CCS in acquired images
x_acq = (A[0:flag,2] - pixCnt_acq[1]/2) * pixSize_acq;
y_acq = -(A[0:flag,3] - pixCnt_acq[0]/2) * pixSize_acq;
A[0:flag,2] = x_acq;
A[0:flag,3] = y_acq;

# Transform from ICS to CCS in reference images
x_ref = (A[flag::,2] - pixCnt_ref[1]/2) * pixSize_ref;
y_ref = -(A[flag::,3] - pixCnt_ref[0]/2) * pixSize_ref;
A[flag::,2] = x_ref
A[flag::,3] = y_ref

## 2. Arranging index of image and IPs
# Change the index of the acquired image
A[0:flag,0] = idxAcqImg
#A1 = A[A[:,1].argsort(),:]
B = A[np.lexsort((A[:,1], A[:,0])),:]

# Make the list of EO corresponding to IP
listEO = np.zeros((idxAcqImg-1,), dtype=np.int)
for n in range(idxAcqImg-1):
    listEO[n] = idxRef[imgList[n+1]]
    
# Change index of images in IP like 1, 2, 3 ...
com = np.unique(B[:,0])
noMatchImg = len(com)   # # of actually matched images with acq. img
for n in range(noMatchImg):
    idxIP = np.nonzero(B[:,0]==com[n])
    B[idxIP,0] = n+1
    
AT.export_IP(B, data_path + 'IP_m.txt')

## 3. Generate EOs corresponding to IP
# Generate EOs
selEO = np.zeros((noMatchImg-1,7))
for n in range(noMatchImg-1):  # Exclude the acquired image(last line)
    idxEO = listEO[int(com[n])-1]
    selEO[n,0] = n+1
    selEO[n,1:7] = EO_ref[idxEO-1,1:7]

'''
selEO[noMatchImg-1,0] = B[B.shape[0]-1,0]
selEO[noMatchImg-1,1:7] = EO_acq[0,1:7]
'''
#AT.export_EO(selEO, data_path + 'EO_m.txt')
AT.export_EO(selEO, 'data_BBA/' + 'EO_c_ref.txt')
AT.export_EO(EO_acq, 'data_BBA/' + 'EO_i_tar.txt')
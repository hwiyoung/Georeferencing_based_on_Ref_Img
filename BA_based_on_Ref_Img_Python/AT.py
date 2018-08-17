# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 00:14:33 2018

@author: Impyeong
"""
import os
import numpy as np
from numpy import linalg as LA

def check_file_exist(file_name):
    if os.path.exists(file_name):
        fExist = True
    else:
        fExist = False
        print("File not Found : " + file_name)
    return fExist

def import_IOP(file_name):
    print("Import IOP from " + file_name)
    if check_file_exist(file_name):
        IO = np.loadtxt(file_name)
        print(" No. IOP: {0:d}".format(len(IO)))
        return IO
    
def import_IOP_with_Name(file_name, image_name):
    print("Import IOP of " + image_name + " from " + file_name)
    if check_file_exist(file_name):
        IO = np.loadtxt(file_name)
        print(" No. IOP: {0:d}".format(len(IO)))
        return IO
        
def import_IP(file_name):
    print("Import IP from " + file_name)
    if check_file_exist(file_name):
        IP = np.loadtxt(file_name, ndmin=2);
        print(" No. IP: {0:d}".format(IP.shape[0]))
        return IP

def import_EOP(file_name):
    print("Import EOP from " + file_name)
    if check_file_exist(file_name):
        EO = np.loadtxt(file_name, ndmin=2)
        print(" No. EOP sets: {0:d}".format(EO.shape[0]))
        return EO

def import_EOP_with_Name(file_name, image_name):
    print("Import EOP of " + image_name + " from " + file_name)
    if check_file_exist(file_name):
        EO = np.loadtxt(file_name, ndmin=2)
        print(" No. EOP sets: {0:d}".format(EO.shape[0]))
        return EO

def import_GP(file_name):
    print("Import GP from " + file_name)
    if check_file_exist(file_name):
        GP = np.loadtxt(file_name, ndmin=2)
        print(" No. GP sets: {0:d}".format(GP.shape[0]))        
        return GP

def import_ImageList(file_name):
    print("Import Image List from " + file_name)
    if check_file_exist(file_name):
        fp = open(file_name, "r")
        ImgList = fp.readlines()
        return ImgList

def import_IndexRef(file_name):
    print("Import Index of Reference Images from " + file_name)
    if check_file_exist(file_name):
        fp = open(file_name, "r")
        strList = fp.readlines()
        idxRef = {}
        for n in range(len(strList)):
            strTemp = strList[n].split()
            if len(strTemp)==2:
                idxRef[strTemp[0]] = int(strTemp[1])
        return idxRef   

def find(A,v):
    for i in range(len(A)):
        if A[i]==v:
            return i
        
def Rot3D(RA):
    om = RA[0]
    ph = RA[1]
    kp = RA[2]

    Rx = np.array([[1, 0, 0], [0, np.cos(om), np.sin(om)], [0, -np.sin(om), np.cos(om)]])
    Ry = np.array([[np.cos(ph), 0, -np.sin(ph)], [0, 1, 0], [np.sin(ph), 0, np.cos(ph)]])
    Rz = np.array([[np.cos(kp), np.sin(kp), 0], [-np.sin(kp), np.cos(kp), 0], [0, 0, 1]])

#    R = Rx * Ry * Rz;
    R = np.matmul(Rz, np.matmul(Ry,Rx))
    return R

def dRot3D(RA):
    om = RA[0]
    ph = RA[1]
    kp = RA[2]

    Rx = np.array([[1, 0, 0], [0, np.cos(om), np.sin(om)], [0, -np.sin(om), np.cos(om)]])
    Ry = np.array([[np.cos(ph), 0, -np.sin(ph)], [0, 1, 0], [np.sin(ph), 0, np.cos(ph)]])
    Rz = np.array([[np.cos(kp), np.sin(kp), 0], [-np.sin(kp), np.cos(kp), 0], [0, 0, 1]])

    dRx_om = np.array([[0, 0, 0], [0, -np.sin(om), np.cos(om)], [0, -np.cos(om), -np.sin(om)]])
    dRy_ph = np.array([[-np.sin(ph), 0, -np.cos(ph)], [0, 0, 0], [np.cos(ph), 0, -np.sin(ph)]])
    dRz_kp = np.array([[-np.sin(kp), np.cos(kp), 0], [-np.cos(kp), -np.sin(kp), 0], [0, 0, 0]])

    dR_om = np.matmul(Rz, np.matmul(Ry,dRx_om))
    dR_ph = np.matmul(Rz, np.matmul(dRy_ph,Rx))
    dR_kp = np.matmul(dRz_kp, np.matmul(Ry,Rx))

    return [dR_om, dR_ph, dR_kp]

def LSE(A,y):
    At = np.transpose(A)
    N = np.matmul(At,A)
    c = np.matmul(At,y)
    iN = LA.inv(N)
    kt = np.matmul(iN,c)
    yt = np.matmul(A,kt)
    et = y - yt
    Om = np.matmul(np.transpose(et), et)
    vct = Om / (A.shape[0]-LA.matrix_rank(A))
    return [kt, vct, et]

def export_IP(IP, file_name):
    print("Export IP to " + file_name)
    outFp = open(file_name, "w")
    for i in range(IP.shape[0]):
        outFp.writelines("{0:d}\t{1:d}\t{2:11.5f}\t{3:11.5f}\n".format(int(IP[i,0]), int(IP[i,1]), IP[i,2], IP[i,3]))
    outFp.close()
	
def export_IP_with_Name(IP, file_name, name):
    print("Export " + name + " IP to " + file_name)
    outFp = open(file_name, "w")
    for i in range(IP.shape[0]):
        outFp.writelines("{0:d}\t{1:d}\t{2:11.5f}\t{3:11.5f}\n".format(int(IP[i,0]), int(IP[i,1]), IP[i,2], IP[i,3]))
    outFp.close()

def export_GP(GP, file_name):
    print("Export GP to " + file_name)
    print(" No. GP: {0:d}".format(GP.shape[0]))        
    outFp = open(file_name, "w")
    for i in range(GP.shape[0]):
        outFp.writelines("{0:d}\t{1:11.3f}\t{2:11.3f}\t{3:11.3f}\n".format(int(GP[i,0]), GP[i,1], GP[i,2], GP[i,3]))
    outFp.close()
    
def export_GP_with_Name(GP, file_name, name):
    print("Export " + name + " GP to " + file_name)
    print(" No. GP: {0:d}".format(GP.shape[0]))        
    outFp = open(file_name, "w")
    for i in range(GP.shape[0]):
        outFp.writelines("{0:d}\t{1:11.3f}\t{2:11.3f}\t{3:11.3f}\n".format(int(GP[i,0]), GP[i,1], GP[i,2], GP[i,3]))
    outFp.close()
        
def export_EO(EO, file_name):
    print("Export EO to " + file_name)
    print(" No. EOP sets: {0:d}".format(EO.shape[0]))        
    outFp = open(file_name, "w")
    for i in range(EO.shape[0]):
        outFp.writelines("{0:d}\t{1:11.3f}\t{2:11.3f}\t{3:11.3f}\t{4:11.6f}\t{5:11.6f}\t{6:11.6f}\n".format(int(EO[i,0]), EO[i,1], EO[i,2], EO[i,3], EO[i,4], EO[i,5], EO[i,6]))
    outFp.close()
    
def export_EO_with_Name(EO, file_name, name):
    print("Export " + name + " EO to " + file_name)
    print(" No. EOP sets: {0:d}".format(EO.shape[0]))        
    outFp = open(file_name, "w")
    for i in range(EO.shape[0]):
        outFp.writelines("{0:d}\t{1:11.3f}\t{2:11.3f}\t{3:11.3f}\t{4:11.6f}\t{5:11.6f}\t{6:11.6f}\n".format(int(EO[i,0]), EO[i,1], EO[i,2], EO[i,3], EO[i,4], EO[i,5], EO[i,6]))
    outFp.close()

def export_IA_Summary(Est_stat, file_name):
    print("Export Summary of Initial Approximation to " + file_name)
    outFp = open(file_name, "w")
    outFp.writelines("Statistics for GP initial appoximation\n")
    outFp.writelines("GP	        vct	         dx	         dy	         dz\n")
    for i in range(Est_stat.shape[0]):
        outFp.writelines("{0:d}\t{1:11.3f}\t{2:11.3f}\t{3:11.3f}\t{4:11.3f}\n".format(int(Est_stat[i,0]), Est_stat[i,1], Est_stat[i,2], Est_stat[i,3], Est_stat[i,4]))
    outFp.close()

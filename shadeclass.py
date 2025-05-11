# -*- coding: utf-8 -*-
"""
Created on Tue May  6 18:11:26 2025

@author: chr
"""
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as scopt
import scipy.fftpack as scfft

import calendar
import time
import datetime
import os

def length(v1):
    return(np.sqrt(np.sum(np.square(v1),axis=0)))
def normalis(v1):
    return (v1/length(v1))
def skalar(v1,v2):
    return(np.sum(v1*v2,axis=0))
def crosspro(v1,v2):
    return(np.cross(v1, v2))

def testsonnenstand(t,geoBreite ,geoLange):
    h=np.sin(t)/5
    azimut=t
    
    e_Sx=np.cos(h)*np.cos(azimut)
    e_Sy=np.cos(h)*np.sin(azimut)
    e_Sz=np.sin(h)
    E_S=np.array((e_Sx,e_Sy,e_Sz))
    return(E_S)
class Sphere: 
    def __init__(self,P_center,radius,transparency=0):
        #an sphere thas blocl light like a tree
        self.p_c=P_center
        self.r=radius
        self.trans=transparency
    def checkrays(self,P0_iab,eS_it):
        #checks if a light ray from P0 in directetion eS go throw this object
        
        # vector d:= P-c-Po
        d_iab=self.p_c[:, np.newaxis,np.newaxis]-P0_iab
        #split d into 2parts parralle an orthogonal to eS
        eSd_tab=skalar(d_iab[:, np.newaxis,:,:],eS_it[:,:, np.newaxis,np.newaxis])
        dParalle_itab=eSd_tab[ np.newaxis,:,:,:] * eS_it[:,:, np.newaxis,np.newaxis]
        dorthogonal_itab=d_iab[:, np.newaxis,:,:]-dParalle_itab
        
        passspher_tab=length(dorthogonal_itab)>self.r
        bevoresphere=eSd_tab[:, np.newaxis,:,:]<0
        notinSphere=length(d_iab[:, np.newaxis,:,:])>self.r
        passed=np.logical_and(notinSphere,np.logical_or(bevoresphere,passspher_tab))
        #print('d',d_tab)
        return(passed)
class Rectangle:
    def __init__(self,P_bottomleft,P_bottomright,P_upleft):
        self.orign=P_bottomleft
        self.l1=length(P_bottomright-self.orign)
        self.e1=normalis(P_bottomright-self.orign)
        
        self.l2=length(P_upleft-self.orign)
        self.e2=normalis(P_upleft-self.orign)
        
        self.e3=crosspro(self.e1, self.e2)
class photovoltaikfläche:
    def __init__(self,rectangle,datadensity=(3,2),bifa=False):
        self.rect=rectangle
        a = np.linspace(0,self.rect.l1,datadensity[0])
        b = np.linspace(0,self.rect.l2,datadensity[1])
        i=np.array((0,1,2))
        A, B = np.meshgrid(a, b,indexing='ij')
        #2Dimensonal
        #ri_iab=.orign[:,,] #https://numpy.org/doc/stable/user/basics.broadcasting.html
        self.ri_iab=self.rect.orign[:, np.newaxis,np.newaxis]\
            +self.rect.e1[:, np.newaxis,np.newaxis]*A\
            +self.rect.e2[:, np.newaxis,np.newaxis]*B

        
    def calcshadow(self,objekt,t):
        eS=testsonnenstand(t, 0, 0)
        return(objekt.checkrays(self.ri_iab,eS))
        
p1=np.array((0,-10,-10))
p2=np.array((0,5,-10))
p3=np.array((0,-10,5))
p4=np.array((-0,0,1))
v1=np.array((0,3,4))
t=(np.array((0,0.01,0.1,1)))

rechteck1=Rectangle(p1,p2,p3)
mesflache1=photovoltaikfläche(rechteck1, (30,30))
tree1=Sphere(p4, 1)


shadow1=mesflache1.calcshadow(tree1,t)


for i in range(len(t)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(shadow1[i,:,:])#,vmin=0, vmax=5*np.pi/180)
    fig.savefig('abstand'+str(i))
    fig.show()
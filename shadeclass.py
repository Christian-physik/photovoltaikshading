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
        directline_iab=self.p_c[:, np.newaxis,np.newaxis]-P0_iab
        direcParale_itab=skalar(directline_iab[:, np.newaxis,:,:],eS_it[:,:, np.newaxis,np.newaxis])
        vecd=directline_iab[:, np.newaxis,:,:]-direcParale_itab
        d_tab=length(vecd)
        print('d',d_tab)
class Rectangle:
    def __init__(self,P_bottomleft,P_bottomright,P_upleft):
        self.orign=P_bottomleft
        self.l1=length(P_bottomright-self.orign)
        self.e1=normalis(P_bottomright-self.orign)
        
        self.l2=length(P_upleft-self.orign)
        self.e2=normalis(P_upleft-self.orign)
        
        self.e3=crosspro(self.e1, self.e2)
class photovoltaikfläche:
    def __init__(self,rectangle,taxis,datadensity=(3,2),bifa=False):
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
        print(self.ri_iab)
        print(length(self.ri_iab))
        
    def calcshadow(self,objekt):
        eS=testsonnenstand(t, 0, 0)
        objekt.checkrays(self.ri_iab,eS)
        
p1=np.array((0,0,0))
p2=np.array((0,5,0))
p3=np.array((0,0,2))
p4=np.array((-5,0,1))
v1=np.array((0,3,4))
t=(np.array((10,11,12,13)))
rechteck1=Rectangle(p1,p2,p3)
mesflache1=photovoltaikfläche(rechteck1, t,(6,10))
tree1=Sphere(p4, 1)
print(mesflache1.calcshadow(tree1))
print('v1length',length(v1))
print('v1*p3',skalar(v1, p3))
print('v1xp3',crosspro(v1, p3))
print('testsonne',testsonnenstand(t, 0, 0))
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
def Sonnenstand(t,geobreitedeg,geolangedeg,azimutcorr=0):
    tJ2000=t-time.mktime((2000,1,1, 12,0,0, 0,0,0))
    n_J2000=tJ2000/24.0/3600.0

    n0_J2000=np.floor(n_J2000)+0.5  #+0.5 da von oh aus gerechnet
    nBruch=n_J2000-n0_J2000
    tBruch=np.mod(t,24*3600)

    geobreite=geobreitedeg*np.pi/180
    geolange=geolangedeg*np.pi/180
    
    #mittlere eliptikische Länge
    Ldeg=(280.460+0.9856474*n_J2000)
    Ldeg=np.mod(Ldeg,360)
    L=Ldeg*np.pi/180
    
    gdeg=(357.528+0.9856003*n_J2000)
    gdeg=np.mod(gdeg,360)

    g=gdeg*np.pi/180
    #e=0.0167    #*0
    #e=edeg*np.pi/180
    #ekliplaenge=L+2*e*np.sin(g)+5/4*np.square(e)*np.sin(2*g)
    ekliplaenge_deg= Ldeg+1.915* np.sin(g)+0.01997* np.sin(2*g)
    ekliplaenge=ekliplaenge_deg*np.pi/180
    
    #tilt of eathaxis
    ee=(23.439-0.0000004*n_J2000)*np.pi/180
    
    #Rektaszension
    rekta=np.arctan(np.cos(ee)*np.tan(ekliplaenge))
#    if np.cos(ekliplaenge)<0:
#        rekta+=np.pi
    rekta[np.cos(ekliplaenge)<0]+=np.pi
    deklin=np.arcsin(np.sin(ee)*np.sin(ekliplaenge))
    T0=n0_J2000/36525
    #LocalHourAngel of Arie
    #teta=((6.697376+12+2400.05134*T0+1.002738*hBruch)*15 +geolange)*np.pi/180
    teta=((6.697376+2400.05134*T0+1.002738*nBruch*24)*15 +geolange)*np.pi/180
    #teta=(6.697376+2400.05134*T0+1.002738*t)*np.pi/12 +geolange
    
    tau=teta-rekta #LHA sun
    teta=np.mod(teta,2*np.pi)
    azimutsouth=np.arctan2(np.sin(tau),(np.cos(tau)*np.sin(geobreite)-np.tan(deklin)*np.cos(geobreite)))
    h=np.arcsin(np.cos(deklin)*np.cos(tau)*np.cos(geobreite)+np.sin(deklin)*np.sin(geobreite))
    azimut=azimutsouth+azimutcorr
    
    e_Sx=np.cos(h)*np.cos(azimut)
    e_Sy=np.cos(h)*np.sin(azimut)
    e_Sz=np.sin(h)
    e_S=np.array((e_Sx,e_Sy,e_Sz))
    fig=plt.figure()
    plt.plot(t,h,label='h')
    plt.plot(t,azimut,label='azimut')
    for k in (0,1,2):
        plt.plot(t,e_S[k,:],label=str(k))
    plt.legend()
    plt.savefig('eS')
        
    print('eS',e_S)
    return(e_S)

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
        print('d',d_iab.shape)
        print('eS_it',eS_it.shape)
        eSd_tab=skalar(d_iab[:, np.newaxis,:,:],eS_it[:,:, np.newaxis,np.newaxis])
        dParalle_itab=eSd_tab[ np.newaxis,:,:,:] * eS_it[:,:, np.newaxis,np.newaxis]
        dorthogonal_itab=d_iab[:, np.newaxis,:,:]-dParalle_itab
        
        passspher_tab=length(dorthogonal_itab)>self.r
        bevoresphere=eSd_tab<0
        notinSphere=length(d_iab[:, np.newaxis,:,:])>self.r
        passed=np.logical_and(notinSphere,np.logical_or(bevoresphere,passspher_tab))
        print('shape')
        print(passspher_tab.shape)
        print(bevoresphere.shape)
        print(notinSphere.shape)
        print(passed.shape)
        shadefactor=passed.astype(float)
        #print('d',d_tab)
        return(shadefactor)
class Rectangle:
    def __init__(self,P_bottomleft,P_bottomright,P_upleft):
        self.orign=P_bottomleft
        self.l1=length(P_bottomright-self.orign)
        self.e1=normalis(P_bottomright-self.orign)
        
        self.l2=length(P_upleft-self.orign)
        self.e2=normalis(P_upleft-self.orign)
        if skalar(self.e1, self.e2)!=0 :
            print('e_1,e_2 are not aorthogonal')
        
        self.e3=crosspro(self.e1, self.e2)
 
class photovoltaikfläche:
    def __init__(self,rectangle,datadensity=(3,2),bifa=False):
        self.rect=rectangle
        self.bifa=bifa
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
        p0=1
        #direction to the sun
        eS_it=Sonnenstand(t, 49, -8) 
        #skalar with normal
        eSd_t=skalar(self.rect.e3[:, np.newaxis],eS_it)#[:,:, np.newaxis,np.newaxis]
        print(' eSd_tabsize', eSd_t)
        p0_t=np.abs(np.copy(eSd_t))*p0#for bifaziat
        p0_t[np.logical_and(eSd_t<0,np.logical_not(self.bifa))]=0
        p0_tab=p0_t[:, np.newaxis,np.newaxis]
        print('posize',p0_tab)
        
        p0_tab=p0_tab*objekt.checkrays(self.ri_iab,eS_it)
        return(p0_tab)
        
        return(objekt.checkrays(self.ri_iab,eS_it))
        
p1=np.array((-10,-10,0))
p2=np.array((-10,10,0))
p3=np.array((10,-10,0))
p4=np.array((-4,0,1))
v1=np.array((0,3,4))
t=(np.array((0,0.01,0.1,1)))
tstart=time.mktime((2000,5,1, 1,0,0, 0,0,0))
tend=time.mktime((2000,5,1, 23,0,0, 0,0,0))
t=np.arange(tstart,tend,3600 )

rechteck1=Rectangle(p1,p3,p2)
mesflache1=photovoltaikfläche(rechteck1, (60,60))
tree1=Sphere(p4, 1)


shadow1=mesflache1.calcshadow(tree1,t)


for i in range(len(t)):
    if i <24:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(shadow1[i,:,:],vmin=0, vmax=0.5)
        fig.savefig('Plots/abstand'+str(i))
        fig.show()
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 18:11:26 2025

@author: chr
"""
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
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
    # fig=plt.figure()
    # plt.plot(t,h,label='h')
    # plt.plot(t,azimut,label='azimut')
    # for k in (0,1,2):
    #     plt.plot(t,e_S[k,:],label=str(k))
    # plt.legend()
    # plt.savefig('eS')
        

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

        eSd_tab=skalar(d_iab[:, np.newaxis,:,:],eS_it[:,:, np.newaxis,np.newaxis])
        dParalle_itab=eSd_tab[ np.newaxis,:,:,:] * eS_it[:,:, np.newaxis,np.newaxis]
        dorthogonal_itab=d_iab[:, np.newaxis,:,:]-dParalle_itab
        
        passspher_tab=length(dorthogonal_itab)>self.r
        bevoresphere=eSd_tab<0
        notinSphere=length(d_iab[:, np.newaxis,:,:])>self.r
        passed=np.logical_and(notinSphere,np.logical_or(bevoresphere,passspher_tab))

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
            print(self.orign)
            print(self.e1)
            print(self.e2)
        
        self.e3=crosspro(self.e1, self.e2)
    def checkrays(self,P0_iab,eS_it):
        #checks if a light ray from P0 in directetion eS go throw this object
        
        # vector d:= P-c-Po
        d_iab=self.orign[:, np.newaxis,np.newaxis]-P0_iab
        
        e3eS_t=skalar(self.e3[:, np.newaxis], eS_it)
        
        #assume e1 and e2 are orthogonal
        #let Phit be the point where the ray hit the rectangle
        #chose lambda a1 a2 to be so that:
        #Phit=P0+eS*lambda
        #    =orign+e1*a1+e2*a2
        #d=eS*lambda-e1*a1-e2*a2
        
        #multiply with e3:
        #    d*e3= e3eS lambda
        de3_ab=skalar(d_iab, self.e3[:, np.newaxis,np.newaxis])
        lambda_tab=de3_ab[np.newaxis,:,:]/e3eS_t[:, np.newaxis,np.newaxis]
        
        #multiply with e1
        #d*e1= lambda*eS*e1-a1
        #a1=-(d-lambda*eS)*e1
        dlambda_itab=d_iab[:, np.newaxis,:,:]-lambda_tab[ np.newaxis,:,:,:]*eS_it[:,:, np.newaxis,np.newaxis]
        a1_tab=-skalar(dlambda_itab,self.e1[:, np.newaxis,np.newaxis,np.newaxis])
        a2_tab=-skalar(dlambda_itab,self.e2[:, np.newaxis,np.newaxis,np.newaxis])
        
        #check if the ray doesnt hit rectangle
        passe1_tab=np.logical_or(a1_tab<0,a1_tab>self.l1)
        passe2_tab=np.logical_or(a2_tab<0,a2_tab>self.l2)
        passrect_tab=np.logical_or(passe1_tab,passe2_tab)
        bevorrect=lambda_tab<0
        passed_tab=np.logical_or( passrect_tab,bevorrect)
        
        return(passed_tab)
            
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

        
    def calcshadow(self,objekts,t,eS_it):
        p0=1
        #direction to the sun
        eS_it=Sonnenstand(t, 49, -8) 
        #skalar with normal
        eSe3_t=skalar(self.rect.e3[:, np.newaxis],eS_it)#[:,:, np.newaxis,np.newaxis]

        p0_t=np.abs(np.copy(eSe3_t))*p0#for bifaziat
        p0_t[np.logical_and(eSe3_t<0,np.logical_not(self.bifa))]=0
        p0_tab=p0_t[:, np.newaxis,np.newaxis]

        
        for objekt in objekts:
            if not objekt is self:
                p0_tab=p0_tab*objekt.checkrays(self.ri_iab,eS_it)
                
        self.lastcalc_tab=p0_tab
        self.lastt=t
        return(p0_tab)
    
class house:
    def __init__(self, geobreitedeg,geolangedeg,azimutcorr):
        self.geobreitedeg=geobreitedeg
        self.geolangedeg=geolangedeg
        self.azimutcorr=azimutcorr
        self.allrect=[]
        self.allsphere=[]
        self.allsolar=[]
    def addsphere(self,P_center,radius):
        newspere=Sphere(np.array(P_center),radius)
        self.allsphere.append(newspere)
        return(newspere)
    def addrect(self,P_bottomleft,P_bottomright,P_upleft):
        newrect=Rectangle(np.array(P_bottomleft), np.array(P_bottomright),np.array( P_upleft))
        self.allrect.append(newrect)
        return(newrect)
    def addsolarrect(self,P_bottomleft,P_bottomright,P_upleft,datadensity=(30,20),bifa=False):
        newrect=Rectangle(np.array(P_bottomleft), np.array(P_bottomright),np.array( P_upleft))
        self.allrect.append(newrect)
        newsolar=photovoltaikfläche(newrect,datadensity,bifa)
        self.allsolar.append(newsolar)
    
    def calcallshadows(self,t):
        allobjekts=self.allrect+self.allsphere
        eS_it=Sonnenstand(t, self.geobreitedeg, self.geolangedeg,self.azimutcorr)
        for solarect in self.allsolar:
            solarect.calcshadow(allobjekts,t,eS_it)
    def plot(self,name='3dplot'):
        #https://matplotlib.org/stable/gallery/mplot3d/box3d.html#sphx-glr-gallery-mplot3d-box3d-py
        #https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
        #https://matplotlib.org/stable/gallery/mplot3d/surface3d_3.html
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111, projection='3d')
        xs, ys, zs, datas=[],[],[],[]#np.array(())
        for solar in self.allsolar:
            x, y, z =solar.ri_iab[0,:,:],solar.ri_iab[1,:,:],solar.ri_iab[2,:,:]
            xs.append(x.flatten())
            ys.append(y.flatten())
            zs.append(z.flatten())
            #sum taxis'
            datapointperh=1
            data=np.mean(solar.lastcalc_tab,axis=0)
            datas.append(data.flatten())
            #colors=cmap=cm.coolwarm(data)
            #ax.scatter(x, y, z, c=data)
            #surf=ax.plot_surface(x, y, z,color=colors,  cmap=cm.coolwarm,linewidth=0, antialiased=False)
            
            # colors = np.empty(X.shape, dtype=str)
            # for y in range(ylen):
            #     for x in range(xlen):
            #         colors[y, x] = colortuple[(x + y) % len(colortuple)]
        # Customize the z axis.
        xflat=np.concatenate(xs)
        yflat=np.concatenate(ys)
        zflat=np.concatenate(zs)
        dataflat=np.concatenate(datas)
        ax.scatter(xflat, yflat, zflat, c=dataflat)
        ax.set_ylim(-0.01, 5.01)
        ax.set_zlim(-0.01, 5.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')

        # Add a color bar which maps values to colors.
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        fig.savefig(name)

    
        
        
p1=np.array((-10,-10,0))
p2=np.array((-10,10,0))
p3=np.array((10,-10,0))
p4=np.array((-4,-3,2))
p4b=np.array((-4,-3,0))

p5=np.array((-0,-0,2))
p6=np.array((8,0,2))
p7=np.array((0,0,5))
v1=np.array((-4,3,4))
t=(np.array((0,0.01,0.1,1)))
i=5
tstart=time.mktime((2000,1,1, 0+i,0,0, 0,0,0))
tend=time.mktime((2000,1,1, 1+i,0,0, 0,0,0))
t=np.arange(tstart,tend,3600 )

tstart=time.mktime((2000,3,1, 8,0,0, 0,0,0))
tend=time.mktime((2000,3,1, 18,0,0, 0,0,0))
t2=np.arange(tstart,tend,3600 )


testhouse=house(0.49, 8, -90*np.pi/180)
#testhouse.addrect((-200,-200,0), (200,-200,0), (-200,200,0))
# testhouse.addsolarrect((0,2,0), (0,0,0), (0,2,4))
# testhouse.addsolarrect((0,0,0), (2,0,0), (0,0,4))
# #testhouse.addsolarrect((0,0,4), (0,0,0), (2,0,4))
# testhouse.addsolarrect((2,0,0), (5,0,0), (2,0,2))
# testhouse.addsolarrect((2,0,2), (5,0,2), (2,1,3))
# testhouse.addsolarrect((5,2,2), (2,2,2), (5,1,3))
# testhouse.addsolarrect((2,0,2), (2,2,2), (2,0,4))
# testhouse.addsolarrect((0,0,4), (2,0,4), (0,2,4))

testhouse.addsolarrect((5,0,0), (5,4,0), (5,0,3))
# testhouse.addsolarrect((5,2,0), (2,2,0), (5,2,2))
# testhouse.addsolarrect((2,2,0), (0,2,0), (2,2,4))

testhouse.addrect((5.5,1,0), (5.5,1.5,0), (5.5,1,6))
testhouse.addsphere((5.5,2,1), 0.3)
testhouse.addsphere((5.5,0.5,1), 0.3)

testhouse.calcallshadows(t)
testhouse.plot()
# for i in range(24):
#     tstart=time.mktime((2000,1,1, 0+i,0,0, 0,0,0))
#     tend=time.mktime((2000,1,1, 1+i,0,0, 0,0,0))
#     t2=np.arange(tstart,tend,3600/6 )
#     testhouse.calcallshadows(t2)
#     testhouse.plot('Plots/3dplot_t'+str(i))

# rechteck1=Rectangle(p1,p3,p2)
# mesflache1=photovoltaikfläche(rechteck1, (60,60))
# tree1=Sphere(p4, 1)
# tree2=Sphere(p4b, 1)
# wall1=Rectangle(p5,p6,p7)


#shadow1=mesflache1.calcshadow(tree1,t)
#shadow1=mesflache1.calcshadow(wall1,t)
# shadow1=mesflache1.calcshadow((tree1,tree2,wall1),t)


# for i in range(len(t)):
#     if i <24:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.imshow(shadow1[i,:,:],vmin=0, vmax=0.5)
#         fig.savefig('Plots/shadow'+str(i))
#         fig.show()


# shadow2=mesflache1.calcshadow((tree1,tree2,wall1),t2)
# for i in range(len(t)):
#     if i <24:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.imshow(shadow2[i,:,:],vmin=0, vmax=0.5)
#         fig.savefig('Plots/shadow'+str(30+i))
#         fig.show()
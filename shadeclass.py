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
    
class Rectangle:
    def __init__(self,P_bottomleft,P_bottomright,P_upleft):
        self.orign=P_bottomleft
        self.l1=length(P_bottomright-self.orign)
        self.e1=normalis(P_bottomright-self.orign)
        
        self.l2=length(P_upleft-self.orign)
        self.e2=normalis(P_upleft-self.orign)
        
        self.e3=crosspro(self.e1, self.e2)
        
p1=np.array((0,0,0))
p2=np.array((0,1,0))
p3=np.array((0,0,2))
v1=np.array((0,3,4))
rechteck1=Rectangle(p1,p2,p3)
print('v1length',length(v1))
print('v1*p3',skalar(v1, p3))
print('v1xp3',crosspro(v1, p3))
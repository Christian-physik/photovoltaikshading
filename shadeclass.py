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

class Rectangle:
    def __init__(self,P_bottomleft,P_bottomright,P_upleft):
        self.orign=P_bottomleft
        self.l1=np.sqrt(np.sum(np.square(P_bottomright-self.orign)))
        self.e1=(P_bottomright-self.orign)/self.l1
        
        self.l2=np.sqrt(np.sum(np.square(P_upleft-self.orign)))
        self.e2=(P_upleft-self.orign)/self.l2
        #self.e3=e1 x e2
        
p1=np.array((0,0,0))
p2=np.array((0,1,0))
p3=np.array((0,0,2))
rechteck1=Rectangle(p1,p2,p3)
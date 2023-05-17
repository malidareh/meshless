# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 08:01:44 2021

@author: com
"""
import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
import math
import scipy.spatial
from scipy.spatial import distance
from scipy.spatial import KDTree
from itertools import combinations
import sys
from scipy import spatial

dmins=[]
dc=[]
def min_ndl_spacing(filterpoints):
    for vorpoints in filterpoints:
       vorpoints1=np.copy(vorpoints) 
       filterpoints1=np.copy(filterpoints) 
       indx=np.argwhere (vorpoints1==filterpoints1) 
       filterpointsremoved=np.delete(filterpoints1, indx[0,0], 0)
       clsest_nde=closest_node(vorpoints1, filterpointsremoved)
       dmin = math.hypot(clsest_nde[0] - vorpoints1[0], clsest_nde[1] - vorpoints1[1])
       dmins.append(dmin)

    return dmins
       

def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]


def avg_ndl_spacing() :
    for dsupport in dmins:
       dc1=3*dsupport
  
       dc.append(dc1)
       
    return dc

def dist(p1, p2):
    (x1, y1), (x2, y2) = p1, p2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

dc=[]
def avg1_ndl_spacing(filterpoints,pts):
    tree = KDTree(filterpoints)
    dd, ii = tree.query([pts], k=8 )
    avg_distance=np.mean(dd[0][1:])
    dc=2.5*avg_distance
    #print(dd,ii,avg_distance1)

    return dc

suportpoints=[]
def pts_subdmn(filterpoints1,dc):
   tree = spatial.KDTree(filterpoints1)
   
   for i in range(0,len(dc)) :
      inexpts1=sorted(tree.query_ball_point(filterpoints1[i], dc[i]))
      suportpoints.append(inexpts1)
   return suportpoints


# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 23:53:46 2017

@author: SETUP COMPUTER
"""
import numpy as np
from math import exp
def  Weight(wtype, para, di, dmi): 
     r = np.abs(di) / dmi 
     if (di >= 0.0): 
          drdx = 1.0/dmi; 
     else: 
          drdx = -1.0/dmi; 
       # EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r 
     if (wtype == 'GAUSS'): 
          w,dwdr,dwdrr,dwdrrr = Gauss(para,r); 
     elif (wtype == 'QUART'): 
          w,dwdr,dwdrr,dwdrrr =Quartic(r); 
     elif (wtype == 'SPLIN'): 
          w,dwdr,dwdrr,dwdrrr=Spline(r); 
     elif (wtype == 'CSRBF'): 
          w,dwdr,dwdrr,dwdrrr=CSRBF2(r); 
     else: 
           print('Invalid type of weight function.'); 
     dwdx  = dwdr * drdx; 
     dwdxx = dwdrr * drdx * drdx; 
     dwxdxx=dwdrrr * drdx * drdx*drdx;
     return w, dwdx, dwdxx,dwxdxx  
  

def  Gauss(para,r): 
      if (r>1.0): 
         w     = 0.0; 
         dwdr  = 0.0; 
         dwdrr = 0.0; 
         dwdrrr=0;
      else: 
        b2 = para*para; 
        r2 = r*r; 
        eb2 = exp(-b2); 
        w = (exp(-b2*r2) - eb2) / (1.0 - eb2); 
        dwdr  = -2*b2*r*exp(-b2*r2) / (1.0 - eb2); 
        dwdrr = -2*b2*exp(-b2*r2)*(1-2*b2*r2) / (1.0 - eb2); 
        dwdrrr=0
      return w,dwdr,dwdrr,dwdrrr        
        
def Quartic(r): 
   if (r>1.0): 
      w     = 0.0; 
      dwdr  = 0.0; 
      dwdrr = 0.0; 
      dwdrrr=0
   else: 
      w     = 1-6*r**2+8*r**3-3*r**4; 
      dwdr  = -12*r+24*r**2-12*r**3; 
      dwdrr = -12+48*r-36*r**2; 
      dwdrrr=0
   return w,dwdr,dwdrr,dwdrrr 
     
def Spline(r): 
  if (r>1.0): 
   w     = 0.0; 
   dwdr  = 0.0; 
   dwdrr = 0.0; 
   dwdrrr=0;
  elif (r<=0.5): 
   w     = 2/3 - 4*r**2 + 4*r**3; 
   dwdr  = -8*r + 12*r**2; 
   dwdrr = -8 + 24*r; 
   dwdrrr=24;
  else: 
   w     = 4/3 - 4*r + 4*r**2 - 4*r**3/3; 
   dwdr  = -4 + 8*r -4*r**2; 
   dwdrr = 8 - 8*r; 
   dwdrrr=-8;
  return w,dwdr,dwdrr,dwdrrr
 
def CSRBF2(r): 
  if (r>1.0): 
   w     = 0.0; 
   dwdr  = 0.0; 
   dwdrr = 0.0; 
   dwdrrr=0
  else: 
    w = (1-r)**6*(6+36*r+82*r**2+72*r**3+30*r**4+5*r**5); 
    dwdr  = 11*r*(r+2)*(5*r**3+15*r**2+18*r+4)*(r-1)**5; 
    dwdrr = 22*(25*r**5+100*r**4+142*r**3+68*r**2-16*r-4)*(r-1)**4; 
    dwdrrr=0
  return w,dwdr,dwdrr,dwdrrr
# """""""""""""""""""""""""""""""""""""
#""""""""""""""""""""""""""""""""""""""
#""""""""""""""""""""""""""""""""""""""
def  Weight2D(X,XI,wtype, para, r, dmi,i,j): 
     if(r==0): 
        drdx=0
        drdy=0
        drdxx=0
        drdyy=0
        drdxxx=0
        drdyyy=0
        drdxy=0
        drdxyy=0
        drdyxx=0
     else: 
        drdx = (X[j][0]-XI[i][0])/(dmi**2*r)
        drdy=(X[j][1]-XI[i][1])/(dmi**2*r)
        drdxx=(r-drdx*(X[j][0]-XI[i][0]))/(dmi**2*r**2)
        drdyy=(r-drdy*(X[j][1]-XI[i][1]))/(dmi**2*r**2)
        drdxxx=((drdx-drdxx*(X[j][0]-XI[i][0])-drdx)*r**2-2*r*drdx*(r-drdx*(X[j][0]-XI[i][0])))/(dmi**2*r**4)
        drdyyy=((drdy-drdyy*(X[j][1]-XI[i][1])-drdy)*r**2-2*r*drdy*(r-drdy*(X[j][1]-XI[i][1])))/(dmi**2*r**4)
        drdxy=(-drdy*(X[j][0]-XI[i][0]))/(dmi**2*r**2)
        drdxyy=((-drdyy*(X[j][0]-XI[i][0])*r**2)+2*drdy*r*drdy*(X[j][0]-XI[i][0]))/(dmi**2*r**4)
        drdyxx=((-drdxx*(X[j][1]-XI[i][1])*r**2)+2*drdx*r*drdx*(X[j][1]-XI[i][1]))/(dmi**2*r**4)
         
       # EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r 
     if (wtype == 'GAUSS'): 
          w,dwdr,dwdrr,dwdrrr = Gauss(para,r); 
     elif (wtype == 'QUART'): 
          w,dwdr,dwdrr,dwdrrr =Quartic(r); 
     elif (wtype == 'SPLIN'): 
          w,dwdr,dwdrr,dwdrrr=Spline(r); 
     elif (wtype == 'CSRBF'): 
          w,dwdr,dwdrr,dwdrrr=CSRBF2(r); 
     else: 
           print('Invalid type of weight function.'); 
     dwdx  = dwdr * drdx
     dwdy= dwdr * drdy
     dwdxx=dwdrr*drdx*drdx+dwdr*drdxx
     dwdyy=dwdrr*drdy*drdy+dwdr*drdyy
     dwdxy=dwdrr*drdx*drdy+dwdr*drdxy
     dwdyx=dwdxy
     dwdxxx=dwdrrr*drdx*drdx*drdx+dwdrr*2*drdxx*drdx+dwdrr*drdxx*drdx+dwdr*drdxxx
     dwdyyy=dwdrrr*drdy*drdy*drdy+dwdrr*2*drdyy*drdy+dwdrr*drdyy*drdy+dwdr*drdyyy
     dwdyxx=dwdrrr*drdx*drdx*drdy+dwdrr*2*drdxy*drdx+dwdrr*drdxx*drdy+dwdr*drdyxx
     dwdxyy=dwdrrr*drdy*drdy*drdx+dwdrr*2*drdxy*drdy+dwdrr*drdyy*drdx+dwdr*drdxyy
   
     return w, dwdx, dwdy, dwdxx, dwdyy,dwdxy,dwdyx,dwdxxx,dwdyyy,dwdyxx,dwdxyy
  

def  Gauss(para,r): 
      if (r>1.0): 
         w     = 0.0; 
         dwdr  = 0.0; 
         dwdrr = 0.0; 
         dwdrrr=0;
      else: 
        b2 = para*para; 
        r2 = r*r; 
        eb2 = exp(-b2); 
        w = (exp(-b2*r2) - eb2) / (1.0 - eb2); 
        dwdr  = -2*b2*r*exp(-b2*r2) / (1.0 - eb2); 
        dwdrr = -2*b2*exp(-b2*r2)*(1-2*b2*r2) / (1.0 - eb2); 
        dwdrrr=8*b2*b2*r*exp(-b2*r2)*(1-b2*r2) / (1.0 - eb2); 
      return w,dwdr,dwdrr,dwdrrr        
        
def Quartic(r): 
   if (r>1.0): 
      w     = 0.0; 
      dwdr  = 0.0; 
      dwdrr = 0.0; 
      dwdrrr=0
   else: 
      w     = 1-6*r**2+8*r**3-3*r**4; 
      dwdr  = -12*r+24*r**2-12*r**3; 
      dwdrr = 48*r-36*r**2; 
      dwdrrr=-72*r
   return w,dwdr,dwdrr,dwdrrr 
     
def Spline(r):     
  if (r>1.0): 
   w     = 0.0; 
   dwdr  = 0.0; 
   dwdrr = 0.0; 
   dwdrrr=0;
  elif (r<=0.5): 
   w     = 2/3 - 4*r**2 + 4*r**3; 
   dwdr  = -8*r + 12*r**2;  
   dwdrr = -8 + 24*r; 
   dwdrrr=24;
  else: 
   w     = 4/3 - 4*r + 4*r**2 - 4*r**3/3; 
   dwdr  = -4 + 8*r -4*r**2; 
   dwdrr = 8 - 8*r; 
   dwdrrr=-8;
  return w,dwdr,dwdrr,dwdrrr
 
def CSRBF2(r): 
  if (r>1.0): 
   w     = 0.0; 
   dwdr  = 0.0; 
   dwdrr = 0.0; 
   dwdrrr=0
  else: 
    w = (1-r)**6*(6+36*r+82*r**2+72*r**3+30*r**4+5*r**5); 
    dwdr  = 11*r*(r+2)*(5*r**3+15*r**2+18*r+4)*(r-1)**5; 
    dwdrr = 22*(25*r**5+100*r**4+142*r**3+68*r**2-16*r-4)*(r-1)**4; 
    dwdrrr=0
  return w,dwdr,dwdrr,dwdrrr
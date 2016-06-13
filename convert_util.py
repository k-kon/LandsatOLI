#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Satellie Image Class for Geometric Correction
	Original data : Landsat 8 OLI : geotif format
	New data : converted image with specified region 	
 6/11/2016 copied from ETM+
"""
import numpy as np
import cv2
#from osgeo import gdal
#from osgeo import osr
import proj_util as pr

class original:
  def __init__(self,fname):
    self.name=fname
    f=open(fname+'MTL.txt','r')
    lines=f.readlines()
    f.close()
    #gain = np.zeros(11,dtype='float32')
    #offset = np.zeros(11,dtype='float32')
    gain=[]
    offset=[]
    for line in lines:
      if line.find('CORNER_UL_PROJECTION_X_PRODUCT') != -1:
        ULX= float(line[36:])
      if line.find('CORNER_UL_PROJECTION_Y_PRODUCT') != -1:
        ULY = float(line[36:])
      if line.find('CORNER_UR_PROJECTION_X_PRODUCT') != -1:
        URX= float(line[36:])
      if line.find('CORNER_UR_PROJECTION_Y_PRODUCT') != -1:
        URY= float(line[36:])
      if line.find('CORNER_LL_PROJECTION_X_PRODUCT') != -1:
        LLX= float(line[36:])
      if line.find('CORNER_LL_PROJECTION_Y_PRODUCT') != -1:
        LLY= float(line[36:])
      if line.find('CORNER_LR_PROJECTION_X_PRODUCT') != -1:
        LRX= float(line[36:])
      if line.find('CORNER_LR_PROJECTION_Y_PRODUCT') != -1:
        LRY= float(line[36:])
      if line.find('REFLECTIVE_LINES') != -1:
        jmax=int(line.split()[2])
      if line.find('REFLECTIVE_SAMPLES') != -1:
        imax=int(line.split()[2])
      if line.find('SUN_AZIMUTH') != -1:
        self.sun_az=float(line.split()[2])
      if line.find('SUN_ELEVATION') != -1:
        self.sun_el=float(line.split()[2])
      if line.find('RADIANCE_MULT_BAND_1 ') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_2') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_3') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_4') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_5') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_6') != -1:
        gain.append(float(line.split()[2]))
      if line.find('RADIANCE_MULT_BAND_7') != -1:
        gain.append(float(line.split()[2]))
      #if line.find('RADIANCE_MULT_BAND_8') != -1:
        #gain[7]=float(line.split()[2])
      if line.find('RADIANCE_MULT_BAND_9') != -1:
        gain.append(float(line.split()[2]))
      #if line.find('RADIANCE_MULT_BAND_10') != -1:
        #gain[9]=float(line.split()[2])
      #if line.find('RADIANCE_MULT_BAND_11') != -1:
        #gain[10]=float(line.split()[2]) 
      if line.find('RADIANCE_ADD_BAND_1 ') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_2') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_3') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_4') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_5') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_6') != -1:
        offset.append(float(line.split()[2]))
      if line.find('RADIANCE_ADD_BAND_7') != -1:
        offset.append(float(line.split()[2]))
      #if line.find('RADIANCE_ADD_BAND_8') != -1:
        #offset[7]=float(line.split()[2])
      if line.find('RADIANCE_ADD_BAND_9') != -1:
        offset.append(float(line.split()[2]))
      #if line.find('RADIANCE_ADD_BAND_10') != -1:
        #offset[9]=float(line.split()[2])
      #if line.find('RADIANCE_ADD_BAND_11') != -1:
        #offset[10]=float(line.split()[2])
    print ULX,ULY,URX,URY,LLX,LLY,LRX,LRY
    print imax,jmax
    print (URX-ULX)/(imax-1) # 30.0
    print (ULY-LLY)/(jmax-1) # 30.0
    self.gain=gain
    self.offset=offset
    fname=self.name+'B1.tif'
    gt,wkt,image=pr.read_tif(fname)
    self.image=image
    self.jmax,self.imax=image.shape
    self.xs=gt[0]
    self.ye=gt[3]
    self.dx=gt[1]
    self.dy=-gt[5]
    self.pc=self.imax/2.0
    self.lc=self.jmax/2.0
    self.xc=self.xs+self.dx*self.pc
    self.yc=self.ye-self.dy*self.lc
    self.t0=0.0
    self.angle=0.0
    self.hsat=691650.00
  def read_band(self,band):
    fname=self.name+'B'+str(band)+'.tif'
    gt,wkt,image=pr.read_tif(fname)
    self.image=image
  def display(self,name,xmax,ymax):
    old=self.image
    min=np.min(old) ; max=np.max(old)
    oldx=np.uint8(255.0*(old-min)/(max-min))
    oldy=cv2.resize(oldx,(xmax,ymax))
    cv2.imshow(name,oldy)
  def xy2pl(self,x,y):
    aa=np.cos(self.t0)/self.dx
    bb=-np.sin(self.t0)/self.dy
    p=aa*(x-self.xs)+bb*(y-self.ye)
    l=bb*(x-self.xs)-aa*(y-self.ye)
    return [p,l]
  def pl2xy(self,p,l):
    aa=np.cos(self.t0)*self.dx
    bb=-np.sin(self.t0)*self.dy
    x=aa*p+bb*l+self.xs
    y=bb*p-aa*l+self.ye
    return [x,y]
  def coverage(self):
    ul=self.pl2xy(0.0,0.0) 
    ur=self.pl2xy(float(self.imax),0.0)
    ll=self.pl2xy(0.0,float(self.jmax))
    lr=self.pl2xy(float(self.imax),float(self.jmax)) 
    left=int(ll[0]/10000.0)
    right=int(np.ceil(ur[0]/10000.0)) 
    lower=int(lr[1]/10000.0) 
    upper=int(np.ceil(ul[1]/10000.0))
    return [left,right,lower,upper]
  def chokka(self):
    mask=np.zeros((self.jmax,self.imax),dtype=np.uint8)
    temp=np.where(self.image > 0) 
    mask[temp]=1
    xtemp = temp[1]
    ytemp = temp[0]
    set1=np.min(ytemp)
    set2=np.max(xtemp)
    set3=np.max(ytemp)
    set4=np.min(xtemp)
    state=0
    position=[0]
    for i in range(self.imax):
      if mask[set1,i] != state:
        position.append(i)
        state=1-state
    state=0
    for j in range(self.jmax):
      if mask[j,set2] != state:
        position.append(j)
        state=1-state
    state=0
    for i in range(self.imax):
      if mask[set3,i] != state:
        position.append(i)
        state=1-state
    state=0
    for j in range(self.jmax):
      if mask[j,set4] != state:
        position.append(j)
        state=1-state
    tq1=float(position[1]-set4)/(position[7]-set1)
    tp1=float(position[3]-set1)/(set2-position[2])
    tq2=float(set2-position[6])/(set3-position[4])
    tp2=float(set3-position[8])/(position[5]-set4)
    tp=np.arctan((tp1+tp2)/2.0)
    tq=np.arctan((tq1+tq2)/2.0)
    self.pv=[np.cos(tp),-np.sin(tp)]
    self.qv=[-np.sin(tq),-np.cos(tq)]
    #pq=transpose([[pv],[qv]])
    #pinv=invert(pq)
    #angle=angle*!dtor
    #angle2=asin((earth+hsat)*sin(angle)/earth)
    #dn=earth*(angle2-angle)


class convert:
  def __init__(self,old,xmax,ymax):
    self.old=old
    self.xs=pr.xs
    self.ye=pr.ye
    self.dx=pr.dx
    self.dy=pr.dy
    self.xmax=xmax
    self.ymax=ymax
  def convert(self,tx,ty):
    x=self.xs+np.arange(self.xmax)*self.dx
    y=self.ye-np.arange(self.ymax)*self.dy
    xso=self.old.xs+tx
    yeo=self.old.ye+ty
    xc=(x-xso)/self.old.dx
    yc=(yeo-y)/self.old.dy
    xx,yy=np.meshgrid(xc,yc)
    xx=np.float32(xx)
    yy=np.float32(yy)
    return cv2.remap(self.old.image,xx,yy,cv2.INTER_LINEAR)
  def hyouka(self,parm):
    new=self.convert(parm[0],parm[1])
    res=np.corrcoef(inc[100:500,100:500].flat,new[100:500,100:500].flat)
    return -res[0,1]

def gwrite(sat,fnew):
  f=open('../'+fnew+'/gparm.txt','w')
  f.write(" * datum : WGS84\n")
  f.write(" image size:\n")  f.write("  pixel= "+str(sat.imax)+"\n")  f.write("  line = "+str(sat.jmax)+"\n")  f.write("image scene center:\n")  f.write("  p0= "+str(sat.pc)+"\n")  f.write("  l0= "+str(sat.lc)+"\n")  f.write("utm scene center:\n")  f.write("  x0= "+str(sat.xc)+"\n")  f.write("  y0= "+str(sat.yc)+"\n")  f.write("spatial resolution:\n")  f.write("  dx0= "+str(sat.dx)+"\n")  f.write("  dy0= "+str(sat.dy)+"\n")  f.write("orientation angle:\n")  f.write("  t0= "+str(sat.t0)+"\n")
  f.write("pointing angle:\n")
  f.write("  angle= "+str(sat.angle)+"\n")  f.write("sun position:\n")  f.write("  el= "+str(sat.sun_el)+"\n")  f.write("  az= "+str(sat.sun_az)+"\n")
  pvx=['{:10.6f}'.format(x) for x in sat.pv]
  qvx=['{:10.6f}'.format(x) for x in sat.qv]
  f.write("scanning direction:\n")
  f.write("  pv= "+" ".join(pvx)+"\n")
  f.write("orbit direction:\n")
  f.write("  qv= "+" ".join(qvx)+"\n")
  f.write("nadir utm point:\n")
  f.write("  xn=  "+str(sat.xc)+"\n")
  f.write("  yn=  "+str(sat.yc)+"\n")
  f.write("satellite height:"+"\n")
  f.write("  hsat= "+str(sat.hsat)+"\n")
  left,right,lower,upper=sat.coverage()
  f.write("coverage:"+"\n")
  f.write("  left = "+str(left)+"\n")
  f.write("  right= "+str(right)+"\n")
  f.write("  lower= "+str(lower)+"\n")
  f.write("  upper= "+str(upper)+"\n")
  offsetx=['{:7.2f}'.format(x) for x in sat.offset]
  gainx=['{:7.5f}'.format(x) for x in sat.gain]
  f.write("calibration:"+"\n")
  f.write("  offset= "+" ".join(offsetx)+"\n")
  f.write("  gain=   "+" ".join(gainx)+"\n")
  f.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Extraction of Landsat Image 
	Original data : Landsat 8 OLI : geotif format
		Specified Region from DEM Geotif Image region 	
"""
# cd /Volumes/Transcend/LandsatOLI
# displacement.py LC81080322015145LGN00 new_utm.tif

import sys
import os
import numpy as np
import cv2
import convert_util as ut
import proj_util as pr

param=sys.argv
if len(param)!=3:
    print 'Usage: displacement.py scene_name dem.tif'
    exit()

fscene=param[1]
fname=param[2]

#fscene='LC81080322015193LGN00'
#fname='new_utm.txt'

gt,wkt,dem=pr.read_tif(fname)
ymax,xmax=dem.shape
pr.xs=gt[0]
pr.dx=gt[1]
pr.xe=gt[0]+xmax*gt[1]
pr.ye=gt[3]
pr.dy=-gt[5]
pr.ys=gt[3]-ymax*gt[5]


os.chdir(fscene)

list=os.listdir('.')
for name in list:
  #print name
  if name.find('MTL.txt')!=-1 and name.find('._')==-1:
    fname=name[:-7]

flag=1
sat=ut.original(fname)
if flag:
  print sat.jmax,sat.imax
  print sat.xs,sat.dx
  print sat.ye,sat.dy
  print sat.offset
  print sat.gain

sat.chokka()
if flag:
  print sat.pv
  print sat.qv

#sat.read_band(40)
#sat.display('sat4',600,600)

#demx=cv2.resize(dem,(600,600))
#cv2.imshow('dem',demx/np.max(demx))


dem[dem<0.0]=0.0
inc=ut.incident(dem,sat)
#incx=cv2.resize(inc,(600,600))
#cv2.imshow('incx',incx/np.max(incx))
#cv2.destroyWindow('incx')

#new32=new.astype(np.float32)
#dd,dt=cv2.phaseCorrelate(inc,new32)
#new=conv.convert(-pr.dx*dd[0],pr.dy*dd[1])
#new32=new.astype(np.float32)
#cv2.phaseCorrelate(inc,new32)

print("*** Phase Only Correlation ***")
for band in [1,2,3,4,5,6,7,9]:
  sat.read_band(band)
  conv=ut.convert(sat,xmax,ymax)
  new=conv.convert(0.0,0.0)
  new32=new.astype(np.float32)
  print ' Band '+str(band)+':'
  #print cv2.phaseCorrelate(inc,new32)
  dd,dt=cv2.phaseCorrelate(inc,new32)
  print -pr.dx*dd[0],pr.dy*dd[1],dt

exit()

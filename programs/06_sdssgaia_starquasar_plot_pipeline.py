import fitsio
import sys
import os
import numpy
import scipy
from scipy import optimize
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db
gaialibdir=os.environ['GAIALIB']
sys.path.append(gaialibdir)
import gaia_db
githubdir=os.environ['GITHUB_DIR']

pngdir=os.environ['SDSS_PNGDIR']

def read_spec(csvfiles,fitsfilenames,plate,mjd):
  #plate=6783 ;  mjd=56284
  #plate=3586  ; mjd=55181
  #plate=3587  ; mjd=55182
  #plate=4221  ; mjd=55443
  #plate=11675 ; mjd=58523

  if(plate<10000): strplate=str("%04i"%(plate))
  if(plate>=10000): strplate=str("%05i"%(plate))
  platedir=pngdir+strplate
  if not (os.path.exists(platedir)): os.mkdir(platedir)
  [fig0,axs0]=plot_ratiosummary_step1()

  for k in range(2):
   # Loading Gaia XP spectra
   gaiaxp=gaia_db.GAIAXP(fitsfilenames[k])
   gaiaxp.readall()

   # Reading SDSS-GAIA csvfile
   df=pd.read_csv(csvfiles[k])
   # A plate example
   dftmp=df[(df['plate']==plate) & (df['mjd']==mjd)]
   dftmp.reset_index()
 
   nobject=len(dftmp)
   if(nobject==0): continue

   image_ratiocurve=numpy.zeros((100,nobject),dtype=numpy.float64)
   image_ratioline =numpy.zeros((100,nobject),dtype=numpy.float64)

   for i in range(len(dftmp)):
   #for i in range(4,5):
      mjd      =dftmp['mjd'].iloc[i]
      fiber    =dftmp['fiber'].iloc[i]
      print('#',plate,mjd,fiber)
      source_id=dftmp['source_id'].iloc[i]
      print('source_id=',source_id)
   
      # Finding GAIA data by source_id
      row_number=gaiaxp.df[gaiaxp.df['source_id']==source_id].index[0]
  
      # Reading GAIA data
      gaiaxp.read(source_id)
      # GAIA data in 20 points
      gaiaxp.photopoints20()

      #print(plate,mjd,fiber,source_id)

      ## Reading SDSS data
      spec=sdss_db.SDSSspec(plate,mjd,fiber)
      spec.read()
      spec.lambdaeff=dftmp['lambdaeff'].iloc[i]
      spec.zoffset  =dftmp['zoffset'].iloc[i]
      spec.xfocal   =dftmp['xfocal'].iloc[i]
      spec.yfocal   =dftmp['yfocal'].iloc[i]
      spec.sdssclass=dftmp['class'].iloc[i]
      spec.subclass =dftmp['subclass'].iloc[i]

# SDSS 20 pionts photometry
      spec.photopoints20()
# Calculate Ratios SDSS / GAIA
      [ratio_ptx,ratio_pty,ratio_err]=sdss_gaira_ratios(gaiaxp,spec)
# Plot
      #gaia_sdss_panels(gaiaxp,spec,ratio_ptx,ratio_pty,ratio_err)
      [curve_x,curve_y,line_x,line_y,a,b]=gaia_sdss_panels(gaiaxp,spec,ratio_ptx,ratio_pty,ratio_err,platedir)
      image_ratiocurve[:,i]=curve_y
      image_ratioline[:,i] =line_y
 
      print('slope=',a,'offset=',b)
      axs0[k,1].plot(curve_x,curve_y,'-',color='salmon',linewidth=0.2)
      axs0[k,1].plot(line_x,line_y,':',color='blue',linewidth=0.2)
      axs0[k,1].errorbar(ratio_ptx,ratio_pty,yerr=ratio_err,fmt='o',color='red',elinewidth=0.2,capsize=0.2,markersize=0.5)

      xmin=3300. ; xmax=10500. ; dx=xmax-xmin
      ymin=0.4   ; ymax=1.6    ; dy=ymax-ymin
      if(k==0):
         square1=axs0[0,0].scatter(spec.xfocal,spec.yfocal,c=a,cmap='rainbow',marker='o',vmin=-1.2,vmax=1.2,s=20.0,alpha=0.9,edgecolors='k')
         square2=axs0[1,0].scatter(spec.xfocal,spec.yfocal,c=b,cmap='rainbow',marker='o',vmin=0.8,vmax=1.2,s=20.0,alpha=0.9,edgecolors='k')
      elif(k==1):
         square1=axs0[0,0].scatter(spec.xfocal,spec.yfocal,c=a,cmap='rainbow',marker='D',vmin=-1.2,vmax=1.2,s=12.0,alpha=0.9,edgecolors='k')
         square2=axs0[1,0].scatter(spec.xfocal,spec.yfocal,c=b,cmap='rainbow',marker='D',vmin=0.8,vmax=1.2,s=12.0,alpha=0.9,edgecolors='k')

   #ax=asx0[0,0]
   median_curve=numpy.median(image_ratiocurve,axis=1)
   median_line =numpy.median(image_ratioline,axis=1)
   axs0[k,1].plot(curve_x,median_curve,'-',color='red',linewidth=3.0)
   axs0[k,1].plot(line_x,median_line,'-',color='blue',linewidth=3.0)

   if(k==0): axs0[0,1].text(xmin+0.7*dx,ymin+dy*0.9,"%3i"%(nobject)+' Stars')
   if(k==1): axs0[1,1].text(xmin+0.7*dx,ymin+dy*0.9,"%3i"%(nobject)+' Quasars')

  platename=spec.strplate+'+'+spec.strmjd
  fig0.suptitle(platename)
  fig0.colorbar(square1,ax=axs0[0,0])
  fig0.colorbar(square2,ax=axs0[1,0])
  fig0.savefig(platedir+'/plate_'+platename+'.png')
  fig0.clf()

def sdss_gaira_ratios(gaiaxp,spec):
   flag_gaia=numpy.where(gaiaxp.ptmask==1,1,0)
   flag_sdss=numpy.where(gaiaxp.ptmask==1,1,0)
   #print('flag_gaia',flag_gaia)
   #print('flag_sdss',flag_sdss)
   flag=flag_gaia*flag_sdss
   #print('total flag',flag)
   flagsum=numpy.sum(flag)
   if(flagsum==0):
      ratio_ptx=numpy.zeros(1)
      ratio_ptx=numpy.ones(1)
      ratio_err=numpy.ones(1)
      ratio_ptx[0]=gaiaxp.ptx[5]
      return [ratio_ptx,ratio_pty,ratio_err]
      
   ratio_ptx=numpy.compress(flag,gaiaxp.ptx)

   gaia_pty=numpy.compress(flag,gaiaxp.pty)
   gaia_err=numpy.compress(flag,gaiaxp.ptyerr)
   gaia_snr=gaia_pty/gaia_err
   sdss_pty=numpy.compress(flag,spec.pty)
   sdss_err=numpy.compress(flag,spec.ptyerr)

   sdss_snr=sdss_pty/sdss_err
   ratio_pty=sdss_pty/gaia_pty

   ratio_snr=numpy.zeros(len(sdss_snr))
   ratio_err=numpy.zeros(len(sdss_snr))
   for i in range(len(sdss_snr)):
      if(sdss_snr[i]<gaia_snr[i]):
         ratio_snr[i]=sdss_snr[i]
      else:
         ratio_snr[i]=gaia_snr[i]
      ratio_err[i]=ratio_pty[i]/ratio_snr[i]
   return [ratio_ptx,ratio_pty,ratio_err]
 
   #print('selected xpt=',gaia_ptx)
   #ratio_snr=min(gaia_snr,sdss_snr)
   #print('snr=',ratio_snr)

def fitting_curve(ratio_ptx,ratio_pty,ratio_err):

   curve_x=3000.0+numpy.arange(100)*75 ; curve_y=numpy.zeros(100)
   if(len(ratio_ptx)<=4):
      curve_y=numpy.ones(100)
      return [curve_x,curve_y]

   fitfunc = lambda p, x: p[0]*((x/10000.0)**p[1])*numpy.exp(p[2]*(x/10000.0))+p[3]
   #errfunc = lambda p, x, y: (fitfunc(p,x) - y)**2
   errfunc = lambda p, x, y, z: (fitfunc(p,x) - y)**2/z**2
   p0 = [1.0,0.1,-1.0,0.1]
   p0 = [1.0,0.0,0.0,0.0]
   #p0 = [1.0,0.1,-1.0]
   #p1,success=scipy.optimize.leastsq(errfunc,p0[:],args=(rwave,rflux),maxfev=10000)
   #p1,success=scipy.optimize.leastsq(errfunc,p0[:],args=(xtmp,ytmp),maxfev=10000)
   #p1,success=scipy.optimize.leastsq(errfunc,p0[:],args=(ratio_ptx,ratio_pty),maxfev=10000)
   p1,success=scipy.optimize.leastsq(errfunc,p0[:],args=(ratio_ptx,ratio_pty,ratio_err),maxfev=10000)

   curve_y=p1[0]*((curve_x/10000.0)**p1[1])*numpy.exp(p1[2]*(curve_x/10000.0))+p1[3]
   
   print(p1)
   return [curve_x,curve_y]

def fitting_line(ratio_ptx,ratio_pty,ratio_err):

   line_x=3000.0+numpy.arange(100)*75 ; line_y=numpy.ones(100)
   if(len(ratio_ptx)<=4):
      a=0.0 ; b=1.0
      return [line_x,line_y,a,b]

   fitfunc = lambda p, x: p[0]*(x/10000.0-0.7)+p[1]
   errfunc = lambda p, x, y, z: (fitfunc(p,x) - y)**2/z**2
   #p0 = [1.0,0.1,-1.0,0.1]
   p0 = [0.0,1.0]
   p1,success=scipy.optimize.leastsq(errfunc,p0[:],args=(ratio_ptx,ratio_pty,ratio_err),maxfev=10000)

   #line_y=p1[0]*(line_x/10000.0-0.7)+p1[1]
   line_y=p1[0]*(line_x/10000.0-0.7)+p1[1]
   a=p1[0] ; b=p1[1]
   #curve_y=p1[0]*((curve_x/10000.0)**p1[1])*numpy.exp(p1[2]*(curve_x/10000.0))+p1[3]
  
   #print(p1)
   return [line_x,line_y,a,b]

def gaia_sdss_panels(gaiaxp,spec,ratio_ptx,ratio_pty,ratio_err,platedir):

# Fit Curve
   [curve_x,curve_y]=fitting_curve(ratio_ptx,ratio_pty,ratio_err)
   [line_x,line_y,a,b]  =fitting_line(ratio_ptx,ratio_pty,ratio_err)

   specname=spec.strplate+'+'+spec.strmjd+'+'+spec.strfiber

   # Limits
   xmin=3300.;  xmax=10500.
   dx=xmax-xmin
   ymin=0.1  ; ymax=max(spec.flux)*1.15 
   dy=ymax-ymin

   matplotlib.use('Agg')
   # Figure Resolution
   plt.rcParams["figure.dpi"] = 144
   # Font to be Times
   plt.rc('font',family='serif')

   fig,axs=plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios': [3, 1]},tight_layout=True)

   # Bottom Panel
   zero_ptx=numpy.array([xmin,xmax]) ; zero_pty=numpy.array([1.0,1.0])
   axs[1].set_xlabel('Wavelength',fontsize=12)
   axs[1].set_ylabel('Ratio',fontsize=12)
   axs[1].set_xlim([xmin,xmax])
   #axs[1].set_ylim([0.6,1.4])
   axs[1].set_ylim([0.4,1.6])
   axs[1].plot(zero_ptx,zero_pty,'--')
   axs[1].plot(curve_x,curve_y,'-',color='salmon')
   axs[1].plot(line_x,line_y,':',color='blue')
   axs[1].errorbar(ratio_ptx,ratio_pty,yerr=ratio_err,fmt='o',color='red',elinewidth=2,capsize=3)

   axs[0].set_ylabel('Flux (10$^{-17}$ erg/Ang/s/cm$^{2}$)',fontsize=12)
   axs[0].set_xlim([xmin,xmax])
   axs[0].set_ylim([ymin,ymax])
   axs[0].tick_params(axis='both',direction='in')
   # SDSS vs GAIA spec
   axs[0].plot(spec.wave,spec.flux,linewidth=1.0,color='b',alpha=0.5)
   axs[0].plot(gaiaxp.wave,gaiaxp.flux,linewidth=2.0,color='r',alpha=0.7)
   # SDSS vs GAIA photometry points
   axs[0].errorbar(gaiaxp.ptx,gaiaxp.pty,yerr=gaiaxp.ptyerr,fmt='o',color='magenta',elinewidth=2,capsize=3)
   axs[0].errorbar(spec.ptx,spec.pty,yerr=spec.ptyerr,fmt='o',color='cyan',elinewidth=2,capsize=3)
   # Spectrum Name
   axs[0].text(xmin+dx*0.48,ymin+0.90*dy,specname,fontsize=18)
   axs[0].text(xmin+dx*0.48,ymin+0.80*dy,'$\lambda_{eff}=$'+"%4.0f"%(spec.lambdaeff),fontsize=12)
   axs[0].text(xmin+dx*0.72,ymin+0.80*dy,'$z_{offset}$='+"%3.0f"%(spec.zoffset),fontsize=12)
   axs[0].text(xmin+dx*0.48,ymin+0.70*dy,"%-10s"%(spec.sdssclass),fontsize=12)
   axs[0].text(xmin+dx*0.72,ymin+0.70*dy,"%-15s"%(spec.subclass),fontsize=12)
   axs[0].text(xmin+dx*0.48,ymin+0.60*dy,'$x_{focal}$='+"%4.1f"%(spec.xfocal),fontsize=12)
   axs[0].text(xmin+dx*0.72,ymin+0.60*dy,'$y_{focal}$='+"%4.1f"%(spec.yfocal),fontsize=12)

   plt.savefig(platedir+'/spec_'+specname+'.png')
   plt.clf()
   plt.close()

   return [curve_x,curve_y,line_x,line_y,a,b]

def gaia_sdss(gaiaxp,spec):
   specname=spec.strplate+'+'+spec.strmjd+'+'+spec.strfiber

   matplotlib.use('Agg')
   # Figure Resolution
   plt.rcParams["figure.dpi"] = 144
   # Font to be Times
   plt.rc('font',family='serif')
   plt.xlabel('Wavelength',fontsize=15)
   plt.ylabel('Flux (10$^{-17}$ erg/Ang/s/cm$^{2}$)',fontsize=15)
   plt.tick_params(axis='both',direction='in')

   xmin=3300.;  xmax=10500.;  ymin=0.1
   plt.xlim([xmin,xmax])
   dx=xmax-xmin
   ymax=max(spec.flux)
   plt.ylim([ymin,ymax])
   dy=ymax-ymin

   plt.plot(spec.wave,spec.flux,linewidth=1.0,color='b',alpha=0.5)
   plt.plot(gaiaxp.wave,gaiaxp.flux,linewidth=2.0,color='r',alpha=0.7)
   plt.errorbar(gaiaxp.ptx,gaiaxp.pty,yerr=gaiaxp.ptyerr,fmt='o',color='magenta',elinewidth=2,capsize=3)
   plt.errorbar(spec.ptx,spec.pty,yerr=spec.ptyerr,fmt='o',color='cyan',elinewidth=2,capsize=3)
   plt.text(xmin+dx*0.48,ymin+0.90*dy,specname,fontsize=18)
   plt.savefig(specname+'.png')
   plt.clf()
   plt.close()


def read_2dspec(fitsfilename):
   gaiaxp=gaia_db.GAIAXP(fitsfilename)
   gaiaxp.readall()
   print(gaiaxp.source_id_list)
   print(gaiaxp.thing_id_list)
   print(gaiaxp.ra_list)
   print(gaiaxp.dec_list)
   print(gaiaxp.df)
   source_id=1544408185954391552
   gaiaxp.read(source_id)
   for i in range(len(gaiaxp.flux)):
      print(gaiaxp.wave[i],gaiaxp.flux[i],gaiaxp.fluxerr[i])


def plot_ratiosummary_step1():
   matplotlib.use('Agg')
   plt.rc('font',family='serif')
   fig0,axs0=plt.subplots(nrows=2,ncols=2,sharex=False,\
           gridspec_kw={'width_ratios': [1.25, 2]},constrained_layout=True)


   xmin=-350.;  xmax=350. ; ymin=-350. ; ymax=350.
   dx=xmax-xmin ; dy=ymax-ymin
   zero_ptx1=numpy.array([xmin,xmax]) ; zero_pty1=numpy.array([0.0,0.0])
   zero_ptx2=numpy.array([0.0,0.0])   ; zero_pty2=numpy.array([ymin,ymax])
   axs0[0,0].text(xmin+0.1*dx,ymin+0.9*dy,'Slope')
   axs0[0,0].set_ylabel('Y Focal',fontsize=12)
   axs0[0,0].set_xlim([-350.,350.])
   axs0[0,0].set_ylim([-350.,350.])
   axs0[0,0].plot(zero_ptx1,zero_pty1,'--',color='k')
   axs0[0,0].plot(zero_ptx2,zero_pty2,'--',color='k')
   circle=Circle((350.,350.), 0.5, color='b', fill=False)
   axs0[0,0].add_patch(circle)

   axs0[1,0].text(xmin+0.1*dx,ymin+0.9*dy,'Offset')
   axs0[1,0].set_xlabel('X Focal',fontsize=12)
   axs0[1,0].set_ylabel('Y Focal',fontsize=12)
   axs0[1,0].set_xlim([-350.,350.])
   axs0[1,0].set_ylim([-350.,350.])
   axs0[1,0].plot(zero_ptx1,zero_pty1,'--',color='k')
   axs0[1,0].plot(zero_ptx2,zero_pty2,'--',color='k')
   circle=Circle((350.,350.), 0.5, color='b', fill=False)
   axs0[1,0].add_patch(circle)

   xmin=3300.;  xmax=10500.
   zero_ptx=numpy.array([xmin,xmax]) ; zero_pty=numpy.array([1.0,1.0])
   axs0[0,1].set_xlim([xmin,xmax])
   axs0[0,1].set_ylim([0.4,1.6])
   axs0[0,1].plot(zero_ptx,zero_pty,'--',color='k')

   axs0[1,1].set_xlim([xmin,xmax])
   axs0[1,1].set_ylim([0.4,1.6])
   axs0[1,1].plot(zero_ptx,zero_pty,'--',color='k')
   
   return [fig0,axs0]

fitsfilenames=[githubdir+'/projects_gaia/data/gaiadr3_xpspec_sdssdr17_star.fits',\
              githubdir+'/projects_gaia/data/gaiadr3_xpspec_sdssdr17_quasar.fits']
csvfiles=[githubdir+'/projects_gaia/csvfiles/gaiadr3_sdssdr17_star_combined.csv',\
         githubdir+'/projects_gaia/csvfiles/gaiadr3_sdssdr17_quasar_combined.csv']

df1=pd.read_csv(csvfiles[0])
df2=pd.read_csv(csvfiles[1])
df3=df1.append(df2,ignore_index=True)
df4=df3.reset_index()
#print(df4)
#df4.to_csv('dfall.csv',index=False)

df5=df4.groupby(['plate','mjd']).size().reset_index()
df5.to_csv('df5.csv',index=False)
#sys.exit(1)

for i in range(28,len(df5)):
  plate=df5['plate'].iloc[i]
  mjd  =df5['mjd'].iloc[i]
  print(i,plate,mjd)
  read_spec(csvfiles,fitsfilenames,plate,mjd)
#spPlate-11675-58523.fits


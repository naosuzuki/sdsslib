import fitsio
import sys
import os
import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Ellipse
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db
gaialibdir=os.environ['GAIALIB']
sys.path.append(gaialibdir)
import gaia_db

def read_spec(csvfile,fitsfilename):
   # Loading Gaia XP spectra
   gaiaxp=gaia_db.GAIAXP(fitsfilename)
   gaiaxp.readall()

   # Reading SDSS-GAIA csvfile
   df=pd.read_csv(csvfile)
   # A plate example
   plate=11675 ; mjd=58523
   plate=6783 ;  mjd=56284
   dftmp=df[(df['plate']==plate) & (df['mjd']==mjd)]
   dftmp.reset_index()

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

   print(len(dftmp))
   for i in range(len(dftmp)):
   #for i in range(4,5):
      mjd=dftmp['mjd'].iloc[i]
      fiber=dftmp['fiber'].iloc[i]
      source_id=dftmp['source_id'].iloc[i]
      print('source_id=',source_id)
   
      # Reading GAIA data
      gaiaxp.read(source_id)
      # GAIA data in 20 points
      gaiaxp.photopoints20()
      #for j in range(20):
        #print(gaiaxp.wave[j],gaiaxp.flux[j],gaiaxp.wid[j])
      #  print(gaiaxp.ptx[j],gaiaxp.pty[j],gaiaxp.ptwid[j])
      #sys.exit(1)

      print(plate,mjd,fiber,source_id)
      # Reading SDSS data
      spec=sdss_db.SDSSspec(plate,mjd,fiber)
      spec.read()
      spec.photopoints20()
      #sys.exit(1)

      specname=spec.strplate+'+'+spec.strmjd+'+'+spec.strfiber
# Plot
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

fitsfilename='../../projects_gaia/data/gaiadr3_xpspec_sdssdr17_quasar.fits'
fitsfilename='../../projects_gaia/data/gaiadr3_xpspec_sdssdr17_star.fits'

#csvfile='../../projects_gaia/data/gaiadr3_sdssdr17_star.csv'
#csvfile='../../projects_gaia/data/gaiadr3_sdssdr17_quasar.csv'
#csvfile='../../projects_gaia/csvfiles/gaiadr3id_sdssdr17_star.csv'
csvfile='../../projects_gaia/csvfiles/gaiadr3_sdssdr17_quasar_combined.csv'
csvfile='../../projects_gaia/csvfiles/gaiadr3_sdssdr17_star_combined.csv'
read_spec(csvfile,fitsfilename)
#spPlate-11675-58523.fits


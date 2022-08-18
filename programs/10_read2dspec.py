import fitsio
import sys
import os
import numpy
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db
import matplotlib.pyplot as plt
from matplotlib import rc

#githubdir=os.environ['GITHUB_DIR']
#gaiacsvdir=githubdir+'projects_gaia/csvfiles/'
#sdsscsvdir=githubdir+'sdsslib/csvfiles/'

def read2dspec(fitsfilename):

   spec2d=sdss_db.SDSSspec2d(fitsfilename)
   print(spec2d.nspec)
   print(spec2d.wave)
   spec2d.read()
   spec2d.read_gaia()

   #df1=pd.DataFrame(list(zip(spec2d.ra_list,spec2d.dec_list,\
   #   spec2d.thing_id_list,spec2d.snr_list,\
   #   spec2d.plate_list,spec2d.mjd_list,spec2d.fiber_list,\
   #   spec2d.parallax_list,spec2d.parallaxsnr_list,\
   #   spec2d.gaia_gmag_list,spec2d.gaia_bprp_list,\
   #   spec2d.gaia_mg_list,spec2d.gaia_pc_list)),\
   #   columns=['ra','dec','thing_id','snall',\
   #   'plate','mjd','fiber','parallax','parallaxsnr',\
   #   'gaia_gmag','gaia_bprp','gaia_mg','gaia_pc'])

   df=pd.DataFrame(list(zip(spec2d.ra_list,spec2d.dec_list,\
      spec2d.thing_id_list,spec2d.snr_list,\
      spec2d.parallax_list,spec2d.parallaxsnr_list,\
      spec2d.gaia_gmag_list,spec2d.gaia_bprp_list,\
      spec2d.gaia_mg_list,spec2d.gaia_pc_list)),\
      columns=['ra','dec','thing_id','snall',\
      'parallax','parallaxsnr',\
      'gaia_gmag','gaia_bprp','gaia_mg','gaia_pc'])

   print('total=',len(df))
   #dfgaia=df[(df['gaia_mg'].notna()) & (df['gaia_bprp'].notna())]
   #dfgaia=df[(df['gaia_mg'].notna()) & (df['gaia_bprp'].notna()) & (df['parallaxsnr']>3.)]
   dfgaia=df[(df['gaia_mg'].notna()) & (df['gaia_bprp'].notna()) & (df['parallaxsnr']>5.)]
   #dfgaia=df[(df['gaia_mg'].notna()) & (df['gaia_bprp'].notna()) & (df['parallaxsnr']>=10.)]
   #dfgaia=df[(df['gaia_mg'].notna()) & (df['gaia_bprp'].notna()) & (df['parallaxsnr']>=20.)]
   ptx=dfgaia['gaia_bprp'].to_numpy()
   pty=dfgaia['gaia_mg'].to_numpy()
   parallaxsnr=dfgaia['parallaxsnr'].to_numpy()

   plot_wd(ptx,pty,parallaxsnr)
   print(dfgaia)
   print(len(dfgaia))
   #print('gaia gmag=',spec2d.gaia_gmag_list)
   #print('gaia bprp=',spec2d.gaia_bprp_list)
   #print('gaia mg=',spec2d.gaia_mg_list)
   #flag=numpy.where(spec2d.gaia_mg_list<0,1,0)
   #print('negative mag=',numpy.sum(flag))

def plot_wd(ptx,pty,parallaxsnr):
   # WD range
   ymin=16.5 ; ymax=8.8
   xmin=-1.0 ; xmax=1.6
   dx=xmax-xmin ; dy=ymax-ymin

   # HR all range
   ymin=16.5 ; ymax=-5.0
   xmin=-1.0 ; xmax=4.0
   dx=xmax-xmin ; dy=ymax-ymin

   # Latex
   # Font to be Times
   plt.rcParams["figure.figsize"] = (6, 8)
   plt.rc('font',family='serif')
   plt.xlabel(r'G$_{BP}$-G$_{RP}$')
   plt.ylabel(r'M$_{G}$')
   plt.tick_params(axis='both',direction='in')
   plt.xlim([xmin,xmax])
   plt.ylim([ymin,ymax])
   #plt.title('SDSS DR8 + GAIA DR3 Parallax S/N : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR8 + GAIA DR3 Parallax S/N > 3 : '+"%6i"%(len(ptx))+" Stars")
   plt.title('SDSS DR8 + GAIA DR3 Parallax S/N > 5 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR8 + GAIA DR3 Parallax S/N > 10 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR8 + GAIA DR3 Parallax S/N > 20 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR17 + GAIA DR3 Parallax S/N : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR17 + GAIA DR3 Parallax S/N > 3 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR17 + GAIA DR3 Parallax S/N > 5 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR17 + GAIA DR3 Parallax S/N > 10 : '+"%6i"%(len(ptx))+" Stars")
   #plt.title('SDSS DR17 + GAIA DR3 Parallax S/N > 20 : '+"%6i"%(len(ptx))+" Stars")
   print(ptx,pty)
   #plt.scatter(ptx,pty,marker='.',c='k',s=3.0,edgecolors='none')
   plt.scatter(ptx,pty,c=parallaxsnr,s=0.5,cmap='rainbow',alpha=1.00,vmin=0.1,vmax=20.0)
   #plt.colorbar(orientation='horizontal',location='top')
   plt.colorbar(orientation='vertical',location='right')
   #plt.savefig('HR_SDSSDR8GAIADR3.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR8GAIADR3SNR3.png',orientation='portrait')
   plt.savefig('HR_SDSSDR8GAIADR3SNR5.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR8GAIADR3SNR10.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR8GAIADR3SNR20.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR17GAIADR3.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR17GAIADR3SNR3.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR17GAIADR3SNR5.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR17GAIADR3SNR10.png',orientation='portrait')
   #plt.savefig('HR_SDSSDR17GAIADR3SNR20.png',orientation='portrait')
   plt.clf()
   plt.close()

fitsfilename='sdssDR17snr10a_wd.fits'
fitsfilename='sdssDR8snr10a_wd.fits'
fitsfilename='sdssDR8_wd.fits'
fitsfilename='sdssDR17_wd.fits'
fitsfilename='../data/sdssDR17_star.fits'
fitsfilename='../data/sdssDR8snr10a_star.fits'
read2dspec(fitsfilename)

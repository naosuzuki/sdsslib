import fitsio
import sys
import os
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db

githubdir=os.environ['GITHUB_DIR']
gaiacsvdir=githubdir+'projects_gaia/csvfiles/'
sdsscsvdir=githubdir+'sdsslib/csvfiles/'
fitsdatadir=githubdir+'sdsslib/fitsdata/'

objtype=sys.argv[1]
dr=sys.argv[2]
snr=int(sys.argv[3])

# Star
#objtype='star'
if(objtype=='star'):
   flag_gaia=True
   flag_restframe=False
   if(dr=='dr8'):
      csvfile=gaiacsvdir+'gaiadr3_sdssdr8_star.csv'
      fitsfilename1='sdssDR8_wd.fits'
      fitsfilename2='sdssDR8_star.fits'
      fitsfilename1=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_wd.fits'
      fitsfilename2=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_star.fits'
   elif(dr=='dr17'):
      csvfile=gaiacsvdir+'gaiadr3_sdssdr17_star.csv'
      fitsfilename1=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'+_wd.fits'
      fitsfilename2=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'_star.fits'

   df=pd.read_csv(csvfile)
   print('df all=',len(df))
   dftmp10=df[df['snall']>=float(snr)]
   print('df sn>10=',len(dftmp10))
   dftmp5=df[df['snall']>5.0]
   print('df sn>5',len(dftmp5))
   dftmp3=df[df['snall']>3.0]
   print('df sn>3',len(dftmp3))
   dftmp1=df[df['snall']>1.0]
   print('df sn>1',len(dftmp1))
   del dftmp5 ; del dftmp3 ; del dftmp1

   dfstar=dftmp10.copy()
# Selecting White Dwarfs
   dfwd=dfstar[dfstar['subclass'].str.contains('WD')]
   dfspec=dfwd.sort_values(by=['teff'],ascending=True)
   dfspec.reset_index()
   sdss_db.create_2dspec(dfspec,fitsfilename1,objtype,flag_gaia,dr)
   del dfspec ; del df ; del dfwd ; dftmp10

# All of Stars
   dfspec2=dfstar.sort_values(by=['teff'],ascending=True)
   dfspec2.reset_index()
   print('dfspec2=',dfspec2)
   sdss_db.create_2dspec(dfspec2,fitsfilename2,objtype,flag_gaia,dr)
   del dfspec2

# Galaxy
#objtype='galaxy'
if(objtype=='galaxy'):
   flag_gaia=False
   flag_restframe=True
   if(dr=='dr8'):
      csvfile=sdsscsvdir+'dr8_spall_galaxy.csv'
      fitsfilename='sdssDR8snr3_galaxy.fits'
      fitsfilename=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_galaxy.fits'
      fitsfilename=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_galaxy_rest.fits'
   elif(dr=='dr17'):
      csvfile=sdsscsvdir+'v5_13_2_spall_galaxy.csv'
      fitsfilename='sdssDR17snr3_galaxy.fits'
      fitsfilename=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'_galaxy.fits'
      fitsfilename=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'_galaxy_rest.fits'
   df=pd.read_csv(csvfile)
   print(df)

   print('dfgalaxy all=',len(df))
   dftmp10=df[df['snall']>=float(snr)]
   print('dfgalaxy sn>10',len(dftmp10))
   dftmp5=df[df['snall']>5.0]
   print('dfgalaxy sn>5',len(dftmp5))
   dftmp3=df[df['snall']>3.0]
   print('dfgalaxy sn>3',len(dftmp3))
   dftmp1=df[df['snall']>1.0]
   print('dfgalaxy sn>1',len(dftmp1))

# S/N > 3 is chosen for DR17
   dfgalaxy=dftmp10.copy()
   del dftmp5 ; del dftmp3 ; del dftmp1

   dfspec=dfgalaxy.sort_values(by=['z'],ascending=True)

   flag_gaia=False
   sdss_db.create_2dspec(dfspec,fitsfilename,objtype,flag_gaia,flag_restframe,dr)
   del dfgalaxy ; del dfspec ; del df

# Quasar
#objtype='quasar'
if(objtype=='quasar'):
   flag_gaia=False
   flag_restframe=True
   if(dr=='dr8'):
      csvfile=sdsscsvdir+'dr8_spall_quasar.csv'
      fitsfilename='sdssDR8_quasar.fits'
      fitsfilename=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_quasar.fits'
      fitsfilename=fitsdatadir+'sdssDR8snr'+"%2i"%(snr)+'_quasar_rest.fits'
   elif(dr=='dr17'):
      csvfile=sdsscsvdir+'v5_13_2_spall_quasar.csv'
      fitsfilename='sdssDR17_quasar.fits'
      fitsfilename=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'_quasar.fits'
      fitsfilename=fitsdatadir+'sdssDR17snr'+"%2i"%(snr)+'_quasar_rest.fits'
   df=pd.read_csv(csvfile)
   print(df)
   print('df all=',len(df))
   dftmp10=df[df['snall']>=float(snr)]
   print('df sn>10',len(dftmp10))
   dftmp5=df[df['snall']>5.0]
   print('df sn>5',len(dftmp5))
   dftmp3=df[df['snall']>3.0]
   print('df sn>3',len(dftmp3))
   dftmp1=df[df['snall']>1.0]
   print('df sn>1',len(dftmp1))

   # S/N>10 for DR17
   dfquasar=dftmp10.copy()
   #dfquasar=df.copy()
   dfspec=dfquasar.sort_values(by=['z'],ascending=True)

   sdss_db.create_2dspec(dfspec,fitsfilename,objtype,flag_gaia,flag_restframe,dr)
   del df ; del dfquasar ;  del dfspec

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

objtype=sys.argv[1]
dr=sys.argv[2]

#spall=sdss_catalog.spall()
#print(spall.dr)
#print(spall.fitstablename)
#spall.read()

# Star
#objtype='star'
if(objtype=='star'):
   flag_gaia=True
   if(dr=='dr8'):
      csvfile=gaiacsvdir+'gaiadr3_sdssdr8_star.csv'
      fitsfilename1='sdssDR8_wd.fits'
      fitsfilename2='sdssDR8_star.fits'
      fitsfilename1='sdssDR8snr10b_wd.fits'
      fitsfilename2='sdssDR8snr10b_star.fits'
   elif(dr=='dr17'):
      csvfile=gaiacsvdir+'gaiadr3_sdssdr17_star.csv'
      fitsfilename1='sdssDR17snr10b_wd.fits'
      fitsfilename2='sdssDR17snr10b_star.fits'

   df=pd.read_csv(csvfile)
   print('df all=',len(df))
   dftmp10=df[df['snall']>=10.0]
   print('df sn>10=',len(dftmp10))
   dftmp5=df[df['snall']>5.0]
   print('df sn>5',len(dftmp5))
   dftmp3=df[df['snall']>3.0]
   print('df sn>3',len(dftmp3))
   dftmp1=df[df['snall']>1.0]
   print('df sn>1',len(dftmp1))
   del dftmp5 ; del dftmp3 ; del dftmp1

   #dfstar=df.copy()
   dfstar=dftmp10.copy()
# Selecting White Dwarfs
   dfwd=dfstar[dfstar['subclass'].str.contains('WD')]
   dfspec=dfwd.sort_values(by=['teff'],ascending=False)
   dfspec.reset_index()
   sdss_db.create_2dspec(dfspec,fitsfilename1,objtype,flag_gaia)
   del dfspec ; del df ; del dfwd ; dftmp10

# All of Stars
   dfspec2=dfstar.sort_values(by=['teff'],ascending=False)
   dfspec2.reset_index()
   print('dfspec2=',dfspec2)
   sdss_db.create_2dspec(dfspec2,fitsfilename2,objtype,flag_gaia)
   del dfspec2

   #for i in range(1000):
#   for i in range(len(dfspec)):
#    plate=dfspec['plate'].iloc[i]
#    mjd  =dfspec['mjd'].iloc[i]
#    fiber=dfspec['fiber'].iloc[i]
    #print(i,plate,mjd,fiber)
#    print(i,plate,mjd,fiber,dfspec['subclass'].iloc[i],dfspec['teff'].iloc[i])

   #for i in range(100,200):
   #for i in range(len(dfspec2)):
   # plate=dfspec2['plate'].iloc[i]
   # mjd  =dfspec2['mjd'].iloc[i]
   # fiber=dfspec2['fiber'].iloc[i]
   # print(i,plate,mjd,fiber)
   # print(i,plate,mjd,fiber,dfspec2['subclass'].iloc[i],dfspec2['teff'].iloc[i])


# Galaxy
#objtype='galaxy'
if(objtype=='galaxy'):
   flag_gaia=False
   if(dr=='dr8'):
      csvfile=sdsscsvdir+'dr8_spall_galaxy.csv'
      fitsfilename='sdssDR8snr3_galaxy.fits'
      fitsfilename='sdssDR8snr10b_galaxy.fits'
   elif(dr=='dr17'):
      csvfile=sdsscsvdir+'v5_13_2_spall_galaxy.csv'
      fitsfilename='sdssDR17snr3_galaxy.fits'
      fitsfilename='sdssDR17snr10b_galaxy.fits'
   df=pd.read_csv(csvfile)
   print(df)
   #df['class']=df['class'].str.strip()
   #df['subclass']=df['subclass'].str.strip()
   #dfgalaxy=df[(df['class']=='GALAXY') & (df['thing_id']!=-1)]
   print('dfgalaxy all=',len(df))
   dftmp10=df[df['snall']>=10.0]
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

   dfspec=dfgalaxy.sort_values(by=['z'],ascending=False)

   flag_gaia=False
   sdss_db.create_2dspec(dfspec,fitsfilename,objtype,flag_gaia)
   del dfgalaxy ; del dfspec ; del df

# Quasar
#objtype='quasar'
if(objtype=='quasar'):
   flag_gaia=False
   if(dr=='dr8'):
      csvfile=sdsscsvdir+'dr8_spall_quasar.csv'
      fitsfilename='sdssDR8_quasar.fits'
      fitsfilename='sdssDR8snr10b_quasar.fits'
   elif(dr=='dr17'):
      csvfile=sdsscsvdir+'v5_13_2_spall_quasar.csv'
      fitsfilename='sdssDR17_quasar.fits'
      fitsfilename='sdssDR17snr10b_quasar.fits'
   df=pd.read_csv(csvfile)
   print(df)
   print('df all=',len(df))
   dftmp10=df[df['snall']>=10.0]
   print('df sn>10',len(dftmp10))
   dftmp5=df[df['snall']>5.0]
   print('df sn>5',len(dftmp5))
   dftmp3=df[df['snall']>3.0]
   print('df sn>3',len(dftmp3))
   dftmp1=df[df['snall']>1.0]
   print('df sn>1',len(dftmp1))
   #del dftmp5 ; del dftmp1
   #del dftmp3  

   # S/N>10 for DR17
   dfquasar=dftmp10.copy()
   #dfquasar=df.copy()
   dfspec=dfquasar.sort_values(by=['z'],ascending=False)

   sdss_db.create_2dspec(dfspec,fitsfilename,objtype,flag_gaia)
   #del products_list ; del df ; del dfquasar
   del df ; del dfquasar ;  del dfspec

import fitsio
import sys
import os
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db

objtype=sys.argv[1]

spall=sdss_catalog.spall()
print(spall.dr)
print(spall.fitstablename)
spall.read()

# Star
#objtype='star'
if(objtype=='star'):
   df=pd.DataFrame(list(zip(spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
               spall.spectromaglist[:,0],spall.spectromagerrlist[:,0],\
               spall.spectromaglist[:,1],spall.spectromagerrlist[:,1],\
               spall.spectromaglist[:,2],spall.spectromagerrlist[:,2],\
               spall.spectromaglist[:,3],spall.spectromagerrlist[:,3],\
               spall.spectromaglist[:,4],spall.spectromagerrlist[:,4],\
               spall.airmasslist[:,0],\
               spall.airmasslist[:,1],\
               spall.airmasslist[:,2],\
               spall.airmasslist[:,3],\
               spall.airmasslist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,\
               spall.elodie_objectlist,\
               spall.elodie_sptypelist,spall.elodie_bvlist,spall.elodie_fehlist,\
               spall.elodie_tefflist,spall.elodie_logglist)),\
               columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'specmag_u','specmagerr_u',\
               'specmag_g','specmagerr_g',\
               'specmag_r','specmagerr_r',\
               'specmag_i','specmagerr_i',\
               'specmag_z','specmagerr_z',\
               'airmass_u',\
               'airmass_g',\
               'airmass_r',\
               'airmass_i',\
               'airmass_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','object','sptype','bv','feh','teff','logg'])
   print(df)
   dftmp=df[(df['class']=='STAR  ') & (df['thing_id']!=-1)]
   dfstar=dftmp.copy()
   update_class  =dftmp[:,'class'].str.strip()
   update_suclass=dftmp[:,'subclass'].str.strip()
   dfstar.loc[:,'class']=update_class
   dfstar.loc[:,'subclass']=update_subclass
   print('dfstar=',dfstar)
   print('subclass=',dfstar['subclass'])
   dfwd=dfstar[((dfstar['subclass']=='WDhotter') | \
               (dfstar['subclass']=='WDcoller') | \
               (dfstar['subclass']=='WDmagnetic') | \
               (dfstar['subclass']=='CalciumWD'))]
   dfspec=dfwd.sort_values(by=['teff'],ascending=False)
   dfspec.reset_index()

   dfstar2=dfstar[((dfstar['subclass']!='WDhotter') | \
               (dfstar['subclass']!='WDcoller') | \
               (dfstar['subclass']!='WDmagnetic') | \
               (dfstar['subclass']!='CalciumWD'))]
   dfspec2=dfstar2.sort_values(by=['teff'],ascending=False)
   dfspec2.reset_index()

   print('dfwd=',dfwd)
   print('dfspec=',dfspec)
   for i in range(100):
    plate=dfspec['plate'].iloc[i]
    mjd  =dfspec['mjd'].iloc[i]
    fiber=dfspec['fiber'].iloc[i]
    print(i,plate,mjd,fiber)
    print(i,plate,mjd,fiber,dfspec['subclass'].iloc[i],dfspec['teff'].iloc[i])

   fitsfilename='sdssDR17_wd.fits'
   sdss_db.create_2dspec(dfspec,fitsfilename,objtype)

   print('dfspec2=',dfspec2)
   for i in range(100,200):
    plate=dfspec2['plate'].iloc[i]
    mjd  =dfspec2['mjd'].iloc[i]
    fiber=dfspec2['fiber'].iloc[i]
    print(i,plate,mjd,fiber)
    print(i,plate,mjd,fiber,dfspec2['subclass'].iloc[i],dfspec2['teff'].iloc[i])

   fitsfilename='sdssDR17_star.fits'
   sdss_db.create_2dspec(dfspec2,fitsfilename,objtype)
   del dfwd ; del dfstar ; del df ; del dfspec ; del dfspec2

# Galaxy
#objtype='galaxy'
if(objtype=='galaxy'):
   df=pd.DataFrame(list(zip(spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.cmodelmaglist[:,0],spall.cmodelmagerrlist[:,0],\
               spall.cmodelmaglist[:,1],spall.cmodelmagerrlist[:,1],\
               spall.cmodelmaglist[:,2],spall.cmodelmagerrlist[:,2],\
               spall.cmodelmaglist[:,3],spall.cmodelmagerrlist[:,3],\
               spall.cmodelmaglist[:,4],spall.cmodelmagerrlist[:,4],\
               spall.spectromaglist[:,0],spall.spectromagerrlist[:,0],\
               spall.spectromaglist[:,1],spall.spectromagerrlist[:,1],\
               spall.spectromaglist[:,2],spall.spectromagerrlist[:,2],\
               spall.spectromaglist[:,3],spall.spectromagerrlist[:,3],\
               spall.spectromaglist[:,4],spall.spectromagerrlist[:,4],\
               spall.airmasslist[:,0],\
               spall.airmasslist[:,1],\
               spall.airmasslist[:,2],\
               spall.airmasslist[:,3],\
               spall.airmasslist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist)),\
               columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'cmodelmag_u','cmodelmagerr_u',\
               'cmodelmag_g','cmodelmagerr_g',\
               'cmodelmag_r','cmodelmagerr_r',\
               'cmodelmag_i','cmodelmagerr_i',\
               'cmodelmag_z','cmodelmagerr_z',\
               'specmag_u','specmagerr_u',\
               'specmag_g','specmagerr_g',\
               'specmag_r','specmagerr_r',\
               'specmag_i','specmagerr_i',\
               'specmag_z','specmagerr_z',\
               'airmass_u',\
               'airmass_g',\
               'airmass_r',\
               'airmass_i',\
               'airmass_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','z','zerr','zwarning'])
   print(df)
   df['class']=df['class'].str.strip()
   df['subclass']=df['subclass'].str.strip()
   dfgalaxy=df[(df['class']=='GALAXY') & (df['thing_id']!=-1)]
   #print(dfgalaxy)
   dfspec=dfgalaxy.sort_values(by=['z'],ascending=False)
   del dfgalaxy 

   fitsfilename='sdssDR17_galaxy.fits'
   sdss_db.create_2dspec(dfspec,fitsfilename,objtype)
   del dfspec

# Quasar
#objtype='quasar'
if(objtype=='quasar'):
   df=pd.DataFrame(list(zip(spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
               spall.spectromaglist[:,0],spall.spectromagerrlist[:,0],\
               spall.spectromaglist[:,1],spall.spectromagerrlist[:,1],\
               spall.spectromaglist[:,2],spall.spectromagerrlist[:,2],\
               spall.spectromaglist[:,3],spall.spectromagerrlist[:,3],\
               spall.spectromaglist[:,4],spall.spectromagerrlist[:,4],\
               spall.airmasslist[:,0],\
               spall.airmasslist[:,1],\
               spall.airmasslist[:,2],\
               spall.airmasslist[:,3],\
               spall.airmasslist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist)),\
               columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'specmag_u','specmagerr_u',\
               'specmag_g','specmagerr_g',\
               'specmag_r','specmagerr_r',\
               'specmag_i','specmagerr_i',\
               'specmag_z','specmagerr_z',\
               'airmass_u',\
               'airmass_g',\
               'airmass_r',\
               'airmass_i',\
               'airmass_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','z','zerr','zwarning'])
               #'lambdaeff','xfocal','yfocal','zoffset',\
               #'snall','object','sptype','bv','feh','teff','logg'])
   print(df)
   df['class']=df['class'].str.strip()
   df['subclass']=df['subclass'].str.strip()
   #dfquasar=df[(df['class']=='QSO   ') & (df['thing_id']!=-1)]
   dfquasar=df[(df['class']=='QSO') & (df['thing_id']!=-1)]
   #print(dfgalaxy)
   dfspec=dfquasar.sort_values(by=['z'],ascending=False)

   fitsfilename='sdssDR17_quasar.fits'
   sdss_db.create_2dspec(dfspec,fitsfilename,objtype)
   #del products_list ; del df ; del dfquasar
   del df ; del dfquasar
   del dfspec

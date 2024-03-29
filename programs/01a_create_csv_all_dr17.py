import fitsio
import sys
import os
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog

spall=sdss_catalog.spall()
print(spall.dr)
print(spall.fitstablename)
spall.read()

# Star
objtype='star'
if(objtype=='star'):
#   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
#               spall.ralist,spall.declist,\
#               spall.thing_idlist,\
#               spall.classlist,spall.subclasslist,\
#               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
#               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
#               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
#               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
#               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
#               spall.lambdaefflist,\
#               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
#               spall.sn_medianalllist,\
#               spall.elodie_objectlist,\
#               spall.elodie_sptypelist,spall.elodie_bvlist,spall.elodie_fehlist,\
#               spall.elodie_tefflist,spall.elodie_logglist]
#   df=pd.DataFrame(products_list).transpose()
#   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
#               'psfmag_u','psfmagerr_u',\
#               'psfmag_g','psfmagerr_g',\
#               'psfmag_r','psfmagerr_r',\
#               'psfmag_i','psfmagerr_i',\
#               'psfmag_z','psfmagerr_z',\
#               'lambdaeff','xfocal','yfocal','zoffset',\
#               'snall','object','sptype','bv','feh','teff','logg']
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
   dfstar=df[(df['class']=='STAR  ') & (df['thing_id']!=-1)]
   print(dfstar)
   dfstar.to_csv('../csvfiles/v5_13_2_spall_star.csv',index=False)
   #del products_list ; del df ; del dfstar
   del df ; del dfstar

# Galaxy
objtype='galaxy'
if(objtype=='galaxy'):
#   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
#               spall.ralist,spall.declist,\
#               spall.thing_idlist,\
#               spall.classlist,spall.subclasslist,\
#               spall.cmodelmaglist[:,0],spall.cmodelmagerrlist[:,0],\
#               spall.cmodelmaglist[:,1],spall.cmodelmagerrlist[:,1],\
#               spall.cmodelmaglist[:,2],spall.cmodelmagerrlist[:,2],\
#               spall.cmodelmaglist[:,3],spall.cmodelmagerrlist[:,3],\
#               spall.cmodelmaglist[:,4],spall.cmodelmagerrlist[:,4],\
#               spall.lambdaefflist,\
#               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
#               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist]
#   df=pd.DataFrame(products_list).transpose()
#
#   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
#               'cmodelmag_u','cmodelmagerr_u',\
#               'cmodelmag_g','cmodelmagerr_g',\
#               'cmodelmag_r','cmodelmagerr_r',\
#               'cmodelmag_i','cmodelmagerr_i',\
#               'cmodelmag_z','cmodelmagerr_z',\
#               'lambdaeff','xfocal','yfocal','zoffset',\
#               'snall','z','zerr','zwarning']
   df=pd.DataFrame(list(zip(spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.cmodelmaglist[:,0],spall.cmodelmagerrlist[:,0],\
               spall.cmodelmaglist[:,1],spall.cmodelmagerrlist[:,1],\
               spall.cmodelmaglist[:,2],spall.cmodelmagerrlist[:,2],\
               spall.cmodelmaglist[:,3],spall.cmodelmagerrlist[:,3],\
               spall.cmodelmaglist[:,4],spall.cmodelmagerrlist[:,4],\
               spall.fibermaglist[:,0],spall.fibermagerrlist[:,0],\
               spall.fibermaglist[:,1],spall.fibermagerrlist[:,1],\
               spall.fibermaglist[:,2],spall.fibermagerrlist[:,2],\
               spall.fibermaglist[:,3],spall.fibermagerrlist[:,3],\
               spall.fibermaglist[:,4],spall.fibermagerrlist[:,4],\
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
               spall.lambdaefflist,spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist)),\
               columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'cmodelmag_u','cmodelmagerr_u',\
               'cmodelmag_g','cmodelmagerr_g',\
               'cmodelmag_r','cmodelmagerr_r',\
               'cmodelmag_i','cmodelmagerr_i',\
               'cmodelmag_z','cmodelmagerr_z',\
               'fibermag_u','fibermagerr_u',\
               'fibermag_g','fibermagerr_g',\
               'fibermag_r','fibermagerr_r',\
               'fibermag_i','fibermagerr_i',\
               'fibermag_z','fibermagerr_z',\
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
   dfgalaxy=df[(df['class']=='GALAXY') & (df['thing_id']!=-1)]
   print(dfgalaxy)
   dfgalaxy.to_csv('../csvfiles/v5_13_2_spall_galaxy.csv',index=False)
   del df ; del dfgalaxy

# Quasar
objtype='quasar'
if(objtype=='quasar'):
#   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
#               spall.ralist,spall.declist,\
#               spall.thing_idlist,\
#               spall.classlist,spall.subclasslist,\
#               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
#               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
#               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
#               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
#               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
#               spall.lambdaefflist,\
#               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
#               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist]
#   df=pd.DataFrame(products_list).transpose()
#   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
#               'psfmag_u','psfmagerr_u',\
#               'psfmag_g','psfmagerr_g',\
#               'psfmag_r','psfmagerr_r',\
#               'psfmag_i','psfmagerr_i',\
#               'psfmag_z','psfmagerr_z',\
#               'lambdaeff','xfocal','yfocal','zoffset',\
#               'snall','z','zerr','zwarning']
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
               spall.lambdaefflist,spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
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
               #'lambdaeff','xfocal','yfocal','zoffset',
               #'snall','object','sptype','bv','feh','teff','logg'])
   print(df)
   dfquasar=df[(df['class']=='QSO   ') & (df['thing_id']!=-1)]
   dfquasar.to_csv('../csvfiles/v5_13_2_spall_quasar.csv',index=False)
   del df ; del dfquasar

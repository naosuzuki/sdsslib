import fitsio
import sys
import os
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog

# Reading DR8 specObj fitsfile
# Missing thing_id, psfmag etc
spobj=sdss_catalog.DR8specobj()
spobj.dr='DR9'
spobj.fitstablename='/Users/suzuki/sdssredux/specObj-SDSS-dr9.fits'
print(spobj.dr)
print(spobj.fitstablename)
spobj.read()

products_list=[spobj.platelist,spobj.mjdlist,spobj.fiberlist,\
               spobj.bestobjidlist,\
               spobj.classlist,spobj.subclasslist,\
               spobj.xfocallist,spobj.yfocallist,\
               spobj.zspzbestlist,spobj.zspzbesterrlist,spobj.zwarninglist,\
               spobj.elodie_objectlist,\
               spobj.elodie_sptypelist,spobj.elodie_bvlist,spobj.elodie_fehlist,\
               spobj.elodie_tefflist,spobj.elodie_logglist,\
               spobj.plug_ralist,spobj.plug_declist]

if(spobj.dr=='DR8'): products_list.append(spobj.sn_medianlist)
if(spobj.dr=='DR9'): products_list.append(spobj.sn_medianalllist)

df0=pd.DataFrame(products_list).transpose()
column_names=['plate','mjd','fiber','objid',\
            'class','subclass','xfocal','yfocal',\
            'z','zerr','zwarning',\
            'object','sptype','bv','feh','teff','logg','plug_ra','plug_dec']

if(spobj.dr=='DR8'): column_names.append('snmedian')
if(spobj.dr=='DR9'): column_names.append('snall')

df0.columns=column_names
print(df0)
# Remove Empty BESTOBJID : not it is in 19A character not an integer
df0a=df0[df0['objid']!='                   ']
df0b = df0a.astype({'objid':'int'})
if(spobj.dr=='DR8'):
   df0b.to_csv('../csvfiles/specObj-dr8.csv',index=False)
elif(spobj.dr=='DR9'):
   df0b.to_csv('../csvfiles/specObj-dr9.csv',index=False)

# Reading DR8 PhotoPlate
# Photometry Data : PSFMAG, CMODELMAG, THING_ID
photoplate=sdss_catalog.DR8photoplate()
print(photoplate.dr)
print(photoplate.fitstablename)
photoplate.read()

products_list1=[photoplate.ralist,photoplate.declist,\
               photoplate.objidlist,\
               photoplate.thing_idlist,\
               photoplate.psfmaglist[:,0],photoplate.psfmagerrlist[:,0],\
               photoplate.psfmaglist[:,1],photoplate.psfmagerrlist[:,1],\
               photoplate.psfmaglist[:,2],photoplate.psfmagerrlist[:,2],\
               photoplate.psfmaglist[:,3],photoplate.psfmagerrlist[:,3],\
               photoplate.psfmaglist[:,4],photoplate.psfmagerrlist[:,4],\
               photoplate.cmodelmaglist[:,0],photoplate.cmodelmagerrlist[:,0],\
               photoplate.cmodelmaglist[:,1],photoplate.cmodelmagerrlist[:,1],\
               photoplate.cmodelmaglist[:,2],photoplate.cmodelmagerrlist[:,2],\
               photoplate.cmodelmaglist[:,3],photoplate.cmodelmagerrlist[:,3],\
               photoplate.cmodelmaglist[:,4],photoplate.cmodelmagerrlist[:,4]]
df1=pd.DataFrame(products_list1).transpose()
df1.columns=['ra','dec','objid','thing_id',\
             'psfmag_u','psfmagerr_u',\
             'psfmag_g','psfmagerr_g',\
             'psfmag_r','psfmagerr_r',\
             'psfmag_i','psfmagerr_i',\
             'psfmag_z','psfmagerr_z',\
             'cmodelmag_u','cmodelmagerr_u',\
             'cmodelmag_g','cmodelmagerr_g',\
             'cmodelmag_r','cmodelmagerr_r',\
             'cmodelmag_i','cmodelmagerr_i',\
             'cmodelmag_z','cmodelmagerr_z']
# Drop Duplicates
df1a=df1.drop_duplicates(subset='thing_id',keep='last')
# Drop No OBJID rows
df1b=df1a[df1a['objid']!='                   ']
# Convert ID to an integer
df1c=df1b.astype({'objid':'int'})
df1d=df1c[df1c['thing_id']!='                   ']
df1e=df1d.astype({'thing_id':'int'})

df=pd.merge(df0b,df1e,how='left',on='objid')
if(spobj.dr=='DR8'):
   df.to_csv('../csvfiles/dr8_spall.csv',index=False)
elif(spobj.dr=='DR9'):
   df.to_csv('../csvfiles/dr9_spall.csv',index=False)

# Star
objtype='star'
if(objtype=='star'):
   output_columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'xfocal','yfocal',\
               'z','zerr','zwarning',\
               'object','sptype','bv','feh','teff','logg']
   if(spobj.dr=='DR8'): output_columns.append('snmedian')
   if(spobj.dr=='DR9'): output_columns.append('snall')
   print(df)
   dftmp=df[(df['class']=='STAR  ') & (df['thing_id']!=-1)]
   dfstar=dftmp.loc[:,output_columns]
   print(dfstar)
   if(spobj.dr=='DR8'):
      dfstar.to_csv('../csvfiles/dr8_spall_star.csv',index=False)
   if(spobj.dr=='DR9'):
      dfstar.to_csv('../csvfiles/dr9_spall_star.csv',index=False)
   del dfstar ;  del dftmp

# Galaxy
objtype='galaxy'
if(objtype=='galaxy'):
   output_columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'cmodelmag_u','cmodelmagerr_u',\
               'cmodelmag_g','cmodelmagerr_g',\
               'cmodelmag_r','cmodelmagerr_r',\
               'cmodelmag_i','cmodelmagerr_i',\
               'cmodelmag_z','cmodelmagerr_z',\
               'xfocal','yfocal',\
               'z','zerr','zwarning']
   if(spobj.dr=='DR8'): output_columns.append('snmedian')
   if(spobj.dr=='DR9'): output_columns.append('snall')
   print(df)
   dftmp=df[(df['class']=='GALAXY') & (df['thing_id']!=-1)]
   dfgalaxy=dftmp.loc[:,output_columns]
   print('galaxy')
   print(dfgalaxy)
   if(spobj.dr=='DR8'):
      dfgalaxy.to_csv('../csvfiles/dr8_spall_galaxy.csv',index=False)
   if(spobj.dr=='DR9'):
      dfgalaxy.to_csv('../csvfiles/dr9_spall_galaxy.csv',index=False)
   del dfgalaxy ; del dftmp

# Quasar
objtype='quasar'
if(objtype=='quasar'):
   output_columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'xfocal','yfocal',\
               'z','zerr','zwarning']
   if(spobj.dr=='DR8'): output_columns.append('snmedian')
   if(spobj.dr=='DR9'): output_columns.append('snall')
   print(df)
   dftmp=df[(df['class']=='QSO   ') & (df['thing_id']!=-1)]
   dfquasar=dftmp.loc[:,output_columns]
   print('quasar')
   print(dfquasar)
   if(spobj.dr=='DR8'):
      dfquasar.to_csv('../csvfiles/dr8_spall_quasar.csv',index=False)
   if(spobj.dr=='DR9'):
      dfquasar.to_csv('../csvfiles/dr9_spall_quasar.csv',index=False)
   del df ; del dfquasar ; del dftmp

sys.exit(1)

# Star
objtype='star'
if(objtype=='star'):
   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,\
               spall.elodie_objectlist,\
               spall.elodie_sptypelist,spall.elodie_bvlist,spall.elodie_fehlist,\
               spall.elodie_tefflist,spall.elodie_logglist]
   df=pd.DataFrame(products_list).transpose()
   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','object','sptype','bv','feh','teff','logg']
   print(df)
   dfstar=df[(df['class']=='STAR  ') & (df['thing_id']!=-1)]
   #dfstar=df[(df['class']=='STAR  ') & (df['thing_id']==-1)]
   print(dfstar)
   dfstar.to_csv('../csvfiles/v5_13_2_spall_star.csv',index=False)
   del products_list ; del df ; del dfstar

# Galaxy
objtype='galaxy'
if(objtype=='galaxy'):
   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.cmodelmaglist[:,0],spall.cmodelmagerrlist[:,0],\
               spall.cmodelmaglist[:,1],spall.cmodelmagerrlist[:,1],\
               spall.cmodelmaglist[:,2],spall.cmodelmagerrlist[:,2],\
               spall.cmodelmaglist[:,3],spall.cmodelmagerrlist[:,3],\
               spall.cmodelmaglist[:,4],spall.cmodelmagerrlist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist]
   df=pd.DataFrame(products_list).transpose()

   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'cmodelmag_u','cmodelmagerr_u',\
               'cmodelmag_g','cmodelmagerr_g',\
               'cmodelmag_r','cmodelmagerr_r',\
               'cmodelmag_i','cmodelmagerr_i',\
               'cmodelmag_z','cmodelmagerr_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','z','zerr','zwarning']
   print(df)
   dfgalaxy=df[(df['class']=='GALAXY') & (df['thing_id']!=-1)]
   print(dfgalaxy)
   dfgalaxy.to_csv('../csvfiles/v5_13_2_spall_galaxy.csv',index=False)
   del dfgalaxy

# Quasar
objtype='quasar'
if(objtype=='quasar'):
   products_list=[spall.platelist,spall.mjdlist,spall.fiberlist,\
               spall.ralist,spall.declist,\
               spall.thing_idlist,\
               spall.classlist,spall.subclasslist,\
               spall.psfmaglist[:,0],spall.psfmagerrlist[:,0],\
               spall.psfmaglist[:,1],spall.psfmagerrlist[:,1],\
               spall.psfmaglist[:,2],spall.psfmagerrlist[:,2],\
               spall.psfmaglist[:,3],spall.psfmagerrlist[:,3],\
               spall.psfmaglist[:,4],spall.psfmagerrlist[:,4],\
               spall.lambdaefflist,\
               spall.xfocallist,spall.yfocallist,spall.zoffsetlist,\
               spall.sn_medianalllist,spall.zspzbestlist,spall.zspzbesterrlist,spall.zwarninglist]
   df=pd.DataFrame(products_list).transpose()
   df.columns=['plate','mjd','fiber','ra','dec','thing_id','class','subclass',\
               'psfmag_u','psfmagerr_u',\
               'psfmag_g','psfmagerr_g',\
               'psfmag_r','psfmagerr_r',\
               'psfmag_i','psfmagerr_i',\
               'psfmag_z','psfmagerr_z',\
               'lambdaeff','xfocal','yfocal','zoffset',\
               'snall','z','zerr','zwarning']
   print(df)
   dfquasar=df[(df['class']=='QSO   ') & (df['thing_id']!=-1)]
   dfquasar.to_csv('../csvfiles/v5_13_2_spall_quasar.csv',index=False)
   del products_list ; del df ; del dfquasar

import fitsio
import sys
import os
import pandas as pd
pylibdir=os.environ['PYLIB']
sys.path.append(pylibdir)
import sdss_catalog
import sdss_db

def read_spec(csvfile):
   df=pd.read_csv(csvfile)
   #plate=11675
   plate=6783
   mjd=56284
   dftmp=df[(df['plate']==6783) & (df['mjd']==56284)]
   dftmp.reset_index()
   print(dftmp)
   print(len(dftmp))
   for i in range(len(dftmp)):
      mjd=dftmp['mjd'].iloc[i]
      fiber=dftmp['fiber'].iloc[i]
      source_id=dftmp['source_id'].iloc[i]
      print(plate,mjd,fiber,source_id)
      #spec=sdss_catalog.SDSSspec(plate,mjd,fiber)
      spec=sdss_db.SDSSspec(plate,mjd,fiber)
      spec.read()

csvfile='../../projects_gaia/data/gaiadr3_sdssdr17_star.csv'
#csvfile='../../projects_gaia/data/gaiadr3_sdssdr17_quasar.csv'
read_spec(csvfile)
#spPlate-11675-58523.fits

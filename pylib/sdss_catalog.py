import sys
import os
import string
import numpy
import scipy
import fitsio

class spall():
      def __init__(self):
          self.platelist=[]
          self.mjdlist=[]
          self.fiberlist=[]
          self.version='v5_13_2'
          self.dr='DR17'
          self.fitstablename=os.environ['SPECTRO_REDUX']+'spAll-'+self.version+'.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def readstar(self):
          columns=['PLATE','MJD','FIBERID','OBJTYPE','CLASS'=='STAR','SUBCLASS']
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']
          self.objtypelist=d['OBJTYPE']
          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']

      def readtest(self):
          columns=['OBJTYPE','CLASS','SUBCLASS','SPECPRIMARY']
          print('Reading ',self.fitstablename)
          print('Number of Spectra is ',len(self.nspec))
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.objtype=d['OBJTYPE']
          self.objclass=d['CLASS']
          self.subclass=d['SUBCLASS']
          self.specprimary=d['SPECPRIMARY']

      def read(self):
          columns=['PLATE','MJD','FIBERID',\
          'PLATEQUALITY','ZOFFSET','ZWARNING',\
          'LAMBDA_EFF','XFOCAL','YFOCAL','PLUG_RA','PLUG_DEC',\
          'OBJTYPE','CLASS','SUBCLASS','SPECPRIMARY',\
          'SPEC1_G','SPEC1_R','SPEC1_I',\
          'SPEC2_G','SPEC2_R','SPEC2_I',\
          'PSFMAG','PSFMAGERR','THING_ID',\
          'SPECTROFLUX','SPECTROFLUX_IVAR',\
          'SPECTROSYNFLUX','SPECTROSYNFLUX_IVAR','EXTINCTION',\
          'AIRMASS','OBJID','Z','Z_ERR']
          print('Reading ',self.fitstablename)
          print('Number of Rows is ',len(self.rows))
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          #d=fitsio.read(fitstablename,columns=columns)
          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']
          self.platequalitylist=d['PLATEQUALITY']
          self.zoffsetlist=d['ZOFFSET']
          self.zwarninglist=d['ZWARNING']

          self.lambdaefflist=d['LAMBDA_EFF']
          self.xfocallist=d['XFOCAL']
          self.yfocallist=d['YFOCAL']

          self.plug_ralist=d['PLUG_RA']
          self.plug_declist=d['PLUG_DEC']

          self.objtypelist=d['OBJTYPE']
          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']
          self.specprimarylist=d['SPECPRIMARY']

          self.spec1glist=d['SPEC1_G']
          self.spec1rlist=d['SPEC1_R']
          self.spec1ilist=d['SPEC1_I']

          self.spec2glist=d['SPEC2_G']
          self.spec2rlist=d['SPEC2_R']
          self.spec2ilist=d['SPEC2_I']

          self.psfmaglist=d['PSFMAG']
          self.psfmagerrlist=d['PSFMAGERR']
          self.thing_idlist=d['THING_ID']

          # DR8 and DR9 definition
          self.spectrofluxlist=d['SPECTROFLUX']
          self.spectroflux_ivarlist=d['SPECTROFLUX_IVAR']
          self.spectrosynfluxlist=d['SPECTROSYNFLUX']
          self.spectrosynflux_ivarlist=d['SPECTROSYNFLUX_IVAR']
          self.extinctionlist=d['EXTINCTION']

          self.airmasslist=d['AIRMASS']
          self.bestobjidlist=d['OBJID']
          self.zspzbestlist=d['Z']
          self.zspzbesterrlist=d['Z_ERR']
          print('Reading spAll FITS file',self.version,' is done')

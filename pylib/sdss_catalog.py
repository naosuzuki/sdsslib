import sys
import os
import string
import numpy
import scipy
import fitsio

class SDSSspec:
      def __init__(self,plate,mjd,fiber):
          self.dr='DR17'
          self.ver='v5_13_2'
          self.plate=int(plate)
          self.mjd=int(mjd)
          self.fiber=int(fiber)
          self.strplate=str(plate)
          self.strmjd=str(mjd)
          self.strfiber=str(fiber)
          self.radeg=0.0
          self.decdeg=0.0
          self.z=0.0
          self.zerr=0.0
          self.zsource='zpipe'
          self.zspzbest=0.0
          self.zspzbesterr=0.0
          self.airmass=1.0
          self.alt=90.0
          self.az=0.0
          self.exptime=0.0
          self.nexp=0
          self.thingid=0

      def read(self):
          # Define Data Directory
          self.fitsfilename=os.environ['SPECTRO_REDUX']+self.ver+'/'+\
          self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
          # Read FITS header
          h=fitsio.read_header(self.fitsfilename,ext=0)
          self.nspec=h['NAXIS2']
          self.coeff0=h['coeff0']
          self.coeff1=h['coeff1']
          self.npix=h['naxis1']
          self.airmass=h['airmass']
          #if(self.plate <=3509) : self.quality=hdulist[0].header['quality']
          self.exptime=h['exptime']
          self.nexp=h['nexp']/4
          self.idlspec2dver=h['VERSCOMB']
          self.idstart=int(h['coeff0']/0.0001)
          self.idend=self.idstart+self.npix-1
          fits=fitsio.FITS(self.fitsfilename)
          #self.flux=numpy.zeros(self.npix,numpy.float32)
          #self.ivar=numpy.zeros(self.npix,numpy.float32)
          #self.flux=fits[1][self.fiber-1,:]
          #self.ivar=fits[2][self.fiber-1,:]
          #self.mask=fits[3][self.fiber-1,:]
          self.flux=fits[1][551,:]
          self.ivar=fits[2][551,:]
          self.mask=fits[3][551,:]
          self.wave=10.0**(self.coeff0+self.coeff1*numpy.arange(self.npix))
          self.err=self.ivar
          self.wid=int(h['coeff0']/0.0001)+numpy.arange(self.npix)

      def write(self):
      # Jun 12, 2012 (Tue) 3pm
      # Write out ascii file
          asciidir=os.environ['SPECTRO_REDUX']+self.ver+'/ascii/'+self.strplate
          if(os.path.exists(asciidir)==False):
            print(asciidir,' directory does not exist')
            os.mkdir(asciidir,775)
            #os.chmod(asciidir,775)

          # Filter Out unnecessary files
          if(self.fiber!=1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-'+self.strfiber+'.dat'
          elif(self.fiber==1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-000.dat'
          ofile=open(asciifilename,'w')

          print(asciifilename)
          for i in range(len(self.wid)):
             ofile.write("%10.4f"%(self.wave[i])+\
                          "%15.5e"%(self.flux[i])+\
                          "%15.5e"%(self.ivar[i])+\
                          "%15i"%(self.mask[i])+\
                          "%15i"%(self.wid[i])+str("%2s"%('\n')))
          ofile.close()

class spall():
      def __init__(self):
          #self.platelist=[]
          #self.mjdlist=[]
          #self.fiberlist=[]
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

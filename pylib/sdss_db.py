import numpy
import fitsio
import os

class SDSSspec:
      def __init__(self,plate,mjd,fiber):
          self.dr='DR17'
          self.plate=plate
          self.mjd=mjd
          self.fiber=fiber
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
          self.bestobjid=0
          self.quality=''
          self.specsn2_g=0.0
          self.specsn2_r=0.0
          self.specsn2_i=0.0
          self.coeff0=0.0
          self.coeff1=0.0
          self.npix=0
          self.idstart=0
          self.idend=0
          self.idlspec2dver='v5_13_2'
          self.specclass=''
          self.specsubclass=''
          self.objtype=''

          self.mwebv=0.0
          self.mwebv2011=0.0

          self.psfmag_u=0.0
          self.psfmag_uerr=0.0
          self.psfmag_g=0.0
          self.psfmag_gerr=0.0
          self.psfmag_r=0.0
          self.psfmag_rerr=0.0
          self.psfmag_i=0.0
          self.psfmag_ierr=0.0
          self.psfmag_z=0.0
          self.psfmag_zerr=0.0

          self.cmodelmag_u=0.0
          self.cmodelmag_uerr=0.0
          self.cmodelmag_g=0.0
          self.cmodelmag_gerr=0.0
          self.cmodelmag_r=0.0
          self.cmodelmag_rerr=0.0
          self.cmodelmag_i=0.0
          self.cmodelmag_ierr=0.0
          self.cmodelmag_z=0.0
          self.cmodelmag_zerr=0.0

          self.spectroflux_g=0.0
          self.spectroflux_ivar_g=0.0
          self.spectroflux_r=0.0
          self.spectroflux_ivar_r=0.0
          self.spectroflux_i=0.0
          self.spectroflux_ivar_i=0.0

          self.spectrosynflux_g=0.0
          self.spectrosynflux_ivar_g=0.0
          self.spectrosynflux_r=0.0
          self.spectrosynflux_ivar_r=0.0
          self.spectrosynflux_i=0.0
          self.spectrosynflux_ivar_i=0.0

          self.specmag_g=0.0
          self.specmag_gerr=0.0
          self.specmag_r=0.0
          self.specmag_rerr=0.0
          self.specmag_i=0.0
          self.specmag_ierr=0.0

          self.modelspecmag_g=0.0
          self.modelspecmag_gerr=0.0
          self.modelspecmag_r=0.0
          self.modelspecmag_rerr=0.0
          self.modelspecmag_i=0.0
          self.modelspecmag_ierr=0.0

          self.specprimary=0
          self.zoffset=0.0
          self.zwarning=0
          self.xfocal=0.0
          self.yfocal=0.0
          self.seeing50=0.0
          self.rmsoff50=0.0
          self.lambdaeff=0.0

      def read(self):
          # Define Data Directory
          #h=fitsio.read_header(self.fitstablename,ext=1)

          self.fitsfilename=os.environ['SPECTRO_REDUX']+'/'+self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
          h=fitsio.read_header(self.fitsfilename,ext=0)
          self.coeff0=h['coeff0']
          self.coeff1=h['coeff1']
          self.npix=h['naxis1']
          self.nspec=h['NAXIS2']
          self.airmass=h['airmass']
          #if(self.plate <=3509) : self.quality=hdulist[0].header['quality']
          self.exptime=h['exptime']
          self.nexp=h['nexp']/4
          self.idlspec2dver=h['VERSCOMB']
          self.idstart=int(h['coeff0']/0.0001)
          self.idend=self.idstart+self.npix-1
          self.wave=10.0**(self.coeff0+self.coeff1*numpy.arange(self.npix))

#          self.flux=hdata[self.fiber-1,]
#          self.ivar=hdulist[1].data[self.fiber-1,]
#          self.mask=hdulist[2].data[self.fiber-1,]
#          self.err=self.ivar
#          self.wid=int(hdulist[0].header['coeff0']/0.0001)+numpy.arange(self.npix)

class spPlate():

      def readimage(self):
          #self.datadir()
          spplatename=self.rootdir+'/'+self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'

          #hdulist=pyfits.open(spplatename)
          h=fitsio.read_header(self.fitstablename,ext=1)
          print(hdulist.info())
          cols=hdulist[5].columns
          print(cols.info())
          self.coeff0=hdulist[0].header['coeff0']
          self.coeff1=hdulist[0].header['coeff1']
          self.npix=hdulist[0].header['naxis1']
          self.exptime=hdulist[0].header['exptime']
          self.dateobs=hdulist[0].header['DATE-OBS']
          self.nexp=hdulist[0].header['nexp']/4
          self.idlspec2dver=hdulist[0].header['VERSCOMB']
          self.idstart=int(hdulist[0].header['coeff0']/0.0001)
          self.idend=hdulist[0].header['coeff0']/0.0001+hdulist[0].header['naxis1']-1
          self.airmass=hdulist[0].header['AIRMASS']
          self.plate=hdulist[0].header['PLATEID']
          self.mjd=hdulist[0].header['MJD']
          self.seeing50=hdulist[0].header['SEEING50']
          self.rmsoff50=hdulist[0].header['RMSOFF50']
          # 2D part
          self.wave=10.0**(self.coeff0+self.coeff1*numpy.arange(self.npix))
          self.imageflux=hdulist[0].data
          self.imageivar=hdulist[1].data
          self.imagemask=hdulist[2].data
          self.wid=int(hdulist[0].header['coeff0']/0.0001)+numpy.arange(self.npix)
          # 5th extension
          tbdata=hdulist[5].data
          self.objtypelist=tbdata.field('OBJTYPE')
          self.zoffsetlist=tbdata.field('ZOFFSET')
          self.lambdaefflist=tbdata.field('LAMBDA_EFF')
          self.fiberidlist=tbdata.field('FIBERID')
          self.xfocallist=tbdata.field('XFOCAL')
          self.yfocallist=tbdata.field('YFOCAL')
          self.objidlist=tbdata.field('OBJID')
          self.ralist=tbdata.field('RA')
          self.declist=tbdata.field('DEC')
          self.xyradiuslist=numpy.sqrt(self.xfocallist**2+self.yfocallist**2)
          self.maglist=tbdata.field('MAG')
          hdulist.close()

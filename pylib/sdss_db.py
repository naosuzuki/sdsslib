import sys
import os
import string
import numpy
import scipy
import fitsio
from astropy.io import fits

## Written by Nao Suzuki

# 2022-08-13 (Sat) Create 2D Fits is added
# 2022-07-16 (Sat) DR8 is added
# 2022-07-03 (Sun) Update

class SDSSspec:
      def __init__(self,plate,mjd,fiber):
          if(plate>3500): 
             self.dr='DR17'
             self.ver='v5_13_2'
          elif(plate<3500): 
             self.dr='DR8'
             self.ver='DR8'
          self.plate=int(plate)
          self.mjd=int(mjd)
          self.fiber=int(fiber)
          self.strmjd=str("%05i"%(self.mjd))
          self.strplate=str("%04i"%(self.plate))
          self.strfiber=str("%03i"%(self.fiber))
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
      #2022-04-22 Fri
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

          # Reading FITS Image
          self.fits=fitsio.FITS(self.fitsfilename)
          self.flux=self.fits[0][self.fiber-1,:][0]
          self.ivar=self.fits[1][self.fiber-1,:][0]
          self.mask=self.fits[2][self.fiber-1,:][0]
          self.err=self.ivar
          # Wavelength Definition
          self.wave=10.0**(self.coeff0+self.coeff1*numpy.arange(self.npix))
          self.wave_plus =10.0**(self.coeff0+self.coeff1*(numpy.arange(self.npix)+0.5))
          self.wave_minus=10.0**(self.coeff0+self.coeff1*(numpy.arange(self.npix)-0.5))
          self.dwave=self.wave_plus-self.wave_minus

          self.wid=int(h['coeff0']/0.0001)+numpy.arange(self.npix)

          # Plate Properties
          self.lambdaeff=5400.
          self.xfocal=0.0
          self.yfocal=0.0
          self.zoffset=0.0

          self.sdssclass='STAR'
          self.subclass=''

      def photopoints20(self):
          # Find 20 photometric points from GAIA XP spectrum
          # We have 251 spectral points = 11pix*1 + 12pix*19 = 20 points
          self.ptx=numpy.zeros(20)
          self.pty=numpy.zeros(20)
          self.ptyerr=numpy.zeros(20)
          self.ptmask=numpy.zeros(20,dtype=numpy.int32)
          self.ptwid=numpy.zeros(20,dtype=numpy.int32)

          # S/N > 5 criterion for good flux value
          spec_mask=numpy.where(self.flux*numpy.sqrt(self.ivar)>5.0,1.0,0)
          # Flux area = flux * dwave
          spec_maskedflux=self.flux*spec_mask*self.dwave
          spec_maskedwave=self.dwave*spec_mask

          wid0=35286        ; wid1=35646
          wid0_end=wid1-120 ; wid1_start=wid1-119 ; wid1_end=wid1+120
          # Calculate the first point
          self.ptx[0]=10.0**(self.coeff1*float(wid0)) ; self.ptwid[0]=wid0

          # Define the first pixel
          if(self.wid[0]<wid0_end):
            nbin=wid0_end-self.wid[0]+1
            if(numpy.sum(spec_mask[0:nbin])>0): 
               self.pty[0]=numpy.sum(spec_maskedflux[0:nbin])/numpy.sum(spec_maskedwave[0:nbin])
               self.ptmask[0]=1

          # Define the pixel 1-19
          for j in range(1,20):
             self.ptwid[j]=wid1+240*(j-1)
             self.ptx[j]=10.0**(self.coeff1*float(self.ptwid[j]))
             wid_start=self.ptwid[j]-119
             wid_end  =self.ptwid[j]+120

             if((self.wid[0] < wid_end) and (self.wid[-1] > wid_start)):
                if(self.wid[0] >= wid_start):
                   nstart=0 ;  nend=wid_end-self.wid[0]+1
                elif((self.wid[0] < wid_start) and (self.wid[-1] >= wid_end)):
                   nstart=wid_start-self.wid[0] ; nend=wid_end-self.wid[0]+1
                elif((self.wid[-1] >= wid_start) and (self.wid[-1] < wid_end)):
                   nstart=wid_start-self.wid[0] ; nend=self.wid[-1]-self.wid[0]+1

                if(numpy.sum(spec_mask[nstart:nend])>0): 
                   self.pty[j]=numpy.sum(spec_maskedflux[nstart:nend])/numpy.sum(spec_maskedwave[nstart:nend])
                   self.ptmask[j]=1
                   self.ptyerr[j]=self.pty[j]/numpy.average(self.flux[nstart:nend]*numpy.sqrt(self.ivar[nstart:nend]))
        
      def write(self):
      # Jun 12, 2012 (Tue) 3pm : Apr 22, 2022 (Fri) 22:52
      # Write out ascii file
          asciidir=os.environ['SPECTRO_REDUX']+self.ver+'/ascii/'+self.strplate
          if(os.path.exists(asciidir)==False):
            os.mkdir(asciidir)

          # Filter Out unnecessary files
          if(self.fiber!=1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-'+self.strfiber+'.dat'
          elif(self.fiber==1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-000.dat'
          ofile=open(asciifilename,'w')

          #print(asciifilename)
          for i in range(len(self.wid)):
             ofile.write("%10.4f"%(self.wave[i])+\
                          "%15.5e"%(self.flux[i])+\
                          "%15.5e"%(self.ivar[i])+\
                          "%15i"%(self.mask[i])+\
                          "%15i"%(self.wid[i])+str("%2s"%('\n')))
          ofile.close()

      def normalize_at7000(self):
          w6900id=38389
          w7000id=38451
          w7100id=38513

          nflux=self.flux[(w6900id-self.wid[0]):(w7100id-self.wid[0])]
          nivar=self.ivar[(w6900id-self.wid[0]):(w7100id-self.wid[0])]
          nflag=numpy.ones(w7100id-w6900id+1,dtype=numpy.int)
          # Flag Good Data Points
          nflag=numpy.where(nivar>0.0,1,0)
          # Flag S/N > 1 :  numpy.where (condition, yes, no)
          nflag=numpy.where(nflux*numpy.sqrt(nivar) > 1.0,1,0)
          # Extract Good Data Points : extract from nflux with nflag
          nflux_good=numpy.compress(nflag,nflux)

          # Taking Median : 
          if(len(nflux_good) > 5):
             nfactor=numpy.median(nflux_good)
             nfactor_std=numpy.std(nflux_good)
          elif(len(nflux) > 1):
             nfactor=numpy.median(nflux)
             nfactor_std=numpy.std(nflux)
          else:
              nfactor=1.0
              nfactor_std=0.0

          self.flux/=nfactor

def create_2dspec(df,fitsfilename,objtype):
# 2022-08-12 LBNL
# Create 2D FITS File from DataFrame

   startID=35486
   endID  =40177
   npixall=endID-startID+1
   # Define Arrays
   imageflux=numpy.zeros((len(df),npixall),dtype=numpy.float32)
   imageivar=numpy.zeros((len(df),npixall),dtype=numpy.float32)
   imagemask=numpy.zeros((len(df),npixall),dtype=numpy.int32)

   plate_list =df['plate'].to_numpy()
   mjd_list   =df['mjd'].to_numpy()
   fiber_list   =df['fiber'].to_numpy()

   ra_list =df['ra'].to_numpy()
   dec_list=df['dec'].to_numpy()
   thing_id_list=df['thing_id'].to_numpy()

   psfmag_u_list=df['psfmag_u'].to_numpy()
   psfmagerr_u_list=df['psfmagerr_u'].to_numpy()
   psfmag_g_list=df['psfmag_g'].to_numpy()
   psfmagerr_g_list=df['psfmagerr_g'].to_numpy()
   psfmag_r_list=df['psfmag_r'].to_numpy()
   psfmagerr_r_list=df['psfmagerr_r'].to_numpy()
   psfmag_i_list=df['psfmag_i'].to_numpy()
   psfmagerr_i_list=df['psfmagerr_i'].to_numpy()
   psfmag_z_list=df['psfmag_z'].to_numpy()
   psfmagerr_z_list=df['psfmagerr_z'].to_numpy()
   snall_list   =df['snall'].to_numpy()

   if(objtype=='star'):
      object_list=df['object'].to_numpy()
      sptype_list=df['sptype'].to_numpy()
      bv_list    =df['bv'].to_numpy()
      feh_list   =df['feh'].to_numpy()
      teff_list  =df['teff'].to_numpy()
      logg_list  =df['logg'].to_numpy()

   #for i in range(len(dfwd)):
   for i in range(100):
      plate=df['plate'].iloc[i]
      mjd  =df['mjd'].iloc[i]
      fiber=df['fiber'].iloc[i]
   
      # Define Spectrum
      spec=SDSSspec(plate,mjd,fiber)
      spec.read()
      #spec.normalize_at7000()
      #spec.read_MWcorrected()

      # Define the spectrum range : first pixel and the last pixel
      if(spec.wid[0] > startID):
         kstart=0
         jstart=spec.wid[0]-startID
      else:
         kstart=spec.wid[0]-startID
         jstart=0

      if(spec.wid[-1] < endID):
         kend=spec.npix-1
         jend=spec.wid[-1]-startID
      else:
         kend=endID-spec.wid[0]
         jend=npixall-1

      #print('reading',i,jstart,jend)
      print('reading',i,plate,mjd,fiber)
      imageflux[i,jstart:jend]=spec.flux[kstart:kend]
      imageivar[i,jstart:jend]=spec.ivar[kstart:kend]
      imagemask[i,jstart:jend]=spec.mask[kstart:kend]


   hdu1=fits.PrimaryHDU(imageflux)
   hdu2=fits.ImageHDU(imageivar)
   hdu3=fits.ImageHDU(imagemask)

   col1=fits.Column(name='RA',format='E',array=ra_list)
   col2=fits.Column(name='DEC',format='E',array=dec_list)
   col3=fits.Column(name='thing_id',format='K',array=thing_id_list)
   col4=fits.Column(name='SNR',format='E',array=snall_list)

   col5=fits.Column(name='psfmag_u',format='E',array=psfmag_u_list)
   col6=fits.Column(name='psfmagerr_u',format='E',array=psfmagerr_u_list)
   col7=fits.Column(name='psfmag_g',format='E',array=psfmag_g_list)
   col8=fits.Column(name='psfmagerr_g',format='E',array=psfmagerr_g_list)
   col9=fits.Column(name='psfmag_r',format='E',array=psfmag_r_list)
   col10=fits.Column(name='psfmagerr_r',format='E',array=psfmagerr_r_list)
   col11=fits.Column(name='psfmag_i',format='E',array=psfmag_i_list)
   col12=fits.Column(name='psfmagerr_i',format='E',array=psfmagerr_i_list)
   col13=fits.Column(name='psfmag_z',format='E',array=psfmag_z_list)
   col14=fits.Column(name='psfmagerr_z',format='E',array=psfmagerr_z_list)

   cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
   tbhdu=fits.BinTableHDU.from_columns(cols)

   fitsio.write(fitsfilename,imageflux)




